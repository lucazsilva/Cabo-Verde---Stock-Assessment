#-------------------------------------------------------------------------------------------#
#   Este script contem toda analise baseada em comprimento de Decapterus macarellus         #
#   Analises possuem foco para o Decapterus macarellus no arquipelago de Cabo Verde         #
# Avaliacao descritiva da estrutura de comprimentos ao longo do tempo, meses, sexos e zonas #
# Analise da variacao temporal e espacial das distribuicoes de frequencia de comprimento    #
# Diagnosticos da cobertura amostral, tamanho de amostra e representatividade dos dados     #
# Estimacao da relacao peso-comprimento e avaliacao de possiveis diferencas entre sexos     #
# Analise da maturidade sexual e estimacao de L50 e L95 a partir de ogivas logisticas       #
# Estagios maturacionais considerados: I = imaturo; II-VII = individuos adultos/maduros     #
# Avaliacao da adequacao dos dados para modelos baseados em comprimento, incluindo LBSPR,   #
# LBB, LBI/aLBI e analises de frequencia de comprimento com TropFishR/ELEFAN                #
# Comparacao entre composicoes de comprimento brutas e padronizadas por evento e por mes    #
# com foco na identificacao de efeitos do desenho amostral sobre a estrutura de comprimentos#
# Codificacao criada por Silva, LVS ; 09/09/2026, Instituto do Mar - IMar, Mindelo          #
#-------------------------------------------------------------------------------------------#
# 0. PACOTES ------------------------------------------------------

pacotes <- c(
  "tidyverse", "lubridate", "janitor", "skimr", "ggridges",
  "patchwork", "broom", "mgcv", "scales"
)

instalar <- pacotes[!pacotes %in% rownames(installed.packages())]
if (length(instalar) > 0) install.packages(instalar)

invisible(lapply(pacotes, library, character.only = TRUE))

# 1. CAMINHOS -----------------------------------------------------
arquivo <- "Length_oficial.csv"
if (!file.exists(arquivo)) arquivo <- file.choose()

pasta_saida <- "resultados_length"
dir.create(pasta_saida, showWarnings = FALSE)

salvar_plot <- function(plot, nome, largura = 10, altura = 7) {
  ggsave(
    filename = file.path(pasta_saida, nome),
    plot = plot,
    width = largura,
    height = altura,
    dpi = 300
  )
}

# 2. IMPORTACAO E LIMPEZA ----------------------------------------
dados_raw <- readr::read_csv(
  arquivo,
  na = c("", "NA", "N/A", "NULL"),
  show_col_types = FALSE,
  trim_ws = TRUE
) %>%
  janitor::remove_empty("cols") %>%
  janitor::clean_names()

# Depois de clean_names(), espera-se:
# data, ano, mes, zona, l_cm, w_g, sexo, em, w_gon, w_fig, obs
print(names(dados_raw))

dados <- dados_raw %>%
  mutate(
    zona_original = zona,
    sexo_original = sexo,

    # Datas do arquivo estao em mes/dia/ano
    data = lubridate::mdy(data),

    # Padronizacao conservadora das zonas: corrige caixa e espacos,
    # mas NAO une automaticamente nomes que possam representar locais diferentes.
    zona = zona %>%
      stringr::str_squish() %>%
      stringr::str_to_upper(),

    # Corrige espacos extras como " M" e " F"
    sexo = sexo %>%
      stringr::str_squish() %>%
      stringr::str_to_upper(),

    # Evita que zona ausente quebre o identificador de amostra
    zona_evento = tidyr::replace_na(zona, "SEM_ZONA"),

    # Um evento = uma combinacao data x zona
    sample_id = interaction(data, zona_evento, drop = TRUE, lex.order = TRUE)
  )

# 3. CONTROLE DE QUALIDADE ---------------------------------------

# 3.1 Dimensoes e periodo
cat("\nNumero de individuos:", nrow(dados), "\n")
cat("Primeira data:", as.character(min(dados$data, na.rm = TRUE)), "\n")
cat("Ultima data:", as.character(max(dados$data, na.rm = TRUE)), "\n")
cat("Anos:", min(dados$ano, na.rm = TRUE), "-", max(dados$ano, na.rm = TRUE), "\n")
cat("Eventos amostrais:", n_distinct(dados$sample_id), "\n")

# 3.2 Missing values
missing <- dados %>%
  summarise(across(everything(), ~ sum(is.na(.x)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variavel",
    values_to = "n_missing"
  ) %>%
  mutate(
    pct_missing = 100 * n_missing / nrow(dados)
  ) %>%
  arrange(desc(pct_missing))

print(missing)
write_csv(missing, file.path(pasta_saida, "01_missing_values.csv"))

# 3.3 Checagem da coerencia DATA x ANO x MES
inconsistencia_data <- dados %>%
  filter(
    !is.na(data) &
      (ano != lubridate::year(data) | mes != lubridate::month(data))
  )

cat("Inconsistencias DATA x ANO/MES:", nrow(inconsistencia_data), "\n")
write_csv(inconsistencia_data, file.path(pasta_saida, "02_inconsistencias_data.csv"))

# 3.4 Duplicatas exatas
# ATENCAO: nao removemos automaticamente. Sem ID individual, dois peixes podem
# legitimamente apresentar as mesmas medidas.
duplicatas_exatas <- dados_raw %>%
  group_by(across(everything())) %>%
  filter(n() > 1) %>%
  ungroup()

cat("Linhas envolvidas em grupos de duplicatas exatas:", nrow(duplicatas_exatas), "\n")
write_csv(duplicatas_exatas, file.path(pasta_saida, "03_duplicatas_para_inspecao.csv"))

# 3.5 Valores impossiveis ou suspeitos sem impor limites biologicos arbitrarios
qc_valores <- tibble(
  teste = c(
    "comprimento <= 0",
    "peso <= 0",
    "peso gonada < 0",
    "peso figado < 0",
    "mes fora de 1-12",
    "ano fora do intervalo da base"
  ),
  n = c(
    sum(dados$l_cm <= 0, na.rm = TRUE),
    sum(dados$w_g <= 0, na.rm = TRUE),
    sum(dados$w_gon < 0, na.rm = TRUE),
    sum(dados$w_fig < 0, na.rm = TRUE),
    sum(!dados$mes %in% 1:12, na.rm = TRUE),
    sum(dados$ano < 2004 | dados$ano > 2024, na.rm = TRUE)
  )
)
print(qc_valores)
write_csv(qc_valores, file.path(pasta_saida, "04_qc_valores.csv"))

# 3.6 Auditoria dos nomes de zona
# Mostra variantes que foram corrigidas apenas por caixa/espacos.
auditoria_zona <- dados %>%
  count(zona_original, zona, sort = TRUE)
print(auditoria_zona, n = Inf)
write_csv(auditoria_zona, file.path(pasta_saida, "05_auditoria_zonas.csv"))

# 3.7 Categorias de sexo e estagio
print(dados %>% count(sexo, sort = TRUE))
print(dados %>% count(em, sort = TRUE))
print(dados %>% count(obs, sort = TRUE))

# 3.8 Resolucao dos comprimentos
comprimentos_unicos <- sort(unique(dados$l_cm[!is.na(dados$l_cm)]))
resolucao_min <- if (length(comprimentos_unicos) > 1) {
  min(diff(comprimentos_unicos))
} else {
  NA_real_
}
cat("Resolucao minima observada do comprimento:", resolucao_min, "cm\n")

# 4. DESCRICAO GERAL ---------------------------------------------

desc_geral <- dados %>%
  summarise(
    n = n(),
    n_eventos = n_distinct(sample_id),
    media_L = mean(l_cm, na.rm = TRUE),
    dp_L = sd(l_cm, na.rm = TRUE),
    cv_L = sd(l_cm, na.rm = TRUE) / mean(l_cm, na.rm = TRUE),
    min_L = min(l_cm, na.rm = TRUE),
    q01_L = quantile(l_cm, 0.01, na.rm = TRUE),
    q05_L = quantile(l_cm, 0.05, na.rm = TRUE),
    q25_L = quantile(l_cm, 0.25, na.rm = TRUE),
    mediana_L = median(l_cm, na.rm = TRUE),
    q75_L = quantile(l_cm, 0.75, na.rm = TRUE),
    q95_L = quantile(l_cm, 0.95, na.rm = TRUE),
    q99_L = quantile(l_cm, 0.99, na.rm = TRUE),
    max_L = max(l_cm, na.rm = TRUE),
    media_W = mean(w_g, na.rm = TRUE),
    mediana_W = median(w_g, na.rm = TRUE),
    min_W = min(w_g, na.rm = TRUE),
    max_W = max(w_g, na.rm = TRUE)
  )

print(desc_geral)
write_csv(desc_geral, file.path(pasta_saida, "06_descritiva_geral.csv"))

# Funcao de resumo para grupos
resumo_comprimento <- function(data, ...) {
  data %>%
    group_by(...) %>%
    summarise(
      n = sum(!is.na(l_cm)),
      n_eventos = n_distinct(sample_id),
      media = mean(l_cm, na.rm = TRUE),
      dp = sd(l_cm, na.rm = TRUE),
      cv = dp / media,
      q05 = quantile(l_cm, 0.05, na.rm = TRUE),
      q25 = quantile(l_cm, 0.25, na.rm = TRUE),
      mediana = median(l_cm, na.rm = TRUE),
      q75 = quantile(l_cm, 0.75, na.rm = TRUE),
      q95 = quantile(l_cm, 0.95, na.rm = TRUE),
      minimo = min(l_cm, na.rm = TRUE),
      maximo = max(l_cm, na.rm = TRUE),
      .groups = "drop"
    )
}

desc_ano <- dados %>%
  group_by(ano) %>%
  summarise(
    n = n(),
    n_eventos = n_distinct(sample_id),
    n_datas = n_distinct(data),
    n_meses = n_distinct(mes),
    n_zonas = n_distinct(zona, na.rm = TRUE),
    media = mean(l_cm, na.rm = TRUE),
    dp = sd(l_cm, na.rm = TRUE),
    q05 = quantile(l_cm, 0.05, na.rm = TRUE),
    q25 = quantile(l_cm, 0.25, na.rm = TRUE),
    mediana = median(l_cm, na.rm = TRUE),
    q75 = quantile(l_cm, 0.75, na.rm = TRUE),
    q95 = quantile(l_cm, 0.95, na.rm = TRUE),
    minimo = min(l_cm, na.rm = TRUE),
    maximo = max(l_cm, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(desc_ano, file.path(pasta_saida, "07_descritiva_por_ano.csv"))
write_csv(resumo_comprimento(dados, mes), file.path(pasta_saida, "08_descritiva_por_mes.csv"))
write_csv(resumo_comprimento(dados, sexo), file.path(pasta_saida, "09_descritiva_por_sexo.csv"))
write_csv(resumo_comprimento(dados, zona), file.path(pasta_saida, "10_descritiva_por_zona.csv"))

# 5. ESTRUTURA DA AMOSTRAGEM -------------------------------------

# 5.1 Tamanho por evento
resumo_eventos <- dados %>%
  group_by(sample_id, data, ano, mes, zona_evento) %>%
  summarise(
    n = n(),
    media_L = mean(l_cm, na.rm = TRUE),
    mediana_L = median(l_cm, na.rm = TRUE),
    sd_L = sd(l_cm, na.rm = TRUE),
    min_L = min(l_cm, na.rm = TRUE),
    max_L = max(l_cm, na.rm = TRUE),
    .groups = "drop"
  )

print(summary(resumo_eventos$n))
write_csv(resumo_eventos, file.path(pasta_saida, "11_resumo_eventos.csv"))

p_event_n <- ggplot(resumo_eventos, aes(n)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
  labs(
    x = "Numero de peixes por evento",
    y = "Numero de eventos",
    title = "Tamanho das amostras por evento"
  ) +
  theme_bw()

salvar_plot(p_event_n, "01_tamanho_amostral_eventos.png", 9, 6)

# 5.2 Cobertura ano x mes
cobertura_mes <- dados %>%
  count(ano, mes, name = "n") %>%
  complete(ano = min(dados$ano):max(dados$ano), mes = 1:12, fill = list(n = 0))

p_cobertura <- ggplot(cobertura_mes, aes(x = factor(mes), y = factor(ano), fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    x = "Mes",
    y = "Ano",
    fill = "n",
    title = "Cobertura temporal da amostragem"
  ) +
  theme_bw()

salvar_plot(p_cobertura, "02_cobertura_ano_mes.png", 10, 8)
write_csv(cobertura_mes, file.path(pasta_saida, "12_cobertura_ano_mes.csv"))

# 5.3 Quantos anos cada mes aparece?
cobertura_meses_serie <- dados %>%
  distinct(ano, mes) %>%
  count(mes, name = "n_anos") %>%
  mutate(prop_anos = n_anos / n_distinct(dados$ano))

write_csv(cobertura_meses_serie, file.path(pasta_saida, "13_cobertura_meses_serie.csv"))
print(cobertura_meses_serie)

# Meses presentes em pelo menos 80% dos anos: diagnostico de sensibilidade
meses_core <- cobertura_meses_serie %>%
  filter(prop_anos >= 0.80) %>%
  pull(mes)
cat("Meses presentes em >=80% dos anos:", paste(meses_core, collapse = ", "), "\n")

# 5.4 Cobertura espacial: zonas mais frequentes
zonas_top <- dados %>%
  filter(!is.na(zona)) %>%
  count(zona, sort = TRUE) %>%
  slice_head(n = 12) %>%
  pull(zona)

cobertura_zona <- dados %>%
  filter(zona %in% zonas_top) %>%
  count(ano, zona, name = "n") %>%
  complete(ano = min(dados$ano):max(dados$ano), zona = zonas_top, fill = list(n = 0))

p_zona_ano <- ggplot(cobertura_zona, aes(x = factor(ano), y = forcats::fct_rev(zona), fill = n)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    x = "Ano",
    y = "Zona",
    fill = "n",
    title = "Cobertura das principais zonas ao longo do tempo"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

salvar_plot(p_zona_ano, "03_cobertura_zona_ano.png", 13, 7)

# 6. DISTRIBUICAO DE COMPRIMENTOS --------------------------------

# 6.1 Distribuicao geral
p_hist <- ggplot(dados, aes(l_cm)) +
  geom_histogram(binwidth = 1, boundary = 0.5) +
  labs(
    x = "Comprimento (cm)",
    y = "Frequencia",
    title = "Distribuicao geral de comprimentos"
  ) +
  theme_bw()

salvar_plot(p_hist, "04_histograma_comprimento_geral.png", 9, 6)

# 6.2 Distribuicao por ano - ridgeline
p_ridge_ano <- ggplot(dados, aes(x = l_cm, y = factor(ano), fill = factor(ano))) +
  ggridges::geom_density_ridges(
    scale = 1.4,
    rel_min_height = 0.01,
    alpha = 0.7,
    show.legend = FALSE
  ) +
  labs(
    x = "Comprimento (cm)",
    y = "Ano",
    title = "Variacao anual da distribuicao de comprimentos"
  ) +
  theme_bw()

salvar_plot(p_ridge_ano, "05_ridgeline_comprimento_ano.png", 10, 10)

# 6.3 Frequencia relativa bruta por ano
freq_ano_raw <- dados %>%
  count(ano, l_cm, name = "n") %>%
  group_by(ano) %>%
  mutate(prop_raw = n / sum(n)) %>%
  ungroup()

p_heat_raw <- ggplot(freq_ano_raw, aes(x = l_cm, y = factor(ano), fill = prop_raw)) +
  geom_tile() +
  scale_fill_viridis_c(labels = scales::percent) +
  labs(
    x = "Comprimento (cm)",
    y = "Ano",
    fill = "Proporcao",
    title = "Composicao anual de comprimentos - dados brutos"
  ) +
  theme_bw()

salvar_plot(p_heat_raw, "06_heatmap_comprimento_ano_bruto.png", 11, 8)

# 6.4 Quantis anuais: mais informativo que olhar apenas a media
p_quantis <- ggplot(desc_ano, aes(x = ano)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.18) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.28) +
  geom_line(aes(y = mediana), linewidth = 0.9) +
  geom_point(aes(y = mediana), size = 2) +
  labs(
    x = "Ano",
    y = "Comprimento (cm)",
    title = "Mudancas na distribuicao: mediana, IQR e intervalo 5-95%",
    subtitle = "Faixa clara = P5-P95; faixa central = P25-P75"
  ) +
  theme_bw()

salvar_plot(p_quantis, "07_quantis_comprimento_ano.png", 10, 6)

# 6.5 Distribuicao por mes
p_mes <- ggplot(dados, aes(x = factor(mes), y = l_cm)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(
    x = "Mes",
    y = "Comprimento (cm)",
    title = "Variacao sazonal do comprimento"
  ) +
  theme_bw()

salvar_plot(p_mes, "08_comprimento_por_mes.png", 10, 6)

# 6.6 Distribuicao por sexo
p_sexo <- dados %>%
  filter(!is.na(sexo)) %>%
  ggplot(aes(x = sexo, y = l_cm)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(
    x = "Sexo",
    y = "Comprimento (cm)",
    title = "Distribuicao de comprimentos por sexo"
  ) +
  theme_bw()

salvar_plot(p_sexo, "09_comprimento_por_sexo.png", 8, 6)

# 6.7 Principais zonas
p_zona <- dados %>%
  filter(zona %in% zonas_top) %>%
  mutate(zona = forcats::fct_reorder(zona, l_cm, median, na.rm = TRUE)) %>%
  ggplot(aes(x = zona, y = l_cm)) +
  geom_boxplot(outlier.alpha = 0.2) +
  coord_flip() +
  labs(
    x = "Zona",
    y = "Comprimento (cm)",
    title = "Variacao de comprimento entre as principais zonas"
  ) +
  theme_bw()

salvar_plot(p_zona, "10_comprimento_principais_zonas.png", 9, 7)

# 7. PADRONIZACAO DA COMPOSICAO DE COMPRIMENTOS ------------------
# Motivo: anos com mais eventos em determinados meses/zonas podem aparentar
# uma mudanca na estrutura de tamanhos que e, na verdade, mudanca no desenho amostral.

# Usa classes inteiras de 1 cm porque essa e a resolucao observada no banco.
grade_comprimento <- seq(
  floor(min(dados$l_cm, na.rm = TRUE)),
  ceiling(max(dados$l_cm, na.rm = TRUE)),
  by = 1
)

meta_evento <- dados %>%
  distinct(sample_id, data, ano, mes, zona_evento)

freq_evento <- dados %>%
  count(sample_id, l_cm, name = "n") %>%
  tidyr::complete(sample_id, l_cm = grade_comprimento, fill = list(n = 0)) %>%
  left_join(meta_evento, by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(prop_evento = n / sum(n)) %>%
  ungroup()

# 7.1 Cada EVENTO recebe o mesmo peso dentro do ano
freq_ano_evento <- freq_evento %>%
  group_by(ano, l_cm) %>%
  summarise(
    prop_event_std = mean(prop_evento),
    .groups = "drop"
  )

p_heat_event <- ggplot(freq_ano_evento, aes(x = l_cm, y = factor(ano), fill = prop_event_std)) +
  geom_tile() +
  scale_fill_viridis_c(labels = scales::percent) +
  labs(
    x = "Comprimento (cm)",
    y = "Ano",
    fill = "Proporcao",
    title = "Composicao anual padronizada por evento"
  ) +
  theme_bw()

salvar_plot(p_heat_event, "11_heatmap_comprimento_evento_padronizado.png", 11, 8)

# 7.2 Cada MES AMOSTRADO recebe o mesmo peso dentro do ano
# OBS.: isso reduz o peso de meses muito amostrados, mas NAO corrige meses totalmente ausentes.
freq_mes_evento <- freq_evento %>%
  group_by(ano, mes, l_cm) %>%
  summarise(
    prop_mes = mean(prop_evento),
    .groups = "drop"
  )

freq_ano_mes_equal <- freq_mes_evento %>%
  group_by(ano, l_cm) %>%
  summarise(
    prop_month_equal = mean(prop_mes),
    .groups = "drop"
  )

p_heat_month <- ggplot(freq_ano_mes_equal, aes(x = l_cm, y = factor(ano), fill = prop_month_equal)) +
  geom_tile() +
  scale_fill_viridis_c(labels = scales::percent) +
  labs(
    x = "Comprimento (cm)",
    y = "Ano",
    fill = "Proporcao",
    title = "Composicao anual com igual peso para cada mes amostrado"
  ) +
  theme_bw()

salvar_plot(p_heat_month, "12_heatmap_comprimento_mes_padronizado.png", 11, 8)

# 7.3 Quantifica quanto a composicao muda apos padronizacao
comparacao_composicao <- freq_ano_raw %>%
  select(ano, l_cm, prop_raw) %>%
  full_join(freq_ano_evento, by = c("ano", "l_cm")) %>%
  full_join(freq_ano_mes_equal, by = c("ano", "l_cm")) %>%
  replace_na(list(prop_raw = 0, prop_event_std = 0, prop_month_equal = 0))

# Distancia de variacao total (0 = identicas; 1 = completamente distintas)
diagnostico_padronizacao <- comparacao_composicao %>%
  group_by(ano) %>%
  summarise(
    TV_raw_vs_event = 0.5 * sum(abs(prop_raw - prop_event_std)),
    TV_raw_vs_month = 0.5 * sum(abs(prop_raw - prop_month_equal)),
    .groups = "drop"
  )

print(diagnostico_padronizacao)
write_csv(diagnostico_padronizacao, file.path(pasta_saida, "14_efeito_padronizacao_composicao.csv"))

p_tv <- diagnostico_padronizacao %>%
  pivot_longer(-ano, names_to = "comparacao", values_to = "TV") %>%
  ggplot(aes(x = ano, y = TV, linetype = comparacao)) +
  geom_line(linewidth = 0.9) +
  geom_point() +
  labs(
    x = "Ano",
    y = "Distancia de variacao total",
    linetype = NULL,
    title = "Quanto o desenho amostral altera a composicao anual?",
    subtitle = "Valores proximos de zero indicam pouca diferenca entre composicoes"
  ) +
  theme_bw()

salvar_plot(p_tv, "13_efeito_padronizacao.png", 10, 6)

# 8. DIAGNOSTICO DAS CAUDAS DA DISTRIBUICAO ----------------------
# Modelos de comprimento dependem de representar bem peixes pequenos e grandes.
q05_global <- quantile(dados$l_cm, 0.05, na.rm = TRUE)
q95_global <- quantile(dados$l_cm, 0.95, na.rm = TRUE)

caudas_ano <- dados %>%
  group_by(ano) %>%
  summarise(
    n = n(),
    prop_abaixo_q05_global = mean(l_cm <= q05_global, na.rm = TRUE),
    prop_acima_q95_global = mean(l_cm >= q95_global, na.rm = TRUE),
    p01 = quantile(l_cm, 0.01, na.rm = TRUE),
    p05 = quantile(l_cm, 0.05, na.rm = TRUE),
    p95 = quantile(l_cm, 0.95, na.rm = TRUE),
    p99 = quantile(l_cm, 0.99, na.rm = TRUE),
    amplitude_90 = p95 - p05,
    .groups = "drop"
  )

write_csv(caudas_ano, file.path(pasta_saida, "15_diagnostico_caudas_por_ano.csv"))

p_caudas <- caudas_ano %>%
  select(ano, prop_abaixo_q05_global, prop_acima_q95_global) %>%
  pivot_longer(-ano, names_to = "cauda", values_to = "proporcao") %>%
  ggplot(aes(x = ano, y = proporcao, linetype = cauda)) +
  geom_line(linewidth = 0.9) +
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Ano",
    y = "Proporcao",
    linetype = NULL,
    title = "Representacao das caudas da distribuicao de comprimento"
  ) +
  theme_bw()

salvar_plot(p_caudas, "14_representacao_caudas.png", 10, 6)

# 9. TENDENCIA TEMPORAL AJUSTADA POR MES E ZONA ------------------
# Este GAM e um DIAGNOSTICO de confusao do desenho amostral.
# Ele pergunta: a tendencia media de comprimento ao longo dos anos permanece
# depois de controlar por sazonalidade e diferencas entre zonas?

dados_gam <- dados %>%
  filter(
    !is.na(l_cm), !is.na(ano), !is.na(mes), !is.na(zona_evento)
  ) %>%
  mutate(zona_evento = factor(zona_evento))

modelo_gam <- mgcv::gam(
  l_cm ~
    s(ano, k = 8) +
    s(mes, bs = "cc", k = 8) +
    s(zona_evento, bs = "re"),
  data = dados_gam,
  method = "REML",
  knots = list(mes = c(0.5, 12.5))
)

capture.output(summary(modelo_gam), file = file.path(pasta_saida, "16_resumo_GAM.txt"))

# Predicao anual padronizada para todos os 12 meses.
# O efeito aleatorio de zona e excluido da predicao.
zona_ref <- levels(dados_gam$zona_evento)[1]
grade_pred <- tidyr::expand_grid(
  ano = sort(unique(dados_gam$ano)),
  mes = 1:12,
  zona_evento = factor(zona_ref, levels = levels(dados_gam$zona_evento))
)

grade_pred$pred <- predict(
  modelo_gam,
  newdata = grade_pred,
  type = "response",
  exclude = "s(zona_evento)"
)

media_ajustada <- grade_pred %>%
  group_by(ano) %>%
  summarise(media_ajustada = mean(pred), .groups = "drop") %>%
  left_join(desc_ano %>% select(ano, media_bruta = media), by = "ano")

p_gam <- media_ajustada %>%
  pivot_longer(
    cols = c(media_bruta, media_ajustada),
    names_to = "serie",
    values_to = "comprimento"
  ) %>%
  ggplot(aes(x = ano, y = comprimento, linetype = serie)) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    x = "Ano",
    y = "Comprimento medio (cm)",
    linetype = NULL,
    title = "Comprimento medio bruto vs. ajustado por mes e zona"
  ) +
  theme_bw()

salvar_plot(p_gam, "15_media_bruta_vs_ajustada_GAM.png", 10, 6)
write_csv(media_ajustada, file.path(pasta_saida, "17_media_bruta_vs_ajustada.csv"))

# 10. COMPRIMENTO-PESO -------------------------------------------
# Util para QC biologico, conversao comprimento-peso e verificacao de outliers.

dados_lw <- dados %>%
  filter(l_cm > 0, w_g > 0)

modelo_lw <- lm(log(w_g) ~ log(l_cm), data = dados_lw)
print(summary(modelo_lw))

coef_lw <- broom::tidy(modelo_lw)
b <- coef(modelo_lw)[["log(l_cm)"]]
a <- exp(coef(modelo_lw)[["(Intercept)"]])
r2 <- summary(modelo_lw)$r.squared

cat("Relacao W = a * L^b\n")
cat("a =", a, "\n")
cat("b =", b, "\n")
cat("R2 log-log =", r2, "\n")

write_csv(
  tibble(a = a, b = b, r2_loglog = r2),
  file.path(pasta_saida, "18_parametros_comprimento_peso.csv")
)

p_lw <- ggplot(dados_lw, aes(x = l_cm, y = w_g)) +
  geom_point(alpha = 0.08, size = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Comprimento (cm, escala log)",
    y = "Peso (g, escala log)",
    title = "Relacao comprimento-peso"
  ) +
  theme_bw()

salvar_plot(p_lw, "16_relacao_comprimento_peso.png", 8, 6)

# 10.1 Comprimento-peso por sexo
modelo_lw_sexo <- dados_lw %>%
  filter(sexo %in% c("F", "M")) %>%
  lm(log(w_g) ~ log(l_cm) * sexo, data = .)

capture.output(summary(modelo_lw_sexo), file = file.path(pasta_saida, "19_comprimento_peso_por_sexo.txt"))

# 11. MATURIDADE E ESTIMATIVA DE L50/L95 --------------------------
# Classificacao utilizada no banco (Fontana & Le Guen, 1969):
# I   = immature
# II  = resting
# III = maturing
# IV  = pre-spawning
# V   = spawning
# VI  = spent
# VII = gonads regression
#
# Para a ogiva principal, adotamos o limiar no estagio III:
#   NAO MADURO/INATIVO = I-II
#   MADURO/ATIVO       = III-VII
#
# ATENCAO BIOLOGICA:
# O estagio II (resting) pode incluir adultos em repouso reprodutivo. Por isso,
# uma ogiva calculada com todos os meses pode superestimar L50 se muitos adultos
# grandes estiverem em repouso. O script:
#   (a) estima III-VII como maduros (analise principal);
#   (b) faz sensibilidade tratando II-VII como maduros;
#   (c) permite restringir a ogiva aos meses da estacao reprodutiva.
#
# Para LBSPR, em geral use a ogiva biologicamente mais defensavel (L50 e L95)
# e documente explicitamente sexo, meses usados e criterio de maturidade.

# 11.1 Rotulos dos estagios
dados <- dados %>%
  mutate(
    em = as.integer(em),
    estagio_maturacional = factor(
      em,
      levels = 1:7,
      labels = c(
        "I - Immature",
        "II - Resting",
        "III - Maturing",
        "IV - Pre-spawning",
        "V - Spawning",
        "VI - Spent",
        "VII - Gonads regression"
      )
    ),
    maduro_principal = case_when(
      em %in% 1:2 ~ 0L,
      em %in% 3:7 ~ 1L,
      TRUE ~ NA_integer_
    ),
    # Apenas para sensibilidade: resting passa a ser considerado individuo
    # que ja atingiu maturidade sexual em algum momento.
    maduro_inclui_resting = case_when(
      em == 1 ~ 0L,
      em %in% 2:7 ~ 1L,
      TRUE ~ NA_integer_
    ),
    gsi = if_else(!is.na(w_gon) & !is.na(w_g) & w_g > 0,
                  100 * w_gon / w_g,
                  NA_real_)
  )

# 11.2 Composicao mensal dos estagios maturacionais
maturidade_mes <- dados %>%
  filter(sexo %in% c("F", "M"), !is.na(estagio_maturacional)) %>%
  count(sexo, mes, estagio_maturacional, name = "n") %>%
  group_by(sexo, mes) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write_csv(maturidade_mes, file.path(pasta_saida, "20_maturidade_por_mes_e_sexo.csv"))

p_estagios_mes <- ggplot(
  maturidade_mes,
  aes(x = factor(mes), y = prop, fill = estagio_maturacional)
) +
  geom_col() +
  facet_wrap(~ sexo, ncol = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Mes",
    y = "Proporcao",
    fill = "Estagio",
    title = "Variacao mensal dos estagios maturacionais por sexo"
  ) +
  theme_bw()

salvar_plot(p_estagios_mes, "17_estagios_maturacionais_por_mes.png", 11, 8)

# 11.3 GSI mensal como diagnostico auxiliar da estacao reprodutiva
# O GSI e usado aqui apenas para identificar a sazonalidade reprodutiva;
# nao e necessario para estimar a ogiva logistica.
gsi_mes <- dados %>%
  filter(sexo %in% c("F", "M"), is.finite(gsi)) %>%
  group_by(sexo, mes) %>%
  summarise(
    n = n(),
    media_gsi = mean(gsi, na.rm = TRUE),
    mediana_gsi = median(gsi, na.rm = TRUE),
    q25_gsi = quantile(gsi, 0.25, na.rm = TRUE),
    q75_gsi = quantile(gsi, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(gsi_mes, file.path(pasta_saida, "21_GSI_mensal_por_sexo.csv"))

p_gsi_mes <- ggplot(gsi_mes, aes(x = mes, y = mediana_gsi, linetype = sexo)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_ribbon(
    aes(ymin = q25_gsi, ymax = q75_gsi, group = sexo),
    alpha = 0.12,
    inherit.aes = TRUE
  ) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    x = "Mes",
    y = "GSI (%)",
    linetype = "Sexo",
    title = "Sazonalidade do indice gonadossomatico",
    subtitle = "Linha = mediana; faixa = intervalo interquartilico"
  ) +
  theme_bw()

salvar_plot(p_gsi_mes, "18_GSI_mensal.png", 10, 6)

# 11.4 Funcao para estimar ogiva, L50 e L95
# A incerteza e obtida por bootstrap de EVENTOS AMOSTRAIS, preservando a
# dependencia entre peixes coletados no mesmo evento. Isso e preferivel a
# tratar os ~23 mil peixes como observacoes totalmente independentes.

estimar_ogiva <- function(data,
                          sexo_alvo,
                          mature_stages = 2:7,
                          meses = NULL,
                          n_boot = 1000,
                          seed = 123,
                          min_por_classe_plot = 10) {

  d <- data %>%
    filter(
      sexo == sexo_alvo,
      !is.na(l_cm),
      !is.na(em),
      !is.na(sample_id)
    )

  if (!is.null(meses)) {
    d <- d %>% filter(mes %in% meses)
  }

  d <- d %>%
    mutate(maduro = as.integer(em %in% mature_stages))

  if (nrow(d) == 0 || n_distinct(d$maduro) < 2) {
    warning("Nao ha dados suficientes para ajustar a ogiva para sexo = ", sexo_alvo)
    return(NULL)
  }

  # Agregacao por classe de comprimento para o ajuste binomial.
  # O ponto estimado e equivalente ao ajuste individual, mas mais eficiente.
  agg <- d %>%
    group_by(l_cm) %>%
    summarise(
      n = n(),
      n_maduros = sum(maduro),
      n_nao_maduros = n - n_maduros,
      prop_maduros = n_maduros / n,
      .groups = "drop"
    )

  modelo <- glm(
    cbind(n_maduros, n_nao_maduros) ~ l_cm,
    data = agg,
    family = binomial(link = "logit")
  )

  b0 <- unname(coef(modelo)[1])
  b1 <- unname(coef(modelo)[2])

  calc_Lp <- function(p, b0, b1) {
    (qlogis(p) - b0) / b1
  }

  L50 <- calc_Lp(0.50, b0, b1)
  L95 <- calc_Lp(0.95, b0, b1)

  # ---------------- BOOTSTRAP POR EVENTO ----------------
  # Resume primeiro cada evento x classe para tornar o bootstrap rapido.
  eventos_len <- d %>%
    group_by(sample_id, l_cm) %>%
    summarise(
      n_maduros = sum(maduro),
      n_nao_maduros = n() - n_maduros,
      .groups = "drop"
    )

  ids <- unique(eventos_len$sample_id)
  set.seed(seed)

  boot_res <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    pesos <- tibble(
      sample_id = sample(ids, size = length(ids), replace = TRUE)
    ) %>%
      count(sample_id, name = "peso_boot")

    boot_agg <- eventos_len %>%
      inner_join(pesos, by = "sample_id") %>%
      group_by(l_cm) %>%
      summarise(
        n_maduros = sum(n_maduros * peso_boot),
        n_nao_maduros = sum(n_nao_maduros * peso_boot),
        .groups = "drop"
      )

    fit_b <- try(
      suppressWarnings(
        glm(
          cbind(n_maduros, n_nao_maduros) ~ l_cm,
          data = boot_agg,
          family = binomial(link = "logit")
        )
      ),
      silent = TRUE
    )

    if (!inherits(fit_b, "try-error") && all(is.finite(coef(fit_b)))) {
      bb0 <- unname(coef(fit_b)[1])
      bb1 <- unname(coef(fit_b)[2])

      # Descarta ajustes biologicamente invertidos/degenerados.
      if (bb1 > 0) {
        boot_res[[i]] <- tibble(
          L50 = calc_Lp(0.50, bb0, bb1),
          L95 = calc_Lp(0.95, bb0, bb1)
        )
      }
    }
  }

  boot <- bind_rows(boot_res) %>%
    filter(is.finite(L50), is.finite(L95))

  if (nrow(boot) >= 50) {
    ci <- boot %>%
      summarise(
        L50_low = quantile(L50, 0.025, na.rm = TRUE),
        L50_high = quantile(L50, 0.975, na.rm = TRUE),
        L95_low = quantile(L95, 0.025, na.rm = TRUE),
        L95_high = quantile(L95, 0.975, na.rm = TRUE)
      )
  } else {
    ci <- tibble(
      L50_low = NA_real_, L50_high = NA_real_,
      L95_low = NA_real_, L95_high = NA_real_
    )
  }

  # Curva ajustada
  grade <- tibble(
    l_cm = seq(min(d$l_cm), max(d$l_cm), by = 0.1)
  )
  grade$prob_maduro <- predict(modelo, newdata = grade, type = "response")

  # Classes observadas: nao mostramos classes com n muito pequeno para evitar
  # dar peso visual excessivo a proporcoes 0/1 baseadas em poucos individuos.
  obs_plot <- agg %>% filter(n >= min_por_classe_plot)

  resultado <- tibble(
    sexo = sexo_alvo,
    mature_stages = paste(mature_stages, collapse = "-"),
    meses = if (is.null(meses)) "todos" else paste(meses, collapse = ","),
    n_individuos = nrow(d),
    n_eventos = n_distinct(d$sample_id),
    n_boot_validos = nrow(boot),
    intercepto = b0,
    slope = b1,
    L50 = L50,
    L50_low = ci$L50_low,
    L50_high = ci$L50_high,
    L95 = L95,
    L95_low = ci$L95_low,
    L95_high = ci$L95_high,
    AIC = AIC(modelo)
  )

  list(
    dados = d,
    observado = agg,
    modelo = modelo,
    grade = grade,
    bootstrap = boot,
    resultado = resultado,
    obs_plot = obs_plot
  )
}

# 11.5 ANALISE PRINCIPAL: III-VII = maduros -----------------------
# Estimamos femeas e machos separadamente. Para parametros de maturidade de
# estoque, normalmente a ogiva de femeas e a mais relevante para spawning output.

ogiva_F <- estimar_ogiva(
  dados,
  sexo_alvo = "F",
  mature_stages = 2:7,
  meses = NULL,
  n_boot = 1000
)

ogiva_M <- estimar_ogiva(
  dados,
  sexo_alvo = "M",
  mature_stages = 2:7,
  meses = NULL,
  n_boot = 1000
)

resultados_L50_principal <- bind_rows(
  ogiva_F$resultado,
  ogiva_M$resultado
)

print(resultados_L50_principal)
write_csv(
  resultados_L50_principal,
  file.path(pasta_saida, "22_L50_L95_principal_II_VII.csv")
)

# Funcao de grafico da ogiva
plot_ogiva <- function(obj, titulo) {
  ggplot() +
    geom_point(
      data = obj$obs_plot,
      aes(x = l_cm, y = prop_maduros, size = n),
      alpha = 0.75
    ) +
    geom_line(
      data = obj$grade,
      aes(x = l_cm, y = prob_maduro),
      linewidth = 1
    ) +
    geom_hline(yintercept = 0.5, linetype = 3) +
    geom_hline(yintercept = 0.95, linetype = 3) +
    geom_vline(xintercept = obj$resultado$L50, linetype = 2) +
    geom_vline(xintercept = obj$resultado$L95, linetype = 2) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent
    ) +
    scale_size_continuous(name = "n na classe") +
    labs(
      x = "Comprimento (cm)",
      y = "Proporcao madura",
      title = titulo,
      subtitle = paste0(
        "L50 = ", round(obj$resultado$L50, 2), " cm (IC95% ",
        round(obj$resultado$L50_low, 2), "-", round(obj$resultado$L50_high, 2),
        "); L95 = ", round(obj$resultado$L95, 2), " cm"
      )
    ) +
    theme_bw()
}

p_L50_F <- plot_ogiva(
  ogiva_F,
  "Ogiva de maturidade - femeas (III-VII = maduras)"
)

p_L50_M <- plot_ogiva(
  ogiva_M,
  "Ogiva de maturidade - machos (III-VII = maduros)"
)

salvar_plot(p_L50_F, "19_ogiva_L50_femeas.png", 9, 6)
salvar_plot(p_L50_M, "20_ogiva_L50_machos.png", 9, 6)

# 11.6 ANALISE DE SENSIBILIDADE: II-VII = maduros -----------------
# NAO e a definicao principal. Serve para mostrar o impacto de tratar individuos
# em resting como adultos que ja atingiram a maturidade sexual.

ogiva_F_resting <- estimar_ogiva(
  dados,
  sexo_alvo = "F",
  mature_stages = 2:7,
  meses = NULL,
  n_boot = 1000
)

ogiva_M_resting <- estimar_ogiva(
  dados,
  sexo_alvo = "M",
  mature_stages = 2:7,
  meses = NULL,
  n_boot = 1000
)

resultados_L50_sens <- bind_rows(
  ogiva_F_resting$resultado,
  ogiva_M_resting$resultado
)

write_csv(
  resultados_L50_sens,
  file.path(pasta_saida, "23_L50_L95_sensibilidade_incluindo_resting.csv")
)

# Compara a decisao sobre o estagio II
comparacao_L50_criterio <- bind_rows(
  resultados_L50_principal %>% mutate(criterio = "III-VII maduros"),
  resultados_L50_sens %>% mutate(criterio = "II-VII maduros")
) %>%
  select(criterio, sexo, n_individuos, n_eventos,
         L50, L50_low, L50_high, L95, L95_low, L95_high)

print(comparacao_L50_criterio)
write_csv(
  comparacao_L50_criterio,
  file.path(pasta_saida, "24_comparacao_criterio_maturidade.csv")
)

# 11.7 OPCIONAL: L50 SOMENTE NA ESTACAO REPRODUTIVA ---------------
# Recomenda-se definir estes meses APOS inspecionar:
#   17_estagios_maturacionais_por_mes.png
#   18_GSI_mensal.png
#
# Exemplo (apenas ilustrativo, NAO ativado automaticamente):
# meses_reproducao_L50 <- c(3, 4, 5, 6)

meses_reproducao_L50 <- NULL

if (!is.null(meses_reproducao_L50)) {

  ogiva_F_repro <- estimar_ogiva(
    dados,
    sexo_alvo = "F",
    mature_stages = 2:7,
    meses = meses_reproducao_L50,
    n_boot = 1000
  )

  ogiva_M_repro <- estimar_ogiva(
    dados,
    sexo_alvo = "M",
    mature_stages = 2:7,
    meses = meses_reproducao_L50,
    n_boot = 1000
  )

  resultados_L50_repro <- bind_rows(
    ogiva_F_repro$resultado,
    ogiva_M_repro$resultado
  )

  print(resultados_L50_repro)
  write_csv(
    resultados_L50_repro,
    file.path(pasta_saida, "25_L50_L95_estacao_reprodutiva.csv")
  )

  salvar_plot(
    plot_ogiva(
      ogiva_F_repro,
      paste0(
        "Ogiva de maturidade - femeas; meses ",
        paste(meses_reproducao_L50, collapse = ", ")
      )
    ),
    "21_ogiva_L50_femeas_estacao_reprodutiva.png",
    9, 6
  )

  salvar_plot(
    plot_ogiva(
      ogiva_M_repro,
      paste0(
        "Ogiva de maturidade - machos; meses ",
        paste(meses_reproducao_L50, collapse = ", ")
      )
    ),
    "22_ogiva_L50_machos_estacao_reprodutiva.png",
    9, 6
  )
}

# 11.8 Estabilidade temporal exploratoria do L50 ------------------
# Esta analise NAO substitui a ogiva principal. Ela ajuda a verificar se o L50
# aparente muda muito entre anos, o que pode sinalizar sazonalidade, mudanca do
# desenho amostral ou mudanca biologica. Exige n razoavel e ambas as categorias.

estimar_L50_rapido <- function(d) {
  d <- d %>%
    filter(!is.na(l_cm), !is.na(em)) %>%
    mutate(maduro = as.integer(em %in% 3:7))

  if (nrow(d) < 150 || n_distinct(d$maduro) < 2 ||
      sum(d$maduro == 0) < 20 || sum(d$maduro == 1) < 20) {
    return(tibble(n = nrow(d), L50 = NA_real_, L95 = NA_real_))
  }

  fit <- try(glm(maduro ~ l_cm, data = d, family = binomial()), silent = TRUE)

  if (inherits(fit, "try-error") || coef(fit)[2] <= 0) {
    return(tibble(n = nrow(d), L50 = NA_real_, L95 = NA_real_))
  }

  b0 <- unname(coef(fit)[1])
  b1 <- unname(coef(fit)[2])

  tibble(
    n = nrow(d),
    L50 = -b0 / b1,
    L95 = (qlogis(0.95) - b0) / b1
  )
}

L50_ano <- dados %>%
  filter(sexo %in% c("F", "M")) %>%
  group_by(sexo, ano) %>%
  group_modify(~ estimar_L50_rapido(.x)) %>%
  ungroup()

write_csv(L50_ano, file.path(pasta_saida, "26_L50_exploratorio_por_ano.csv"))

p_L50_ano <- L50_ano %>%
  filter(!is.na(L50)) %>%
  ggplot(aes(x = ano, y = L50, linetype = sexo)) +
  geom_line(linewidth = 0.8) +
  geom_point() +
  labs(
    x = "Ano",
    y = "L50 aparente (cm)",
    linetype = "Sexo",
    title = "Estabilidade temporal exploratoria do L50",
    subtitle = "Interpretar com cautela: cobertura mensal e espacial varia entre anos"
  ) +
  theme_bw()

salvar_plot(p_L50_ano, "23_L50_exploratorio_por_ano.png", 10, 6)

# 12. MATRIZES PARA MODELOS BASEADOS EM COMPRIMENTO ---------------

# 12.1 Contagens anuais brutas: ano x classe de comprimento
matriz_ano_contagem <- dados %>%
  count(ano, l_cm, name = "n") %>%
  complete(ano, l_cm = grade_comprimento, fill = list(n = 0)) %>%
  pivot_wider(names_from = l_cm, values_from = n, names_prefix = "L_")

write_csv(matriz_ano_contagem, file.path(pasta_saida, "21_matriz_comprimento_anual_contagens.csv"))

# 12.2 Proporcoes anuais brutas
matriz_ano_prop <- freq_ano_raw %>%
  select(ano, l_cm, prop_raw) %>%
  complete(ano, l_cm = grade_comprimento, fill = list(prop_raw = 0)) %>%
  pivot_wider(names_from = l_cm, values_from = prop_raw, names_prefix = "L_")

write_csv(matriz_ano_prop, file.path(pasta_saida, "22_matriz_comprimento_anual_proporcoes.csv"))

# 12.3 Proporcoes padronizadas por evento
matriz_ano_event_std <- freq_ano_evento %>%
  pivot_wider(names_from = l_cm, values_from = prop_event_std, names_prefix = "L_")

write_csv(matriz_ano_event_std, file.path(pasta_saida, "23_matriz_comprimento_evento_padronizado.csv"))

# 12.4 LFQ mensal para exploracao de coortes/crescimento
lfq_mensal <- dados %>%
  count(ano, mes, l_cm, name = "n") %>%
  mutate(periodo = sprintf("%04d-%02d", ano, mes)) %>%
  select(periodo, l_cm, n) %>%
  complete(periodo, l_cm = grade_comprimento, fill = list(n = 0))

write_csv(lfq_mensal, file.path(pasta_saida, "24_LFQ_mensal_long.csv"))

# Heatmap mensal de LFQ: cada coluna = um periodo amostrado
p_lfq <- lfq_mensal %>%
  group_by(periodo) %>%
  mutate(prop = n / sum(n, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = periodo, y = l_cm, fill = prop)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    x = "Periodo",
    y = "Comprimento (cm)",
    fill = "Proporcao",
    title = "Serie mensal de frequencias de comprimento"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)
  )

salvar_plot(p_lfq, "18_LFQ_mensal_heatmap.png", 18, 8)

# 13. TABELA DE DIAGNOSTICO POR ANO -------------------------------
# Os limiares abaixo sao APENAS sinalizadores exploratorios, nao criterios
# universais de aceitacao para um modelo.
limiar_n <- 300
limiar_meses <- 8
limiar_eventos <- 8

# Numero de eventos e meses por ano ja esta em desc_ano
diagnostico_ano <- desc_ano %>%
  left_join(caudas_ano %>% select(ano, amplitude_90), by = "ano") %>%
  mutate(
    flag_n_baixo = n < limiar_n,
    flag_poucos_meses = n_meses < limiar_meses,
    flag_poucos_eventos = n_eventos < limiar_eventos,
    n_flags = flag_n_baixo + flag_poucos_meses + flag_poucos_eventos
  ) %>%
  arrange(desc(n_flags), ano)

print(diagnostico_ano)
write_csv(diagnostico_ano, file.path(pasta_saida, "25_diagnostico_modelos_length_based.csv"))

# 14. RELATORIO AUTOMATICO NO CONSOLE -----------------------------
cat("\n========================================================\n")
cat("DIAGNOSTICOS PRINCIPAIS PARA MODELOS LENGTH-BASED\n")
cat("========================================================\n")
cat("N individuos:", nrow(dados), "\n")
cat("N eventos:", n_distinct(dados$sample_id), "\n")
cat("Periodo:", min(dados$ano), "-", max(dados$ano), "\n")
cat("Comprimento:", min(dados$l_cm), "-", max(dados$l_cm), "cm\n")
cat("Resolucao observada:", resolucao_min, "cm\n")
cat("Mediana de peixes/evento:", median(resumo_eventos$n), "\n")
cat("Faixa de meses amostrados/ano:", min(desc_ano$n_meses), "-", max(desc_ano$n_meses), "\n")
cat("Faixa de zonas amostradas/ano:", min(desc_ano$n_zonas), "-", max(desc_ano$n_zonas), "\n")
cat("\nAnos com <", limiar_meses, "meses amostrados:\n")
print(desc_ano %>% filter(n_meses < limiar_meses) %>% select(ano, n, n_eventos, n_meses, n_zonas))

cat("\nATENCAO PARA INTERPRETACAO:\n")
cat("1) Nao use n de 75 peixes/evento como indice de abundancia.\n")
cat("2) Antes de agregar por ano, compare composicao bruta e padronizada.\n")
cat("3) Verifique se os peixes foram selecionados aleatoriamente da captura.\n")
cat("4) Confirme o tipo de comprimento (total, furcal etc.) e se mudou no tempo.\n")
cat("5) Confirme artes de pesca, seletividade e se o desenho espacial mudou.\n")
cat("6) Para LBSPR sera necessario, alem da LF, L_inf, M/K e L50/L95 de maturidade.\n")
cat("7) Para LBB/LBI/aLBI, confirme os parametros biologicos e os pressupostos de cada metodo.\n")
cat("8) Se a amostragem for dependente da pesca, a ausencia de peixes pequenos pode ser seletividade, nao ausencia populacional.\n")
cat("========================================================\n")

# 15. SKIM FINAL --------------------------------------------------
skimr::skim(dados %>% select(data, ano, mes, zona, l_cm, w_g, sexo, em, w_gon, w_fig))


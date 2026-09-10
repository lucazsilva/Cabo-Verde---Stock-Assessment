#-------------------------------------------------------------------------------------------#
#     Este script contem a avaliacao LBSPR de Decapterus macarellus em Cabo Verde          #
#   Analises possuem foco na serie de comprimentos amostrada entre 2004 e 2024              #
#   Modelo principal: Length-Based Spawning Potential Ratio - LBSPR                         #
#   Estrutura do modelo: Growth-Type-Group (GTG), configuracao padrao do pacote LBSPR       #
#   Parametros biologicos principais: Linf, M/K, L50 e L95                                 #
#   Maturidade principal: estadio I = imaturo; estadios II-VII = adultos/maduros            #
#   Modelo-base: Linf = 40.6 cm FL; M/K = 0.92/0.45; classes de comprimento de 1 cm         #
#   Sensibilidades one-at-a-time para crescimento, M/K e largura de classe                  #
#   Sensibilidades expandidas podem ser ativadas sem cruzar todos os parametros entre si    #
#   Codificacao criada por Silva, MLS ; 10/09/2026, Instituto do Mar - IMar, Mindelo       #
#-------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------#
# 0. PACOTES E CONFIGURACOES
#-------------------------------------------------------------------------------------------#

pacotes <- c("LBSPR", "dplyr", "tidyr", "readr", "ggplot2", "stringr", "purrr")

instalar <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(instalar) > 0) install.packages(instalar)

library(LBSPR)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(purrr)

# Caminho do banco de dados
arquivo_length <- "Length_oficial.csv"

# FALSE = roda apenas modelo-base + sensibilidades principais
# TRUE  = acrescenta as sensibilidades expandidas
RUN_EXTENDED <- FALSE

# O LBSPR usa por padrao o modelo GTG. Deixamos explicitamente definido.
MODTYPE <- "GTG"

# Coeficiente de variacao de Linf entre individuos.
# 0.10 e o valor padrao usual do pacote; aqui fica explicito para evitar default oculto.
CV_LINF <- 0.10

# Referencias de SPR apenas para visualizacao. Ajuste se houver regra de manejo especifica.
SPR_TARGET <- 0.40
SPR_LIMIT  <- 0.20

# Pasta para resultados
out_dir <- "LBSPR_results"
dir.create(out_dir, showWarnings = FALSE)

#-------------------------------------------------------------------------------------------#
# 1. LEITURA E LIMPEZA DOS DADOS
#-------------------------------------------------------------------------------------------#

dados_raw <- read_csv(
  arquivo_length,
  show_col_types = FALSE,
  na = c("", "NA", "N/A", "-", ".")
)

dados <- dados_raw %>%
  transmute(
    ano  = as.integer(ANO),
    mes  = as.integer(MES),
    zona = str_squish(as.character(ZONA)),
    l_cm = as.numeric(`L(cm)`),
    sexo = str_to_upper(str_squish(as.character(SEXO))),
    em   = as.numeric(EM)
  ) %>%
  filter(
    !is.na(ano),
    !is.na(l_cm),
    l_cm > 0
  )

# Mantem apenas sexos padronizados para analises sexuais/maturidade
# Os individuos de sexo I ou sem sexo continuam no banco para a composicao de comprimento.

cat("\nNumero total de comprimentos utilizados:", nrow(dados), "\n")
cat("Anos:", min(dados$ano), "a", max(dados$ano), "\n")
cat("Amplitude de comprimento:", min(dados$l_cm), "a", max(dados$l_cm), "cm\n\n")

#-------------------------------------------------------------------------------------------#
# 2. FUNCAO PARA ESTIMAR L50 E L95
#-------------------------------------------------------------------------------------------#

# Classificacao utilizada em TODAS as analises:
# I      = imaturo/jovem
# II-VII = adulto/maduro
#
# Nao e realizada sensibilidade para uma definicao alternativa de maturidade.

estimar_maturidade <- function(dat, sexo_alvo = "F") {
  
  d <- dat %>%
    filter(
      sexo == sexo_alvo,
      !is.na(em),
      em >= 1,
      em <= 7,
      !is.na(l_cm)
    ) %>%
    mutate(
      adulto = if_else(em >= 2, 1L, 0L)
    )
  
  if (nrow(d) < 50) stop("Poucos individuos para estimar a ogiva de maturidade.")
  
  mod <- glm(adulto ~ l_cm, family = binomial(link = "logit"), data = d)
  
  a <- unname(coef(mod)[1])
  b <- unname(coef(mod)[2])
  
  if (!is.finite(b) || b <= 0) {
    stop("A ogiva estimada nao possui inclinacao positiva. Verifique os dados.")
  }
  
  L50 <- -a / b
  L95 <- (log(19) - a) / b
  
  tibble(
    sexo = sexo_alvo,
    definicao = "I = jovem; II-VII = adultos",
    n = nrow(d),
    L50 = L50,
    L95 = L95
  )
}

# Ogiva utilizada no modelo-base e em TODAS as sensibilidades
mat_base <- estimar_maturidade(dados, sexo_alvo = "F")

print(mat_base)

write_csv(
  mat_base,
  file.path(out_dir, "maturity_L50_L95.csv")
)

#-------------------------------------------------------------------------------------------#
# 3. FUNCAO PARA CONSTRUIR FREQUENCIAS DE COMPRIMENTO
#-------------------------------------------------------------------------------------------#

# A origem das classes e definida a partir da precisao original dos dados.
# Para 1 cm, comprimentos inteiros 13, 14, 15... sao os pontos medios.
# Para 2 cm, as classes agrupam pares consecutivos e os pontos medios tornam-se 13.5, 15.5...

construir_lfq <- function(dat, bin_width = 1) {
  
  min_len <- floor(min(dat$l_cm, na.rm = TRUE))
  max_len <- ceiling(max(dat$l_cm, na.rm = TRUE))
  
  origin <- min_len - 0.5
  n_bins <- ceiling(((max_len + 0.5) - origin) / bin_width)
  breaks <- origin + (0:n_bins) * bin_width
  mids <- breaks[-length(breaks)] + bin_width / 2
  
  anos <- sort(unique(dat$ano))
  
  aux <- dat %>%
    mutate(
      bin_id = cut(
        l_cm,
        breaks = breaks,
        include.lowest = TRUE,
        right = FALSE,
        labels = FALSE
      )
    ) %>%
    filter(!is.na(bin_id)) %>%
    count(ano, bin_id, name = "n")
  
  grade <- tidyr::expand_grid(
    ano = anos,
    bin_id = seq_along(mids)
  ) %>%
    left_join(aux, by = c("ano", "bin_id")) %>%
    mutate(
      n = replace_na(n, 0L),
      LMids = mids[bin_id]
    )
  
  wide <- grade %>%
    select(LMids, ano, n) %>%
    pivot_wider(
      names_from = ano,
      values_from = n,
      values_fill = 0
    ) %>%
    arrange(LMids)
  
  list(
    data = wide,
    years = anos,
    bin_width = bin_width
  )
}

# Salva as LFQs de 1, 2 e 3 cm para auditoria
for (bw in c(1, 2, 3)) {
  lfq_tmp <- construir_lfq(dados, bw)$data
  write_csv(lfq_tmp, file.path(out_dir, paste0("LFQ_", bw, "cm.csv")))
}

#-------------------------------------------------------------------------------------------#
# 4. CENARIOS DE MODELO E SENSIBILIDADE
#-------------------------------------------------------------------------------------------#

L50_base <- mat_base$L50
L95_base <- mat_base$L95

MK_VIEIRA <- 0.92 / 0.45
MK_ALMADA_PAULY <- 0.64 / 0.34
MK_ALMADA_TANAKA <- 0.43 / 0.34

# Sensibilidades PRINCIPAIS: uma fonte de incerteza alterada por vez.
cenarios_core <- tibble::tribble(
  ~scenario,         ~set,    ~Linf,  ~MK,               ~L50,      ~L95,      ~BinWidth, ~maturity_definition, ~source,                                      ~rationale,
  "BASE",            "core",  40.6,   MK_VIEIRA,         L50_base,  L95_base,  1,         "II-VII adults",      "Vieira (2019) + current maturity data",     "Reference model using local growth/mortality and observed maturity",
  "S_Linf_high",     "core",  46.0,   MK_VIEIRA,         L50_base,  L95_base,  1,         "II-VII adults",      "da Cruz Delgado et al. (2024), 2017-2021", "Tests sensitivity to a higher recent estimate of asymptotic length",
  "S_MK_1.59",       "core",  40.6,   1.59,              L50_base,  L95_base,  1,         "II-VII adults",      "da Luz & Vieira (2020)",                    "Tests sensitivity to a lower published M/K estimate",
  "S_bin_2cm",       "core",  40.6,   MK_VIEIRA,         L50_base,  L95_base,  2,         "II-VII adults",      "Current length data",                      "Tests sensitivity to aggregation of the original 1-cm length classes"
)

# Sensibilidades EXPANDIDAS: uteis caso os cenarios principais indiquem forte sensibilidade.
# Nao sao combinacoes fatoriais; cada uma altera apenas um elemento em relacao ao BASE.
cenarios_extended <- tibble::tribble(
  ~scenario,         ~set,        ~Linf,  ~MK,               ~L50,      ~L95,      ~BinWidth, ~maturity_definition, ~source,                                      ~rationale,
  "E_Linf_41.48",    "extended",  41.48,  MK_VIEIRA,         L50_base,  L95_base,  1,         "II-VII adults",      "da Cruz Delgado et al. (2024), 2003-2007", "Alternative local Linf estimate close to the reference value",
  "E_MK_1.88",       "extended",  40.6,   MK_ALMADA_PAULY,   L50_base,  L95_base,  1,         "II-VII adults",      "Almada (1997), M by Pauly",                "Historical alternative M/K based on Pauly natural mortality",
  "E_MK_1.26",       "extended",  40.6,   MK_ALMADA_TANAKA,  L50_base,  L95_base,  1,         "II-VII adults",      "Almada (1997), M by Tanaka",               "Lower historical M/K estimate; represents an extreme alternative",
  "E_bin_3cm",       "extended",  40.6,   MK_VIEIRA,         L50_base,  L95_base,  3,         "II-VII adults",      "Current length data",                      "More aggressive aggregation to test sensitivity to class width"
)

cenarios_all <- bind_rows(cenarios_core, cenarios_extended)

# Salva a tabela completa de configuracao, mesmo quando RUN_EXTENDED = FALSE
write_csv(cenarios_all, file.path(out_dir, "LBSPR_sensitivity_scenarios.csv"))

if (RUN_EXTENDED) {
  cenarios_run <- cenarios_all
} else {
  cenarios_run <- cenarios_core
}

cat("\nCenarios que serao executados:\n")
print(cenarios_run %>% select(scenario, set, Linf, MK, L50, L95, BinWidth, maturity_definition))

#-------------------------------------------------------------------------------------------#
# 5. FUNCAO PARA RODAR UM CENARIO LBSPR
#-------------------------------------------------------------------------------------------#

rodar_lbspr <- function(sc, dat) {
  
  sc <- as.list(sc)
  
  lfq <- construir_lfq(dat, bin_width = sc$BinWidth)
  
  # O construtor LB_lengths le frequencias em CSV:
  # primeira coluna = pontos medios das classes; demais colunas = anos.
  arq_temp <- tempfile(fileext = ".csv")
  write_csv(lfq$data, arq_temp)
  
  pars <- new("LB_pars", verbose = FALSE)
  pars@Species  <- "Decapterus macarellus"
  pars@L_units  <- "cm"
  pars@Linf     <- sc$Linf
  pars@MK       <- sc$MK
  pars@L50      <- sc$L50
  pars@L95      <- sc$L95
  pars@CVLinf   <- CV_LINF
  pars@BinWidth <- sc$BinWidth
  
  lens <- new(
    "LB_lengths",
    LB_pars = pars,
    file = arq_temp,
    dataType = "freq",
    header = TRUE,
    verbose = FALSE
  )
  
  fit <- LBSPRfit(
    pars,
    lens,
    Control = list(modtype = MODTYPE),
    verbose = FALSE
  )
  
  # Estimativas suavizadas pelo filtro/smoother multianual do LBSPR
  sm <- as.data.frame(fit@Ests)
  if (!all(c("SL50", "SL95", "FM", "SPR") %in% names(sm))) {
    names(sm)[seq_len(min(4, ncol(sm)))] <- c("SL50", "SL95", "FM", "SPR")[seq_len(min(4, ncol(sm)))]
  }
  
  n_est <- nrow(sm)
  anos_use <- lfq$years[seq_len(min(length(lfq$years), n_est))]
  sm <- sm[seq_len(length(anos_use)), , drop = FALSE]
  
  est_smoothed <- sm %>%
    transmute(
      year = anos_use,
      SL50 = SL50,
      SL95 = SL95,
      FM = FM,
      SPR = SPR,
      estimate_type = "smoothed"
    )
  
  # Estimativas anuais brutas, antes da suavizacao temporal
  n_raw <- min(length(fit@SPR), length(lfq$years))
  est_raw <- tibble(
    year = lfq$years[seq_len(n_raw)],
    SL50 = fit@SL50[seq_len(n_raw)],
    SL95 = fit@SL95[seq_len(n_raw)],
    FM = fit@FM[seq_len(n_raw)],
    SPR = fit@SPR[seq_len(n_raw)],
    estimate_type = "raw"
  )
  
  estimates <- bind_rows(est_smoothed, est_raw) %>%
    mutate(
      scenario = sc$scenario,
      set = sc$set,
      Linf = sc$Linf,
      MK = sc$MK,
      L50_input = sc$L50,
      L95_input = sc$L95,
      BinWidth = sc$BinWidth
    )
  
  list(
    fit = fit,
    pars = pars,
    lengths = lens,
    lfq = lfq$data,
    estimates = estimates
  )
}

#-------------------------------------------------------------------------------------------#
# 6. RODAR TODOS OS CENARIOS
#-------------------------------------------------------------------------------------------#

fits <- vector("list", nrow(cenarios_run))
names(fits) <- cenarios_run$scenario

erros <- list()

for (i in seq_len(nrow(cenarios_run))) {
  
  sc <- cenarios_run[i, ]
  cat("\nRodando:", sc$scenario, "...\n")
  
  res <- tryCatch(
    rodar_lbspr(sc, dados),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    warning(paste("Falha no cenario", sc$scenario, ":", res$message))
    erros[[sc$scenario]] <- res$message
  } else {
    fits[[sc$scenario]] <- res
  }
}

fits <- fits[!vapply(fits, is.null, logical(1))]

if (length(fits) == 0) {
  stop("Nenhum cenario LBSPR foi ajustado com sucesso.")
}

resultados <- bind_rows(lapply(fits, `[[`, "estimates"))

write_csv(resultados, file.path(out_dir, "LBSPR_all_estimates.csv"))

if (length(erros) > 0) {
  tabela_erros <- tibble(
    scenario = names(erros),
    error = unlist(erros)
  )
  write_csv(tabela_erros, file.path(out_dir, "LBSPR_errors.csv"))
  print(tabela_erros)
}

#-------------------------------------------------------------------------------------------#
# 7. RESUMOS DO MODELO
#-------------------------------------------------------------------------------------------#

# Usaremos as estimativas suavizadas para os graficos principais.
res_sm <- resultados %>%
  filter(estimate_type == "smoothed")

# Ultimo ano de cada cenario
terminal <- res_sm %>%
  group_by(scenario) %>%
  filter(year == max(year, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(
    cenarios_all %>% select(scenario, source, rationale, maturity_definition),
    by = "scenario"
  )

write_csv(terminal, file.path(out_dir, "LBSPR_terminal_estimates.csv"))

# Media dos tres anos finais de cada cenario
ultimos3 <- res_sm %>%
  group_by(scenario) %>%
  arrange(year) %>%
  slice_tail(n = 3) %>%
  summarise(
    year_start = min(year),
    year_end = max(year),
    SPR_mean3 = mean(SPR, na.rm = TRUE),
    FM_mean3 = mean(FM, na.rm = TRUE),
    SL50_mean3 = mean(SL50, na.rm = TRUE),
    SL95_mean3 = mean(SL95, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(ultimos3, file.path(out_dir, "LBSPR_last3yr_summary.csv"))

print(terminal %>% select(scenario, year, SPR, FM, SL50, SL95))

#-------------------------------------------------------------------------------------------#
# 8. GRAFICOS PRINCIPAIS - SENSIBILIDADES CORE
#-------------------------------------------------------------------------------------------#

res_core <- res_sm %>% filter(set == "core")

# SPR ao longo do tempo
p_spr <- ggplot(res_core, aes(x = year, y = SPR, color = scenario)) +
  geom_hline(yintercept = SPR_TARGET, linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = SPR_LIMIT, linetype = "dotted", linewidth = 0.7) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = sort(unique(res_core$year))) +
  labs(
    x = "Year",
    y = "SPR",
    color = "Scenario",
    title = "LBSPR - annual spawning potential ratio"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_spr)
ggsave(file.path(out_dir, "LBSPR_SPR_core.png"), p_spr, width = 11, height = 6, dpi = 300)

# F/M ao longo do tempo
p_fm <- ggplot(res_core, aes(x = year, y = FM, color = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = sort(unique(res_core$year))) +
  labs(
    x = "Year",
    y = "F/M",
    color = "Scenario",
    title = "LBSPR - annual relative fishing mortality"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_fm)
ggsave(file.path(out_dir, "LBSPR_FM_core.png"), p_fm, width = 11, height = 6, dpi = 300)

# SL50 e SL95
sel_long <- res_core %>%
  select(year, scenario, SL50, SL95) %>%
  pivot_longer(
    cols = c(SL50, SL95),
    names_to = "parameter",
    values_to = "length"
  )

p_sel <- ggplot(sel_long, aes(x = year, y = length, color = scenario)) +
  geom_hline(yintercept = 20, color = "red", linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.3) +
  facet_wrap(~parameter, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(sel_long$year))) +
  labs(
    x = "Year",
    y = "Length (cm FL)",
    color = "Scenario",
    title = "LBSPR - estimated selectivity lengths"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p_sel)
ggsave(file.path(out_dir, "LBSPR_selectivity_core.png"), p_sel, width = 11, height = 8, dpi = 300)

# Comparacao do SPR terminal entre todos os cenarios efetivamente rodados
p_terminal <- terminal %>%
  mutate(scenario = reorder(scenario, SPR)) %>%
  ggplot(aes(x = scenario, y = SPR)) +
  geom_hline(yintercept = SPR_TARGET, linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = SPR_LIMIT, linetype = "dotted", linewidth = 0.7) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Terminal SPR",
    title = "LBSPR sensitivity analysis - terminal SPR"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

print(p_terminal)
ggsave(file.path(out_dir, "LBSPR_terminal_SPR_sensitivity.png"), p_terminal, width = 8, height = 6, dpi = 300)

#-------------------------------------------------------------------------------------------#
# 9. DIAGNOSTICOS DO MODELO-BASE
#-------------------------------------------------------------------------------------------#

# Os graficos nativos do pacote sao uteis para verificar ajuste das frequencias,
# seletividade e estimativas temporais.

if ("BASE" %in% names(fits)) {
  fit_base <- fits[["BASE"]]$fit
  
  cat("\n--- Diagnosticos do modelo-base ---\n")
  cat("Execute individualmente no painel Plots se desejar inspecionar:\n")
  cat("plotSize(fit_base)\n")
  cat("plotMat(fit_base)\n")
  cat("plotEsts(fit_base)\n\n")
  
  # Descomente para abrir diretamente ao final da execucao:
  # plotSize(fit_base)
  # plotMat(fit_base)
  # plotEsts(fit_base)
}

#-------------------------------------------------------------------------------------------#
# 10. COMPARACAO BASE VS SENSIBILIDADES EM TERMOS RELATIVOS
#-------------------------------------------------------------------------------------------#

# Percentual de diferenca no SPR terminal em relacao ao modelo-base.
if ("BASE" %in% terminal$scenario) {
  
  spr_base_terminal <- terminal %>%
    filter(scenario == "BASE") %>%
    pull(SPR) %>%
    first()
  
  comparacao_terminal <- terminal %>%
    mutate(
      delta_SPR = SPR - spr_base_terminal,
      delta_SPR_pct = 100 * (SPR - spr_base_terminal) / spr_base_terminal
    ) %>%
    select(
      scenario, year, SPR, delta_SPR, delta_SPR_pct,
      FM, SL50, SL95, Linf, MK, L50_input, L95_input, BinWidth,
      maturity_definition, source, rationale
    )
  
  write_csv(
    comparacao_terminal,
    file.path(out_dir, "LBSPR_terminal_sensitivity_comparison.csv")
  )
  
  print(comparacao_terminal)
}

#-------------------------------------------------------------------------------------------#
# 11. INFORMACOES IMPORTANTES PARA INTERPRETACAO
#-------------------------------------------------------------------------------------------#

cat("\n=====================================================================\n")
cat("LBSPR finalizado.\n")
cat("Resultados salvos em:", out_dir, "\n")
cat("\nInterpretacao recomendada:\n")
cat("1. Comece pelos cenarios CORE (RUN_EXTENDED = FALSE).\n")
cat("2. Nao combine automaticamente todos os valores de Linf e M/K.\n")
cat("3. Rode os cenarios expandidos apenas se a conclusao for sensivel nos CORE.\n")
cat("4. Inspecione plotSize(fit_base), plotMat(fit_base) e plotEsts(fit_base).\n")
cat("5. Lembre que o LBSPR assume uma composicao representativa da populacao explorada\n")
cat("   e uma aproximacao de equilibrio; interprete mudancas anuais junto ao desenho amostral.\n")
cat("=====================================================================\n")

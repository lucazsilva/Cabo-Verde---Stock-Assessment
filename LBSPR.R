#-------------------------------------------------------------------------------------------#
#     Este script contem a avaliacao LBSPR de Decapterus macarellus em Cabo Verde          #
#   Analises possuem foco na serie de comprimentos amostrada entre 2004 e 2024              #
#   Modelo principal: Length-Based Spawning Potential Ratio - LBSPR                         #
#   Parametros biologicos avaliados: Linf, M/K, L50 e L95                                  #
#   Maturidade dos dados atuais: somente estadio I = imaturo; II-VII = adultos/maduros      #
#   Modelo BASE: Vieira (2019) + maturidade das femeas dos dados atuais + bin de 1 cm        #
#   Sensibilidades one-at-a-time: crescimento, maturidade, largura de classe e M/K          #
#   O objetivo e representar alternativas biologicamente plausiveis sem ampliar             #
#   artificialmente a incerteza por combinacoes fatoriais extremas                          #
#   Codificacao criada por Silva, LVS ; 10/09/2026, Instituto do Mar - IMar, Mindelo       #
#-------------------------------------------------------------------------------------------#

# ---- 1. Pacotes -------------------------------------------------------------
required_pkgs <- c("LBSPR", "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "purrr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---- 2. Arquivos / pastas ----------------------------------------------------
length_file <- "Base_Comprimentos_Combinada_1988_2024.xlsx"
lh_file     <- "Parametros_Historia_de_vida.xlsx"
output_dir  <- "LBSPR_output"
dir.create(output_dir, showWarnings = FALSE)

ANO_MIN <- 2004
ANO_MAX <- 2024
sex_filter    <- "F"   # foco do LBSPR: fracao desovante (femeas)
BinWidth_base <- 1      # cm — largura de classe do cenario BASE
sat_incluido  <- FALSE  # estadio "SAT" nao faz parte da escala I-VII; excluido do ajuste da ogiva

# =============================================================================
# ---- 3. DADOS DE COMPRIMENTO (2004-2024) -----------------------------------
# =============================================================================
raw_len <- read_excel(length_file)

df_len <- raw_len %>%
  filter(!is.na(`L(cm)`), `L(cm)` > 0, Ano >= ANO_MIN, Ano <= ANO_MAX)

if (sex_filter != "ALL") df_len <- df_len %>% filter(Sexo == sex_filter)

cat("Medicoes de comprimento utilizadas (", ANO_MIN, "-", ANO_MAX, ", sexo=", sex_filter, "): ",
    nrow(df_len), "\n", sep = "")

# =============================================================================
# ---- 4. OGIVA DE MATURIDADE A PARTIR DOS DADOS ATUAIS (cenario BASE) -------
# Estadio I = imatura (0); estadio II-VII = madura (1); "SAT" excluido por
# nao pertencer a escala numerica I-VII usada nos dados.
# =============================================================================
mat_data <- df_len %>%
  filter(!is.na(Maturidade)) %>%
  mutate(estagio = as.character(Maturidade)) %>%
  filter(estagio %in% c(as.character(1:7), if (sat_incluido) "SAT")) %>%
  mutate(Madura = if_else(estagio == "1", 0, 1))  # SAT (se incluido) conta como madura

cat("Registos usados na ogiva de maturidade:", nrow(mat_data),
    "(imaturas:", sum(mat_data$Madura == 0), "| maduras:", sum(mat_data$Madura == 1), ")\n\n")

mat_model <- glm(Madura ~ `L(cm)`, data = mat_data, family = binomial)
a_mat <- unname(coef(mat_model)[1]); b_mat <- unname(coef(mat_model)[2])

base_L50 <- -a_mat / b_mat
base_L95 <- (log(0.95 / 0.05) - a_mat) / b_mat
ratio_L95_L50 <- base_L95 / base_L50   # usado para projetar bounds de L50 da literatura em L95 equivalente

cat("Ogiva de maturidade (dados atuais, femeas):\n")
cat("  L50 =", round(base_L50, 2), "cm | L95 =", round(base_L95, 2), "cm\n\n")

# =============================================================================
# ---- 5. PARAMETROS DE CRESCIMENTO E M/K DO CENARIO BASE (Vieira, 2019) -----
# =============================================================================
hv <- read_excel(lh_file, sheet = "Historia_de_vida")
pr <- read_excel(lh_file, sheet = "Parametros_relativos")
mf <- read_excel(lh_file, sheet = "Medidas_FL")
cf <- read_excel(lh_file, sheet = "Confiabilidade_Fontes") %>%
  select(fonte, pontuacao_media, nivel_confiabilidade) %>%
  filter(!is.na(pontuacao_media))

vieira_row <- hv %>% filter(str_detect(fonte, "^Vieira"), !is.na(linf_fl), !is.na(k), !is.na(m))
base_Linf <- vieira_row$linf_fl[1] / 10
base_K    <- vieira_row$k[1]
base_M    <- vieira_row$m[1]
base_MK   <- base_M / base_K

cat("Cenario BASE (crescimento/mortalidade — Vieira 2019):\n")
cat("  Linf =", round(base_Linf, 2), "cm | K =", base_K, "| M =", base_M,
    "| M/K =", round(base_MK, 3), "\n\n")

## -- Filtro de confiabilidade para as SENSIBILIDADES (mesmo criterio da
##    propria planilha: pontuacao_media >= 2,5) --------------------------------
strip_year <- function(x) str_trim(str_replace(x, "\\s*\\(.*", ""))
cf_keys <- cf %>% mutate(chave = strip_year(fonte)) %>%
  filter(pontuacao_media >= 2.5) %>% pull(chave) %>% unique()
is_reliable <- function(fonte_vec) strip_year(fonte_vec) %in% cf_keys

cat("Fontes confiaveis (pontuacao >= 2,5) usadas nas sensibilidades:\n")
cat(paste(" -", cf_keys), sep = "\n"); cat("\n\n")

## -- Intervalo de Linf (crescimento) -----------------------------------------
linf_cand <- hv %>%
  filter(!is.na(linf_fl), sex == "A", is_reliable(fonte)) %>%
  distinct(fonte, linf_fl) %>%
  transmute(fonte, Linf_cm = linf_fl / 10)
range_Linf <- range(linf_cand$Linf_cm)

## -- Intervalo de M/K ---------------------------------------------------------
mk_reportado <- pr %>% filter(parametro == "M/K", is_reliable(fonte)) %>%
  transmute(fonte, MK = valor)
mk_derivado <- hv %>% filter(!is.na(m), !is.na(k), sex == "A", is_reliable(fonte)) %>%
  transmute(fonte, MK = m / k)
mk_cand <- bind_rows(mk_reportado, mk_derivado) %>% distinct()
range_MK <- range(mk_cand$MK)

## -- Intervalo de L50 (para projetar o bound de maturidade) -------------------
l50_cand <- mf %>% filter(medida == "l50_fl", is_reliable(fonte)) %>%
  transmute(fonte, sex, L50_cm = valor_mm / 10)
range_L50_lit <- range(l50_cand$L50_cm)
# L95 equivalente projetado usando a razao L95/L50 da PROPRIA ogiva ajustada
# aos dados atuais (Secao 4), para manter a mesma forma de ogiva ao testar
# apenas um deslocamento de tamanho de maturacao.
range_L95_lit <- range_L50_lit * ratio_L95_L50

## -- Tabela-resumo para auditoria --------------------------------------------
resumo_parametros <- tibble(
  parametro = c("Linf_cm (base=Vieira)", "M_K (base=Vieira)",
                "L50_cm (base=dados atuais)", "L95_cm (base=dados atuais)"),
  base      = c(base_Linf, base_MK, base_L50, base_L95),
  minimo    = c(range_Linf[1], range_MK[1], range_L50_lit[1], range_L95_lit[1]),
  maximo    = c(range_Linf[2], range_MK[2], range_L50_lit[2], range_L95_lit[2]),
  observacao = c("intervalo: fontes confiaveis (Historia_de_vida)",
                 "intervalo: fontes confiaveis (Parametros_relativos + m/k derivado)",
                 "intervalo: L50 de fontes confiaveis (Medidas_FL)",
                 "estimado: L50 da literatura x razao L95/L50 da ogiva atual")
)
write.csv(resumo_parametros, file.path(output_dir, "parametros_base_sensibilidade.csv"), row.names = FALSE)
print(resumo_parametros); cat("\n")

# =============================================================================
# ---- 6. CLASSES DE COMPRIMENTO (largura base = 1 cm) ------------------------
# Estende ate cobrir 1.25 x o maior Linf testado (bom-senso do LBSPR: o maior
# bin de comprimento deve superar o comprimento assintotico).
# =============================================================================
maxL_mult <- 1.25

build_lengths <- function(df, bin_width) {
  min_len <- floor(min(df$`L(cm)`) / bin_width) * bin_width
  max_len <- ceiling(max(max(df$`L(cm)`), maxL_mult * range_Linf[2]) / bin_width) * bin_width
  breaks  <- seq(min_len, max_len + bin_width, by = bin_width)
  
  d <- df %>%
    mutate(LenBin = cut(`L(cm)`, breaks = breaks, include.lowest = TRUE,
                        right = FALSE, labels = breaks[-length(breaks)]))
  
  freq_table <- d %>%
    count(Ano, LenBin, .drop = FALSE) %>%
    tidyr::complete(Ano, LenBin, fill = list(n = 0)) %>%
    arrange(Ano, as.numeric(as.character(LenBin))) %>%
    tidyr::pivot_wider(names_from = Ano, values_from = n) %>%
    arrange(as.numeric(as.character(LenBin)))
  
  LMids    <- as.numeric(as.character(freq_table$LenBin))
  Years    <- as.numeric(colnames(freq_table)[-1])
  LFreqMat <- as.matrix(freq_table[, -1, drop = FALSE])
  storage.mode(LFreqMat) <- "numeric"
  
  Lengths <- new("LB_lengths")
  Lengths@LMids  <- LMids
  Lengths@LData  <- LFreqMat
  Lengths@Years  <- Years
  Lengths@NYears <- length(Years)
  Lengths
}

MyLengths_base <- build_lengths(df_len, BinWidth_base)

# =============================================================================
# ---- 7. FUNCAO AUXILIAR: roda o LBSPR para um cenario -----------------------
# =============================================================================
run_lbspr <- function(cenario, eixo_sensibilidade, valor_testado,
                      Linf, MK, L50, L95, BinWidth, Lengths) {
  Pars <- new("LB_pars")
  Pars@Species  <- "Decapterus macarellus"
  Pars@Linf     <- Linf
  Pars@L50      <- L50
  Pars@L95      <- L95
  Pars@MK       <- MK
  Pars@L_units  <- "cm"
  Pars@BinWidth <- BinWidth
  # NOTA: esta versao do LBSPR nao tem o slot @maxL em LB_pars. O comprimento
  # maximo efetivo e definido pelo maior valor em Lengths@LMids, que ja cobre
  # 1.25*Linf para todos os cenarios (Secao 6).
  
  out <- tryCatch({
    fit <- LBSPRfit(Pars, Lengths, verbose = FALSE)
    df <- data.frame(
      cenario = cenario, eixo_sensibilidade = eixo_sensibilidade,
      valor_testado = valor_testado, Ano = Lengths@Years,
      SPR = fit@Ests[, "SPR"], FM = fit@Ests[, "FM"],
      SL50 = fit@Ests[, "SL50"], SL95 = fit@Ests[, "SL95"],
      erro = NA_character_
    )
    list(df = df, fit = fit)
  }, error = function(e) {
    df <- data.frame(cenario = cenario, eixo_sensibilidade = eixo_sensibilidade,
                     valor_testado = valor_testado, Ano = NA, SPR = NA, FM = NA,
                     SL50 = NA, SL95 = NA, erro = conditionMessage(e))
    list(df = df, fit = NULL)
  })
  out
}

# Atalho para as sensibilidades, que so precisam da tabela de resultados
# (nao dos objetos LB_obj completos, usados apenas no cenario BASE abaixo).
run_lbspr_df <- function(...) run_lbspr(...)$df

# =============================================================================
# ---- 8. CENARIO BASE ---------------------------------------------------------
# =============================================================================
resultado_base <- run_lbspr("Base", "base", NA,
                            Linf = base_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
                            BinWidth = BinWidth_base, Lengths = MyLengths_base)
res_base <- resultado_base$df
fit_base <- resultado_base$fit   # objeto LB_obj — usado nos graficos nativos do LBSPR (Secao 11d)

if (is.null(fit_base)) {
  stop("O ajuste do cenario BASE falhou (", res_base$erro[1], "). Corrija antes de prosseguir.")
}

## -- 8a. GRAFICOS NATIVOS DO PACOTE LBSPR (Hordyk et al.) para o cenario BASE
## plotSize: estrutura de tamanho observada vs. ajustada
## plotMat:  ogivas de maturidade e de selectividade da frota
## plotEsts: series temporais dos parametros estimados (SL50, SL95, F/M, SPR)
##           com intervalos de confianca do proprio pacote
png(file.path(output_dir, "LBSPR_base_estrutura_tamanho.png"), width = 1200, height = 850, res = 120)
print(plotSize(fit_base))
dev.off()

png(file.path(output_dir, "LBSPR_base_maturidade_selectividade.png"), width = 1200, height = 850, res = 120)
print(plotMat(fit_base))
dev.off()

png(file.path(output_dir, "LBSPR_base_series_temporais.png"), width = 1200, height = 850, res = 120)
print(plotEsts(fit_base))
dev.off()

cat("Graficos nativos do LBSPR (cenario BASE) salvos em:", output_dir, "\n")
cat("  - LBSPR_base_estrutura_tamanho.png (plotSize)\n")
cat("  - LBSPR_base_maturidade_selectividade.png (plotMat)\n")
cat("  - LBSPR_base_series_temporais.png (plotEsts)\n\n")

# =============================================================================
# ---- 9. SENSIBILIDADES ONE-AT-A-TIME (OAT) -----------------------------------
# Cada eixo varia sozinho, mantendo os demais no cenario BASE — evita a
# explosao combinatoria de cenarios factoriais extremos.
# =============================================================================

## 9a. Crescimento (Linf) -------------------------------------------------------
res_crescimento <- bind_rows(
  run_lbspr_df("Crescimento_min", "crescimento", range_Linf[1],
               Linf = range_Linf[1], MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_max", "crescimento", range_Linf[2],
               Linf = range_Linf[2], MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

## 9b. Maturidade (L50 e L95 variam juntos, mesma forma de ogiva) --------------
res_maturidade <- bind_rows(
  run_lbspr_df("Maturidade_min", "maturidade", range_L50_lit[1],
               Linf = base_Linf, MK = base_MK, L50 = range_L50_lit[1], L95 = range_L95_lit[1],
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Maturidade_max", "maturidade", range_L50_lit[2],
               Linf = base_Linf, MK = base_MK, L50 = range_L50_lit[2], L95 = range_L95_lit[2],
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

## 9c. Largura de classe (bin width) -------------------------------------------
BinWidths_alt <- c(2, 3)  # cm — alternativas a 1 cm do cenario BASE
res_binwidth <- purrr::map_dfr(BinWidths_alt, function(bw) {
  Lengths_bw <- build_lengths(df_len, bw)
  run_lbspr_df(paste0("LarguraClasse_", bw, "cm"), "largura_classe", bw,
               Linf = base_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = bw, Lengths = Lengths_bw)
})

## 9d. M/K ----------------------------------------------------------------------
res_mk <- bind_rows(
  run_lbspr_df("MK_min", "M_K", range_MK[1],
               Linf = base_Linf, MK = range_MK[1], L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_max", "M_K", range_MK[2],
               Linf = base_Linf, MK = range_MK[2], L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

# =============================================================================
# ---- 10. CONSOLIDAR E EXPORTAR -----------------------------------------------
# =============================================================================
todos_resultados <- bind_rows(res_base, res_crescimento, res_maturidade, res_binwidth, res_mk)

if (any(!is.na(todos_resultados$erro))) {
  cat("\n\u26a0 Cenarios com erro no ajuste do LBSPR:\n")
  print(todos_resultados %>% filter(!is.na(erro)) %>%
          distinct(cenario, eixo_sensibilidade, valor_testado, erro))
}

write.csv(todos_resultados, file.path(output_dir, "LBSPR_sensibilidade_resultados.csv"), row.names = FALSE)
cat("\nResultados (base + sensibilidades) salvos em:",
    file.path(output_dir, "LBSPR_sensibilidade_resultados.csv"), "\n\n")

# =============================================================================
# ---- 11. GRAFICOS -------------------------------------------------------------
# =============================================================================

## 11a. SPR ao longo do tempo, por eixo de sensibilidade ----------------------
eixos <- unique(todos_resultados$eixo_sensibilidade[todos_resultados$eixo_sensibilidade != "base"])

base_repetido <- todos_resultados %>% filter(cenario == "Base") %>%
  tidyr::crossing(eixo_facet = eixos) %>%
  mutate(eixo_sensibilidade = eixo_facet) %>% select(-eixo_facet)

plot_dat <- bind_rows(
  todos_resultados %>% filter(is.na(erro), cenario != "Base"),
  base_repetido
)

p_time <- ggplot(plot_dat, aes(x = Ano, y = SPR, color = cenario)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0.40, linetype = "dashed", color = "darkgreen", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.20, linetype = "dashed", color = "firebrick", inherit.aes = FALSE) +
  facet_wrap(~eixo_sensibilidade, scales = "free_y") +
  labs(title = "Sensibilidade do SPR (LBSPR) — crescimento, maturidade, largura de classe e M/K",
       subtitle = "Linhas tracejadas: 0.20 (critico) e 0.40 (alvo) — convencoes comuns",
       y = "Spawning Potential Ratio (SPR)", x = "Ano", color = "Cenario") +
  theme_minimal(base_size = 12)
p_time
ggsave(file.path(output_dir, "sensibilidade_SPR_por_ano.png"), p_time, width = 11, height = 7, dpi = 150)

## 11b. Grafico "tornado" — impacto no SPR medio -------------------------------
spr_medio_base <- mean(res_base$SPR, na.rm = TRUE)

tornado_dat <- todos_resultados %>%
  filter(is.na(erro), eixo_sensibilidade != "base") %>%
  group_by(eixo_sensibilidade, cenario) %>%
  summarise(SPR_medio = mean(SPR, na.rm = TRUE), .groups = "drop") %>%
  mutate(delta = SPR_medio - spr_medio_base)

p_tornado <- ggplot(tornado_dat, aes(x = eixo_sensibilidade, y = delta, fill = cenario)) +
  geom_col(position = "identity", alpha = 0.8) +
  geom_hline(yintercept = 0, color = "black") +
  coord_flip() +
  labs(title = "Grafico tornado — sensibilidade do SPR medio",
       subtitle = paste0("Referencia (cenario BASE): SPR medio = ", round(spr_medio_base, 3)),
       x = NULL, y = "Variacao do SPR medio em relacao ao cenario BASE", fill = "Cenario") +
  theme_minimal(base_size = 12)
p_tornado
ggsave(file.path(output_dir, "sensibilidade_tornado_SPR_medio.png"), p_tornado, width = 9, height = 5.5, dpi = 150)

cat("Graficos salvos em:", output_dir, "\n")
cat("  - sensibilidade_SPR_por_ano.png\n")
cat("  - sensibilidade_tornado_SPR_medio.png\n\n")

# =============================================================================
# NOTAS FINAIS 
# -----------------------------------------------------------------------------
# 1. O cenario BASE combina crescimento/M/K de Vieira (2019) com uma ogiva de
#    maturidade ajustada diretamente aos dados atuais (femeas, 2004-2024,
#    estadios I a VII). "SAT" foi excluido do ajuste por nao pertencer a essa
#    escala — revise 'sat_incluido' na Secao 2 se quiser tratar esse estadio
#    de outra forma.
# 2. As sensibilidades sao one-at-a-time (um parametro varia, os demais ficam
#    no cenario BASE), conforme pedido no cabecalho, para nao inflar a
#    incerteza com combinacoes factoriais extremas.
# 3. O intervalo de L95 na sensibilidade de maturidade e ESTIMADO: aplica a
#    razao L95/L50 da propria ogiva ajustada aos dados atuais sobre os
#    extremos de L50 reportados na literatura.
# 4. Repita a extracao das Secoes 4-5 se as planilhas de entrada mudarem — o
#    script e totalmente reprodutivel a partir delas.
# =============================================================================
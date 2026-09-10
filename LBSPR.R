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

#-------------------------------------------------------------------------------------------#
# 0. PACOTES E CONFIGURACOES
#-------------------------------------------------------------------------------------------#

pacotes <- c(
  "LBSPR", "dplyr", "tidyr", "readr", "ggplot2", "stringr", "purrr", "tibble"
)

instalar <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(instalar) > 0) install.packages(instalar)

library(LBSPR)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)
library(purrr)
library(tibble)

arquivo_length <- "Length_oficial.csv"
out_dir <- "LBSPR_reduced_sensitivity_MK"
cache_dir <- file.path(out_dir, "cache_models")

dir.create(out_dir, showWarnings = FALSE)
dir.create(cache_dir, showWarnings = FALSE)

# Configuracao do modelo LBSPR
MODTYPE <- "GTG"
CV_LINF <- 0.10

# Linhas de referencia apenas para visualizacao
SPR_TARGET <- 0.40
SPR_LIMIT  <- 0.20

# Tamanho minimo legal de captura (cm FL), usado apenas em figuras quando necessario
MLS <- 20

# Temperatura usada para derivar M de Delgado 2003-2007 pela equacao de Pauly
T_PAULY <- 24.7

# TRUE = ignora modelos salvos em cache e recalcula tudo
RECALCULATE_ALL <- FALSE

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

cat("\nNumero total de comprimentos:", nrow(dados), "\n")
cat("Anos:", min(dados$ano), "a", max(dados$ano), "\n")
cat("Amplitude:", min(dados$l_cm), "a", max(dados$l_cm), "cm FL\n")

#-------------------------------------------------------------------------------------------#
# 2. ESTIMACAO DE L50 E L95 COM OS DADOS ATUAIS
#-------------------------------------------------------------------------------------------#

# Classificacao maturacional fixa:
# I       = imaturo/jovem
# II-VII  = adulto/maduro
#
# Para o modelo principal, a ogiva e estimada utilizando somente femeas.

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
  
  if (nrow(d) < 50) {
    stop("Poucos individuos para estimar a ogiva de maturidade.")
  }
  
  mod <- glm(
    adulto ~ l_cm,
    family = binomial(link = "logit"),
    data = d
  )
  
  a <- unname(coef(mod)[1])
  b <- unname(coef(mod)[2])
  
  if (!is.finite(b) || b <= 0) {
    stop("A ogiva estimada nao possui inclinacao positiva. Verifique os dados.")
  }
  
  L50 <- -a / b
  L95 <- (log(19) - a) / b
  
  tibble(
    maturity_id = paste0("OWN_", sexo_alvo),
    sex = sexo_alvo,
    maturity_source = "Current data 2004-2024",
    L50 = L50,
    L95 = L95,
    slope_b = b,
    L50_origin = "estimated from current individual maturity data",
    L95_origin = "estimated from current individual maturity data",
    maturity_note = "Stage I immature; stages II-VII adult/mature"
  )
}

mat_own_f <- estimar_maturidade(dados, sexo_alvo = "F")

cat("\n--- Maturidade estimada com os dados atuais ---\n")
print(mat_own_f)

write_csv(
  mat_own_f,
  file.path(out_dir, "maturity_current_data.csv")
)

#-------------------------------------------------------------------------------------------#
# 3. PARAMETROS BIOLOGICOS UTILIZADOS
#-------------------------------------------------------------------------------------------#

# Equacao de Pauly usada apenas para derivar M no cenario de crescimento de
# da Cruz Delgado et al. (2003-2007), pois esse conjunto fornece Linf e K,
# mas M/K nao esta disponivel diretamente na parametrizacao utilizada aqui.

pauly_M <- function(Linf_cm, K, T = T_PAULY) {
  10^(
    -0.0066 -
      0.2790 * log10(Linf_cm) +
      0.6543 * log10(K) +
      0.4630 * log10(T)
  )
}

# Modelo BASE - Vieira (2019)
Linf_base <- 40.60
K_base    <- 0.45
M_base    <- 0.92
MK_base   <- M_base / K_base

# Sensibilidade de crescimento - da Cruz Delgado et al. (2003-2007)
Linf_growth <- 41.48
K_growth    <- 0.39
M_growth    <- pauly_M(Linf_growth, K_growth, T_PAULY)
MK_growth   <- M_growth / K_growth

# Sensibilidades isoladas de M/K
# Mantem Linf, maturidade e bin width iguais ao BASE e altera somente M/K.
# O limite inferior vem de da Luz & Vieira (2020).
# O limite superior usa o M/K derivado para Delgado 2003-2007.
MK_low  <- 1.590
MK_high <- MK_growth

# Maturidade dos dados atuais
L50_own <- mat_own_f$L50
L95_own <- mat_own_f$L95

# Sensibilidade de maturidade - Costa et al. (2020), femeas
L50_costa_f <- 24.1
L95_costa_f <- 27.8

#-------------------------------------------------------------------------------------------#
# 4. CENARIOS REDUZIDOS DE SENSIBILIDADE
#-------------------------------------------------------------------------------------------#

# IMPORTANTE:
# Cada sensibilidade altera apenas UM componente em relacao ao BASE.
# Isso evita o cruzamento de extremos de crescimento, maturidade, bin width e M/K.
#
# BASE       = Vieira 2019 + maturidade atual F + bin 1 cm
# S_GROWTH   = muda o conjunto de crescimento para Delgado 2003-2007
# S_MATURITY = muda somente L50/L95 para Costa 2020 F
# S_BIN2     = muda somente a largura de classe para 2 cm
# S_MK_LOW   = muda somente M/K para 1.59
# S_MK_HIGH  = muda somente M/K para o valor derivado de Delgado 2003-2007

cenarios_run <- tribble(
  ~scenario,     ~scenario_label,                       ~sensitivity_component,
  ~growth_id,    ~growth_source,                        ~Linf,        ~K,          ~M,          ~MK,          ~MK_origin,
  ~maturity_id,  ~maturity_source,                      ~sex,         ~L50,        ~L95,        ~L50_origin,  ~L95_origin,
  ~BinWidth,     ~is_base,
  
  "BASE",        "BASE",                               "Reference model",
  "VIEIRA_2019", "Vieira (2019)",                      Linf_base,    K_base,       M_base,       MK_base,       "same-study M estimated by Pauly",
  "OWN_F",       "Current data 2004-2024",             "F",          L50_own,      L95_own,      "estimated",  "estimated",
  1,              TRUE,
  
  "S_GROWTH",    "Growth - Delgado 2003-2007",         "Growth and M/K",
  "DELGADO_2003_2007", "da Cruz Delgado et al. (2024)", Linf_growth, K_growth,     M_growth,     MK_growth,     paste0("derived with Pauly equation; T = ", T_PAULY, " C"),
  "OWN_F",       "Current data 2004-2024",             "F",          L50_own,      L95_own,      "estimated",  "estimated",
  1,              FALSE,
  
  "S_MATURITY",  "Maturity - Costa 2020 F",            "Maturity",
  "VIEIRA_2019", "Vieira (2019)",                      Linf_base,    K_base,       M_base,       MK_base,       "same-study M estimated by Pauly",
  "COSTA_2020_F", "Costa et al. (2020)",               "F",          L50_costa_f,  L95_costa_f,  "published",  "published",
  1,              FALSE,
  
  "S_BIN2",      "Bin width - 2 cm",                   "Length-class width",
  "VIEIRA_2019", "Vieira (2019)",                      Linf_base,    K_base,       M_base,       MK_base,       "same-study M estimated by Pauly",
  "OWN_F",       "Current data 2004-2024",             "F",          L50_own,      L95_own,      "estimated",  "estimated",
  2,              FALSE,
  
  "S_MK_LOW",    "M/K - 1.59",                         "M/K",
  "VIEIRA_2019", "Vieira (2019) Linf; M/K from da Luz & Vieira (2020)", Linf_base, K_base, NA_real_, MK_low, "M/K reported by da Luz & Vieira (2020); only M/K changed from BASE",
  "OWN_F",       "Current data 2004-2024",             "F",          L50_own,      L95_own,      "estimated",  "estimated",
  1,              FALSE,
  
  "S_MK_HIGH",   "M/K - 2.13",                         "M/K",
  "VIEIRA_2019", "Vieira (2019) Linf; M/K derived from Delgado 2003-2007", Linf_base, K_base, NA_real_, MK_high, paste0("M/K from Delgado 2003-2007 derived with Pauly equation; T = ", T_PAULY, " C; only M/K changed from BASE"),
  "OWN_F",       "Current data 2004-2024",             "F",          L50_own,      L95_own,      "estimated",  "estimated",
  1,              FALSE
)

# Verificacao automatica dos cenarios
cenarios_run <- cenarios_run %>%
  mutate(
    valid_biology =
      is.finite(Linf) & is.finite(MK) & is.finite(L50) & is.finite(L95) &
      MK > 0 & L50 > 0 & L95 > L50 & L95 < Linf
  )

if (any(!cenarios_run$valid_biology)) {
  print(cenarios_run %>% filter(!valid_biology))
  stop("Existe pelo menos um cenario biologicamente invalido.")
}

cat("\n--- Cenarios LBSPR utilizados ---\n")
print(
  cenarios_run %>%
    select(
      scenario, scenario_label, growth_id, Linf, K, M, MK,
      maturity_id, L50, L95, BinWidth
    )
)

write_csv(
  cenarios_run,
  file.path(out_dir, "LBSPR_reduced_sensitivity_scenarios.csv")
)

# Tabela simplificada para manuscrito/relatorio
config_table <- cenarios_run %>%
  transmute(
    Scenario = scenario,
    Description = scenario_label,
    `Sensitivity component` = sensitivity_component,
    `Growth source` = growth_source,
    `Linf (cm FL)` = Linf,
    K = K,
    M = M,
    `M/K` = MK,
    `Maturity source` = maturity_source,
    L50 = L50,
    L95 = L95,
    `Bin width (cm)` = BinWidth
  )

write_csv(
  config_table,
  file.path(out_dir, "LBSPR_sensitivity_configuration_table.csv")
)

#-------------------------------------------------------------------------------------------#
# 5. CONSTRUIR FREQUENCIAS DE COMPRIMENTO
#-------------------------------------------------------------------------------------------#

construir_lfq <- function(dat, bin_width = 1) {
  
  min_len <- floor(min(dat$l_cm, na.rm = TRUE))
  max_len <- ceiling(max(dat$l_cm, na.rm = TRUE))
  
  # Dados foram registrados em centimetros inteiros.
  # Para bin = 1 cm, os valores inteiros sao os pontos medios das classes.
  origin <- min_len - 0.5
  
  n_bins <- ceiling(
    ((max_len + 0.5) - origin) / bin_width
  )
  
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
  
  grade <- expand_grid(
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

# Nesta analise reduzida precisamos apenas de bins de 1 e 2 cm.
lfq_cache <- list(
  `1` = construir_lfq(dados, 1),
  `2` = construir_lfq(dados, 2)
)

for (bw in c(1, 2)) {
  write_csv(
    lfq_cache[[as.character(bw)]]$data,
    file.path(out_dir, paste0("LFQ_", bw, "cm.csv"))
  )
}

#-------------------------------------------------------------------------------------------#
# 6. FUNCAO PARA RODAR UM CENARIO LBSPR
#-------------------------------------------------------------------------------------------#

rodar_lbspr <- function(sc) {
  
  sc <- as.list(sc)
  lfq <- lfq_cache[[as.character(sc$BinWidth)]]
  
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
  
  # Estimativas suavizadas do ajuste multianual
  sm <- as.data.frame(fit@Ests)
  
  if (!all(c("SL50", "SL95", "FM", "SPR") %in% names(sm))) {
    names(sm)[seq_len(min(4, ncol(sm)))] <-
      c("SL50", "SL95", "FM", "SPR")[seq_len(min(4, ncol(sm)))]
  }
  
  n_est <- min(nrow(sm), length(lfq$years))
  
  est_smoothed <- sm[seq_len(n_est), , drop = FALSE] %>%
    transmute(
      year = lfq$years[seq_len(n_est)],
      SL50 = SL50,
      SL95 = SL95,
      FM = FM,
      SPR = SPR,
      estimate_type = "smoothed"
    )
  
  # Estimativas anuais nao suavizadas
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
      scenario_label = sc$scenario_label,
      sensitivity_component = sc$sensitivity_component,
      is_base = sc$is_base,
      growth_id = sc$growth_id,
      growth_source = sc$growth_source,
      Linf = sc$Linf,
      K = sc$K,
      M = sc$M,
      MK = sc$MK,
      MK_origin = sc$MK_origin,
      maturity_id = sc$maturity_id,
      maturity_source = sc$maturity_source,
      maturity_sex = sc$sex,
      L50_input = sc$L50,
      L95_input = sc$L95,
      L50_origin = sc$L50_origin,
      L95_origin = sc$L95_origin,
      BinWidth = sc$BinWidth
    )
  
  list(
    fit = fit,
    pars = pars,
    lengths = lens,
    estimates = estimates
  )
}

#-------------------------------------------------------------------------------------------#
# 7. RODAR OS SEIS CENARIOS
#-------------------------------------------------------------------------------------------#

fits <- vector("list", nrow(cenarios_run))
names(fits) <- cenarios_run$scenario

erros <- list()

for (i in seq_len(nrow(cenarios_run))) {
  
  sc <- cenarios_run[i, ]
  nome <- sc$scenario
  cache_file <- file.path(cache_dir, paste0(nome, ".rds"))
  
  cat(
    "\n[", i, "/", nrow(cenarios_run), "] ",
    "Rodando ", nome,
    " | Bin=", sc$BinWidth,
    " | Linf=", round(sc$Linf, 2),
    " | M/K=", round(sc$MK, 3),
    " | L50=", round(sc$L50, 2),
    " | L95=", round(sc$L95, 2),
    "\n",
    sep = ""
  )
  
  if (file.exists(cache_file) && !RECALCULATE_ALL) {
    
    res <- readRDS(cache_file)
    cat("  -> carregado do cache\n")
    
  } else {
    
    res <- tryCatch(
      rodar_lbspr(sc),
      error = function(e) e
    )
    
    if (!inherits(res, "error")) {
      saveRDS(res, cache_file)
    }
  }
  
  if (inherits(res, "error")) {
    
    warning(paste("Falha no cenario", nome, ":", res$message))
    erros[[nome]] <- res$message
    
  } else {
    
    fits[[nome]] <- res
  }
}

fits <- fits[!vapply(fits, is.null, logical(1))]

if (length(fits) == 0) {
  stop("Nenhum cenario LBSPR foi ajustado com sucesso.")
}

resultados <- bind_rows(lapply(fits, `[[`, "estimates"))

write_csv(
  resultados,
  file.path(out_dir, "LBSPR_all_estimates.csv")
)

if (length(erros) > 0) {
  
  tabela_erros <- tibble(
    scenario = names(erros),
    error = unlist(erros)
  )
  
  write_csv(
    tabela_erros,
    file.path(out_dir, "LBSPR_fit_errors.csv")
  )
  
  print(tabela_erros)
}

#-------------------------------------------------------------------------------------------#
# 8. RESUMOS DOS RESULTADOS
#-------------------------------------------------------------------------------------------#

res_sm <- resultados %>%
  filter(estimate_type == "smoothed")

# Estimativa terminal de cada cenario
terminal <- res_sm %>%
  group_by(scenario) %>%
  filter(year == max(year, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()

# Comparacao com o BASE
if ("BASE" %in% terminal$scenario) {
  
  spr_base <- terminal %>%
    filter(scenario == "BASE") %>%
    pull(SPR) %>%
    first()
  
  fm_base <- terminal %>%
    filter(scenario == "BASE") %>%
    pull(FM) %>%
    first()
  
  terminal <- terminal %>%
    mutate(
      delta_SPR = SPR - spr_base,
      delta_SPR_pct = 100 * delta_SPR / spr_base,
      delta_FM = FM - fm_base,
      delta_FM_pct = 100 * delta_FM / fm_base
    )
}

write_csv(
  terminal,
  file.path(out_dir, "LBSPR_terminal_estimates.csv")
)

# Media dos tres ultimos anos
ultimos3 <- res_sm %>%
  group_by(scenario, scenario_label) %>%
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
  ) %>%
  left_join(
    cenarios_run %>%
      select(
        scenario, sensitivity_component,
        growth_id, growth_source, Linf, K, M, MK,
        maturity_id, maturity_source, L50, L95,
        BinWidth, is_base
      ),
    by = "scenario"
  )

write_csv(
  ultimos3,
  file.path(out_dir, "LBSPR_last3yr_summary.csv")
)

#-------------------------------------------------------------------------------------------#
# 9. ENVELOPE REDUZIDO
#-------------------------------------------------------------------------------------------#

# Este envelope representa somente a amplitude entre os seis cenarios selecionados.
# NAO corresponde a intervalo de confianca estatistico.

annual_envelope <- res_sm %>%
  group_by(year) %>%
  summarise(
    n_scenarios = n(),
    SPR_min = min(SPR, na.rm = TRUE),
    SPR_median = median(SPR, na.rm = TRUE),
    SPR_max = max(SPR, na.rm = TRUE),
    FM_min = min(FM, na.rm = TRUE),
    FM_median = median(FM, na.rm = TRUE),
    FM_max = max(FM, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  annual_envelope,
  file.path(out_dir, "LBSPR_reduced_sensitivity_envelope.csv")
)

#-------------------------------------------------------------------------------------------#
# 10. GRAFICOS
#-------------------------------------------------------------------------------------------#

base_est <- res_sm %>%
  filter(scenario == "BASE")

# 10.1 Envelope reduzido de SPR
p_envelope <- ggplot(annual_envelope, aes(x = year)) +
  geom_ribbon(
    aes(ymin = SPR_min, ymax = SPR_max),
    alpha = 0.20
  ) +
  geom_line(
    aes(y = SPR_median),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_line(
    data = base_est,
    aes(x = year, y = SPR),
    linewidth = 1.1
  ) +
  geom_hline(
    yintercept = SPR_TARGET,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = SPR_LIMIT,
    linetype = "dotted",
    linewidth = 0.6
  ) +
  scale_x_continuous(
    breaks = sort(unique(annual_envelope$year))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "SPR",
    title = "LBSPR reduced sensitivity envelope",
    subtitle = "Ribbon = range among selected plausible scenarios; solid line = BASE; dashed line = scenario median"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.minor = element_blank()
  )

print(p_envelope)

ggsave(
  file.path(out_dir, "01_SPR_reduced_sensitivity_envelope.png"),
  p_envelope,
  width = 12,
  height = 6.5,
  dpi = 300
)

# 10.2 Trajetorias de SPR dos seis cenarios
p_spr_scenarios <- ggplot(
  res_sm,
  aes(
    x = year,
    y = SPR,
    color = scenario_label,
    linetype = scenario_label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = SPR_TARGET,
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = SPR_LIMIT,
    linetype = "dotted",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = sort(unique(res_sm$year))
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "SPR",
    color = "Scenario",
    linetype = "Scenario",
    title = "LBSPR sensitivity scenarios"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_spr_scenarios)

ggsave(
  file.path(out_dir, "02_SPR_selected_scenarios.png"),
  p_spr_scenarios,
  width = 12,
  height = 7,
  dpi = 300
)

# 10.3 SPR terminal dos seis cenarios
terminal <- terminal %>%
  mutate(
    model_type = if_else(
      scenario == "BASE",
      "BASE",
      "Sensitivity"
    )
  )

p_terminal <- ggplot(
  terminal,
  aes(
    x = scenario_label,
    y = SPR,
    shape = model_type
  )
) +
  geom_point(
    size = 4,
    stroke = 1.2
  ) +
  
  geom_hline(
    yintercept = SPR_TARGET,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_hline(
    yintercept = SPR_LIMIT,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  
  scale_shape_manual(
    values = c(
      "BASE" = 18,
      "Sensitivity" = 16
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  
  labs(
    x = NULL,
    y = "Terminal SPR",
    shape = NULL
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

p_terminal

ggsave(
  file.path(out_dir, "03_terminal_SPR_selected_scenarios.png"),
  p_terminal,
  width = 9,
  height = 6,
  dpi = 300
)

# 10.4 Trajetorias de F/M
p_fm <- ggplot(
  res_sm,
  aes(
    x = year,
    y = FM,
    color = scenario_label,
    linetype = scenario_label
  )
) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(
    breaks = sort(unique(res_sm$year))
  ) +
  labs(
    x = "Year",
    y = "F/M",
    color = "Scenario",
    linetype = "Scenario",
    title = "Fishing mortality relative to natural mortality"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_fm)

ggsave(
  file.path(out_dir, "04_FM_selected_scenarios.png"),
  p_fm,
  width = 12,
  height = 7,
  dpi = 300
)

# 10.5 SL50 estimado
p_sl50 <- ggplot(
  res_sm,
  aes(
    x = year,
    y = SL50,
    color = scenario_label,
    linetype = scenario_label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = MLS,
    color = "red",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = sort(unique(res_sm$year))
  ) +
  labs(
    x = "Year",
    y = expression(SL[50]~"(cm FL)"),
    color = "Scenario",
    linetype = "Scenario",
    title = expression("Estimated fishery selectivity "*SL[50]),
    subtitle = "Red dashed line = minimum legal catch size (20 cm FL)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_sl50)

ggsave(
  file.path(out_dir, "05_SL50_selected_scenarios.png"),
  p_sl50,
  width = 12,
  height = 7,
  dpi = 300
)

# 10.6 SL95 estimado
p_sl95 <- ggplot(
  res_sm,
  aes(
    x = year,
    y = SL95,
    color = scenario_label,
    linetype = scenario_label
  )
) +
  geom_line(linewidth = 0.9) +
  geom_hline(
    yintercept = MLS,
    color = "red",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_x_continuous(
    breaks = sort(unique(res_sm$year))
  ) +
  labs(
    x = "Year",
    y = expression(SL[95]~"(cm FL)"),
    color = "Scenario",
    linetype = "Scenario",
    title = expression("Estimated fishery selectivity "*SL[95]),
    subtitle = "Red dashed line = minimum legal catch size (20 cm FL)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p_sl95)

ggsave(
  file.path(out_dir, "06_SL95_selected_scenarios.png"),
  p_sl95,
  width = 12,
  height = 7,
  dpi = 300
)

#-------------------------------------------------------------------------------------------#
# 11. DIAGNOSTICOS DO MODELO BASE
#-------------------------------------------------------------------------------------------#

if ("BASE" %in% names(fits)) {
  
  fit_base <- fits[["BASE"]]$fit
  
  cat("\n--- Modelo BASE ---\n")
  cat(
    "Growth: Vieira (2019): Linf=40.6 cm, K=0.45, M=0.92, M/K=",
    round(MK_base, 4),
    "\n",
    sep = ""
  )
  cat(
    "Maturity: current female data; Stage I immature, II-VII adult; L50=",
    round(L50_own, 2),
    "; L95=",
    round(L95_own, 2),
    " cm FL\n",
    sep = ""
  )
  cat("Bin width: 1 cm\n\n")
  
  cat("Execute para diagnosticos graficos do BASE:\n")
  cat("plotSize(fit_base)\n")
  cat("plotMat(fit_base)\n")
  cat("plotEsts(fit_base)\n\n")
  
  # Descomente caso queira abrir automaticamente:
  # plotSize(fit_base)
  # plotMat(fit_base)
  # plotEsts(fit_base)
}

#-------------------------------------------------------------------------------------------#
# 12. RESUMO FINAL NO CONSOLE
#-------------------------------------------------------------------------------------------#

cat("\n=====================================================================\n")
cat("LBSPR REDUCED SENSITIVITY FINALIZADO\n")
cat("Resultados em:", out_dir, "\n\n")

cat("Cenarios avaliados:\n")
cat("1. BASE       - Vieira 2019 + own female maturity + bin 1 cm\n")
cat("2. S_GROWTH   - Delgado 2003-2007 growth/MK; demais parametros = BASE\n")
cat("3. S_MATURITY - Costa 2020 female L50/L95; demais parametros = BASE\n")
cat("4. S_BIN2     - bin width 2 cm; demais parametros = BASE\n")
cat("5. S_MK_LOW   - M/K = 1.59; todos os demais parametros = BASE\n")
cat("6. S_MK_HIGH  - M/K = ", round(MK_high, 3), "; todos os demais parametros = BASE\n\n", sep = "")

cat("Arquivos principais:\n")
cat(" - LBSPR_reduced_sensitivity_scenarios.csv\n")
cat(" - LBSPR_sensitivity_configuration_table.csv\n")
cat(" - LBSPR_all_estimates.csv\n")
cat(" - LBSPR_terminal_estimates.csv\n")
cat(" - LBSPR_last3yr_summary.csv\n")
cat(" - LBSPR_reduced_sensitivity_envelope.csv\n")
cat(" - 01_SPR_reduced_sensitivity_envelope.png\n")
cat(" - 02_SPR_selected_scenarios.png\n")
cat(" - 03_terminal_SPR_selected_scenarios.png\n")
cat(" - 04_FM_selected_scenarios.png\n")
cat(" - 05_SL50_selected_scenarios.png\n")
cat(" - 06_SL95_selected_scenarios.png\n\n")

cat("IMPORTANTE:\n")
cat("1. As sensibilidades sao one-at-a-time: apenas um componente muda em cada cenario.\n")
cat("2. O envelope representa a amplitude entre cenarios deterministicos selecionados; nao e IC.\n")
cat("3. L50 e L95 do BASE sao estimados diretamente das femeas do banco atual.\n")
cat("4. O cenario de maturidade usa a ogiva feminina completa de Costa et al. (2020).\n")
cat("5. O cenario de crescimento usa Delgado 2003-2007; M e M/K sao derivados por Pauly.\n")
cat("6. S_MK_LOW e S_MK_HIGH alteram somente M/K, mantendo Linf, L50, L95 e bin do BASE.\n")
cat("=====================================================================\n")

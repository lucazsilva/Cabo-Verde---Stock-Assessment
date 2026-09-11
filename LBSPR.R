#-------------------------------------------------------------------------------------------#
#     Este script contem a avaliacao LBSPR de Decapterus macarellus em Cabo Verde          #
#   Analises possuem foco na serie de comprimentos amostrada entre 2004 e 2024              #
#   Modelo principal: Length-Based Spawning Potential Ratio - LBSPR                         #
#   Parametros biologicos avaliados: Linf, M/K, L50 e L95                                  #
#   Maturidade dos dados atuais: somente estadio I = imaturo; II-VII = adultos/maduros      #
#   Modelo BASE: da Cruz Delgado et al. (2024), periodo pos-politica (2017-2021)            #
#     + L95 ajustado aos dados atuais + bin de 1 cm                                         #
#   Sensibilidades one-at-a-time: crescimento, maturidade, largura de classe e M/K          #
#   O objetivo e representar alternativas biologicamente plausiveis sem ampliar             #
#   artificialmente a incerteza por combinacoes fatoriais extremas                          #
#   Codificacao criada por Silva, LVS ; 10/09/2026, Instituto do Mar - IMar, Mindelo       #
#-------------------------------------------------------------------------------------------#
#
# BASE DE DADOS DE COMPRIMENTOS
#   Arquivo:      Base_Comprimentos_Combinada_1988_2024.xlsx
#   Serie principal (matriz de comprimento para o LBSPR): Ano entre 2004 e 2024
#   Subconjuntos adicionais (so para as sensibilidades de crescimento/M-K/maturidade,
#     replicando os dois periodos do artigo-base): 2003-2007 (pre-politica) e
#     2017-2021 (pos-politica) — usam a base completa (inclui 2003, presente no
#     arquivo original Base_Final_Cavala_1988_2003)
#   Colunas:      Ano, L(cm), Sexo, Maturidade (EM)
#   Sexo:         somente femeas (Sexo == "F")
#   Classe base:  largura de 1 cm (ver sensibilidade "largura de classe")
#
# ARTIGO-BASE DO CENARIO PRINCIPAL
#   da Cruz Delgado, K., Osemwegie, I., Medina, A. D., Nascimento da Luz, A.,
#   Kubik, Z., & Kouamelan, E. P. (2024). Ex-post evaluation of fishery
#   management policies on wild fisheries production in northern Cabo Verde:
#   An example of mackerel scad (Decapterus macarellus, Carangidae).
#   Journal of Fish Biology, 105(4), 1212-1226. https://doi.org/10.1111/jfb.15861
#
#   O artigo reporta, via FiSAT II (VBGF de von Bertalanffy/Beverton & Holt) e
#   modelo logistico de Haddon (L50), dois periodos de amostragem no porto de
#   Cova d'Inglesa, Sao Vicente (n = 7512 especimes):
#     PRE-POLITICA  (2003-2007, n=4424): Linf = 41,48 cm FL | K = 0,39 /ano | L50 = 24,2 cm
#     POS-POLITICA  (2017-2021, n=3039): Linf = 46,00 cm FL | K = 0,47 /ano | L50 = 28,2 cm
#   O artigo NAO reporta M nem M/K diretamente. Por isso, M e estimado aqui pela
#   equacao empirica de Pauly (1980) a partir do proprio Linf/K de cada periodo,
#   usando as temperaturas de superficie do mar (SST) citadas na propria
#   Discussao do artigo: SST = 24,1 C (Almada, 1997b) para o periodo pre-politica
#   e SST = 24,7 C (Vieira, 2018) para o periodo pos-politica.
#   O artigo tambem NAO reporta L95 (so L50). Por isso, L95 e estimado aqui
#   ajustando uma ogiva de maturidade aos dados atuais (femeas, mesmo periodo do
#   artigo) e aplicando a razao L95/L50 dessa ogiva sobre o L50 REPORTADO PELO
#   ARTIGO — ou seja, a escala (L50) vem do artigo, e a forma da ogiva (L95/L50)
#   vem dos dados. Isso e necessario porque a ogiva ajustada aos dados atuais,
#   sozinha, resulta em L50 sistematicamente menor que o do artigo (ver log do
#   script), provavelmente por diferenca de criterio de estadiamento — por isso
#   NAO usamos o L50 dos dados atuais como valor absoluto, so a forma da curva.
#
#-------------------------------------------------------------------------------------------#

# =============================================================================
# AVALIACAO DE ESTOQUE BASEADA EM COMPRIMENTO — PACOTE LBSPR (Hordyk et al. 2015)
# ICES J. Mar. Sci. 72(1): 217-231.
#
# IMPORTANTE:
#   - Este script nao pode ser executado no ambiente onde foi escrito (sem R
#     instalado). Rode-o localmente e revise as mensagens de log.
#   - LBSPR usa apenas a razao M/K (slot @MK), nao K isoladamente — por isso
#     "crescimento" e testado via Linf, e M/K e testado separadamente.
# =============================================================================

# ---- 1. Pacotes -------------------------------------------------------------
required_pkgs <- c("LBSPR", "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "purrr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---- 2. Arquivos / pastas ----------------------------------------------------
length_file <- "Base_Comprimentos_Combinada_1988_2024.xlsx"
output_dir  <- "LBSPR_output"
dir.create(output_dir, showWarnings = FALSE)

ANO_MIN <- 1988   # janela principal (serie de comprimento para o LBSPR)
ANO_MAX <- 2024
sex_filter    <- "F"
BinWidth_base <- 1   # cm
sat_incluido  <- FALSE  # "SAT" nao pertence a escala I-VII; excluido do ajuste da ogiva

# =============================================================================
# ---- 3. DADOS DE COMPRIMENTO -------------------------------------------------
# =============================================================================
raw_len <- read_excel(length_file)   # base completa (1988-2024), usada tambem
# para os subconjuntos pre/pos-politica

df_len <- raw_len %>%
  filter(!is.na(`L(cm)`), `L(cm)` > 0, Ano >= ANO_MIN, Ano <= ANO_MAX)
if (sex_filter != "ALL") df_len <- df_len %>% filter(Sexo == sex_filter)

cat("Medicoes de comprimento (serie principal, ", ANO_MIN, "-", ANO_MAX, ", sexo=", sex_filter, "): ",
    nrow(df_len), "\n", sep = "")

## Diagnostico: tamanho amostral por ano. Anos com poucas amostras (ou
## composicao de tamanho atipica, por exemplo dominada por peixes pequenos
## antes da vigencia do MLS/BRP a partir de 2008) podem gerar estimativas de
## SPR pouco confiaveis ou proximas de 0 no LBSPR — ver NOTA 6 no final do
## script para interpretacao.
n_por_ano <- df_len %>% count(Ano, name = "n_medicoes")
write.csv(n_por_ano, file.path(output_dir, "diagnostico_n_por_ano.csv"), row.names = FALSE)
cat("Tamanho amostral por ano (serie principal):\n")
print(n_por_ano, n = Inf)
cat("\n")

MIN_N_ANO <- 150  # limiar minimo de medicoes por ano para entrar na analise

anos_excluidos <- n_por_ano %>% filter(n_medicoes < MIN_N_ANO)
if (nrow(anos_excluidos) > 0) {
  cat("Anos excluidos por amostra insuficiente (<", MIN_N_ANO, "medicoes):\n")
  print(anos_excluidos)
  df_len <- df_len %>% filter(!(Ano %in% anos_excluidos$Ano))
}
# =============================================================================
# ---- 4. PARAMETROS DE CRESCIMENTO (da Cruz Delgado et al., 2024) -----------
# =============================================================================
# Valores reportados diretamente no artigo (Secao 3.2 e Figura 2)
crescimento_pre <- list(periodo = "2003-2007 (pre-politica)", Linf = 41.48, K = 0.39, SST = 24.1)
crescimento_pos <- list(periodo = "2017-2021 (pos-politica)", Linf = 46.00, K = 0.47, SST = 24.7)

# ---- 4a. M via equacao empirica de Pauly (1980) ------------------------------
# log10(M) = -0.0066 - 0.279*log10(Linf) + 0.6543*log10(K) + 0.4634*log10(T)
pauly_M <- function(Linf, K, T) {
  10^(-0.0066 - 0.279 * log10(Linf) + 0.6543 * log10(K) + 0.4634 * log10(T))
}

M_pre <- pauly_M(crescimento_pre$Linf, crescimento_pre$K, crescimento_pre$SST)
M_pos <- pauly_M(crescimento_pos$Linf, crescimento_pos$K, crescimento_pos$SST)
MK_pre <- M_pre / crescimento_pre$K
MK_pos <- M_pos / crescimento_pos$K

cat("\nCrescimento e M/K (da Cruz Delgado et al. 2024 + M via Pauly 1980):\n")
cat("  Pre-politica  (2003-2007): Linf=", crescimento_pre$Linf, "cm | K=", crescimento_pre$K,
    "| M=", round(M_pre, 3), "| M/K=", round(MK_pre, 3), "\n")
cat("  Pos-politica  (2017-2021): Linf=", crescimento_pos$Linf, "cm | K=", crescimento_pos$K,
    "| M=", round(M_pos, 3), "| M/K=", round(MK_pos, 3), "\n\n")

# =============================================================================
# ---- 5. MATURIDADE — L50 do artigo + L95 estimado a partir dos dados atuais -
# =============================================================================
L50_artigo_pre <- 24.2   # cm, Haddon logistic model, da Cruz Delgado et al. (2024)
L50_artigo_pos <- 28.2   # cm

fit_ogiva <- function(df_full, ano_min, ano_max) {
  d <- df_full %>%
    filter(Sexo == "F", !is.na(Maturidade), Ano >= ano_min, Ano <= ano_max) %>%
    mutate(estagio = as.character(Maturidade)) %>%
    filter(estagio %in% c(as.character(1:7), if (sat_incluido) "SAT")) %>%
    mutate(Madura = if_else(estagio == "1", 0, 1))
  
  modelo <- glm(Madura ~ `L(cm)`, data = d, family = binomial)
  a <- unname(coef(modelo)[1]); b <- unname(coef(modelo)[2])
  L50 <- -a / b
  L95 <- (log(0.95 / 0.05) - a) / b
  list(L50 = L50, L95 = L95, n = nrow(d))
}

ogiva_pre <- fit_ogiva(raw_len, 2003, 2007)
ogiva_pos <- fit_ogiva(raw_len, 2017, 2021)

cat("Ogivas de maturidade ajustadas aos dados atuais (femeas, validacao vs. artigo):\n")
cat("  Pre-politica : L50_dados=", round(ogiva_pre$L50, 2), "| L95_dados=", round(ogiva_pre$L95, 2),
    "| n=", ogiva_pre$n, "| (artigo reporta L50=", L50_artigo_pre, ")\n")
cat("  Pos-politica : L50_dados=", round(ogiva_pos$L50, 2), "| L95_dados=", round(ogiva_pos$L95, 2),
    "| n=", ogiva_pos$n, "| (artigo reporta L50=", L50_artigo_pos, ")\n")
cat("  NOTA: o L50 ajustado aos dados atuais tende a ficar abaixo do L50 do artigo\n")
cat("  (possivel diferenca no criterio de estadiamento). Por isso, usamos o L50 do\n")
cat("  ARTIGO como valor absoluto e so a RAZAO L95/L50 dos dados para estimar L95.\n\n")

# Razao L95/L50 (forma da ogiva) aplicada ao L50 reportado pelo artigo
ratio_L95_L50_pre <- ogiva_pre$L95 / ogiva_pre$L50
ratio_L95_L50_pos <- ogiva_pos$L95 / ogiva_pos$L50
L95_pre <- L50_artigo_pre * ratio_L95_L50_pre
L95_pos <- L50_artigo_pos * ratio_L95_L50_pos

cat("Maturidade final usada no LBSPR (L50 do artigo + L95 estimado):\n")
cat("  Pre-politica : L50=", L50_artigo_pre, "| L95=", round(L95_pre, 2), "\n")
cat("  Pos-politica : L50=", L50_artigo_pos, "| L95=", round(L95_pos, 2), "\n\n")

# =============================================================================
# ---- 6. CENARIO BASE = POS-POLITICA (2017-2021) ------------------------------
# =============================================================================
base_Linf <- crescimento_pos$Linf
base_MK   <- MK_pos
base_L50  <- L50_artigo_pos
base_L95  <- L95_pos

# Valores do periodo PRE-POLITICA — usados como alternativa em cada
# sensibilidade OAT (comparando com o artigo, nao com uma faixa arbitraria)
alt_Linf <- crescimento_pre$Linf
alt_MK   <- MK_pre
alt_L50  <- L50_artigo_pre
alt_L95  <- L95_pre

resumo_parametros <- tibble(
  parametro = c("Linf_cm", "M_K", "L50_cm", "L95_cm"),
  base_pos_politica = c(base_Linf, round(base_MK, 3), base_L50, round(base_L95, 2)),
  alt_pre_politica  = c(alt_Linf, round(alt_MK, 3), alt_L50, round(alt_L95, 2)),
  fonte = c("da Cruz Delgado et al. (2024)",
            "M via Pauly (1980), a partir de Linf/K do artigo + SST citada no artigo",
            "da Cruz Delgado et al. (2024), modelo de Haddon",
            "razao L95/L50 de ogiva ajustada aos dados atuais x L50 do artigo")
)
write.csv(resumo_parametros, file.path(output_dir, "parametros_base_sensibilidade.csv"), row.names = FALSE)
print(resumo_parametros); cat("\n")

# =============================================================================
# ---- 7. CLASSES DE COMPRIMENTO (largura base = 1 cm) ------------------------
# Estende ate cobrir 1.25 x o maior Linf testado (bom-senso do LBSPR).
# =============================================================================
maxL_mult <- 1.25
maior_Linf_testado <- max(base_Linf, alt_Linf)

build_lengths <- function(df, bin_width) {
  min_len <- floor(min(df$`L(cm)`) / bin_width) * bin_width
  max_len <- ceiling(max(max(df$`L(cm)`), maxL_mult * maior_Linf_testado) / bin_width) * bin_width
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
# ---- 8. FUNCAO AUXILIAR: RMSE entre comprimento observado e previsto --------
# O LBSPR guarda a composicao de comprimento OBSERVADA em @LData (contagens)
# e a composicao PREVISTA pelo modelo ajustado em @pLCatch (proporcoes). Aqui
# normalizamos @LData para proporcoes por ano e comparamos com @pLCatch.
# NOTA: nomes de slots ja variaram entre versoes do LBSPR nesta conversa
# (ex.: @maxL). Se este calculo falhar, rode slotNames(fit_base) e confirme
# se o slot de comprimento previsto se chama mesmo "pLCatch" na sua versao.
# =============================================================================
compute_rmse <- function(fit) {
  obs_prop <- tryCatch(
    sweep(fit@LData, 2, colSums(fit@LData), "/"),
    error = function(e) NULL
  )
  pred_prop <- tryCatch(fit@pLCatch, error = function(e) NULL)
  
  if (is.null(obs_prop) || is.null(pred_prop) || !identical(dim(obs_prop), dim(pred_prop))) {
    return(data.frame(Ano = fit@Years, RMSE = NA_real_))
  }
  
  rmse_por_ano <- sqrt(colMeans((obs_prop - pred_prop)^2, na.rm = TRUE))
  data.frame(Ano = fit@Years, RMSE = as.numeric(rmse_por_ano))
}

# =============================================================================
# ---- 9. FUNCAO AUXILIAR: roda o LBSPR para um cenario -----------------------
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
  # maximo efetivo e definido pelo maior valor em Lengths@LMids (Secao 7).
  
  out <- tryCatch({
    fit <- LBSPRfit(Pars, Lengths, verbose = FALSE)
    rmse_df <- compute_rmse(fit)
    df <- data.frame(
      cenario = cenario, eixo_sensibilidade = eixo_sensibilidade,
      valor_testado = valor_testado, Ano = Lengths@Years,
      SPR = fit@Ests[, "SPR"], FM = fit@Ests[, "FM"],
      SL50 = fit@Ests[, "SL50"], SL95 = fit@Ests[, "SL95"],
      erro = NA_character_
    )
    df <- left_join(df, rmse_df, by = "Ano")
    list(df = df, fit = fit)
  }, error = function(e) {
    df <- data.frame(cenario = cenario, eixo_sensibilidade = eixo_sensibilidade,
                     valor_testado = valor_testado, Ano = NA, SPR = NA, FM = NA,
                     SL50 = NA, SL95 = NA, RMSE = NA, erro = conditionMessage(e))
    list(df = df, fit = NULL)
  })
  out
}
run_lbspr_df <- function(...) run_lbspr(...)$df

# =============================================================================
# ---- 10. CENARIO BASE (pos-politica, 2017-2021) -------------------------------
# =============================================================================
resultado_base <- run_lbspr("Base_PosPolitica", "base", NA,
                            Linf = base_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
                            BinWidth = BinWidth_base, Lengths = MyLengths_base)
res_base <- resultado_base$df
fit_base <- resultado_base$fit

if (is.null(fit_base)) {
  stop("O ajuste do cenario BASE falhou (", res_base$erro[1], "). Corrija antes de prosseguir.")
}

## -- 10a. GRAFICOS NATIVOS DO PACOTE LBSPR (Hordyk et al.) para o cenario BASE
png(file.path(output_dir, "LBSPR_base_estrutura_tamanho.png"), width = 1200, height = 850, res = 120)
print(plotSize(fit_base)); dev.off()

png(file.path(output_dir, "LBSPR_base_maturidade_selectividade.png"), width = 1200, height = 850, res = 120)
print(plotMat(fit_base)); dev.off()

png(file.path(output_dir, "LBSPR_base_series_temporais.png"), width = 1200, height = 850, res = 120)
print(plotEsts(fit_base)); dev.off()

cat("Graficos nativos do LBSPR (cenario BASE) salvos em:", output_dir, "\n")
cat("  - LBSPR_base_estrutura_tamanho.png (plotSize)\n")
cat("  - LBSPR_base_maturidade_selectividade.png (plotMat)\n")
cat("  - LBSPR_base_series_temporais.png (plotEsts)\n\n")

# =============================================================================
# ---- 11. SENSIBILIDADES ONE-AT-A-TIME (OAT) ----------------------------------
# Cada eixo compara o cenario BASE (pos-politica) com alternativas
# biologicamente plausiveis: (i) o periodo PRE-POLITICA do MESMO artigo,
# (ii) fontes de maior confiabilidade ("Alta") da planilha de historia de
# vida — Jardim (1996/1999) e Costa et al. (2020) — e (iii) as fontes de
# confiabilidade "Média" da mesma planilha (aba Confiabilidade_Fontes):
# Carvalho & Caramelo (1996/1999), Santos (2018), Almada (1997),
# da Luz & Vieira (2020) e Vieira (2019). Fontes "Baixa" (Stobberup & Erzini,
# Tariche & Martins) continuam fora, como ja discutido antes.
# =============================================================================

## -- Jardim (1996/1999) [Alta]: crescimento e M/K --------------------------------
## Linf=31,5cm / K=0,43 (Historia_de_vida). M/K reportado por dois metodos:
## Tanaka (M/K=1,00) e Pauly (M/K=1,86) — Parametros_relativos.
jardim_Linf <- 31.5
jardim_MK_tanaka <- 1.00
jardim_MK_pauly  <- 1.86

## -- Costa et al. (2020) [Alta]: maturidade (femeas) -----------------------------
## L50=24,1cm / L95=27,8cm — unica fonte com L95 reportado diretamente
## (nao precisa de razao estimada, ao contrario do L95 do cenario BASE).
costa_L50 <- 24.1
costa_L95 <- 27.8

## -- Fontes de confiabilidade "Média" ---------------------------------------------
## Crescimento (Linf, cm): valor unico reportado por fonte.
carvalho_caramelo_Linf <- 41.0    # Carvalho & Caramelo (1996/1999)
santos_Linf            <- 41.83  # Santos (2018)
almada_Linf            <- 30.1    # Almada (1997)
da_luz_vieira_Linf     <- 40.9    # da Luz & Vieira (2020)
vieira_Linf            <- 40.6    # Vieira (2019)

## M/K: alguns autores reportam mais de um metodo — mantidos separados por
## serem estimativas genuinamente distintas da literatura, nao combinacoes
## artificiais.
carvalho_caramelo_MK <- 1.00                 # M via Pauly / K=0,30
santos_MK_tanaka     <- 1.54                 # M via Tanaka / K=0,39
santos_MK_pauly      <- 1.71                 # M via Pauly / K=0,39
almada_MK_pauly      <- 0.64 / 0.34          # M via Pauly / K=0,34
almada_MK_tanaka     <- 0.43 / 0.34          # M via Tanaka / K=0,34
da_luz_vieira_MK      <- 1.59                # mediana entre anos, metodo LBB
vieira_MK             <- 0.92 / 0.45         # M/K derivado (m e k na mesma linha)

## Maturidade: so as fontes que reportam L50 podem entrar neste eixo.
## Carvalho & Caramelo e da Luz & Vieira NAO reportam L50 (esta ultima so
## reporta comprimento de SELETIVIDADE, nao de maturidade) — ficam de fora.
## L95 nao reportado por nenhuma delas -> estimado com a mesma razao
## L95/L50 da ogiva ajustada aos dados atuais (pos-politica), como no
## cenario BASE.
santos_L50 <- 22.9                                  # Santos (2018), sexos combinados
santos_L95 <- santos_L50 * ratio_L95_L50_pos
almada_L50 <- 21.44                                 # Almada (1997), femeas
almada_L95 <- almada_L50 * ratio_L95_L50_pos
vieira_L50 <- 20.3                                  # Vieira (2019), sexos combinados
vieira_L95 <- vieira_L50 * ratio_L95_L50_pos

## 11a. Crescimento (Linf) ----------------------------------------------------------
res_crescimento <- bind_rows(
  run_lbspr_df("Crescimento_PrePolitica", "crescimento", alt_Linf,
               Linf = alt_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_Jardim_1996_1999", "crescimento", jardim_Linf,
               Linf = jardim_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_CarvalhoCaramelo", "crescimento", carvalho_caramelo_Linf,
               Linf = carvalho_caramelo_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_Santos_2018", "crescimento", santos_Linf,
               Linf = santos_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_Almada_1997", "crescimento", almada_Linf,
               Linf = almada_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_daLuzVieira_2020", "crescimento", da_luz_vieira_Linf,
               Linf = da_luz_vieira_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_Vieira_2019", "crescimento", vieira_Linf,
               Linf = vieira_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

## 11b. Maturidade (L50 e L95 variam juntos, mesma forma de ogiva) -----------------
res_maturidade <- bind_rows(
  run_lbspr_df("Maturidade_PrePolitica", "maturidade", alt_L50,
               Linf = base_Linf, MK = base_MK, L50 = alt_L50, L95 = alt_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Maturidade_Costa_2020_Femea", "maturidade", costa_L50,
               Linf = base_Linf, MK = base_MK, L50 = costa_L50, L95 = costa_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Maturidade_Santos_2018", "maturidade", santos_L50,
               Linf = base_Linf, MK = base_MK, L50 = santos_L50, L95 = santos_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Maturidade_Almada_1997_Femea", "maturidade", almada_L50,
               Linf = base_Linf, MK = base_MK, L50 = almada_L50, L95 = almada_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Maturidade_Vieira_2019", "maturidade", vieira_L50,
               Linf = base_Linf, MK = base_MK, L50 = vieira_L50, L95 = vieira_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

## 11c. Largura de classe (bin width) -----------------------------------------------
BinWidths_alt <- c(2, 3)  # cm — alternativas a 1 cm do cenario BASE
res_binwidth <- purrr::map_dfr(BinWidths_alt, function(bw) {
  Lengths_bw <- build_lengths(df_len, bw)
  run_lbspr_df(paste0("LarguraClasse_", bw, "cm"), "largura_classe", bw,
               Linf = base_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = bw, Lengths = Lengths_bw)
})

## 11d. M/K ---------------------------------------------------------------------------
res_mk <- bind_rows(
  run_lbspr_df("MK_PrePolitica", "M_K", alt_MK,
               Linf = base_Linf, MK = alt_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_Jardim_Tanaka", "M_K", jardim_MK_tanaka,
               Linf = base_Linf, MK = jardim_MK_tanaka, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_Jardim_Pauly", "M_K", jardim_MK_pauly,
               Linf = base_Linf, MK = jardim_MK_pauly, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_CarvalhoCaramelo", "M_K", carvalho_caramelo_MK,
               Linf = base_Linf, MK = carvalho_caramelo_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_Santos_Tanaka", "M_K", santos_MK_tanaka,
               Linf = base_Linf, MK = santos_MK_tanaka, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_Santos_Pauly", "M_K", santos_MK_pauly,
               Linf = base_Linf, MK = santos_MK_pauly, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_Almada_Pauly", "M_K", almada_MK_pauly,
               Linf = base_Linf, MK = almada_MK_pauly, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_Almada_Tanaka", "M_K", almada_MK_tanaka,
               Linf = base_Linf, MK = almada_MK_tanaka, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_daLuzVieira_2020", "M_K", da_luz_vieira_MK,
               Linf = base_Linf, MK = da_luz_vieira_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_Vieira_2019", "M_K", vieira_MK,
               Linf = base_Linf, MK = vieira_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

# =============================================================================
# ---- 12. CONSOLIDAR E EXPORTAR -----------------------------------------------
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

## -- 12a. RESUMO DO RMSE POR CENARIO -------------------------------------------
## RMSE medio (entre todos os anos) da composicao de comprimento prevista vs.
## observada — quanto MENOR o RMSE, melhor o conjunto de parametros descreve
## os dados observados. Util para comparar objetivamente se o cenario BASE
## (ou alguma sensibilidade) se ajusta melhor aos dados.
rmse_resumo <- todos_resultados %>%
  filter(is.na(erro), !is.na(RMSE)) %>%
  group_by(cenario, eixo_sensibilidade, valor_testado) %>%
  summarise(RMSE_medio = mean(RMSE, na.rm = TRUE),
            RMSE_min   = min(RMSE, na.rm = TRUE),
            RMSE_max   = max(RMSE, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(RMSE_medio)

write.csv(rmse_resumo, file.path(output_dir, "LBSPR_RMSE_por_cenario.csv"), row.names = FALSE)
cat("RMSE medio por cenario (composicao de comprimento observada vs. prevista):\n")
print(rmse_resumo)
cat("\n")

if (all(is.na(todos_resultados$RMSE))) {
  cat("\u26a0 RMSE nao pode ser calculado — confira se o slot '@pLCatch' existe na\n")
  cat("  sua versao do LBSPR (rode slotNames(fit_base) para verificar o nome\n")
  cat("  correto do slot de composicao de comprimento PREVISTA).\n\n")
}

# =============================================================================
# ---- 13. GRAFICOS -------------------------------------------------------------
# =============================================================================

## 13a-preview. RMSE por cenario (grafico de barras) ----------------------------
p_rmse <- ggplot(rmse_resumo, aes(x = reorder(cenario, RMSE_medio), y = RMSE_medio, fill = eixo_sensibilidade)) +
  geom_col() +
  coord_flip() +
  labs(title = "RMSE medio por cenario — ajuste da composicao de comprimento",
       subtitle = "Menor RMSE = melhor ajuste do modelo aos dados observados",
       x = NULL, y = "RMSE medio (proporcao por classe de comprimento)", fill = "Eixo") +
  theme_minimal(base_size = 12)
p_rmse
ggsave(file.path(output_dir, "RMSE_por_cenario.png"), p_rmse, width = 10, height = 9, dpi = 150)
cat("Grafico de RMSE salvo em:", file.path(output_dir, "RMSE_por_cenario.png"), "\n\n")

## 13a. SPR ao longo do tempo, por eixo de sensibilidade -----------------------
eixos <- unique(todos_resultados$eixo_sensibilidade[todos_resultados$eixo_sensibilidade != "base"])

base_repetido <- todos_resultados %>% filter(cenario == "Base_PosPolitica") %>%
  tidyr::crossing(eixo_facet = eixos) %>%
  mutate(eixo_sensibilidade = eixo_facet) %>% select(-eixo_facet)

plot_dat <- bind_rows(
  todos_resultados %>% filter(is.na(erro), cenario != "Base_PosPolitica"),
  base_repetido
)
plot_dat
p_time <- ggplot(plot_dat, aes(x = Ano, y = SPR, color = cenario)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0.40, linetype = "dashed", color = "darkgreen", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.20, linetype = "dashed", color = "firebrick", inherit.aes = FALSE) +
  facet_wrap(~eixo_sensibilidade, scales = "free_y") +
  labs(title = "Sensibilidade do SPR (LBSPR) — crescimento, maturidade, largura de classe e M/K",
       subtitle = "Base = pos-politica (2017-2021); alternativa = pre-politica (2003-2007), da Cruz Delgado et al. (2024)",
       y = "Spawning Potential Ratio (SPR)", x = "Ano", color = "Cenario") +
  theme_minimal(base_size = 12)
p_time

ggsave(file.path(output_dir, "sensibilidade_SPR_por_ano.png"), p_time, width = 14, height = 9, dpi = 150)

## 13b. Grafico "tornado" — impacto no SPR medio --------------------------------
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
       subtitle = paste0("Referencia (cenario BASE, pos-politica): SPR medio = ", round(spr_medio_base, 3)),
       x = NULL, y = "Variacao do SPR medio em relacao ao cenario BASE", fill = "Cenario") +
  theme_minimal(base_size = 12)
p_tornado

ggsave(file.path(output_dir, "sensibilidade_tornado_SPR_medio.png"), p_tornado, width = 10, height = 8, dpi = 150)

cat("Graficos salvos em:", output_dir, "\n")
cat("  - sensibilidade_SPR_por_ano.png\n")
cat("  - sensibilidade_tornado_SPR_medio.png\n\n")

# =============================================================================
# NOTAS FINAIS
# -----------------------------------------------------------------------------
# 1. O cenario BASE usa o periodo POS-POLITICA (2017-2021) de da Cruz Delgado
#    et al. (2024) para Linf, K (via M/K) e L50; L95 e estimado combinando o
#    L50 do artigo com a forma (razao L95/L50) de uma ogiva ajustada aos
#    dados atuais no mesmo periodo.
# 2. Cada sensibilidade OAT compara o cenario BASE com o periodo PRE-POLITICA
#    (2003-2007) do MESMO artigo — nao com uma faixa arbitraria de literatura.
#    Isso reflete diretamente o antes/depois documentado pelo artigo-base.
# 3. M e M/K nao sao reportados no artigo; foram estimados pela equacao de
#    Pauly (1980) a partir do proprio Linf/K de cada periodo e das SST citadas
#    na Discussao do artigo (24,1 C e 24,7 C). Revise essas SST se dispuser de
#    valores mais especificos para o periodo/area de amostragem.
# 4. A largura de classe (bin width) e a unica sensibilidade que nao vem do
#    artigo — e puramente metodologica (1 cm vs. 2 cm vs. 3 cm).
# 5. Repita as Secoes 4-5 se os valores do artigo forem revisados, ou se
#    quiser trocar as SST assumidas na equacao de Pauly.
# 6. SPR baixo/proximo de 0 nos primeiros anos da serie NAO e necessariamente
#    um erro. O LBSPR infere a mortalidade por pesca a partir de quanto a
#    distribuicao de comprimentos observada esta "truncada" (poucos peixes
#    grandes) em relacao ao esperado sem pesca. O MLS (tamanho minimo de
#    captura) e o BRP (periodo de defeso) só passaram a vigorar a partir de
#    2008 (Resolucao 11/2007) — antes disso, nao havia tamanho minimo
#    legal, entao e esperado que os desembarques incluíssem mais peixes
#    pequenos/imaturos, o que o LBSPR interpreta como alta pressao de pesca
#    sobre a fracao desovante (SPR baixo). Isso e, na verdade, consistente
#    com o proprio achado do artigo-base: o tamanho de desembarque aumentou
#    apos a intervencao de politica. Vale, ainda assim, checar
#    'diagnostico_n_por_ano.csv' — anos com poucas amostras podem produzir
#    estimativas de SPR instaveis por baixo tamanho amostral, nao so por
#    sinal biologico real. Compare os primeiros anos da serie (nos graficos
#    de saida) com o numero de medicoes disponiveis antes de interpretar.
# 7. O RMSE (Secao 12a) compara a composicao de comprimento PREVISTA por cada
#    cenario com a OBSERVADA nos dados — nao deve ser usado sozinho para
#    "escolher" o cenario mais correto biologicamente (um M/K muito alto,
#    por exemplo, pode ajustar melhor o formato da distribuicao observada
#    sem que isso signifique que o valor e biologicamente realista). Use o
#    RMSE como um diagnostico complementar, nao como criterio unico.
# =============================================================================

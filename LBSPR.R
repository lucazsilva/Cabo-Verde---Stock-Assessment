#-------------------------------------------------------------------------------------------#
#     Este script contem a avaliacao LBSPR de Decapterus macarellus em Cabo Verde          #
#   Analises usam a serie de comprimentos completa, desde 1988 ate 2024                    #
#   Modelo principal: Length-Based Spawning Potential Ratio - LBSPR                         #
#   Parametros biologicos avaliados: Linf, M/K, L50 e L95                                  #
#   Maturidade dos dados atuais: somente estadio I = imaturo; II-VII = adultos/maduros      #
#   Modelo BASE: Vieira (2019, crescimento/M-K) + Costa et al. (2020, maturidade)           #
#   Anos com menos de 150 individuos amostrados sao excluidos da serie principal            #
#   Sensibilidades one-at-a-time: crescimento, maturidade, largura de classe e M/K          #
#   O objetivo e representar alternativas biologicamente plausiveis sem ampliar             #
#   artificialmente a incerteza por combinacoes fatoriais extremas                          #
#   Codificacao criada por Silva, LVS ; 10/09/2026, Instituto do Mar - IMar, Mindelo       #
#-------------------------------------------------------------------------------------------#
#
# BASE DE DADOS DE COMPRIMENTOS
#   Arquivo:      Base_Comprimentos_Combinada_1988_2024.xlsx
#   Serie principal (matriz de comprimento para o LBSPR): Ano entre 1988 e 2024,
#     excluindo anos com menos de MIN_N_ANO individuos amostrados
#   Subconjuntos adicionais (so para as sensibilidades de crescimento/M-K/maturidade,
#     replicando os dois periodos do artigo da Cruz Delgado et al. 2024):
#     2003-2007 (pre-politica) e 2017-2021 (pos-politica) — usam a base completa
#   Colunas:      Ano, L(cm), Sexo, Maturidade (EM)
#   Sexo:         somente femeas (Sexo == "F")
#   Classe base:  largura de 1 cm (ver sensibilidade "largura de classe")
#
# CENARIO BASE
#   Crescimento e M/K: Vieira, N. (2019). Stock assessment and the influence
#     of environmental parameters on the distribution of mackerel scad
#     (Decapterus macarellus) in Cabo Verde waters. United Nations University
#     Fisheries Training Programme.
#     Linf = 40,6 cm FL | K = 0,45 /ano | M = 0,92 /ano | M/K = 2,04
#   Maturidade: Costa, M.P.V., Cruz, D.R.S., Monteiro, L.S., Evora, K.S.M., &
#     Cardoso, L.G. (2020). Reproductive biology of the mackerel scad
#     Decapterus macarellus from Cabo Verde and the implications for its
#     fishery management. African Journal of Marine Science, 42(1), 35-42.
#     Femeas: L50 = 24,1 cm | L95 = 27,8 cm
#
# FONTES DAS SENSIBILIDADES (planilha Parametros_Historia_de_vida.xlsx)
#   - da Cruz Delgado et al. (2024): periodos pre-politica (2003-2007) e
#     pos-politica (2017-2021) do mesmo artigo de avaliacao de politicas —
#     usados nos eixos de crescimento, M/K e maturidade. M/K desses periodos
#     nao e reportado no artigo; foi estimado pela equacao de Pauly (1980) a
#     partir do Linf/K de cada periodo e das SST citadas na Discussao do
#     artigo (24,1 C e 24,7 C, respectivamente).
#   - Fontes de confiabilidade "Alta" (aba Confiabilidade_Fontes): Jardim
#     (1996/1999) — crescimento e M/K.
#   - Fontes de confiabilidade "Média": Carvalho & Caramelo (1996/1999),
#     Santos (2018), Almada (1997), da Luz & Vieira (2020) — crescimento e/ou
#     M/K e/ou maturidade, conforme disponibilidade de cada fonte. Vieira
#     (2019) tambem contribui um L50 alternativo (maturidade), distinto do
#     seu proprio papel como fonte de crescimento/M-K do cenario BASE.
#   - Fontes de confiabilidade "Baixa" (Stobberup & Erzini 2006, Tariche &
#     Martins 2009) NAO sao usadas.
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
#   - Nomes de slots ja variaram entre versoes do LBSPR (ex.: @maxL nao existe
#     em todas as versoes). Se algo falhar, rode slotNames() no objeto
#     envolvido e confira contra a documentacao da sua versao instalada.
# =============================================================================
setwd("C:/Users/lucas/OneDrive/Desktop/Consultorias/Brasil - Sustentamares/Cabo Verde/Análises/Cabo-Verde---Stock-Assessment/LBSPR_Hordyk2015")   # onde está o .xlsx
# ---- 1. Pacotes -------------------------------------------------------------
required_pkgs <- c("LBSPR", "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "purrr")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---- 2. Arquivos / parametros gerais -----------------------------------------
length_file <- "Base_Comprimentos_Combinada_1988_2024.xlsx"
output_dir  <- "LBSPR_output"
dir.create(output_dir, showWarnings = FALSE)

ANO_MIN       <- 1988   # serie completa
ANO_MAX       <- 2024
sex_filter    <- "F"
BinWidth_base <- 1        # cm
MIN_N_ANO     <- 150      # individuos minimos por ano para entrar na serie principal
sat_incluido  <- FALSE    # "SAT" nao pertence a escala I-VII; excluido do ajuste da ogiva

# =============================================================================
# ---- 3. DADOS DE COMPRIMENTO -------------------------------------------------
# =============================================================================
raw_len <- read_excel(length_file)   # base completa (1988-2024); tambem usada
# para os subconjuntos pre/pos-politica

df_len_bruto <- raw_len %>%
  filter(!is.na(`L(cm)`), `L(cm)` > 0, Ano >= ANO_MIN, Ano <= ANO_MAX)
if (sex_filter != "ALL") df_len_bruto <- df_len_bruto %>% filter(Sexo == sex_filter)

## -- 3a. Diagnostico e exclusao de anos com amostra insuficiente ---------------
n_por_ano <- df_len_bruto %>% count(Ano, name = "n_medicoes")
write.csv(n_por_ano, file.path(output_dir, "diagnostico_n_por_ano.csv"), row.names = FALSE)

anos_excluidos <- n_por_ano %>% filter(n_medicoes < MIN_N_ANO)
cat("Medicoes de comprimento (", ANO_MIN, "-", ANO_MAX, ", sexo=", sex_filter, "), antes do filtro: ",
    nrow(df_len_bruto), "\n", sep = "")

if (nrow(anos_excluidos) > 0) {
  cat("Anos excluidos por amostra insuficiente (< ", MIN_N_ANO, " individuos):\n", sep = "")
  print(anos_excluidos)
  df_len <- df_len_bruto %>% filter(!(Ano %in% anos_excluidos$Ano))
} else {
  cat("Nenhum ano abaixo do limiar de", MIN_N_ANO, "individuos.\n")
  df_len <- df_len_bruto
}

cat("Medicoes de comprimento apos o filtro: ", nrow(df_len),
    " | Anos incluidos: ", paste(range(df_len$Ano), collapse = "-"), "\n\n", sep = "")

# =============================================================================
# ---- 4. CRESCIMENTO E M/K — CENARIO BASE (Vieira, 2019) ---------------------
# =============================================================================
base_Linf <- 40.6           # Vieira (2019): Linf (cm)
base_MK   <- 0.92 / 0.45    # Vieira (2019): M/K (M=0,92, K=0,45)

# =============================================================================
# ---- 5. MATURIDADE — CENARIO BASE (Costa et al., 2020, femeas) --------------
# =============================================================================
base_L50 <- 24.1   # cm, Costa et al. (2020), femeas
base_L95 <- 27.8   # cm, Costa et al. (2020), femeas

# =============================================================================
# ---- 6. PARAMETROS DAS SENSIBILIDADES ----------------------------------------
# =============================================================================

## -- 6a. da Cruz Delgado et al. (2024): pre e pos-politica --------------------
## Linf e K reportados diretamente no artigo. M nao e reportado; estimado
## pela equacao empirica de Pauly (1980):
##   log10(M) = -0.0066 - 0.279*log10(Linf) + 0.6543*log10(K) + 0.4634*log10(T)
pauly_M <- function(Linf, K, T) {
  10^(-0.0066 - 0.279 * log10(Linf) + 0.6543 * log10(K) + 0.4634 * log10(T))
}

crescimento_pre <- list(periodo = "2003-2007 (pre-politica)", Linf = 41.48, K = 0.39, SST = 24.1)
crescimento_pos <- list(periodo = "2017-2021 (pos-politica)", Linf = 46.00, K = 0.47, SST = 24.7)

M_pre <- pauly_M(crescimento_pre$Linf, crescimento_pre$K, crescimento_pre$SST)
M_pos <- pauly_M(crescimento_pos$Linf, crescimento_pos$K, crescimento_pos$SST)
MK_pre <- M_pre / crescimento_pre$K
MK_pos <- M_pos / crescimento_pos$K

L50_pre <- 24.2   # cm, Haddon logistic model, da Cruz Delgado et al. (2024)
L50_pos <- 28.2   # cm

## L95 nao e reportado no artigo (so L50). Estimado ajustando uma ogiva de
## maturidade aos dados atuais (femeas, mesmo periodo do artigo) e aplicando
## a razao L95/L50 dessa ogiva sobre o L50 REPORTADO PELO ARTIGO — a escala
## (L50) vem do artigo, a forma da ogiva (razao L95/L50) vem dos dados.
fit_ogiva <- function(df_full, ano_min, ano_max) {
  d <- df_full %>%
    filter(Sexo == "F", !is.na(Maturidade), Ano >= ano_min, Ano <= ano_max) %>%
    mutate(estagio = as.character(Maturidade)) %>%
    filter(estagio %in% c(as.character(1:7), if (sat_incluido) "SAT")) %>%
    mutate(Madura = if_else(estagio == "1", 0, 1))
  modelo <- glm(Madura ~ `L(cm)`, data = d, family = binomial)
  a <- unname(coef(modelo)[1]); b <- unname(coef(modelo)[2])
  list(L50 = -a / b, L95 = (log(0.95 / 0.05) - a) / b, n = nrow(d))
}

ogiva_pre <- fit_ogiva(raw_len, 2003, 2007)
ogiva_pos <- fit_ogiva(raw_len, 2017, 2021)

cat("Ogivas de maturidade ajustadas aos dados atuais (femeas, validacao vs. artigo):\n")
cat("  Pre-politica : L50_dados=", round(ogiva_pre$L50, 2), "| n=", ogiva_pre$n,
    "| (artigo reporta L50=", L50_pre, ")\n")
cat("  Pos-politica : L50_dados=", round(ogiva_pos$L50, 2), "| n=", ogiva_pos$n,
    "| (artigo reporta L50=", L50_pos, ")\n")
cat("  NOTA: o L50 ajustado aos dados atuais fica sistematicamente abaixo do L50\n")
cat("  do artigo (possivel diferenca de criterio de estadiamento) — por isso usamos\n")
cat("  o L50 do ARTIGO como valor absoluto e so a razao L95/L50 dos dados.\n\n")

ratio_L95_L50_pre <- ogiva_pre$L95 / ogiva_pre$L50
ratio_L95_L50_pos <- ogiva_pos$L95 / ogiva_pos$L50
L95_pre <- L50_pre * ratio_L95_L50_pre
L95_pos <- L50_pos * ratio_L95_L50_pos

## -- 6b. Jardim (1996/1999) [Alta]: crescimento e M/K --------------------------
jardim_Linf       <- 31.5
jardim_MK_tanaka  <- 1.00
jardim_MK_pauly   <- 1.86

## -- 6c. Fontes de confiabilidade "Média" --------------------------------------
carvalho_caramelo_Linf <- 41.0    # Carvalho & Caramelo (1996/1999)
santos_Linf            <- 41.83  # Santos (2018)
almada_Linf            <- 30.1    # Almada (1997)
da_luz_vieira_Linf     <- 40.9    # da Luz & Vieira (2020)

carvalho_caramelo_MK <- 1.00                 # M via Pauly / K=0,30
santos_MK_tanaka     <- 1.54                 # M via Tanaka / K=0,39
santos_MK_pauly      <- 1.71                 # M via Pauly / K=0,39
almada_MK_pauly      <- 0.64 / 0.34          # M via Pauly / K=0,34
almada_MK_tanaka     <- 0.43 / 0.34          # M via Tanaka / K=0,34
da_luz_vieira_MK     <- 1.59                 # mediana entre anos, metodo LBB

## Maturidade: so as fontes que reportam L50 podem entrar neste eixo.
## Carvalho & Caramelo e da Luz & Vieira NAO reportam L50 (esta ultima so
## reporta comprimento de SELETIVIDADE, nao de maturidade) — ficam de fora.
## L95 nao reportado por nenhuma delas -> estimado com a razao L95/L50 da
## ogiva pos-politica ajustada aos dados atuais (mesma logica da Secao 6a).
santos_L50 <- 22.9;  santos_L95 <- santos_L50 * ratio_L95_L50_pos
almada_L50 <- 21.44; almada_L95 <- almada_L50 * ratio_L95_L50_pos   # femeas
vieira_L50 <- 20.3;  vieira_L95 <- vieira_L50 * ratio_L95_L50_pos   # sexos combinados

# =============================================================================
# ---- 7. CLASSES DE COMPRIMENTO (largura base = 1 cm) -------------------------
# Estende ate cobrir 1.25 x o maior Linf testado (bom-senso do LBSPR).
# =============================================================================
maxL_mult <- 1.25
todos_Linf_testados <- c(base_Linf, crescimento_pre$Linf, crescimento_pos$Linf,
                         jardim_Linf, carvalho_caramelo_Linf, santos_Linf,
                         almada_Linf, da_luz_vieira_Linf)
maior_Linf_testado <- max(todos_Linf_testados)

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
# ---- 8. FUNCOES AUXILIARES: RMSE + execucao do LBSPR -------------------------
# O LBSPR guarda a composicao de comprimento OBSERVADA em @LData (contagens) e
# a composicao PREVISTA pelo modelo ajustado em @pLCatch (proporcoes). Aqui
# normalizamos @LData para proporcoes por ano e comparamos com @pLCatch.
# =============================================================================
compute_rmse <- function(fit) {
  obs_prop <- tryCatch(sweep(fit@LData, 2, colSums(fit@LData), "/"), error = function(e) NULL)
  pred_prop <- tryCatch(fit@pLCatch, error = function(e) NULL)
  if (is.null(obs_prop) || is.null(pred_prop) || !identical(dim(obs_prop), dim(pred_prop))) {
    return(data.frame(Ano = fit@Years, RMSE = NA_real_))
  }
  data.frame(Ano = fit@Years, RMSE = as.numeric(sqrt(colMeans((obs_prop - pred_prop)^2, na.rm = TRUE))))
}

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
# ---- 9. CENARIO BASE (Vieira 2019 + Costa et al. 2020) -----------------------
# =============================================================================
resultado_base <- run_lbspr("Base_VieiraCosta", "base", NA,
                            Linf = base_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
                            BinWidth = BinWidth_base, Lengths = MyLengths_base)
res_base <- resultado_base$df
fit_base <- resultado_base$fit

if (is.null(fit_base)) {
  stop("O ajuste do cenario BASE falhou (", res_base$erro[1], "). Corrija antes de prosseguir.")
}

## -- 9a. GRAFICOS NATIVOS DO PACOTE LBSPR (Hordyk et al.) para o cenario BASE --
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
# ---- 10. SENSIBILIDADES ONE-AT-A-TIME (OAT) ----------------------------------
# Cada eixo varia UM parametro por vez, mantendo os demais no cenario BASE.
# =============================================================================

## 10a. Crescimento (Linf) -------------------------------------------------------
res_crescimento <- bind_rows(
  run_lbspr_df("Crescimento_PrePolitica", "crescimento", crescimento_pre$Linf,
               Linf = crescimento_pre$Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Crescimento_PosPolitica", "crescimento", crescimento_pos$Linf,
               Linf = crescimento_pos$Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
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
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

## 10b. Maturidade (L50 e L95 variam juntos, mesma forma de ogiva) ---------------
res_maturidade <- bind_rows(
  run_lbspr_df("Maturidade_PrePolitica", "maturidade", L50_pre,
               Linf = base_Linf, MK = base_MK, L50 = L50_pre, L95 = L95_pre,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("Maturidade_PosPolitica", "maturidade", L50_pos,
               Linf = base_Linf, MK = base_MK, L50 = L50_pos, L95 = L95_pos,
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

## 10c. Largura de classe (bin width) --------------------------------------------
BinWidths_alt <- c(2, 3)  # cm — alternativas a 1 cm do cenario BASE
res_binwidth <- purrr::map_dfr(BinWidths_alt, function(bw) {
  Lengths_bw <- build_lengths(df_len, bw)
  run_lbspr_df(paste0("LarguraClasse_", bw, "cm"), "largura_classe", bw,
               Linf = base_Linf, MK = base_MK, L50 = base_L50, L95 = base_L95,
               BinWidth = bw, Lengths = Lengths_bw)
})

## 10d. M/K -----------------------------------------------------------------------
res_mk <- bind_rows(
  run_lbspr_df("MK_PrePolitica", "M_K", MK_pre,
               Linf = base_Linf, MK = MK_pre, L50 = base_L50, L95 = base_L95,
               BinWidth = BinWidth_base, Lengths = MyLengths_base),
  run_lbspr_df("MK_PosPolitica", "M_K", MK_pos,
               Linf = base_Linf, MK = MK_pos, L50 = base_L50, L95 = base_L95,
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
               BinWidth = BinWidth_base, Lengths = MyLengths_base)
)

## 10e. Tabela unica com TODOS os parametros usados em cada sensibilidade -------
tabela_MK <- tibble(
  cenario = c("Base_VieiraCosta", "MK_PrePolitica", "MK_PosPolitica",
              "MK_Jardim_Tanaka", "MK_Jardim_Pauly", "MK_CarvalhoCaramelo",
              "MK_Santos_Tanaka", "MK_Santos_Pauly", "MK_Almada_Pauly",
              "MK_Almada_Tanaka", "MK_daLuzVieira_2020"),
  fonte = c("Vieira (2019)",
            "da Cruz Delgado et al. (2024), pre-politica",
            "da Cruz Delgado et al. (2024), pos-politica",
            "Jardim (1996/1999)", "Jardim (1996/1999)",
            "Carvalho & Caramelo (1996/1999)",
            "Santos (2018)", "Santos (2018)",
            "Almada (1997)", "Almada (1997)",
            "da Luz & Vieira (2020)"),
  metodo_M = c("reportado (m e k na mesma linha)",
               "Pauly (1980), via Linf/K do artigo + SST=24,1C",
               "Pauly (1980), via Linf/K do artigo + SST=24,7C",
               "Tanaka", "Pauly", "Pauly", "Tanaka", "Pauly", "Pauly", "Tanaka",
               "LBB (mediana entre anos)"),
  confiabilidade = c("Média (base)", "Alta (artigo-base)", "Alta (artigo-base)",
                     "Alta", "Alta", "Média", "Média", "Média", "Média", "Média", "Média"),
  MK = c(base_MK, MK_pre, MK_pos, jardim_MK_tanaka, jardim_MK_pauly,
         carvalho_caramelo_MK, santos_MK_tanaka, santos_MK_pauly,
         almada_MK_pauly, almada_MK_tanaka, da_luz_vieira_MK)
) %>% arrange(MK)

write.csv(tabela_MK, file.path(output_dir, "MK_parametros_utilizados.csv"), row.names = FALSE)
cat("Tabela de todos os valores de M/K utilizados (base + sensibilidades):\n")
print(tabela_MK); cat("\n")

tabela_sensibilidades <- bind_rows(
  tibble(eixo = "crescimento", cenario = "Crescimento_PrePolitica",
         fonte = "da Cruz Delgado et al. (2024), pre-politica",
         parametro = "Linf_cm", valor = crescimento_pre$Linf, confiabilidade = "Alta (artigo-base)"),
  tibble(eixo = "crescimento", cenario = "Crescimento_PosPolitica",
         fonte = "da Cruz Delgado et al. (2024), pos-politica",
         parametro = "Linf_cm", valor = crescimento_pos$Linf, confiabilidade = "Alta (artigo-base)"),
  tibble(eixo = "crescimento", cenario = "Crescimento_Jardim_1996_1999",
         fonte = "Jardim (1996/1999)", parametro = "Linf_cm", valor = jardim_Linf, confiabilidade = "Alta"),
  tibble(eixo = "crescimento", cenario = "Crescimento_CarvalhoCaramelo",
         fonte = "Carvalho & Caramelo (1996/1999)", parametro = "Linf_cm", valor = carvalho_caramelo_Linf, confiabilidade = "Média"),
  tibble(eixo = "crescimento", cenario = "Crescimento_Santos_2018",
         fonte = "Santos (2018)", parametro = "Linf_cm", valor = santos_Linf, confiabilidade = "Média"),
  tibble(eixo = "crescimento", cenario = "Crescimento_Almada_1997",
         fonte = "Almada (1997)", parametro = "Linf_cm", valor = almada_Linf, confiabilidade = "Média"),
  tibble(eixo = "crescimento", cenario = "Crescimento_daLuzVieira_2020",
         fonte = "da Luz & Vieira (2020)", parametro = "Linf_cm", valor = da_luz_vieira_Linf, confiabilidade = "Média"),
  
  tabela_MK %>% filter(cenario != "Base_VieiraCosta") %>%
    transmute(eixo = "M_K", cenario, fonte, parametro = "M_K", valor = MK, confiabilidade),
  
  tibble(eixo = "maturidade", cenario = "Maturidade_PrePolitica",
         fonte = "da Cruz Delgado et al. (2024), pre-politica",
         parametro = c("L50_cm", "L95_cm"), valor = c(L50_pre, L95_pre), confiabilidade = "Alta (artigo-base)"),
  tibble(eixo = "maturidade", cenario = "Maturidade_PosPolitica",
         fonte = "da Cruz Delgado et al. (2024), pos-politica",
         parametro = c("L50_cm", "L95_cm"), valor = c(L50_pos, L95_pos), confiabilidade = "Alta (artigo-base)"),
  tibble(eixo = "maturidade", cenario = "Maturidade_Santos_2018",
         fonte = "Santos (2018)", parametro = c("L50_cm", "L95_cm"),
         valor = c(santos_L50, santos_L95), confiabilidade = "Média"),
  tibble(eixo = "maturidade", cenario = "Maturidade_Almada_1997_Femea",
         fonte = "Almada (1997)", parametro = c("L50_cm", "L95_cm"),
         valor = c(almada_L50, almada_L95), confiabilidade = "Média"),
  tibble(eixo = "maturidade", cenario = "Maturidade_Vieira_2019",
         fonte = "Vieira (2019)", parametro = c("L50_cm", "L95_cm"),
         valor = c(vieira_L50, vieira_L95), confiabilidade = "Média"),
  
  tibble(eixo = "largura_classe", cenario = paste0("LarguraClasse_", BinWidths_alt, "cm"),
         fonte = "metodologico (nao vem da literatura)",
         parametro = "BinWidth_cm", valor = BinWidths_alt, confiabilidade = NA_character_),
  
  tibble(eixo = "base", cenario = "Base_VieiraCosta",
         fonte = "Vieira (2019) [Linf, M/K] + Costa et al. (2020), femeas [L50, L95]",
         parametro = c("Linf_cm", "M_K", "L50_cm", "L95_cm"),
         valor = c(base_Linf, base_MK, base_L50, base_L95),
         confiabilidade = "Alta (artigo-base)")
)

write.csv(tabela_sensibilidades, file.path(output_dir, "tabela_todos_parametros_sensibilidade.csv"), row.names = FALSE)
cat("Tabela consolidada de todos os parametros de sensibilidade salva em:\n")
cat(" ", file.path(output_dir, "tabela_todos_parametros_sensibilidade.csv"), "\n\n")
print(tabela_sensibilidades, n = Inf); cat("\n")

# =============================================================================
# ---- 11. CONSOLIDAR E EXPORTAR -----------------------------------------------
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

## -- 11a. Resumo do RMSE por cenario -------------------------------------------
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
print(rmse_resumo); cat("\n")

if (all(is.na(todos_resultados$RMSE))) {
  cat("\u26a0 RMSE nao pode ser calculado — confira se o slot '@pLCatch' existe na\n")
  cat("  sua versao do LBSPR (rode slotNames(fit_base) para verificar).\n\n")
}

# =============================================================================
# ---- 12. GRAFICOS -------------------------------------------------------------
# =============================================================================

## 12a. RMSE por cenario (grafico de barras) -------------------------------------
p_rmse <- ggplot(rmse_resumo, aes(x = reorder(cenario, RMSE_medio), y = RMSE_medio, fill = eixo_sensibilidade)) +
  geom_col() +
  coord_flip() +
  labs(title = "RMSE medio por cenario — ajuste da composicao de comprimento",
       subtitle = "Menor RMSE = melhor ajuste do modelo aos dados observados",
       x = NULL, y = "RMSE medio (proporcao por classe de comprimento)", fill = "Eixo") +
  theme_minimal(base_size = 12)
p_rmse
ggsave(file.path(output_dir, "RMSE_por_cenario.png"), p_rmse, width = 10, height = 9, dpi = 150)

## 12b. SPR ao longo do tempo, por eixo de sensibilidade --------------------------
eixos <- unique(todos_resultados$eixo_sensibilidade[todos_resultados$eixo_sensibilidade != "base"])

base_repetido <- todos_resultados %>% filter(cenario == "Base_VieiraCosta") %>%
  tidyr::crossing(eixo_facet = eixos) %>%
  mutate(eixo_sensibilidade = eixo_facet) %>% select(-eixo_facet)

plot_dat <- bind_rows(
  todos_resultados %>% filter(is.na(erro), cenario != "Base_VieiraCosta"),
  base_repetido
)
plot_dat

p_time <- ggplot(plot_dat, aes(x = Ano, y = SPR, color = cenario)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = 0.40, linetype = "dashed", color = "darkgreen", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.20, linetype = "dashed", color = "firebrick", inherit.aes = FALSE) +
  facet_wrap(~eixo_sensibilidade) +
  labs(title = "Sensibilidade do SPR (LBSPR) — crescimento, maturidade, largura de classe e M/K",
       subtitle = "Base = Vieira (2019, crescimento/M-K) + Costa et al. (2020, maturidade)",
       y = "Spawning Potential Ratio (SPR)", x = "Ano", color = "Cenario") +
  theme_minimal(base_size = 12)
p_time
ggsave(file.path(output_dir, "sensibilidade_SPR_por_ano.png"), p_time, width = 14, height = 9, dpi = 150)

## 12c. Grafico "tornado" — impacto no SPR medio ----------------------------------
spr_medio_base <- mean(res_base$SPR, na.rm = TRUE)

tornado_dat <- todos_resultados %>%
  filter(is.na(erro), eixo_sensibilidade != "base") %>%
  group_by(eixo_sensibilidade, cenario) %>%
  summarise(SPR_medio = mean(SPR, na.rm = TRUE), .groups = "drop") %>%
  mutate(delta = SPR_medio - spr_medio_base)

p_tornado <- ggplot(tornado_dat, aes(x = cenario, y = delta, fill = eixo_sensibilidade)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "black") +
  coord_flip() +
  labs(title = "Grafico tornado — sensibilidade do SPR medio",
       subtitle = paste0("Referencia (cenario BASE, Vieira+Costa): SPR medio = ", round(spr_medio_base, 3)),
       x = NULL, y = "Variacao do SPR medio em relacao ao cenario BASE", fill = "Eixo") +
  theme_minimal(base_size = 12)
p_tornado
ggsave(file.path(output_dir, "sensibilidade_tornado_SPR_medio.png"), p_tornado, width = 10, height = 8, dpi = 150)

cat("Graficos salvos em:", output_dir, "\n")
cat("  - RMSE_por_cenario.png\n")
cat("  - sensibilidade_SPR_por_ano.png\n")
cat("  - sensibilidade_tornado_SPR_medio.png\n\n")

# =============================================================================
# NOTAS FINAIS
# -----------------------------------------------------------------------------
# 1. O cenario BASE usa Vieira (2019) para Linf e M/K (M=0,92, K=0,45,
#    M/K=2,04) e Costa et al. (2020), femeas, para L50/L95 (24,1cm / 27,8cm)
#    — ambos reportados diretamente nas fontes.
# 2. Anos com menos de MIN_N_ANO (150) individuos amostrados sao excluidos da
#    serie principal (Secao 3a) — confira 'diagnostico_n_por_ano.csv' para
#    ver exatamente quais anos entraram/sairam.
# 3. As sensibilidades OAT comparam o cenario BASE com: os periodos
#    PRE-POLITICA e POS-POLITICA de da Cruz Delgado et al. (2024), Jardim
#    (1996/1999) [Alta], e Carvalho & Caramelo, Santos, Almada e
#    da Luz & Vieira [Média]. Fontes "Baixa" nao sao usadas.
# 4. M/K do periodo pre/pos-politica nao e reportado no artigo da Cruz
#    Delgado et al. (2024); foi estimado pela equacao de Pauly (1980) a
#    partir do Linf/K de cada periodo e da SST citada na Discussao do artigo.
# 5. A largura de classe (bin width) e a unica sensibilidade que nao vem da
#    literatura — e puramente metodologica (1 cm vs. 2 cm vs. 3 cm).
# 6. O RMSE compara a composicao de comprimento PREVISTA por cada cenario com
#    a OBSERVADA nos dados — e um diagnostico de ajuste, nao um criterio de
#    escolha do cenario biologicamente "correto".
# 7. SPR baixo/proximo de 0 em anos com poucos peixes grandes na amostra nao
#    e necessariamente erro — pode refletir composicao de captura real antes
#    da vigencia do MLS/BRP (a partir de 2008). Compare com
#    'diagnostico_n_por_ano.csv' antes de interpretar anos especificos.
# =============================================================================

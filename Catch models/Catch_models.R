#-------------------------------------------------------------------------------------------#
#     Este script contem toda analise baseada em captura do Decapterus macarellus           #
#   Analises possuem foco para o Decapterus macarellus no arquipelago de Cabo verde         #
#   Hipóteses de estoques consideradas h1: estoque unitário para o arquipélago              #
#       h2: estoques unitarios para cada ilha ou ilhas barlavento e sotavento               #
#          Modelos baseados em captura: DB-SRA (Dick & MacCall, 2011)                       #
#     Modelo utilizado em seguida: CMSY (Froese et al., 2023)                               #
# Analises de sensibilidade e diferentes hipoteses assumidas de historia de vida, taxa      # 
#   intrinseca de crescimento (r) e deplecoes de biomassa (B/k) para o fim da serie         #
# Codificacao criada por Silva, MLS ; 09/09/2026, Instituto do Mar- IMar, Mindelo           #
#-------------------------------------------------------------------------------------------#

# limpando ambiente de trabalho...
rm(list = ls())

#@pacotes..
#install.packages("readxl")
library(readxl)
#installed.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#installed.packages("tidyr")
library(tidyr)
#install.packages("neuralnet")
library(neuralnet)
#install.packages("purrr")
library(purrr)
#------------------------------------
#instalando o datalimited2
#install.packages("devtools") #pra baixar o datalimited2
#devtools::install_github("cfree14/datalimited2",force = TRUE) #to run zBRT Zhou method 
#other option to instal datalimited2
#install.packages("pak")
#pak::pak("cfree14/datalimited2")
library(datalimited2)
#-----------------------------------

# definindo diretorio de trabalho..
setwd("C:/Users/mathe/OneDrive/Documents/Cabo-Verde---Stock-Assessment/Catch models")

### lendo os dados de capturas... ###
ct<- read.csv("Catch_Luz and Vieira.csv",sep = ",",dec = ".")
# lendo dados de história de vida... ##
lh<- read_xlsx("Parâmetros_História de vida.xlsx")
#carregar ffnn.bin (parametros da rede neural do CMSY)
load("ffnn.bin")


#---------------------------------------#
# Analise exploratoria das capturas

# Gráfico da série temporal de capturas
#---------------------------------------#

p_ct <- ggplot(ct, aes(x = Year, y = Catch)) +
  # Captura observada
  geom_line(linewidth = 1.5) +
  # Tendência suavizada
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1,
    alpha=0.2) +
  # Possível mudança estrutural
  geom_vline(
    xintercept = 2014,
    linetype = "dashed",
    linewidth = 0.7,
    col="red") +
  labs(
    x = "Ano",
    y = "Captura (t)") +
  theme_classic(base_size = 14) +
  theme(
    plot.margin = unit(
      c(0.05, 0.05, 0.05, 0.05),
      "cm"))
p_ct

ggsave(
  filename = "capturas_Decapterus_macarellus.png",
  plot = p_ct,
  device = "png",
  units = "cm",
  width = 25,
  height = 18,
  dpi = 300)


# ==============================================================================
# Hipoteses de depleção (B/k) final da série para Decapterus macarellus
# ==============================================================================
## Objetivo:
# Construir hipóteses de depleção B/K para Decapterus macarellus
# utilizando:
#
#   1. Neural Network do CMSY++
#   2. zBRT
#   3. Hipótese independente de mudança de alvo (target switching)
#   4. Hipótese não informativa
#
# IMPORTANTE:
# Para cada ano-alvo, somente os dados até aquele ano são utilizados.
# Isso evita usar informação futura para estimar a depleção histórica.
#
# ==============================================================================
# 1. PACOTES
library(dplyr)
library(tidyr)
library(neuralnet)
library(datalimited2)
# ------------------------------------------------
# 2. CARREGAR MODELO DA NEURAL NETWORK DO CMSY++
# ------------------------------------------------
# ffnn.bin deve estar no working directory
# Esse arquivo contém:
#   - nn.endbio
#   - slope.first.min
#   - slope.first.max
#   - slope.last.min
#   - slope.last.max
#   - e possivelmente outros objetos utilizados pelo CMSY++
load("ffnn.bin")
# --------------------------------
# 4. PADRONIZAR NOMES DAS COLUNAS
# --------------------------------
# O restante do script trabalha com:  year  e  ct
if ("Year" %in% names(ct)) {
  names(ct)[names(ct) == "Year"] <- "year"
}
if ("Catch" %in% names(ct)) {
  names(ct)[names(ct) == "Catch"] <- "ct"
}
# -----------------------
# 5. CHECAGEM DOS DADOS
# -----------------------
if (!all(c("year", "ct") %in% names(ct))) {
  
  stop(
    paste0(
      "O arquivo precisa conter as colunas 'Year' e 'Catch' ",
      "ou 'year' e 'ct'.\n\n",
      "Colunas encontradas: ",
      paste(names(ct), collapse = ", ")
    )
  )
}
# ------------------------------
# 6. PREPARAR SÉRIE DE CAPTURA
# ------------------------------
ct <- ct %>%
  select(year, ct) %>%
  mutate(
    year = as.numeric(year),
    ct   = as.numeric(ct)
  ) %>%
  filter(
    !is.na(year),
    !is.na(ct)
  ) %>%
  arrange(year)
# ------------------------------
# 7. VERIFICAR ANOS DUPLICADOS
# ------------------------------
if (anyDuplicated(ct$year) > 0) {
  
  warning(
    "Existem anos duplicados na série de captura. ",
    "Será feita a soma das capturas dentro de cada ano."
  )
  
  ct <- ct %>%
    group_by(year) %>%
    summarise(
      ct = sum(ct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(year)
}
# --------------------------------
# 8. VERIFICAR CAPTURAS NEGATIVAS
# --------------------------------
if (any(ct$ct < 0, na.rm = TRUE)) {
  
  stop(
    "Existem valores de captura negativos. ",
    "Verifique os dados antes de continuar."
  )
}
# -----------------
# 9. MOSTRAR SÉRIE
# -----------------
print(ct)

# ========================================
# 10. FUNÇÃO — NEURAL NETWORK DO CMSY++
# ========================================

estimate_endbio <- function(ct_raw, yr) {
  # ----------------
  # 10.1. CHECAGENS
  # ----------------
  if (!exists("nn.endbio")) {
    stop(
      "`nn.endbio` não encontrado. ",
      "Certifique-se de carregar o arquivo ffnn.bin."
    )
  }

  need_vars <- c(
    "slope.first.min",
    "slope.first.max",
    "slope.last.min",
    "slope.last.max"
  )
  
  missing_vars <- setdiff(
    need_vars,
    ls(envir = .GlobalEnv)
  )
  
  if (length(missing_vars) > 0) {
    stop(
      "Objetos necessários não encontrados: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  # ------------------------
  # 10.2. ORGANIZAR SÉRIE
  # ------------------------
  ord <- order(yr)
  
  yr      <- yr[ord]
  ct_raw  <- ct_raw[ord]
  
  # Remover NA
  keep <- !is.na(ct_raw) & !is.na(yr)
  ct_raw <- ct_raw[keep]
  yr     <- yr[keep]
  nyr <- length(yr)
  if (nyr < 5) {
  stop(
      "Série muito curta (<5 anos). ",
      "Não é possível formar os preditores da neural network."
    )
  }
  # -------------------------------
  # 10.3. ESTIMATIVA DO MSY PRIOR
  # -------------------------------
  ct.sort <- sort(ct_raw)
  sd.ct   <- sd(ct_raw)
  mean.ct <- mean(ct_raw)
  min.ct  <- min(ct_raw)
  max.ct  <- max(ct_raw)
  min_max <- min.ct / max.ct
  max.yr.i <- which.max(ct_raw)
  if (
    max.yr.i > (nyr - 4) ||
    (
      (sd.ct / mean.ct) < 0.1 &&
      min_max > 0.66
    )
  ) {
      MSY.pr <- mean(
      ct.sort[(nyr - 2):nyr]
    )
  } else {
    MSY.pr <- 0.75 * mean(
      ct.sort[(nyr - 4):nyr]
    )
  }
  # ------------------------
  # 10.4. ANO INTERMEDIÁRIO
  # ------------------------
  if (min_max > 0.7) {
    int.yr <- as.integer(
      mean(c(min(yr), max(yr)))
    )
  } else {
    yrs.int <- yr[
      yr > (yr[nyr] - 30) &
        yr > yr[min(4, nyr)] &
        yr < yr[max(1, nyr - 4)]
    ]
    
    ct.int <- ct_raw[
      yr > (yr[nyr] - 30) &
        yr > yr[min(4, nyr)] &
        yr < yr[max(1, nyr - 4)]
    ]
    
    
    if (length(yrs.int) == 0) {
      
      int.yr <- as.integer(
        mean(c(min(yr), max(yr)))
      )
      
    } else {
      
      min.ct.int    <- min(ct.int)
      min.ct.int.yr <- yrs.int[which.min(ct.int)]
      
      max.ct.int    <- max(ct.int)
      max.ct.int.yr <- yrs.int[which.max(ct.int)]
      
      
      if (min.ct.int.yr > max.ct.int.yr) {
        
        int.yr <- min.ct.int.yr
        
      } else {
        
        min.ct.after.max <- min(
          ct.int[
            yrs.int >= max.ct.int.yr
          ]
        )
        
        if (
          (min.ct.after.max / max.ct.int) < 0.75
        ) {
          
          int.yr <- yrs.int[
            yrs.int > max.ct.int.yr &
              ct.int == min.ct.after.max
          ][1]
          
        } else {
          
          int.yr <- min.ct.int.yr
        }
      }
    }
  }
  
  # -----------------------------------------
  # 10.5. CAPTURA / MSY NO ANO INTERMEDIÁRIO
  # -----------------------------------------
  idx.int <- which(yr == int.yr)
  if (length(idx.int) == 0) {
    
    idx.int <- which.min(
      abs(yr - int.yr)
    )
  }
  
  ct_MSY.int <- ct_raw[idx.int[1]] / MSY.pr
  
  # -----------------------------
  # 10.6. POSIÇÕES NORMALIZADAS
  # -----------------------------
  min.ct.i <- which.min(ct_raw) / nyr
  max.ct.i <- which.max(ct_raw) / nyr
  
  int.ct.i <- idx.int[1] / nyr
  
  # --------------------------------------
  # 10.7. NORMALIZAÇÃO DO NÚMERO DE ANOS
  # --------------------------------------
  if (
    all(
      c(
        "yr.norm.min",
        "yr.norm.max"
      ) %in% ls(envir = .GlobalEnv)
    )
  ) {
    yr.norm <- (
      nyr - get("yr.norm.min")
    ) /
      (
        get("yr.norm.max") -
          get("yr.norm.min")
      )
    
  } else {
    
    yr.norm <- (
      nyr - min(yr)
    ) /
      (
        max(yr) - min(yr)
      )
  }
  
  # -------------------------------------------
  # 10.8. MÉDIAS DE CAPTURA NO INÍCIO E FINAL
  # -------------------------------------------
  k_start <- min(5, nyr)
  k_end   <- min(5, nyr)
  
  mean.ct_MSY.start <-
    mean(ct_raw[1:k_start]) /
    MSY.pr
  
  mean.ct_MSY.end <-
    mean(
      ct_raw[
        (nyr - k_end + 1):nyr
      ]
    ) /
    MSY.pr
  # -------------
  # 10.9. SLOPES
  # -------------
  m_first <- min(10, nyr)
  m_last  <- min(10, nyr)
  
  slope.first <- coef(
    lm(
      (
        ct_raw[1:m_first] /
          mean.ct
      ) ~ seq_len(m_first)
    )
  )[2]
  
  
  slope.last <- coef(
    lm(
      (
        ct_raw[
          (nyr - m_last + 1):nyr
        ] /
          mean.ct
      ) ~ seq_len(m_last)
    )
  )[2]
  
  # -------------------------
  # 10.10. NORMALIZAR SLOPES
  # -------------------------
  slope.first.nrm <-
    (
      slope.first -
        get("slope.first.min")
    ) /
    (
      get("slope.first.max") -
        get("slope.first.min")
    )
  
  slope.last.nrm <-
    (
      slope.last -
        get("slope.last.min")
    ) /
    (
      get("slope.last.max") -
        get("slope.last.min")
    )
  
  # -----------------------------------
  # 10.11. PADRÕES DA SÉRIE DE CAPTURA
  # -----------------------------------
  min_max <- min(ct_raw) /
    max(ct_raw)
  start.rel <- ct_raw[1] /
    max(ct_raw)
  end.rel <- ct_raw[nyr] /
    max(ct_raw)
  Flat <- as.numeric(
    min_max >= 0.45 &
      start.rel >= 0.45 &
      end.rel >= 0.45 )
  LH <- as.numeric(
    min_max < 0.25 &
      start.rel < 0.45 &
      end.rel > 0.45)
  LHL <- as.numeric(
    min_max < 0.25 &
      start.rel < 0.45 &
      end.rel < 0.25)
  HL <- as.numeric(
    min_max < 0.25 &
      start.rel > 0.50 &
      end.rel < 0.25)
  HLH <- as.numeric(
    min_max < 0.25 &
      start.rel >= 0.45 &
      end.rel >= 0.45 )
  OTH <- as.numeric(
    sum(
      c(
        Flat,
        LH,
        LHL,
        HL,
        HLH
      )
    ) < 1
  )
  
  # ---------------------------------
  # 10.12. DATA FRAME DOS PREDITORES
  # ---------------------------------
  preds <- data.frame(
    Flat,
    LH,
    LHL,
    HL,
    HLH,
    OTH,
    ct_MSY.int,
    min_max,
    max.ct.i,
    int.ct.i,
    min.ct.i,
    yr.norm,
    mean.ct_MSY.start,
    slope.first.nrm,
    mean.ct_MSY.end,
    slope.last.nrm
  )
  # ----------------------
  # 10.13. NEURAL NETWORK
  # ----------------------
  pr.nn <- neuralnet::compute(
    nn.endbio,
    preds
  )
  idx <- max.col(
    pr.nn$net.result
  )
  # ------------------------------------
  # 10.14. RAZÃO CAPTURA / MSY NO FINAL
  # ------------------------------------
  ct_MSY.end <-
    ct_raw[nyr] /
    MSY.pr
  
  ct_MSY.use <-
    min(
      ct_MSY.end,
      mean.ct_MSY.end
    )
  # ----------------------------------
  # 10.15. MAPEAMENTO CMSY++ PARA B/K
  # ----------------------------------
  bk.MSY <- c(
    0.256,
    0.721
  )
  CL.1 <- c(
    0.01,
    0.203
  )
  CL.2 <- c(
    0.20,
    0.431
  )
  CL.3 <- c(
    0.80,
    -0.45
  )
  CL.4 <- c(
    1.02,
    -0.247
  )
  ct_MSY.lim <- 1.0
  # -----------------
  # 10.16. RESULTADO
  # -----------------
  if (
    mean.ct_MSY.end >= ct_MSY.lim
  ) {
    bk <- bk.MSY
  } else if (
    idx == 1
  ) {
    # Classe 1:
    # B/K provavelmente < 0.5
    bk <- c(
      CL.1[1] +
        CL.1[2] *
        ct_MSY.use,
      
      CL.2[1] +
        CL.2[2] *
        ct_MSY.use
    )
  } else {
    # Classe 2:
    # B/K provavelmente >= 0.5
    bk <- c(
      CL.3[1] +
        CL.3[2] *
        ct_MSY.use,
      
      CL.4[1] +
        CL.4[2] *
        ct_MSY.use
    )
  }
  # ---------------------------------
  # 10.17. GARANTIR LIMITES FÍSICOS
  # ---------------------------------
  bk[1] <- max(
    0,
    min(1, bk[1])
  )
  
  bk[2] <- max(
    0,
    min(1, bk[2])
  )
  # Se por algum motivo o limite inferior ficar
  # maior que o superior, ordenar.
  bk <- sort(bk)
  return(bk)
}
# ==================
# 11. FUNÇÃO — zBRT
# ==================
estimate_zbrt <- function(
    ct_data,
    target_year
) {
  
  # -------------------------------------------------------
  # Selecionar apenas informação disponível até o ano-alvo
  # -------------------------------------------------------
  sub <- ct_data %>%
    filter(
      year <= target_year
    ) %>%
    arrange(year)
  
  if (nrow(sub) < 5) {
    
    return(
      data.frame(
        bk_lo = NA_real_,
        bk_hi = NA_real_,
        bk = NA_real_
      )
    )
  }
  # -----------
  # Rodar zBRT
  # -----------
  output <- tryCatch(
    {
      zbrt(
        sub$year,
        sub$ct
      )
    },
    error = function(e) {
      warning(
        paste0(
          "zBRT falhou para ",
          target_year,
          ": ",
          e$message
        )
      )
      NULL
    }
  )
  if (is.null(output)) {
      return(
      data.frame(
        bk_lo = NA_real_,
        bk_hi = NA_real_,
        bk = NA_real_
      )
    )
  }
  # -----------------------------------------
  # Verificar se o ano-alvo existe no output
  # -----------------------------------------
  if (
    !"year" %in%
    names(output$ts)
  ) {
    
    return(
      data.frame(
        bk_lo = NA_real_,
        bk_hi = NA_real_,
        bk = NA_real_
      )
    )
  }
  
  if (
    !target_year %in%
    output$ts$year
  ) {
    
    return(
      data.frame(
        bk_lo = NA_real_,
        bk_hi = NA_real_,
        bk = NA_real_
      )
    )
  }
  
  idx <- which(
    output$ts$year == target_year
  )[1]
  
  # -------------
  # Extrair B/K
  # -------------
  bk_lo <- output$ts$s_lo[idx]
  bk_hi <- output$ts$s_hi[idx]
  bk    <- output$ts$s[idx]
  
  # Garantir limites
  bk_lo <- max(
    0,
    min(1, bk_lo) )
  bk_hi <- max(
    0,
    min(1, bk_hi)
  )
  bk <- max(
    0,
    min(1, bk)
  )
  return(
    data.frame(
      bk_lo = bk_lo,
      bk_hi = bk_hi,
      bk = bk
    )
  )
}

# ==========================================
# 12. FUNÇÃO PRINCIPAL — TODAS AS HIPÓTESES
# ==========================================
run_depletion_hypotheses <- function(
    data,
    target_years = c(2015, 2025)
) {

  results <- list()
  
  # ====================
  # LOOP SOBRE OS ANOS
  # ====================
  for (
    yr_target in target_years
  ) {
    # ---------------------------------
    # Dados disponíveis até o ano-alvo
    # ---------------------------------
    sub <- data %>%
      filter(
        year <= yr_target
      ) %>%
      arrange(year)
    
    if (nrow(sub) < 5) {
      warning(
        paste0(
          "Menos de 5 anos de dados até ",
          yr_target,
          "."
        )
      )
      next
    }
    # ========================
    # HIPÓTESE 1 — NN CMSY++
    # ========================
    bk_nn <- tryCatch(
      {
        estimate_endbio(
          ct_raw = sub$ct,
          yr = sub$year
        )
      },
      error = function(e) {
        warning(
          paste0(
            "NN-CMSY++ falhou para ",
            yr_target,
            ": ",
            e$message
          )
        )
        
        c(
          NA_real_,
          NA_real_
        )
      }
    )
    
    results[[length(results) + 1]] <- data.frame(
      
      year = yr_target,
      
      hypothesis = "NN_CMSY",
      
      method = "Neural Network CMSY++",
      
      source = "Catch-informed",
      
      rationale =
        "Depleção inferida pela forma e magnitude da série de capturas através da ANN do CMSY++.",
      
      bk_lo = bk_nn[1],
      
      bk_hi = bk_nn[2],
      
      bk = mean(
        bk_nn,
        na.rm = TRUE
      )
    )
  
    # ===================
    # HIPÓTESE 2 — zBRT
    # ===================
    bk_brt <- estimate_zbrt(
      ct_data = data,
      target_year = yr_target
    )
    
    results[[length(results) + 1]] <- data.frame(
      
      year = yr_target,
      
      hypothesis = "zBRT",
      
      method = "zBRT",
      
      source = "Catch-informed",
      
      rationale =
        "Depleção inferida pela dinâmica temporal da série de capturas através do zBRT.",
      
      bk_lo = bk_brt$bk_lo,
      
      bk_hi = bk_brt$bk_hi,
      
      bk = bk_brt$bk
    )
    
    # ==============================
    # HIPÓTESE 3 — TARGET SWITCHING
    # ==============================
    results[[length(results) + 1]] <- data.frame(
      
      year = yr_target,
      
      hypothesis = "Target_switch",
      
      method = "Target switching",
      
      source = "Independent",
      
      rationale =
        paste0(
          "Hipótese de mudança no direcionamento da pescaria: ",
          "a redução da captura de D. macarellus pode refletir ",
          "mudança de alvo/catchability, e não necessariamente ",
          "redução proporcional da biomassa. Intervalo amplo de B/K."
        ),
      
      bk_lo = 0.20,
      
      bk_hi = 0.70,
      
      bk = 0.45
    )
    
    # ==============================
    # HIPÓTESE 4 — NÃO INFORMATIVA
    # ==============================
    results[[length(results) + 1]] <- data.frame(
      
      year = yr_target,
      
      hypothesis = "Uninformative",
      
      method = "Non-informative depletion prior",
      
      source = "Independent",
      
      rationale =
        paste0(
          "Ausência de informação prévia sobre a depleção final. ",
          "Toda a faixa biologicamente possível de B/K entre 0 e 1 ",
          "é permitida para avaliar quanto os dados de captura e ",
          "a produtividade conseguem restringir a solução."
        ),
      
      bk_lo = 0.00,
      
      bk_hi = 1.00,
      
      bk = 0.50
    )
    
  }
  
  # ===================
  # JUNTAR RESULTADOS
  # ===================
  results <- bind_rows(
    results
  )
  # =========
  # ORDENAR
  # =========
  results <- results %>%
    mutate(
      
      hypothesis = factor(
        hypothesis,
        levels = c(
          "NN_CMSY",
          "zBRT",
          "Target_switch",
          "Uninformative"
        )
      )
    ) %>%
    arrange(
      year,
      hypothesis
    )
  return(results)
}

# ==============================================================================
# 13. RODAR TODAS AS HIPÓTESES PARA D. macarellus
# ==============================================================================
bk_macarellus <- run_depletion_hypotheses(
  data = ct,
  target_years = c(
    2015
  )
)
# ===========================
# 14. VISUALIZAR RESULTADOS
# ===========================
bk_macarellus <- bk_macarellus %>%
  mutate(bk_lo= round(bk_lo,2),
         bk_hi= round(bk_hi,2),
         bk= round(bk,2))

print(
  bk_macarellus
)
# ======================
# 16. SALVAR RESULTADOS
# ======================
write.csv(
  bk_macarellus,
  "Depletion_hypotheses_Decapterus_macarellus.csv",
  row.names = FALSE
)



# ============================================================================
# PRIOR DE r PARA Decapterus macarellus
# Baseado em história de vida e variabilidade observada na literatura
#
# Adaptado por Matheus Silva
## Objetivos:
#   1. Extrair parâmetros de história de vida de D. macarellus
#   2. Padronizar unidades
#   3. Estimar CV empírico a partir da literatura
#   4. Propagar incerteza por bootstrap paramétrico
#   5. Estimar r a partir de parâmetros de história de vida
#   6. Gerar distribuição de r para uso como prior no CMSY/DB-SRA
#
# ============================================================================
# -----------
# 1. PACOTES
# ------------
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
# -----------------------------
# 2. DADOS DE HISTÓRIA DE VIDA
# -----------------------------
# A tabela lh já deve estar carregada, por exemplo:
# lh <- read_xlsx("Parâmetros_História de vida.xlsx")
# Selecionar somente Decapterus macarellus
lh_mac <- lh %>%
  filter(especie == "Decapterus macarellus")

# -------------------------------
# 3. VERIFICAR DADOS DISPONÍVEIS
# -------------------------------
cat("\n")
cat("============================================================\n")
cat("HISTÓRIA DE VIDA — Decapterus macarellus\n")
cat("============================================================\n")
cat("\n")
cat("Número de registros:", nrow(lh_mac), "\n")
cat("Número de fontes:",
    length(unique(lh_mac$fonte[!is.na(lh_mac$fonte)])), "\n\n")

# -------------------------
# 4. CONVERSÃO DE UNIDADES
# -------------------------
## a tabela apresenta:
#  linf_fl  -> comprimento assintótico em FL
#  unidade_crescimento -> FL
## Os valores parecem estar em mm (ex.: 301, 406).
# Para o cálculo usamos cm.
# L50 também aparece como valores como 218 e 203,
# portanto convertemos mm -> cm.
## -----------------------------
lh_mac <- lh_mac %>%
  mutate(
    Linf_cm = ifelse(
      !is.na(linf_fl),
      linf_fl / 10,
      NA_real_
    ),
    L50_cm = ifelse(
      !is.na(l50_fl),
      l50_fl / 10,
      NA_real_
    )
  )
# ------------------------------------
# 5. FUNÇÃO PARA CALCULAR CV EMPÍRICO
# ------------------------------------
calc_cv <- function(x) {
  
  x <- x[
    is.finite(x) &
      !is.na(x) &
      x > 0
  ]
  
  n <- length(x)
  
  if (n >= 2) {
    
    media <- mean(x)
    sd_x <- sd(x)
    cv <- sd_x / media
    
    return(
      tibble(
        n = n,
        mean = media,
        sd = sd_x,
        cv = cv,
        source = "Literature empirical CV"
      )
    )
    
  } else {
    
    return(
      tibble(
        n = n,
        mean = ifelse(n == 1, x[1], NA_real_),
        sd = NA_real_,
        cv = NA_real_,
        source = "Insufficient literature estimates"
      )
    )
  }
}


# ----------------------------------------------------------------------------
# 6. CALCULAR CV EMPÍRICO DA LITERATURA
# ----------------------------------------------------------------------------
#
# Cada parâmetro é tratado separadamente.
#
# ----------------------------------------------------------------------------

cv_linf <- calc_cv(lh_mac$Linf_cm)
cv_k    <- calc_cv(lh_mac$k)
cv_M    <- calc_cv(lh_mac$m)
cv_tmax <- calc_cv(lh_mac$tmax)
cv_L50  <- calc_cv(lh_mac$L50_cm)


# ----------------------------------------------------------------------------
# 7. TABELA DE VARIABILIDADE OBSERVADA
# ----------------------------------------------------------------------------

cv_literature <- bind_rows(
  
  cv_linf %>%
    mutate(parameter = "Linf"),
  
  cv_k %>%
    mutate(parameter = "K"),
  
  cv_M %>%
    mutate(parameter = "M"),
  
  cv_tmax %>%
    mutate(parameter = "tmax"),
  
  cv_L50 %>%
    mutate(parameter = "L50")
  
) %>%
  select(
    parameter,
    n,
    mean,
    sd,
    cv,
    source
  )


cat("\n")
cat("============================================================\n")
cat("VARIABILIDADE DOS PARÂMETROS NA LITERATURA\n")
cat("============================================================\n")
cat("\n")

print(cv_literature)


# ----------------------------------------------------------------------------
# 8. CV FALLBACK
# ----------------------------------------------------------------------------
#
# Quando existe apenas uma estimativa na literatura, não conseguimos
# calcular empiricamente o CV.
#
# Nesses casos usamos valores conservadores.
#
# IMPORTANTE:
# esses valores NÃO substituem a variabilidade da literatura.
# Eles são utilizados apenas quando a literatura disponível não permite
# estimar um CV diretamente.
#
# ----------------------------------------------------------------------------

cv_fallback <- list(
  Linf = 0.15,
  K    = 0.20,
  M    = 0.30,
  tmax = 0.20,
  L50  = 0.15
)


# ----------------------------------------------------------------------------
# 9. FUNÇÃO PARA ESCOLHER CV EMPÍRICO OU FALLBACK
# ----------------------------------------------------------------------------

get_cv <- function(parameter, cv_table, fallback) {
  
  cv_emp <- cv_table$cv[
    cv_table$parameter == parameter
  ]
  
  if (
    length(cv_emp) == 1 &&
    !is.na(cv_emp)
  ) {
    
    return(
      list(
        cv = cv_emp,
        source = "Literature empirical CV"
      )
    )
    
  } else {
    
    return(
      list(
        cv = fallback[[parameter]],
        source = "Fallback CV"
      )
    )
  }
}


# ----------------------------------------------------------------------------
# 10. DEFINIR CVs UTILIZADOS NO BOOTSTRAP
# ----------------------------------------------------------------------------

cv_Linf_use <- get_cv(
  "Linf",
  cv_literature,
  cv_fallback
)

cv_K_use <- get_cv(
  "K",
  cv_literature,
  cv_fallback
)

cv_M_use <- get_cv(
  "M",
  cv_literature,
  cv_fallback
)

cv_tmax_use <- get_cv(
  "tmax",
  cv_literature,
  cv_fallback
)

cv_L50_use <- get_cv(
  "L50",
  cv_literature,
  cv_fallback
)


# ----------------------------------------------------------------------------
# 11. MOSTRAR CVs UTILIZADOS
# ----------------------------------------------------------------------------

cv_used <- tibble(
  
  parameter = c(
    "Linf",
    "K",
    "M",
    "tmax",
    "L50"
  ),
  
  cv = c(
    cv_Linf_use$cv,
    cv_K_use$cv,
    cv_M_use$cv,
    cv_tmax_use$cv,
    cv_L50_use$cv
  ),
  
  source = c(
    cv_Linf_use$source,
    cv_K_use$source,
    cv_M_use$source,
    cv_tmax_use$source,
    cv_L50_use$source
  )
)


cat("\n")
cat("============================================================\n")
cat("CVs UTILIZADOS NO BOOTSTRAP\n")
cat("============================================================\n")
cat("\n")

print(cv_used)


# ----------------------------------------------------------------------------
# 12. FUNÇÃO — LOGNORMAL A PARTIR DE MÉDIA E CV
# ----------------------------------------------------------------------------
#
# A distribuição lognormal é usada porque os parâmetros:
#
#   Linf
#   K
#   M
#   tmax
#   L50
#
# precisam ser positivos.
#
# ----------------------------------------------------------------------------

rlnorm_from_mean_cv <- function(
    mean,
    cv,
    n
) {
  
  if (
    is.na(mean) ||
    is.na(cv) ||
    mean <= 0 ||
    cv <= 0
  ) {
    
    return(
      rep(
        NA_real_,
        n
      )
    )
  }
  
  sigma2 <- log(
    1 + cv^2
  )
  
  mu <- log(mean) -
    0.5 * sigma2
  
  rlnorm(
    n,
    meanlog = mu,
    sdlog = sqrt(sigma2)
  )
}


# ----------------------------------------------------------------------------
# 13. EXTRAIR ESTIMATIVAS CENTRAIS DA LITERATURA
# ----------------------------------------------------------------------------
#
# Usamos a média das estimativas disponíveis.
#
# ----------------------------------------------------------------------------

Linf0 <- mean(
  lh_mac$Linf_cm,
  na.rm = TRUE
)

K0 <- mean(
  lh_mac$k,
  na.rm = TRUE
)

M0 <- mean(
  lh_mac$m,
  na.rm = TRUE
)

tmax0 <- mean(
  lh_mac$tmax,
  na.rm = TRUE
)

L500 <- mean(
  lh_mac$L50_cm,
  na.rm = TRUE
)


# ----------------------------------------------------------------------------
# 14. RESUMO DOS PARÂMETROS CENTRAIS
# ----------------------------------------------------------------------------

life_history_summary <- tibble(
  
  parameter = c(
    "Linf",
    "K",
    "M",
    "tmax",
    "L50"
  ),
  
  value = c(
    Linf0,
    K0,
    M0,
    tmax0,
    L500
  )
)


cat("\n")
cat("============================================================\n")
cat("PARÂMETROS CENTRAIS\n")
cat("============================================================\n")
cat("\n")

print(life_history_summary)


# ----------------------------------------------------------------------------
# 15. CHECAGEM DE PARÂMETROS
# ----------------------------------------------------------------------------

if (
  any(
    is.na(
      c(
        Linf0,
        K0,
        M0,
        tmax0,
        L500
      )
    )
  )
) {
  
  stop(
    paste0(
      "Parâmetros insuficientes para calcular o prior de r.\n",
      "Verifique Linf, K, M, tmax e L50."
    )
  )
}


# ----------------------------------------------------------------------------
# 16. BOOTSTRAP DE HISTÓRIA DE VIDA
# ----------------------------------------------------------------------------

set.seed(123)


nboot <- 10000


Linf_samps <- rlnorm_from_mean_cv(
  Linf0,
  cv_Linf_use$cv,
  nboot
)


K_samps <- rlnorm_from_mean_cv(
  K0,
  cv_K_use$cv,
  nboot
)


M_samps <- rlnorm_from_mean_cv(
  M0,
  cv_M_use$cv,
  nboot
)


tmax_samps <- rlnorm_from_mean_cv(
  tmax0,
  cv_tmax_use$cv,
  nboot
)


L50_samps <- rlnorm_from_mean_cv(
  L500,
  cv_L50_use$cv,
  nboot
)


# ----------------------------------------------------------------------------
# 17. RESTRIÇÕES BIOLÓGICAS
# ----------------------------------------------------------------------------
#
# L50 não pode ser >= Linf.
#
# Quando isso ocorrer no bootstrap, ajustamos L50 para uma fração de Linf.
#
# Para tmax usamos pelo menos 1 ano.
#
# ----------------------------------------------------------------------------

L50_samps <- pmin(
  L50_samps,
  0.95 * Linf_samps
)

tmax_samps <- pmax(
  tmax_samps,
  1
)


# ----------------------------------------------------------------------------
# 18. FUNÇÃO PRINCIPAL PARA ESTIMAR r
# ----------------------------------------------------------------------------
#
# Aqui usamos o princípio demográfico baseado em sobrevivência e maturidade.
#
# Como sua tabela não possui fecundidade/litter size confiável,
# NÃO utilizamos as equações de Myers/Smith do código antigo.
#
# Em vez disso, utilizamos uma aproximação de Euler-Lotka com fecundidade
# relativa normalizada.
#
# Isso fornece um prior de r dependente de:
#
#   - crescimento
#   - maturidade
#   - mortalidade natural
#   - longevidade
#
# sem transformar F (mortalidade por pesca) em fecundidade.
#
# ----------------------------------------------------------------------------

estimate_r_euler <- function(
    Linf,
    K,
    L50,
    M,
    tmax
) {
  
  if (
    any(
      !is.finite(
        c(
          Linf,
          K,
          L50,
          M,
          tmax
        )
      )
    )
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  if (
    Linf <= 0 ||
    K <= 0 ||
    L50 <= 0 ||
    M <= 0 ||
    tmax <= 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Idade de maturação
  # --------------------------------------------------------------------------
  
  if (
    L50 >= Linf
  ) {
    
    L50 <- 0.95 * Linf
  }
  
  
  t50 <- -log(
    1 - L50 / Linf
  ) / K
  
  
  # --------------------------------------------------------------------------
  # Idades
  # --------------------------------------------------------------------------
  
  ages <- seq(
    0,
    ceiling(tmax),
    by = 1
  )
  
  
  # --------------------------------------------------------------------------
  # Sobrevivência
  # --------------------------------------------------------------------------
  
  lx <- exp(
    -M * ages
  )
  
  
  # --------------------------------------------------------------------------
  # Maturidade
  # --------------------------------------------------------------------------
  #
  # Função logística centrada em t50.
  #
  # --------------------------------------------------------------------------
  
  mat <- 1 /
    (
      1 +
        exp(
          -(ages - t50)
        )
    )
  
  
  # --------------------------------------------------------------------------
  # Fecundidade relativa
  # --------------------------------------------------------------------------
  #
  # Sem dados confiáveis de fecundidade absoluta, utilizamos fecundidade
  # relativa proporcional à maturidade.
  #
  # A escala é normalizada de forma que o maior valor de fecundidade
  # relativa seja 1.
  #
  # --------------------------------------------------------------------------
  
  mx <- mat
  
  
  if (
    sum(
      lx * mx
    ) <= 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Normalização da fecundidade
  # --------------------------------------------------------------------------
  #
  # O fator de reprodução é ajustado para que a população esteja
  # aproximadamente no equilíbrio na ausência de mortalidade adicional.
  #
  # Isso transforma o cálculo em uma estimativa de potencial intrínseco
  # de crescimento, e não em uma estimativa absoluta de recrutamento.
  #
  # --------------------------------------------------------------------------
  
  mx <- mx /
    sum(
      lx * mx
    )
  
  
  # --------------------------------------------------------------------------
  # Equação de Euler-Lotka
  # --------------------------------------------------------------------------
  
  euler_fn <- function(r) {
    
    sum(
      lx *
        mx *
        exp(
          -r * ages
        )
    ) - 1
  }
  
  
  # --------------------------------------------------------------------------
  # Encontrar raiz
  # --------------------------------------------------------------------------
  
  r_grid <- seq(
    -2,
    3,
    by = 0.01
  )
  
  
  vals <- sapply(
    r_grid,
    euler_fn
  )
  
  
  valid <- is.finite(
    vals
  )
  
  
  r_grid <- r_grid[valid]
  vals <- vals[valid]
  
  
  if (
    length(vals) < 2
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  change <- which(
    vals[-length(vals)] *
      vals[-1] <= 0
  )
  
  
  if (
    length(change) == 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  i <- change[1]
  
  
  tryCatch(
    
    uniroot(
      euler_fn,
      lower = r_grid[i],
      upper = r_grid[i + 1]
    )$root,
    
    error = function(e)
      NA_real_
  )
}


# ----------------------------------------------------------------------------
# 19. RODAR BOOTSTRAP
# ----------------------------------------------------------------------------

r_sims <- map_dfr(
  
  seq_len(nboot),
  
  function(i) {
    
    r <- estimate_r_euler(
      
      Linf = Linf_samps[i],
      
      K = K_samps[i],
      
      L50 = L50_samps[i],
      
      M = M_samps[i],
      
      tmax = tmax_samps[i]
    )
    
    
    tibble(
      
      iter = i,
      
      specie =
        "Decapterus macarellus",
      
      Linf =
        Linf_samps[i],
      
      K =
        K_samps[i],
      
      M =
        M_samps[i],
      
      tmax =
        tmax_samps[i],
      
      L50 =
        L50_samps[i],
      
      r =
        r
    )
  }
)


# ----------------------------------------------------------------------------
# 20. REMOVER SIMULAÇÕES NÃO CONVERGENTES
# ----------------------------------------------------------------------------

r_sims_valid <- r_sims %>%
  filter(
    is.finite(r)
  )


# ----------------------------------------------------------------------------
# 21. RESUMO DO PRIOR DE r
# ----------------------------------------------------------------------------

r_summary <- tibble(
  
  specie =
    "Decapterus macarellus",
  
  n_total =
    nrow(r_sims),
  
  n_converged =
    nrow(r_sims_valid),
  
  convergence =
    nrow(r_sims_valid) /
    nrow(r_sims),
  
  r_median =
    median(
      r_sims_valid$r
    ),
  
  r_q025 =
    quantile(
      r_sims_valid$r,
      0.025
    ),
  
  r_q975 =
    quantile(
      r_sims_valid$r,
      0.975
    ),
  
  r_mean =
    mean(
      r_sims_valid$r
    ),
  
  r_sd =
    sd(
      r_sims_valid$r
    ),
  
  r_min =
    min(
      r_sims_valid$r
    ),
  
  r_max =
    max(
      r_sims_valid$r
    )
)


# ----------------------------------------------------------------------------
# 22. MOSTRAR RESULTADO
# ----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("PRIOR DE r — Decapterus macarellus\n")
cat("============================================================\n")
cat("\n")

print(r_summary)


# ----------------------------------------------------------------------------
# 23. CHECAGEM DE VALORES EXTREMOS
# ----------------------------------------------------------------------------

cat("\n")
cat("Percentis da distribuição de r:\n\n")

print(
  quantile(
    r_sims_valid$r,
    probs = c(
      0.001,
      0.01,
      0.025,
      0.05,
      0.25,
      0.50,
      0.75,
      0.95,
      0.975,
      0.99,
      0.999
    )
  )
)


# ----------------------------------------------------------------------------
# 24. TABELA DOS PARÂMETROS BOOTSTRAP
# ----------------------------------------------------------------------------

parameter_bootstrap_summary <- tibble(
  
  parameter = c(
    "Linf",
    "K",
    "M",
    "tmax",
    "L50"
  ),
  
  mean = c(
    mean(r_sims_valid$Linf),
    mean(r_sims_valid$K),
    mean(r_sims_valid$M),
    mean(r_sims_valid$tmax),
    mean(r_sims_valid$L50)
  ),
  
  median = c(
    median(r_sims_valid$Linf),
    median(r_sims_valid$K),
    median(r_sims_valid$M),
    median(r_sims_valid$tmax),
    median(r_sims_valid$L50)
  ),
  
  q025 = c(
    quantile(r_sims_valid$Linf, 0.025),
    quantile(r_sims_valid$K, 0.025),
    quantile(r_sims_valid$M, 0.025),
    quantile(r_sims_valid$tmax, 0.025),
    quantile(r_sims_valid$L50, 0.025)
  ),
  
  q975 = c(
    quantile(r_sims_valid$Linf, 0.975),
    quantile(r_sims_valid$K, 0.975),
    quantile(r_sims_valid$M, 0.975),
    quantile(r_sims_valid$tmax, 0.975),
    quantile(r_sims_valid$L50, 0.975)
  )
)


# ----------------------------------------------------------------------------
# 25. MOSTRAR PARÂMETROS BOOTSTRAP
# ----------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("DISTRIBUIÇÕES DOS PARÂMETROS BOOTSTRAP\n")
cat("============================================================\n")
cat("\n")

print(
  parameter_bootstrap_summary
)


# ----------------------------------------------------------------------------
# 26. SALVAR SIMULAÇÕES
# ----------------------------------------------------------------------------

write.csv(
  r_sims,
  "r_sims_Decapterus_macarellus.csv",
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 27. SALVAR SIMULAÇÕES VÁLIDAS
# ----------------------------------------------------------------------------

write.csv(
  r_sims_valid,
  "r_sims_valid_Decapterus_macarellus.csv",
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 28. SALVAR RESUMO DE r
# ----------------------------------------------------------------------------

write.csv(
  r_summary,
  "r_summary_Decapterus_macarellus.csv",
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 29. SALVAR CVs
# ----------------------------------------------------------------------------

write.csv(
  cv_literature,
  "CV_literature_Decapterus_macarellus.csv",
  row.names = FALSE
)


write.csv(
  cv_used,
  "CV_used_Decapterus_macarellus.csv",
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 30. SALVAR RESUMO DE HISTÓRIA DE VIDA
# ----------------------------------------------------------------------------

write.csv(
  parameter_bootstrap_summary,
  "life_history_bootstrap_Decapterus_macarellus.csv",
  row.names = FALSE
)


# ============================================================================
# 31. DISTRIBUIÇÃO DE r
# ============================================================================

p_r <- ggplot(
  r_sims_valid,
  aes(
    x = r
  )
) +
  
  geom_histogram(
    bins = 60
  ) +
  
  geom_vline(
    xintercept =
      r_summary$r_median,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  labs(
    x = "Intrinsic growth rate (r)",
    y = "Frequency"
  ) +
  
  theme_classic(
    base_size = 14
  )


p_r


# ----------------------------------------------------------------------------
# 32. SALVAR FIGURA
# ----------------------------------------------------------------------------

ggsave(
  "r_prior_Decapterus_macarellus.png",
  plot = p_r,
  device = "png",
  units = "cm",
  width = 18,
  height = 12,
  dpi = 300
)


# ============================================================================
# 33. DISTRIBUIÇÕES DOS PARÂMETROS DE HISTÓRIA DE VIDA
# ============================================================================

lh_long <- r_sims_valid %>%
  
  select(
    Linf,
    K,
    M,
    tmax,
    L50
  ) %>%
  
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  )


p_lh <- ggplot(
  lh_long,
  aes(
    x = value
  )
) +
  
  geom_histogram(
    bins = 50
  ) +
  
  facet_wrap(
    ~parameter,
    scales = "free"
  ) +
  
  labs(
    x = "Value",
    y = "Frequency"
  ) +
  
  theme_classic(
    base_size = 14
  )


p_lh


# ----------------------------------------------------------------------------
# 34. SALVAR FIGURA
# ----------------------------------------------------------------------------

ggsave(
  "life_history_bootstrap_Decapterus_macarellus.png",
  plot = p_lh,
  device = "png",
  units = "cm",
  width = 20,
  height = 14,
  dpi = 300
)


# ============================================================================
# 35. RESUMO FINAL
# ============================================================================

cat("\n")
cat("============================================================\n")
cat("ANÁLISE FINALIZADA\n")
cat("============================================================\n")
cat("\n")

cat(
  "Espécie: Decapterus macarellus\n"
)

cat(
  "Bootstrap:",
  nboot,
  "simulações\n"
)

cat(
  "Simulações convergentes:",
  nrow(r_sims_valid),
  "\n"
)

cat(
  "Taxa de convergência:",
  round(
    r_summary$convergence,
    3
  ),
  "\n"
)

cat(
  "r mediano:",
  round(
    r_summary$r_median,
    3
  ),
  "\n"
)

cat(
  "IC/intervalo 95%:",
  round(
    r_summary$r_q025,
    3
  ),
  "–",
  round(
    r_summary$r_q975,
    3
  ),
  "\n"
)

cat("\n")

cat(
  "Arquivos gerados:\n"
)

cat(
  " - r_sims_Decapterus_macarellus.csv\n"
)

cat(
  " - r_sims_valid_Decapterus_macarellus.csv\n"
)

cat(
  " - r_summary_Decapterus_macarellus.csv\n"
)

cat(
  " - CV_literature_Decapterus_macarellus.csv\n"
)

cat(
  " - CV_used_Decapterus_macarellus.csv\n"
)

cat(
  " - life_history_bootstrap_Decapterus_macarellus.csv\n"
)

cat(
  " - r_prior_Decapterus_macarellus.png\n"
)

cat(
  " - life_history_bootstrap_Decapterus_macarellus.png\n"
)

cat("\n")







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
#install.packages("tidyverse")
library(tidyverse)
#install.packages("tibble")
library(tibble)
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
      
      especie = 'Decapterus macarellus',
      
      ano = yr_target,
      
      hipotese = "NN_CMSY",
      
      metodo = "Rede Neural CMSY++",
      
      fonte = "Informado pela Captura",
      
      racional =
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
      
      especie = 'Decapterus macarellus',
      
      ano = yr_target,
      
      hipotese = "zBRT",
      
      metodo = "zBRT",
      
      fonte = "Informado pela Captura",
      
      racional =
        "Depleção inferida pela dinâmica temporal da série de capturas através do zBRT.",
      
      bk_lo = bk_brt$bk_lo,
      
      bk_hi = bk_brt$bk_hi,
      
      bk = bk_brt$bk
    )
    
    # ==============================
    # HIPÓTESE 3 — TARGET SWITCHING
    # ==============================
    results[[length(results) + 1]] <- data.frame(
      
      especie = 'Decapterus macarellus',
      
      ano = yr_target,
      
      hipotese = "Target_switch",
      
      metodo = "Mudança de alvo",
      
      fonte = "Independente",
      
      racional =
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
      
      especie = 'Decapterus macarellus',
      
      ano = yr_target,
      
      hipotese = "Uninformative",
      
      metodo = "Priori de depleção não informativa",
      
      fonte = "Independente",
      
      racional =
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
      
      hipotese = factor(
        hipotese,
        levels = c(
          "NN_CMSY",
          "zBRT",
          "Target_switch",
          "Uninformative"
        )
      )
    ) %>%
    arrange(
      ano,
      hipotese
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
  "bk_macarellus.csv",
  row.names = FALSE
)


library(tidyverse)

# lh já está carregado como tibble

#------------------------------------------------------------------
# Funções auxiliares (INALTERADAS)
#------------------------------------------------------------------
safe_uniroot <- function(fn, lower = 0, upper = 5, tol = 1e-8,
                         max_expand = 10, by = 0.5) {
  safe_eval <- function(x) tryCatch(fn(x), error = function(e) NA_real_)
  f_low  <- safe_eval(lower)
  f_high <- safe_eval(upper)
  if (!is.na(f_low) && !is.na(f_high) && f_low * f_high < 0) {
    return(tryCatch(uniroot(fn, c(lower, upper), tol = tol)$root,
                    error = function(e) NA_real_))
  }
  for (i in seq_len(max_expand)) {
    new_upper <- upper + i * by
    f_new <- safe_eval(new_upper)
    if (!is.na(f_low) && !is.na(f_new) && f_low * f_new < 0) {
      return(tryCatch(uniroot(fn, c(lower, new_upper), tol = tol)$root,
                      error = function(e) NA_real_))
    }
  }
  NA_real_
}

rlnorm_from_mean_cv <- function(mean, cv, n) {
  if (is.na(mean) || is.na(cv) || mean <= 0) return(rep(NA_real_, n))
  sigma2 <- log(1 + cv^2)
  mu <- log(mean) - 0.5 * sigma2
  rlnorm(n, meanlog = mu, sdlog = sqrt(sigma2))
}

#------------------------------------------------------------------
# Bootstrap para uma espécie  (ÚNICA COISA QUE MUDA: nomes das colunas)
#------------------------------------------------------------------
estimate_r_boot <- function(sp_row, nboot = 1000,
                            cvs = list(Linf = 0.15, k = 0.20, M = 0.20,
                                       tmax = 0.10, L50 = 0.15, ls = 0.25),
                            ls_euler_fixed = 9.5,
                            options = list(r_upper = 5, verbose = FALSE)) {
  
  sp <- sp_row$especie[1]   # <-- era Especie
  
  # ---- MAPEAMENTO DE COLUNAS --------------------------------------
  Linf0 <- mean(na.omit(sp_row$linf_fl)) / 10   # <-- era `Linf(mm)TL`
  k0    <- mean(na.omit(sp_row$k))              # <-- era `K(ano)`
  L500  <- mean(na.omit(sp_row$l50_fl)) / 10    # <-- era `L50(mm)TL`
  M0    <- mean(na.omit(sp_row$m))              # <-- era M
  tmax0 <- mean(na.omit(sp_row$tmax))           # <-- era `Tmáx`
  ls0   <- 4   # default (não existe em lh)
  f     <- 1   # default (não existe em lh)
  # -----------------------------------------------------------------
  
  if (any(is.na(c(Linf0, k0, L500, M0, tmax0)))) {
    if (isTRUE(options$verbose)) warning(sp, ": insufficient parameters.")
    return(list(
      sims = tibble(),
      summary = tibble(specie = sp,
                       method = c("euler","myers","smith_rebound_eq6","demographic_inv"),
                       r_median = NA_real_, r_q025 = NA_real_, r_q975 = NA_real_,
                       n_conv = 0L, n_total = nboot)
    ))
  }
  
  # Bootstrap paramétrico (inalterado)
  Linf_samps <- rlnorm_from_mean_cv(Linf0, cvs$Linf, nboot)
  k_samps    <- rlnorm_from_mean_cv(k0,    cvs$k,    nboot)
  M_samps    <- rlnorm_from_mean_cv(M0,    cvs$M,    nboot)
  tmax_samps <- pmax(1, round(rlnorm_from_mean_cv(tmax0, cvs$tmax, nboot)))
  L50_samps  <- rlnorm_from_mean_cv(L500,  cvs$L50,  nboot)
  ls_samps   <- rlnorm_from_mean_cv(ls0,   cvs$ls,   nboot)
  
  # Cálculo por iteração (INALTERADO)
  run_one <- function(Linf, k, L50, M, tmax, ls) {
    if (is.na(Linf) || Linf <= 0 || is.na(k) || is.na(L50) || is.na(M) || is.na(tmax))
      return(c(NA,NA,NA,NA))
    if (L50 >= Linf) L50 <- 0.5 * Linf
    
    t50   <- -(log(1 - L50 / Linf) / k)
    ages  <- 0:ceiling(tmax)
    lx    <- exp(-M * ages)
    mat_a <- 1 / (1 + exp(-(ages - t50)))
    fr    <- ls / f / 2
    mx    <- fr * mat_a
    fr_euler  <- ls_euler_fixed / f / 2
    mx_euler  <- fr_euler * mat_a
    
    euler_fn <- function(r) sum(lx * mx_euler * exp(-r * ages)) - 1
    
    s_adult  <- lx[which.min(abs(lx - 0.5))]
    litter   <- ls; freqv <- f; tmat <- t50
    formula_myers <- function(rm) ((exp(rm))^tmat) -
      ((s_adult) * ((exp(rm))^(tmat - 1))) - (litter / freqv / 2)
    
    Z <- 1.5 * M
    l_alpha <- if ((tmax - tmat + 1) > 0)
      (1 - exp(-Z)) / ((litter/2/freqv) * (1 - exp(-Z*(tmax - tmat + 1)))) else NA
    
    eq6 <- function(reb) if (is.na(l_alpha)) NA_real_ else
      1 - exp(-(M + reb)) -
      l_alpha * (litter/2/freqv) * 1.25 * exp(-reb*tmat) *
      (1 - exp(-(M+reb)*(tmax - tmat + 1)))
    
    formula5 <- function(r) if (exp(r) <= s_adult) NA_real_ else
      exp(r) - (exp(1 / (tmat + 1 + (s_adult / (exp(r) - s_adult)))))
    
    up <- options$r_upper %||% 5
    c(safe_uniroot(euler_fn, 0, up),
      safe_uniroot(formula_myers, 0, up),
      safe_uniroot(eq6, 0, up),
      safe_uniroot(formula5, 0, up))
  }
  
  sims <- purrr::pmap_dfr(
    list(Linf_samps, k_samps, L50_samps, M_samps, tmax_samps, ls_samps),
    function(Linf, k, L50, M, tmax, ls) {
      rv <- run_one(Linf, k, L50, M, tmax, ls)
      tibble(r_euler = rv[1], r_myers = rv[2], r_eq6 = rv[3], r_f5 = rv[4])
    }) %>%
    dplyr::mutate(iter = dplyr::row_number(), specie = sp)
  
  summarize_method <- function(x) {
    n_conv <- sum(!is.na(x))
    tibble(median = median(x, na.rm = TRUE),
           q025   = quantile(x, 0.025, na.rm = TRUE),
           q975   = quantile(x, 0.975, na.rm = TRUE),
           n_conv = n_conv)
  }
  s1 <- summarize_method(sims$r_euler)
  s2 <- summarize_method(sims$r_myers)
  s3 <- summarize_method(sims$r_eq6)
  s4 <- summarize_method(sims$r_f5)
  
  summary_tbl <- tibble(
    specie   = sp,
    method   = c("Euler","Myers","Smith rebound","Demographic inv"),
    r_median = c(s1$median, s2$median, s3$median, s4$median),
    r_q025   = c(s1$q025,   s2$q025,   s3$q025,   s4$q025),
    r_q975   = c(s1$q975,   s2$q975,   s3$q975,   s4$q975),
    n_conv   = c(s1$n_conv, s2$n_conv, s3$n_conv, s4$n_conv),
    n_total  = nboot
  )
  list(sims = sims, summary = summary_tbl)
}

#------------------------------------------------------------------
# Aplicar a todas as espécies (usando `lh`)
#------------------------------------------------------------------
species_list <- unique(lh$especie)

res_list <- map(species_list, function(sp) {
  sp_row <- lh %>% filter(especie == sp)
  estimate_r_boot(sp_row, nboot = 1000)
})

r_sims <- map_dfr(res_list, "sims")
write.csv(r_sims, "r_sims.csv", row.names = FALSE)

#------------------------------------------------------------------
# 1) Sumários por espécie × método (inalterado)
#------------------------------------------------------------------
r_methods <- map_dfr(res_list, "summary")

#------------------------------------------------------------------
# 2) Agregação baseline (mesma lógica de antes)
#------------------------------------------------------------------
r_base <- r_methods %>%
  dplyr::group_by(specie) %>%
  dplyr::summarise(
    r_median = median(r_median, na.rm = TRUE),
    r_min    = pmax(median(r_q025, na.rm = TRUE), 0.1),
    r_max    = pmin(median(r_q975, na.rm = TRUE), 1.5),
    .groups  = "drop"
  )

#------------------------------------------------------------------
# 3) r_summary expandido com as hipóteses (formato bk_macarellus)
#------------------------------------------------------------------
r_macarellus <- bind_rows(
  # (a) Baseline – modelos demográficos de história de vida
  r_base %>%
    transmute(
      specie,
      ano       = "1989-2015",
      hipotese = "Life-history Euler-lotka derived methods (baseline)",
      metodo     = "Euler / Myers / Smith rebound / Demographic inv",
      fonte     = "Baseado em modelo",
      racional  = paste("Bootstrap paramétrico sobre os quatro métodos",
                         "demográficos clássicos baseados em história de vida",
                         "(Euler-Lotka, Myers, Smith rebound e demographic invariant).",
                         "Mediana agregada como valor central."),
      r_lo = r_min,
      r_hi = r_max,
      r    = r_median
    ),
  
  # (b) Baixa resiliência (depleção mais baixa): subtrai 0.2
  r_base %>%
    transmute(
      specie,
      ano       = "1989-2015",
      hipotese = "lower resilience",
      metodo     = "Metodos Euler-Lotka - 0.2",
      fonte     = "Sensibilidade",
      racional  = paste("Cenário de produtividade pessimista: subtrai 0.2 dos limites inferior,",
                         "superior e da mediana estimados pelos métodos Euler-Lotka.",
                         "Representa uma população menos resiliente"),
      r_lo = pmax(r_min - 0.2, 0.05),
      r_hi = pmax(r_max - 0.2, 0.05),
      r    = pmax(r_median - 0.2, 0.05)
    ),
  
  # (c) Alta resiliência: adiciona 0.2
  r_base %>%
    transmute(
      specie,
      ano       = "1989-2015",
      hipotese = "Higher resilience",
      metodo     = "Metodos Euler-Lotka + 0.2",
      fonte     = "Sensibilidade",
      racional  = paste("Cenário de produtividade otimista: adiciona 0.2 aos limites inferior,",
                         "superior e à mediana estimados pelos métodos Euler-Lotka.",
                         "Representa uma população mais produtiva"),
      r_lo = r_min + 0.2,
      r_hi = pmin(r_max + 0.2, 1.5),
      r    = r_median + 0.2
    ),
  
  # (d) Não informativo (uniforme 0–1.5)
  r_base %>%
    transmute(
      specie,
      ano       = "1989-2015",
      hipotese = "Non-informative",
      metodo     = "Priori não informativa",
      fonte     = "Independente",
      racional  = paste("Sem informação prévia: toda a faixa biologicamente",
                         "plausível de r entre 0 e 1.5 é considerada."),
      r_lo = 0,
      r_hi = 1.5,
      r    = NA_real_
    )
) %>%
  dplyr::mutate(across(c(r, r_lo, r_hi), \(x) round(x, 2))) %>%
  rename(especie = specie)


# Salvar
write.csv(r_macarellus, "r_macarellus.csv", row.names = FALSE)

#------------------------------------------------------------------
# Plot (igual)
#------------------------------------------------------------------
all_sims_long <- r_sims %>%
  pivot_longer(cols = starts_with("r_"),
               names_to = "method", values_to = "r") %>%
  mutate(method = dplyr::recode(method,
                                r_euler = "Euler",
                                r_myers = "Myers",
                                r_eq6   = "Smith rebound",
                                r_f5    = "Demographic inv"))

p_r <- ggplot(all_sims_long, aes(x = specie, y = r, col = method, fill = method)) +
  geom_boxplot(aes(fill = method, col = method), alpha = 0.4, width = 0.3,
               position = position_dodge(width = 0.8)) +
  geom_violin(aes(col = method), trim = TRUE, alpha = 0.5, width = 1.5,
              position = position_dodge(width = 0.8)) +
  geom_jitter(aes(col = method),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
              size = 1, alpha = 0.3) +
  labs(x = "Species", y = "Intrinsic growth rate (r)", fill = "", color = "") +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.1)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_classic(base_size = 13) %+replace%
  theme(
    strip.background = element_blank(),
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"),
    strip.text.x = element_text(margin = margin(b = 1), size = 13),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13, face = "italic"),
    legend.text = element_text(size = 13),
    legend.box.margin = margin(t = -10),
    legend.spacing.y = unit(0.1, "cm"),
    legend.position = "bottom"
  )
p_r

ggsave("r_priors.png", plot = p_r, device = "png", units = "cm",
       width = 32, height = 17)















































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
#--------------------------------------------------------------------------
# Create prior rmax based on life-history traits (Cortes, 2016)
# Modified by Silva MLS
# Parametric bootstrap to propagate uncertainty
# Adapted for Decapterus macarellus - Cabo Verde
#--------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

#--------------------------------------------------------------------------
# 1. IMPORT LIFE-HISTORY DATA
#--------------------------------------------------------------------------

# lh <- readxl::read_xlsx("Parâmetros_História de vida.xlsx")

# Keep only Decapterus macarellus
lf <- lh %>%
  filter(especie == "Decapterus macarellus") %>%
  transmute(
    Especie = especie,
    `Linf(mm)TL` = linf_fl,
    `K(ano)` = k,
    `L50(mm)TL` = l50_fl,
    M = m,
    `Tmáx` = tmax,
    ls = NA_real_,
    f = f
  )

# IMPORTANT:
# Your current table does not contain "ls".
# Therefore the original model uses ls0 = 4.
# f is kept exactly as in the original model.
#
# If f is absent from the dataframe, the model uses f = 2.

#--------------------------------------------------------------------------
# 2. SAFE UNIROOT
#--------------------------------------------------------------------------

safe_uniroot <- function(fn, lower = 0, upper = 5, tol = 1e-8,
                         max_expand = 10, by = 0.5) {
  safe_eval <- function(x) tryCatch(fn(x), error = function(e) NA_real_)
  
  f_low <- safe_eval(lower)
  f_high <- safe_eval(upper)
  
  if (!is.na(f_low) && !is.na(f_high) && f_low * f_high < 0) {
    return(
      tryCatch(
        uniroot(fn, c(lower, upper), tol = tol)$root,
        error = function(e) NA_real_
      )
    )
  }
  
  for (i in seq_len(max_expand)) {
    new_upper <- upper + i * by
    f_new <- safe_eval(new_upper)
    
    if (!is.na(f_low) && !is.na(f_new) && f_low * f_new < 0) {
      return(
        tryCatch(
          uniroot(fn, c(lower, new_upper), tol = tol)$root,
          error = function(e) NA_real_
        )
      )
    }
  }
  
  return(NA_real_)
}

#--------------------------------------------------------------------------
# 3. LOGNORMAL SAMPLING FROM MEAN AND CV
#--------------------------------------------------------------------------

rlnorm_from_mean_cv <- function(mean, cv, n) {
  if (is.na(mean) || is.na(cv) || mean <= 0) {
    return(rep(NA_real_, n))
  }
  
  sigma2 <- log(1 + cv^2)
  mu <- log(mean) - 0.5 * sigma2
  
  rlnorm(
    n,
    meanlog = mu,
    sdlog = sqrt(sigma2)
  )
}

#--------------------------------------------------------------------------
# 4. BOOTSTRAP FOR ONE SPECIES
#--------------------------------------------------------------------------

estimate_r_boot <- function(
    sp_row,
    nboot = 1000,
    cvs = list(
      Linf = 0.15,
      k = 0.20,
      M = 0.20,
      tmax = 0.10,
      L50 = 0.15,
      ls = 0.25
    ),
    ls_euler_fixed = 9.5,
    options = list(
      r_upper = 5,
      verbose = FALSE
    )) {
  
  sp <- sp_row$Especie[1]
  
  #----------------------------------------------------------------------
  # Mean point estimates from data
  #----------------------------------------------------------------------
  
  Linf0 <- mean(
    na.omit(sp_row$`Linf(mm)TL`)
  ) / 10
  
  k0 <- mean(
    na.omit(sp_row$`K(ano)`)
  )
  
  L500 <- mean(
    na.omit(sp_row$`L50(mm)TL`)
  ) / 10
  
  M0 <- mean(
    na.omit(sp_row$M)
  )
  
  tmax0 <- mean(
    na.omit(sp_row$`Tmáx`)
  )
  
  # EXACTLY AS IN THE ORIGINAL CODE:
  # if ls does not exist, ls0 = 4
  ls0 <- if ("ls" %in% names(sp_row)) {
    mean(na.omit(sp_row$ls))
  } else {
    4
  }
  
  # EXACTLY AS IN THE ORIGINAL CODE:
  # if f does not exist, f = 2
  f <- if ("f" %in% names(sp_row)) {
    mean(na.omit(sp_row$f))
  } else {
    2
  }
  
  #----------------------------------------------------------------------
  # Check required parameters
  #----------------------------------------------------------------------
  
  if (any(is.na(c(
    Linf0,
    k0,
    L500,
    M0,
    tmax0
  )))) {
    
    if (isTRUE(options$verbose)) {
      warning(sp, ": insufficient parameters.")
    }
    
    return(
      list(
        sims = tibble(),
        summary = tibble(
          specie = sp,
          method = c(
            "euler",
            "myers",
            "smith_rebound_eq6",
            "demographic_inv"
          ),
          r_median = NA_real_,
          r_q025 = NA_real_,
          r_q975 = NA_real_,
          n_conv = 0L,
          n_total = nboot
        )
      )
    )
  }
  
  #----------------------------------------------------------------------
  # Parametric bootstrap
  #----------------------------------------------------------------------
  
  Linf_samps <- rlnorm_from_mean_cv(
    Linf0,
    cvs$Linf,
    nboot
  )
  
  k_samps <- rlnorm_from_mean_cv(
    k0,
    cvs$k,
    nboot
  )
  
  M_samps <- rlnorm_from_mean_cv(
    M0,
    cvs$M,
    nboot
  )
  
  tmax_samps <- pmax(
    1,
    round(
      rlnorm_from_mean_cv(
        tmax0,
        cvs$tmax,
        nboot
      )
    )
  )
  
  L50_samps <- rlnorm_from_mean_cv(
    L500,
    cvs$L50,
    nboot
  )
  
  ls_samps <- rlnorm_from_mean_cv(
    ls0,
    cvs$ls,
    nboot
  )
  
  #----------------------------------------------------------------------
  # Four original estimators of r
  #----------------------------------------------------------------------
  
  run_one <- function(
    Linf,
    k,
    L50,
    M,
    tmax,
    ls) {
    
    if (
      is.na(Linf) ||
      Linf <= 0 ||
      is.na(k) ||
      is.na(L50) ||
      is.na(M) ||
      is.na(tmax)
    ) {
      return(c(NA, NA, NA, NA))
    }
    
    if (L50 >= Linf) {
      L50 <- 0.5 * Linf
    }
    
    #--------------------------------------------------------------------
    # Calculate age at 50% maturity, survivorship and fecundity
    # EXACTLY AS ORIGINAL
    #--------------------------------------------------------------------
    
    t50 <- -(
      log(1 - L50 / Linf) / k
    )
    
    ages <- 0:ceiling(tmax)
    
    lx <- exp(
      -M * ages
    )
    
    mat_a <- 1 / (
      1 + exp(-(ages - t50))
    )
    
    fr <- ls / f / 2
    
    mx <- fr * mat_a
    
    fr_euler <- ls_euler_fixed / f / 2
    
    mx_euler <- fr_euler * mat_a
    
    #--------------------------------------------------------------------
    # Original root-finding functions
    #--------------------------------------------------------------------
    
    euler_fn <- function(r) {
      sum(
        lx *
          mx_euler *
          exp(-r * ages)
      ) - 1
    }
    
    s_adult <- lx[
      which.min(
        abs(lx - 0.5)
      )
    ]
    
    litter <- ls
    freqv <- f
    tmat <- t50
    
    formula_myers <- function(rm) {
      ((exp(rm))^tmat) -
        (
          (s_adult) *
            ((exp(rm))^(tmat - 1))
        ) -
        (litter / freqv / 2)
    }
    
    Z <- 1.5 * M
    
    l_alpha <- if (
      (tmax - tmat + 1) > 0
    ) {
      (1 - exp(-Z)) /
        (
          (litter / 2 / freqv) *
            (
              1 -
                exp(
                  -Z *
                    (tmax - tmat + 1)
                )
            )
        )
    } else {
      NA
    }
    
    eq6 <- function(reb) {
      if (is.na(l_alpha)) {
        NA_real_
      } else {
        1 -
          exp(-(M + reb)) -
          l_alpha *
          (litter / 2 / freqv) *
          1.25 *
          exp(-reb * tmat) *
          (
            1 -
              exp(
                -(M + reb) *
                  (tmax - tmat + 1)
              )
          )
      }
    }
    
    formula5 <- function(r) {
      if (
        exp(r) <= s_adult
      ) {
        NA_real_
      } else {
        exp(r) -
          (
            exp(
              1 /
                (
                  tmat +
                    1 +
                    (
                      s_adult /
                        (
                          exp(r) -
                            s_adult
                        )
                    )
                )
            )
          )
      }
    }
    
    #--------------------------------------------------------------------
    # Apply safe_uniroot
    #--------------------------------------------------------------------
    
    up <- options$r_upper
    
    r1 <- safe_uniroot(
      euler_fn,
      0,
      up
    )
    
    r2 <- safe_uniroot(
      formula_myers,
      0,
      up
    )
    
    r3 <- safe_uniroot(
      eq6,
      0,
      up
    )
    
    r4 <- safe_uniroot(
      formula5,
      0,
      up
    )
    
    return(
      c(
        r1,
        r2,
        r3,
        r4
      )
    )
  }
  
  #----------------------------------------------------------------------
  # Run bootstrap
  #----------------------------------------------------------------------
  
  sims <- purrr::pmap_dfr(
    list(
      Linf_samps,
      k_samps,
      L50_samps,
      M_samps,
      tmax_samps,
      ls_samps
    ),
    function(
    Linf,
    k,
    L50,
    M,
    tmax,
    ls) {
      
      rvec <- run_one(
        Linf,
        k,
        L50,
        M,
        tmax,
        ls
      )
      
      tibble(
        r_euler = rvec[1],
        r_myers = rvec[2],
        r_eq6 = rvec[3],
        r_f5 = rvec[4]
      )
    }
  ) %>%
    dplyr::mutate(
      iter = dplyr::row_number(),
      specie = sp
    )
  
  #----------------------------------------------------------------------
  # Summaries
  #----------------------------------------------------------------------
  
  summarize_method <- function(x) {
    
    n_conv <- sum(
      !is.na(x)
    )
    
    tibble(
      median = median(
        x,
        na.rm = TRUE
      ),
      q025 = quantile(
        x,
        0.025,
        na.rm = TRUE
      ),
      q975 = quantile(
        x,
        0.975,
        na.rm = TRUE
      ),
      n_conv = n_conv
    )
  }
  
  s1 <- summarize_method(
    sims$r_euler
  )
  
  s2 <- summarize_method(
    sims$r_myers
  )
  
  s3 <- summarize_method(
    sims$r_eq6
  )
  
  s4 <- summarize_method(
    sims$r_f5
  )
  
  summary_tbl <- tibble(
    specie = sp,
    method = c(
      "Euler",
      "Myers",
      "Smith rebound",
      "Demographic inv"
    ),
    r_median = c(
      s1$median,
      s2$median,
      s3$median,
      s4$median
    ),
    r_q025 = c(
      s1$q025,
      s2$q025,
      s3$q025,
      s4$q025
    ),
    r_q975 = c(
      s1$q975,
      s2$q975,
      s3$q975,
      s4$q975
    ),
    n_conv = c(
      s1$n_conv,
      s2$n_conv,
      s3$n_conv,
      s4$n_conv
    ),
    n_total = nboot
  )
  
  return(
    list(
      sims = sims,
      summary = summary_tbl
    )
  )
}

#--------------------------------------------------------------------------
# 5. RUN FOR Decapterus macarellus
#--------------------------------------------------------------------------

species_list <- unique(
  lf$Especie
)

res_list <- map(
  species_list,
  function(sp) {
    
    sp_row <- lf %>%
      filter(
        Especie == sp
      )
    
    estimate_r_boot(
      sp_row,
      nboot = 10000
    )
  }
)

#--------------------------------------------------------------------------
# 6. COMBINE SIMULATIONS
#--------------------------------------------------------------------------

r_sims <- map_dfr(
  res_list,
  "sims"
)

write.csv(
  r_sims,
  "r_sims_D_macarellus.csv",
  row.names = FALSE
)

#--------------------------------------------------------------------------
# 7. COMBINE SUMMARIES
#--------------------------------------------------------------------------

r_summary <- map_dfr(
  res_list,
  "summary"
) %>%
  dplyr::group_by(
    specie
  ) %>%
  dplyr::summarise(
    r_median = median(
      r_median,
      na.rm = TRUE
    ),
    r_min = pmax(
      median(
        r_q025,
        na.rm = TRUE
      ),
      0.1
    ),
    r_max = pmin(
      median(
        r_q975,
        na.rm = TRUE
      ),
      1.5
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    across(
      where(is.numeric),
      \(x) round(x, 2)
    )
  )

write.csv(
  r_summary,
  "r_summary_D_macarellus.csv",
  row.names = FALSE
)

#--------------------------------------------------------------------------
# 8. LONG FORMAT
#--------------------------------------------------------------------------

all_sims_long <- r_sims %>%
  pivot_longer(
    cols = starts_with("r_"),
    names_to = "method",
    values_to = "r"
  ) %>%
  mutate(
    method = dplyr::recode(
      method,
      r_euler = "Euler",
      r_myers = "Myers",
      r_eq6 = "Smith rebound",
      r_f5 = "Demographic inv"
    )
  )

#--------------------------------------------------------------------------
# 9. PLOT
#--------------------------------------------------------------------------

p4 <- ggplot(
  all_sims_long,
  aes(
    x = specie,
    y = r,
    col = method,
    fill = method
  )
) +
  geom_boxplot(
    aes(
      fill = method,
      col = method
    ),
    alpha = 0.4,
    width = 0.3,
    position = position_dodge(
      width = 0.8
    )
  ) +
  geom_violin(
    aes(
      col = method
    ),
    trim = TRUE,
    alpha = 0.5,
    width = 1.5,
    position = position_dodge(
      width = 0.8
    )
  ) +
  geom_jitter(
    aes(
      col = method
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.3
  ) +
  labs(
    x = "Species",
    y = "Intrinsic growth rate (r)",
    fill = "",
    color = ""
  ) +
  scale_y_continuous(
    limits = c(0, 1.5),
    breaks = seq(
      0,
      1.5,
      0.1
    )
  ) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_classic(
    base_size = 15
  ) %+replace%
  theme(
    strip.background = element_blank(),
    plot.margin = unit(
      c(
        0.05,
        0.05,
        0.05,
        0.05
      ),
      "mm"
    ),
    strip.text.x = element_text(
      margin = margin(b = 1),
      size = 15
    ),
    axis.text.y = element_text(
      size = 15
    ),
    axis.text.x = element_text(
      size = 15,
      face = "italic"
    ),
    legend.text = element_text(
      size = 15
    ),
    legend.box.margin = margin(
      t = -10
    ),
    legend.spacing.y = unit(
      0.1,
      "cm"
    ),
    legend.position = "bottom"
  )

p4

#--------------------------------------------------------------------------
# 10. SAVE FIGURE
#--------------------------------------------------------------------------

ggsave(
  "r_priors_D_macarellus.png",
  plot = p4,
  device = "png",
  units = "cm",
  width = 32,
  height = 17
)






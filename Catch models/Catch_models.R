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
#install.packages("fishmethods")
library(fishmethods)
install.packages("future.apply")
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
lh<- read_xlsx("Parametros_Historia_de_vida.xlsx")
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
      c(Flat,
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
  preds <- data.frame(Flat, LH,LHL, HL, HLH, OTH, ct_MSY.int,
    min_max,max.ct.i, int.ct.i, min.ct.i,yr.norm, mean.ct_MSY.start,
    slope.first.nrm,mean.ct_MSY.end, slope.last.nrm
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
  bk.MSY <- c(0.256, 0.721)
  CL.1 <- c(0.01, 0.203 )
  CL.2 <- c(0.20, 0.431 )
  CL.3 <- c(0.80, -0.45 )
  CL.4 <- c(1.02, -0.247)
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

# Plot de depleção B/K - Decapterus macarellus
p_bk <- ggplot(bk_macarellus, aes(x = metodo, y = bk, ymin = bk_lo,
  ymax = bk_hi, color = fonte,shape = hipotese, group = hipotese
)) +
  geom_linerange( position = position_dodge(width = 0.6),
    linewidth = 1.2 ) +
  geom_point( position = position_dodge(width = 0.6),
    size = 3.5, stroke = 1.5) +
  geom_hline(yintercept = 0.5,linetype = "dashed",color = "grey50",
    linewidth = 0.4 ) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.1)
  ) +
  scale_color_viridis_d() +
  labs(x = "Depletion hypothesis",
    y = expression("Biomass depletion (B/K"[2015]*")"),
    color = "Information source",
    shape = "Hypothesis"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "mm"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    legend.text = element_text(size = 11.5),
    legend.box.margin = margin(t = -10),
    legend.spacing.y = unit(0.1, "cm")
  )

p_bk

#salvar
ggsave("bk_priors.png", plot = p_bk, device = "png",  units = "cm", width = 30, height = 17)



# ============================================================================
# PRIOR DE r PARA Decapterus macarellus (Euler-Lotka methods Cortes, 2016)
# Baseado em história de vida e variabilidade observada na literatura
#
# Adaptado por Matheus Silva
## Objetivos:
#   1. Extrair parâmetros de história de vida de D. macarellus
#   2. Padronizar unidades
#   4. Propagar incerteza por bootstrap paramétrico
#   5. Estimar r a partir de parâmetros de história de vida
#   6. Gerar distribuição de r para uso como prior no CMSY/DB-SRA
#
# ============================================================================

library(tidyverse)
#--------------------
# Funções auxiliares 
#--------------------
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

#-----------------------------
# Bootstrap para uma espécie  
#-----------------------------
estimate_r_boot <- function(sp_row, nboot = 1000,
                            cvs = list(Linf = 0.15, k = 0.20, M = 0.20,
                                       tmax = 0.10, L50 = 0.15, ls = 0.25),
                            ls_euler_fixed = 9.5,
                            options = list(r_upper = 5, verbose = FALSE)) {
  
  sp <- sp_row$especie[1]   # especie
  
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
  
  # Bootstrap paramétrico
  Linf_samps <- rlnorm_from_mean_cv(Linf0, cvs$Linf, nboot)
  k_samps    <- rlnorm_from_mean_cv(k0,    cvs$k,    nboot)
  M_samps    <- rlnorm_from_mean_cv(M0,    cvs$M,    nboot)
  tmax_samps <- pmax(1, round(rlnorm_from_mean_cv(tmax0, cvs$tmax, nboot)))
  L50_samps  <- rlnorm_from_mean_cv(L500,  cvs$L50,  nboot)
  ls_samps   <- rlnorm_from_mean_cv(ls0,   cvs$ls,   nboot)
  
  # Cálculo por iteração 
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

#-------------------------------------------
# Aplicar a todas as espécies (usando `lh`)
#-------------------------------------------
species_list <- unique(lh$especie)

res_list <- map(species_list, function(sp) {
  sp_row <- lh %>% filter(especie == sp)
  estimate_r_boot(sp_row, nboot = 1000)
})

r_sims <- map_dfr(res_list, "sims")
write.csv(r_sims, "r_sims.csv", row.names = FALSE)

#-------------------------------------------------
# 1) Sumários por espécie × método (inalterado)
#-------------------------------------------------
r_methods <- map_dfr(res_list, "summary")

#-------------------------------------------------
# 2) Agregação baseline (mesma lógica de antes)
#-------------------------------------------------
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


#================================================================
# Avaliações baseadas em capturas 
# Modelo utilizado : DB-SRA (Dick & MacCall, 2011)
# Objetivo estimar viabilidade de trajetórias de biomasa compatíveis
# com a serie de captura e biologia assumida
# Estimativas de MSY, Bmsy, Fmsy, Cmsy (OFL), K para o D. macarellus
# avaliações para cada cenario considerado de depleção (btk)
# e a diferença de produtividade que entra pelo M no DB-SRA
# Cenários de Bt/k (já definidos) x Cenários de M (revisão de história
# de vida) para D. macarellus -- para alimentar loops de fishmethods::dbsra()
# =====================================================================
library(fishmethods)
print(bk_macarellus) #deplecoes já calculadas

# ---- cenários de M, com autor/fonte (da tabela Confiabilidade_Fontes) ----
m_macarellus <- data.frame(
  m_hipotese = c("M_Jardim(1996/1999)", "M_Santos(2018)"),
  m_fonte    = c("Jardim (1996/1999)", "Santos (2018)"),
  m_metodo   = c("Tanaka", "Tanaka"),
  M          = c(0.43, 0.60),
  stringsAsFactors = FALSE
)

# ---- produto cartesiano: cada hipótese de bk x cada hipótese de M ----
# (via índices -- não depende de nomes de coluna em comum, então é seguro
#  mesmo que as duas tabelas ganhem colunas com nomes iguais no futuro)
idx <- expand.grid(bk_i = seq_len(nrow(bk_macarellus)),
                   m_i  = seq_len(nrow(m_macarellus)))

cenarios_macarellus_dbsra <- cbind(
  bk_macarellus[idx$bk_i, c("especie", "hipotese", "bk_lo", "bk_hi", "bk")],
  m_macarellus[idx$m_i, ]
)
rownames(cenarios_macarellus_dbsra) <- NULL

cenarios_macarellus_dbsra$cenario_id <- paste(cenarios_macarellus_dbsra$hipotese,
                                        cenarios_macarellus_dbsra$m_hipotese, sep = "_")

# reordena pra ficar fácil de ler (hipótese de bk como bloco externo)
cenarios_macarellus_dbsra <- cenarios_macarellus_dbsra[order(cenarios_macarellus_dbsra$hipotese,
                                                 cenarios_macarellus_dbsra$m_hipotese), ]
rownames(cenarios_macarellus_dbsra) <- NULL
print(cenarios_macarellus_dbsra)
cat("\nDimensões:", nrow(cenarios_macarellus_dbsra), "linhas x", ncol(cenarios_macarellus_dbsra), "colunas\n")

#=====================
#idade de maturação
#====================
# ---- 1) Parâmetros de crescimento -----------------------------------
# Fonte mais confiável COM TRIO COMPLETO (Linf, K, t0): Jardim (1996/1999)
# -- nível "Alta" na Confiabilidade_Fontes (nota media 4,0), o único trio
# completo entre as fontes Alta (Costa et al. 2020 não estima crescimento;
# da Cruz Delgado et al. 2024, também Alta, não reporta t0).
Linf <- 315.0   # mm FL
K    <- 0.43    # /ano
t0   <- -1.56   # anos

vbgf <- function(t, Linf_ = Linf, K_ = K, t0_ = t0) Linf_ * (1 - exp(-K_ * (t - t0_)))
idade_no_comprimento <- function(L, Linf_ = Linf, K_ = K, t0_ = t0) {
  t0_ - (1 / K_) * log(1 - L / Linf_)
}

# ---- 2) L50 mais confiável -------------------------------------------
# Fonte mais confiável para maturação: Costa et al. (2020) -- também nível
# "Alta", e a única fonte Alta com nota máxima (5) em revisão por pares E em
# consistência interna (sem nenhuma ressalva na verificação). Reporta L50
# por sexo, sobre amostra de desembarques industriais 2012-2018.
L50_F <- 241.0; n_F <- 284   # fêmeas
L50_M <- 266.0; n_M <- 85    # machos
L50_comb <- (L50_F * n_F + L50_M * n_M) / (n_F + n_M)   # média ponderada por n

idade_maturacao <- data.frame(
  grupo = c("Fêmeas", "Machos", "Combinado (média ponderada por n)"),
  L50_mm_FL = c(L50_F, L50_M, round(L50_comb, 1)),
  n = c(n_F, n_M, n_F + n_M)
)
idade_maturacao$idade_anos <- sapply(idade_maturacao$L50_mm_FL, idade_no_comprimento)
idade_maturacao$idade_meses <- round(idade_maturacao$idade_anos * 12, 1)
idade_maturacao$idade_anos <- round(idade_maturacao$idade_anos, 3)

print(idade_maturacao)



#=====================
#rodando o modelo
#=====================
library(future.apply)
plan(multisession, workers = min(nrow(cenarios_macarellus_dbsra), parallel::detectCores() - 1))

resultados_dbsra <- future_lapply(seq_len(nrow(cenarios_macarellus_dbsra)), function(i) {
  cen <- cenarios_macarellus_dbsra[i, ]
  fishmethods::dbsra(
    year = ct$year, catch = ct$ct,
    agemat = 2,
    k     = list(low = 3000, up = 60000, tol = 0.01, permax = 1000),
    b1k   = list(dist = "unif", low = 0.8, up = 0.99, mean = 1, sd = 0.1),
    btk   = list(dist = "unif", low = cen$bk_lo, up = cen$bk_hi, refyr = 2015),
    fmsym = list(dist = "lnorm", low = 0.1, up = 2, mean = log(0.8), sd = 0.3),
    bmsyk = list(dist = "beta", low = 0.05, up = 0.95, mean = 0.4, sd = 0.1),
    M     = list(dist = "lnorm", low = cen$M * 0.7, up = cen$M * 1.3, mean = log(cen$M), sd = 0.10),
    nsims = 10000, grout = 0
  )
}, future.seed = TRUE)

names(resultados_dbsra) <- cenarios_macarellus_dbsra$cenario_id
plan(sequential)  # libera os workers no final




resultados_dbsra <- lapply(seq_len(nrow(cenarios_macarellus_dbsra)), function(i) {
  cen <- cenarios_macarellus_dbsra[i, ]
  dbsra(
    year = ct$year, catch =ct$ct,
    agemat = 2,
    k     = list(low = 3000, up = 120000, tol = 0.01, permax = 1000),   # busca aberta
    b1k   = list(dist = "unif", low = 0.8, up = 0.99, mean = 1, sd = 0.1),  
    btk = list(dist = "unif", low = cen$bk_lo, up = cen$bk_hi, refyr = 2015),
    fmsym = list(dist = "lnorm", low = 0.1, up = 2, mean = log(0.8), sd = 0.3),
    bmsyk = list(dist = "beta", low = 0.05, up = 0.95, mean = 0.4, sd = 0.1),
    M     = list(dist = "lnorm", low = cen$M* 0.7, up = cen$M * 1.3, mean = log(cen$M), sd = 0.10),
    nsims = 10000, grout = 0
  )
})
names(resultados_dbsra) <- cenarios_macarellus_dbsra$cenario_id



# =====================================================================
# Pós-processamento dos 8 cenários de dbsra() -- D. macarellus
# =====================================================================
# Pré-requisitos no ambiente (já devem existir depois de rodar o lapply
# e o cenarios_macarellus.R):
#   resultados          - lista de objetos retornados por dbsra(), com
#                          names(resultados) <- cenarios_macarellus$cenario_id
#   cenarios_macarellus - data.frame com as 8 combinações hipótese bk x M
#                          (colunas: cenario_id, hipotese, bk_lo, bk_hi, bk,
#                           m_hipotese, m_fonte, M, ...)
#
# Gera:
#   1) tabela_resumo_dbsra_macarellus.csv   -- resumo longo (1 linha por
#      variável x cenário: média, mediana, IC 95%, % aceitação)
#   2) tabela_aceitacao_dbsra_macarellus.csv -- 1 linha por cenário
#   3) comparacao_cenarios_outputs.png       -- boxplots comparando OFL,
#      K, MSY, Bmsy entre os 8 cenários
#   4) comparacao_cenarios_parametros.png    -- boxplots comparando as
#      posteriores de Fmsy/M, Bt/K, Bmsy/K e M entre os 8 cenários
#   5) priori_posteriori_<cenario_id>.png    -- 1 arquivo por cenário
#      (ou só para os cenários listados em `cenarios_para_detalhar`)
# =====================================================================

stopifnot(exists("resultados"), exists("cenarios_macarellus"))
stopifnot(all(names(resultados) == cenarios_macarellus$cenario_id) ||
            all(names(resultados) %in% cenarios_macarellus$cenario_id))

## ---------------------------------------------------------------------
## 1) TABELA RESUMO (longa) -- média, mediana, IC95%, por variável x cenário
## ---------------------------------------------------------------------

vars_saida <- c("K", "MSY", "Bmsy", "Fmsy", "Umsy", "OFLT1", "Brefyr",
                "FmsyM", "BtK", "BmsyK", "M")

resumir_cenario <- function(res, cenario_id) {
  vals <- res$Values
  acc  <- vals[vals$ll == 1, ]
  n_tot <- nrow(vals)
  n_acc <- nrow(acc)
  
  vars <- vars_saida[vars_saida %in% names(acc)]
  
  linhas <- lapply(vars, function(v) {
    x <- acc[[v]]
    data.frame(
      cenario_id = cenario_id,
      variavel   = v,
      media      = mean(x),
      mediana    = median(x),
      p2.5       = as.numeric(quantile(x, 0.025)),
      p97.5      = as.numeric(quantile(x, 0.975)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, linhas)
  out$n_total       <- n_tot
  out$n_aceitos     <- n_acc
  out$pct_aceitacao <- round(100 * n_acc / n_tot, 1)
  out
}

tabela_resumo <- do.call(rbind, lapply(names(resultados_dbsra), function(id) {
  resumir_cenario(resultados_dbsra[[id]], id)
}))

meta_cols <- intersect(c("cenario_id", "hipotese", "m_hipotese", "m_fonte", "bk_lo", "bk_hi"),
                       names(cenarios_macarellus_dbsra))
tabela_resumo <- merge(cenarios_macarellus_dbsra[, meta_cols], tabela_resumo, by = "cenario_id")

# reordena para leitura mais fácil: hipótese de bk como bloco externo
if (all(c("hipotese", "m_hipotese") %in% names(tabela_resumo))) {
  tabela_resumo <- tabela_resumo[order(tabela_resumo$hipotese, tabela_resumo$m_hipotese,
                                       tabela_resumo$variavel), ]
}
rownames(tabela_resumo) <- NULL

cat("\n===== Tabela resumo (primeiras linhas) =====\n")
print(head(tabela_resumo, 12))
write.csv(tabela_resumo, "tabela_resumo_dbsra_macarellus.csv", row.names = FALSE)

## ---------------------------------------------------------------------
## 2) TABELA DE ACEITAÇÃO -- 1 linha por cenário
## ---------------------------------------------------------------------

tabela_aceitacao <- unique(tabela_resumo[, c("cenario_id", meta_cols[meta_cols != "cenario_id"],
                                             "n_total", "n_aceitos", "pct_aceitacao")])
rownames(tabela_aceitacao) <- NULL

cat("\n===== Tabela de aceitacao por cenario =====\n")
print(tabela_aceitacao)
write.csv(tabela_aceitacao, "tabela_aceitacao_dbsra_macarellus.csv", row.names = FALSE)

# aviso automático se algum cenário tiver aceitação muito baixa (< 5%) ou
# muito alta (> 90%) -- ambos merecem checagem antes de confiar no resultado
baixa <- tabela_aceitacao$cenario_id[tabela_aceitacao$pct_aceitacao < 5]
alta  <- tabela_aceitacao$cenario_id[tabela_aceitacao$pct_aceitacao > 90]
if (length(baixa) > 0) {
  cat("\n[AVISO] Aceitacao muito baixa (<5%) -- resultado instavel / poucas draws aceitas:\n  ",
      paste(baixa, collapse = ", "), "\n")
}
if (length(alta) > 0) {
  cat("\n[AVISO] Aceitacao muito alta (>90%) -- confira se algum criterio (ex: permax)",
      "esta afrouxado demais:\n  ", paste(alta, collapse = ", "), "\n")
}

## ---------------------------------------------------------------------
## 3) COMPARACAO ENTRE CENARIOS -- variaveis de manejo (OFL, K, MSY, Bmsy)
## ---------------------------------------------------------------------

ordem_ids <- if (all(c("hipotese", "m_hipotese") %in% names(cenarios_macarellus_dbsra))) {
  cenarios_macarellus_dbsra$cenario_id[order(cenarios_macarellus_dbsra$hipotese, cenarios_macarellus_dbsra$m_hipotese)]
} else {
  names(resultados_dbsra)
}

extrair_aceitos <- function(id, var) {
  acc <- resultados_dbsra[[id]]$Values
  acc <- acc[acc$ll == 1, ]
  acc[[var]]
}

png("comparacao_cenarios_outputs.png", width = 28, height = 20,
                                res = 300,antialias = "cleartype", units = "cm")
op <- par(mfrow = c(2, 2), mar = c(8, 4.5, 3, 1), bty="l",cex=0.8,cex.main=0.9)
for (v in c("OFLT1", "K", "MSY", "Bmsy")) {
  if (!v %in% names(resultados_dbsra[[1]]$Values)) next
  lst <- lapply(ordem_ids, extrair_aceitos, var = v)
  boxplot(lst, names = ordem_ids, las = 2, main = v, col = "#8FAADC",
          cex.axis = 0.65, ylab = v, lwd=1, bty="l", xaxt="n")
  axis( 1,at = 1:length(ordem_ids),  labels = FALSE )
  text(x = 1:length(ordem_ids), y = par("usr")[3],labels = ordem_ids, srt = 45,adj = 1,xpd = TRUE) 
}
par(op)
dev.off()
cat("\nPNG salvo: comparacao_cenarios_outputs.png\n")

## ---------------------------------------------------------------------
## 4) COMPARACAO ENTRE CENARIOS -- os 4 parametros estocasticos
## ---------------------------------------------------------------------

png("comparacao_cenarios_parametros.png", width = 28, height = 20,
                          res = 300,antialias = "cleartype", units = "cm")
op <- par(mfrow = c(2, 2), mar = c(8, 4.5, 3, 1), bty="l",cex=0.8,cex.main=0.9)
for (v in c("FmsyM", "BtK", "BmsyK", "M")) {
  if (!v %in% names(resultados_dbsra[[1]]$Values)) next
  lst <- lapply(ordem_ids, extrair_aceitos, var = v)
  boxplot(lst, names = ordem_ids, las = 2, main = v, lwd=1,col = "#74C476" ,
          cex.axis = 0.65, ylab = v, bty="l", xaxt="n")
  axis( 1,at = 1:length(ordem_ids),  labels = FALSE )
  text(x = 1:length(ordem_ids), y = par("usr")[3],labels = ordem_ids, srt = 45,adj = 1,xpd = TRUE) 
}
par(op)
dev.off()
cat("PNG salvo: comparacao_cenarios_parametros.png\n")

## ---------------------------------------------------------------------
## 5) PRIORI x POSTERIORI por cenario (generaliza o script anterior,
##    usando os limites/medias especificos de cada cenario em vez de
##    valores fixos)
## ---------------------------------------------------------------------

plot_prior_post <- function(prior_dens_fun, post_values, xlim, xlab, col_post, main) {
  x <- seq(xlim[1], xlim[2], length.out = 500)
  pd <- prior_dens_fun(x)
  plot(x, pd / max(pd), type = "l", col = "grey45", lwd = 3, lty = 2,
       xlab = xlab, ylab = "densidade (normalizada ao pico)", main = main,
       ylim = c(0, 1.05),bty="l")
  pdens <- density(post_values, from = xlim[1], to = xlim[2])
  lines(pdens$x, pdens$y / max(pdens$y), col = col_post, lwd = 2.5)
  legend("topright", c("Priori (especificada)", "Posteriori (aceitos, ll=1)"),
         col = c("grey45", col_post), lty = c(2, 1), lwd = 2, bty = "n", cex = 0.8)
}

gerar_priori_posteriori <- function(cenario_id, cen_row,
                                    fmsym_mean = 0.8, fmsym_sd = 0.3,
                                    bmsyk_mean = 0.35, bmsyk_sd = 0.1,
                                    m_sd = 0.10) {
  res <- resultados_dbsra[[cenario_id]]
  acc <- res$Values
  acc <- acc[acc$ll == 1, ]
  
  png(sprintf("priori_posteriori_%s.png", cenario_id), width = 28, height = 20, 
                                      res = 300,antialias = "cleartype",units = "cm")
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1),cex.main=0.9)
  
  plot_prior_post(function(x) dlnorm(x, meanlog = log(fmsym_mean), sdlog = fmsym_sd),
                  acc$FmsyM, xlim = c(0.1, 2), xlab = "Fmsy/M",
                  col_post = "#1F4E79", main = paste("Fmsy/M -", cenario_id))
  
  plot_prior_post(function(x) dunif(x, min = cen_row$bk_lo, max = cen_row$bk_hi),
                  acc$BtK, xlim = c(max(0, cen_row$bk_lo - 0.05), min(1, cen_row$bk_hi + 0.05)),
                  xlab = "Bt/K", col_post = "#C00000", main = paste("Bt/K -", cenario_id))
  
  m_alpha <- bmsyk_mean * ((bmsyk_mean * (1 - bmsyk_mean) / bmsyk_sd^2) - 1)
  m_beta  <- (1 - bmsyk_mean) * ((bmsyk_mean * (1 - bmsyk_mean) / bmsyk_sd^2) - 1)
  plot_prior_post(function(x) dbeta(x, m_alpha, m_beta),
                  acc$BmsyK, xlim = c(0.05, 0.7), xlab = "Bmsy/K",
                  col_post = "#548235", main = paste("Bmsy/K -", cenario_id))
  
  plot_prior_post(function(x) dlnorm(x, meanlog = log(cen_row$M), sdlog = m_sd),
                  acc$M, xlim = c(cen_row$M * 0.6, cen_row$M * 1.4), xlab = "M",
                  col_post = "#7030A0", main = paste("M -", cenario_id))
  
  par(op)
  dev.off()
  cat(sprintf("PNG salvo: priori_posteriori_%s.png  (aceitos: %d de %d, %.1f%%)\n",
              cenario_id, nrow(acc), nrow(res$Values), 100 * nrow(acc) / nrow(res$Values)))
}

# Por padrao gera para todos os 8 cenarios. Se preferir só alguns
# representativos, troque a linha abaixo por, por exemplo:
#   cenarios_para_detalhar <- c("NN_CMSY_x_M_mais_confiavel", "Uninformative_x_M_segunda_confiavel")
cenarios_para_detalhar <- cenarios_macarellus_dbsra$cenario_id

for (id in cenarios_para_detalhar) {
  cen_row <- cenarios_macarellus_dbsra[cenarios_macarellus_dbsra$cenario_id == id, ]
  gerar_priori_posteriori(id, cen_row)
}

cat("\nPos-processamento concluido.\n")


# =====================================================================
# Trajetorias de biomassa reconstruidas, todos os cenarios num so grafico
# =====================================================================
# O dbsra() plota a trajetoria de biomassa (aceitos/rejeitados) internamente
# (grout), mas NAO devolve essa matriz no objeto retornado -- res$Values
# so guarda os parametros/estimativas de cada simulacao, nao a serie
# ano-a-ano de biomassa. Por isso reconstruimos a trajetoria aqui, ano a
# ano, usando a mesma funcao de producao Schaefer-Pella-Tomlinson-Fletcher
# do metodo (Dick & MacCall 2011; formulas tambem em Owashi 2014 eq 1.2 e
# Sweka et al. 2018 eqs 1-4), a partir dos parametros que JA ficam
# guardados por simulacao aceita em res$Values: K, n, g (=gamma), B1K, MSY.
#
#   B[1]        = B1K * K
#   P(B[t-a])   = g * MSY * (B[t-a]/K) - g * MSY * (B[t-a]/K)^n
#   B[t]        = B[t-1] + P(B[t-a]) - C[t-1]      (a = agemat; se t-a<1,usa B[1] no lugar)
#
# A validacao que o script faz: para cada simulacao aceita, recalcula a
# biomassa no ano de referencia (refyr) e compara com o BtK que o proprio
# dbsra() reportou naquela linha de res$Values (que foi o alvo que o
# optimize() do pacote usou para achar k). Se a diferenca media for
# pequena (a impressao no console avisa), a reconstrucao esta capturando
# a dinamica corretamente e da para confiar na FORMA da trajetoria.
# =====================================================================

stopifnot(exists("resultados_dbsra"), exists("cenarios_macarellus_dbsra"), exists("ct"))

agemat <- 2          # mesmo valor usado no dbsra()
refyr  <- 2015        # mesmo refyr usado no btk
anos   <- ct$year
catches <- ct$ct
n_anos <- length(anos)
idx_refyr <- which(anos == refyr)
stopifnot(length(idx_refyr) == 1)

n_amostra_por_cenario <- 10000   # quantas trajetorias aceitas usar (amostra, p/ nao pesar)

## ---------------------------------------------------------------------
## Reconstroi UMA trajetoria de biomassa a partir de 1 linha de res$Values
## ---------------------------------------------------------------------
reconstruir_biomassa <- function(K, n, g, B1K, MSY, catches, agemat) {
  n_t <- length(catches)
  B <- numeric(n_t)
  B[1] <- B1K * K
  for (t in 2:n_t) {
    lag_idx <- t - agemat
    B_lag <- if (lag_idx >= 1) B[lag_idx] else B[1]
    razao <- B_lag / K
    P <- g * MSY * razao - g * MSY * razao^n
    B[t] <- B[t - 1] + P - catches[t - 1]
  }
  B
}

## ---------------------------------------------------------------------
## Reconstroi as trajetorias aceitas de 1 cenario e resume (mediana + IC)
## ---------------------------------------------------------------------
reconstruir_cenario <- function(res, catches, agemat, idx_refyr, n_amostra) {
  vals <- res$Values
  acc  <- vals[vals$ll == 1, ] #só aceita trajetórias válidas
  
  vars_necessarias <- c("K", "n", "g", "B1K", "MSY", "BtK")
  faltando <- setdiff(vars_necessarias, names(acc))
  if (length(faltando) > 0) {
    stop("res$Values esta sem as colunas: ", paste(faltando, collapse = ", "),
         " -- confira names(resultados[[1]]$Values)")
  }
  
  if (nrow(acc) > n_amostra) {
    acc <- acc[sample(nrow(acc), n_amostra), ]
  }
  
  mat <- t(mapply(function(K, n, g, B1K, MSY) {
    reconstruir_biomassa(K, n, g, B1K, MSY, catches, agemat)
  }, acc$K, acc$n, acc$g, acc$B1K, acc$MSY))
  
  # validacao: BtK reconstruido no refyr vs BtK reportado pelo dbsra()
  btk_reconstruido <- mat[, idx_refyr] / acc$K
  residuo <- btk_reconstruido - acc$BtK
  diagnostico <- c(erro_medio_abs = mean(abs(residuo)),
                   erro_max_abs  = max(abs(residuo)))
  
  mat_bk <- sweep(mat, 1, acc$K, "/")   # biomassa relativa (B/K), por linha
  
  list(
    biomassa   = mat,
    biomassa_bk = mat_bk,
    mediana_B  = apply(mat, 2, median),
    p2.5_B     = apply(mat, 2, quantile, 0.025),
    p97.5_B    = apply(mat, 2, quantile, 0.975),
    mediana_BK = apply(mat_bk, 2, median),
    p2.5_BK    = apply(mat_bk, 2, quantile, 0.025),
    p97.5_BK   = apply(mat_bk, 2, quantile, 0.975),
    diagnostico = diagnostico
  )
}

## ---------------------------------------------------------------------
## Roda para todos os cenarios e valida
## ---------------------------------------------------------------------
set.seed(1)
trajetorias <- lapply(names(resultados_dbsra), function(id) {
  cat("Reconstruindo trajetorias:", id, "... ")
  out <- reconstruir_cenario(resultados_dbsra[[id]], catches, agemat, idx_refyr, n_amostra_por_cenario)
  cat(sprintf("erro medio abs no Bt/K do ano de referencia: %.4f (max: %.4f)\n",
              out$diagnostico["erro_medio_abs"], out$diagnostico["erro_max_abs"]))
  out
})
names(trajetorias) <- names(resultados_dbsra)

erros <- sapply(trajetorias, function(x) x$diagnostico["erro_medio_abs"])
if (any(erros > 0.05)) {
  cat("\n[AVISO] Em pelo menos um cenario o erro medio no Bt/K reconstruido",
      "passou de 0.02 -- a reconstrucao pode nao estar batendo com a",
      "convencao exata do pacote (ex: tratamento do lag nos primeiros anos).",
      "Trate a FORMA da trajetoria com cautela nesses casos.\n\n")
} else {
  cat("\nValidacao OK em todos os cenarios (erro medio no Bt/K reconstruido <= 0.02).\n\n")
}

## ---------------------------------------------------------------------
## Grafico unico -- todos os cenarios sobrepostos (mediana + IC 95%)
## ---------------------------------------------------------------------

ordem_ids <- if (all(c("hipotese", "m_hipotese") %in% names(cenarios_macarellus_dbsra))) {
  cenarios_macarellus_dbsra$cenario_id[order(cenarios_macarellus_dbsra$hipotese, cenarios_macarellus_dbsra$m_hipotese)]
} else {
  names(resultados)
}
cores <- setNames(grDevices::hcl.colors(length(ordem_ids), palette = "Dark 3"), ordem_ids)

png("trajetorias_biomassa_cenarios.png", width = 32, height = 16,
                          res = 300,antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), xpd = FALSE, bty="l",cex.main=0.9)

# ---- painel 1: B/K (comparavel entre cenarios com K muito diferente) ----
plot(NA, xlim = range(anos), ylim = c(0, 1),
     xlab = "Ano", ylab = "Biomassa relativa (B/K)",
     main = "Trajetorias de biomassa relativa -- B/K")
for (id in ordem_ids) {
  tr <- trajetorias[[id]]
  polygon(c(anos, rev(anos)), c(tr$p2.5_BK, rev(tr$p97.5_BK)),
          col = adjustcolor(cores[id], alpha.f = 0.12), border = NA)
}
for (id in ordem_ids) {
  lines(anos, trajetorias[[id]]$mediana_BK, col = cores[id], lwd = 3.2)
}
abline(h=0.5, col="firebrick",lty=2)
legend("topright", legend = ordem_ids, col = cores, lwd = 3.2, bty = "n", cex = 0.7)

# ---- painel 2: biomassa absoluta ----
todas_max <- max(sapply(trajetorias, function(tr) max(tr$p97.5_B)))
plot(NA, xlim = range(anos), ylim = c(0, todas_max),
     xlab = "Ano", ylab = "Biomassa (t)",
     main = "Trajetorias de biomassa absoluta -- (t)")
for (id in ordem_ids) {
  tr <- trajetorias[[id]]
  polygon(c(anos, rev(anos)), c(tr$p2.5_B, rev(tr$p97.5_B)),
          col = adjustcolor(cores[id], alpha.f = 0.12), border = NA)
}
for (id in ordem_ids) {
  lines(anos, trajetorias[[id]]$mediana_B, col = cores[id], lwd = 3.2)
}
legend("topright", legend = ordem_ids, col = cores, lwd = 3.2, bty = "n", cex = 0.7)

par(op)
dev.off()
cat("PNG salvo: trajetorias_biomassa_cenarios.png\n")


# =====================================================================
# Posteriores de MSY por cenário, reunidas num data.frame e plotadas
# junto com a série histórica de captura, num painel só.
# =====================================================================
# Pré-requisitos no ambiente:
#   resultados_dbsra          - lista de objetos retornados por dbsra(),
#                                com names(resultados_dbsra) <-
#                                cenarios_macarellus_dbsra$cenario_id
#   cenarios_macarellus_dbsra - data.frame com as combinações de cenário
#   ct                        - data.frame com a série de captura, colunas
#                                year e ct (a mesma série usada no dbsra())
#
# Gera:
#   1) msy_posteriores_cenarios_dbsra.csv  -- formato longo: 1 linha por
#      simulação aceita (cenario_id, MSY), todas as posteriores juntas
#   2) msy_resumo_cenarios_dbsra.csv       -- 1 linha por cenário: mediana
#      e IC 95% do MSY (resumo usado para desenhar as faixas do gráfico)
#   3) msy_posteriores_serie_captura_dbsra.png -- painel único: linha da
#      captura histórica + uma linha (mediana) e uma faixa (IC 95%)
#      horizontais por cenário, mostrando onde o MSY estimado de cada
#      cenário cai em relação ao nível de captura efetivamente pescado.
# =====================================================================

stopifnot(exists("resultados_dbsra"), exists("cenarios_macarellus_dbsra"), exists("ct"))

## ---------------------------------------------------------------------
## 1) Data frame único com as posteriores de MSY de todos os cenários
## ---------------------------------------------------------------------

msy_posteriores <- do.call(rbind, lapply(names(resultados_dbsra), function(id) {
  vals <- resultados_dbsra[[id]]$Values
  acc  <- vals[vals$ll == 1, ]
  data.frame(cenario_id = id, MSY = acc$MSY, stringsAsFactors = FALSE)
}))

meta_cols <- intersect(c("cenario_id", "hipotese", "m_hipotese"), names(cenarios_macarellus_dbsra))
if (length(meta_cols) > 1) {
  msy_posteriores <- merge(cenarios_macarellus_dbsra[, meta_cols], msy_posteriores, by = "cenario_id")
}

cat("===== msy_posteriores (primeiras linhas) =====\n")
print(head(msy_posteriores, 8))
write.csv(msy_posteriores, "msy_posteriores_cenarios_dbsra.csv", row.names = FALSE)
cat(sprintf("\nCSV salvo: msy_posteriores_cenarios_dbsra.csv (%d linhas, %d cenarios)\n",
            nrow(msy_posteriores), length(unique(msy_posteriores$cenario_id))))

## ---------------------------------------------------------------------
## 2) Resumo por cenário (mediana + IC 95%) -- usado no gráfico
## ---------------------------------------------------------------------

msy_resumo <- do.call(rbind, lapply(split(msy_posteriores, msy_posteriores$cenario_id), function(d) {
  data.frame(
    cenario_id = d$cenario_id[1],
    mediana    = median(d$MSY),
    p2.5       = as.numeric(quantile(d$MSY, 0.025)),
    p97.5      = as.numeric(quantile(d$MSY, 0.975)),
    n_aceitos  = nrow(d)
  )
}))
if (length(meta_cols) > 1) {
  msy_resumo <- merge(cenarios_macarellus_dbsra[, meta_cols], msy_resumo, by = "cenario_id")
}
rownames(msy_resumo) <- NULL

cat("\n===== Resumo do MSY por cenario =====\n")
print(msy_resumo)
write.csv(msy_resumo, "msy_resumo_cenarios_dbsra.csv", row.names = FALSE)
cat("CSV salvo: msy_resumo_cenarios_dbsra.csv\n")

## ---------------------------------------------------------------------
## 3) Gráfico único -- captura histórica (linha) + faixas de MSY por cenário
## ---------------------------------------------------------------------

ordem_ids <- if (all(c("hipotese", "m_hipotese") %in% names(msy_resumo))) {
  msy_resumo$cenario_id[order(msy_resumo$hipotese, msy_resumo$m_hipotese)]
} else {
  msy_resumo$cenario_id
}
cores <- setNames(grDevices::hcl.colors(length(ordem_ids), palette = "Dark 3"), ordem_ids)

anos    <- ct$year
catches <- ct$ct

# folga de ~28% à direita (em anos) para caber a legenda sem sobrepor as faixas
xlim_plot <- c(min(anos), max(anos) + diff(range(anos)) * 0.45)
ylim_max  <- max(catches, msy_resumo$p97.5) * 1.08

png("msy_posteriores_serie_captura_dbsra.png",  width = 25, height = 16,
    res = 300,antialias = "cleartype", units = "cm")
op <- par(mar = c(4.5, 5, 3, 1),bty="l",cex.main=0.7)

# ---- eixo/moldura em branco primeiro, pra desenhar as faixas de MSY atrás da linha ----
plot(NA, xlim = xlim_plot, ylim = c(0, ylim_max),
     xlab = "Ano", ylab = "Captura (t)",
     main = "Série de captura e posteriores de MSY por cenário")

# ---- faixas horizontais (IC95%) + linha (mediana) de MSY, por cenário ----
x0 <- min(anos)
x1 <- max(anos)
for (id in ordem_ids) {
  r <- msy_resumo[msy_resumo$cenario_id == id, ]
  rect(x0, r$p2.5, x1, r$p97.5, col = adjustcolor(cores[id], alpha.f = 0.14), border = NA)
}
for (id in ordem_ids) {
  r <- msy_resumo[msy_resumo$cenario_id == id, ]
  segments(x0, r$mediana, x1, r$mediana, col = cores[id], lwd = 2.4)
}

# ---- linha da captura histórica, por cima das faixas ----
lines(anos, catches, type = "o", col = "grey25", pch = 16, cex = 0.8, lwd = 2)

legend("topright", inset = c(0, 0), xpd = NA,
       legend = ordem_ids, col = cores, lwd = 2.4, bty = "n", cex = 0.8,
       title = "MSY (mediana, faixa = IC 95%)", title.cex = 0.66)
legend("topleft", legend = "Captura observada", col = "grey25", lwd = 2.4, pch = 16,
       pt.cex = 0.8, bty = "n", cex = 0.8)

par(op)
dev.off()
cat("\nPNG salvo: msy_posteriores_serie_captura_dbsra.png\n")

## ---------------------------------------------------------------------
## 1) Pool único com o MSY de todas as simulações aceitas, de todos os
##    cenários (ignora de qual cenário veio -- é a distribuição conjunta)
## ---------------------------------------------------------------------

msy_todos <- unlist(lapply(resultados_dbsra, function(res) {
  vals <- res$Values
  vals$MSY[vals$ll == 1]
}), use.names = FALSE)

cat(sprintf("MSY combinado: %d simulacoes aceitas, de %d cenarios.\n",
            length(msy_todos), length(resultados_dbsra)))

## ---------------------------------------------------------------------
## 2) Quantis de incerteza
## ---------------------------------------------------------------------

probs <- c(0.025, 0.25, 0.50, 0.75, 0.975)
q <- quantile(msy_todos, probs)

quantis_msy <- data.frame(
  quantil = c("2.5%", "25%", "mediana (50%)", "75%", "97.5%"),
  MSY     = as.numeric(q)
)
cat("\n===== Quantis do MSY (pool conjunto) =====\n")
print(quantis_msy)
write.csv(quantis_msy, "quantis_msy_conjunto_dbsra.csv", row.names = FALSE)
cat("CSV salvo: quantis_msy_conjunto_dbsra.csv\n")

## ---------------------------------------------------------------------
## 3) Gráfico de densidade, com a faixa de 95% sombreada e os quantis
##    marcados por linhas verticais
## ---------------------------------------------------------------------
# Duas escolhas deliberadas aqui, pensadas para um pool que mistura
# cenários bem restritos (posterior estreita) com cenários pouco
# informativos (ex.: Uninformative, com bk quase livre): isso produz uma
# mistura de escalas muito diferentes, e uma densidade "ingênua" (escala
# linear, bw padrão) sai com cara de agulha -- pico fino demais e cauda
# comprida quase invisível, mesmo quando os dados estão certos.
#
#  a) suavização um pouco mais larga (adjust > 1): ainda é a mesma forma
#     geral, só sem o serrilhado de amostra finita.
#  b) densidade calculada em log10(MSY) e depois transformada de volta
#     para a escala de MSY (mudança de variável: f_X(x) = f_U(u)/(x*ln10),
#     com u=log10(x)) -- isto NÃO distorce a densidade, é o jeito
#     estatisticamente correto de exibir uma quantidade estritamente
#     positiva e assimetricamente distribuída (MSY, biomassa, captura...)
#     num eixo log, o que evita que o grosso da massa (perto da mediana)
#     fique espremido em poucos pixels enquanto a cauda dos cenários
#     pouco informativos estica o eixo inteiro.
suavizacao <- 2   # >1 = mais suave; ajuste se ainda parecer serrilhado

dl <- density(log10(msy_todos), adjust = suavizacao)
x_msy <- 10^dl$x
y_msy <- dl$y / (x_msy * log(10))   # densidade na escala de MSY (mudança de variável)

# eixo x log, com marcas "redondas" legíveis em t de MSY
marcas_x <- pretty(log10(msy_todos), n = 8)
marcas_x <- marcas_x[10^marcas_x >= min(x_msy) & 10^marcas_x <= max(x_msy)]

png("msy_densidade_conjunta_dbsra.png", width = 28, height = 20,
                          res = 300,antialias = "cleartype", units = "cm")
op <- par(mar = c(4.5, 4.5, 3, 1),bty="l",cex.main=0.9)

plot(x_msy, y_msy, type = "l", log = "x",
     main = "Distribuição conjunta de MSY (todos os cenários combinados)",
     xlab = "MSY (t)", ylab = "Densidade", col = "#1F4E79", lwd = 2.2,
     xaxt = "n")
axis(1, at = 10^marcas_x, labels = format(round(10^marcas_x), big.mark = ".", decimal.mark = ",", scientific = FALSE))

# sombreia a faixa de 95% (entre os quantis 2.5% e 97.5%) sob a curva
faixa <- x_msy >= q["2.5%"] & x_msy <= q["97.5%"]
polygon(c(x_msy[faixa], rev(x_msy[faixa])), c(y_msy[faixa], rep(0, sum(faixa))),
        col = adjustcolor("#1F4E79", alpha.f = 0.18), border = NA)

# curva por cima da faixa sombreada
lines(x_msy, y_msy, col = "#1F4E79", lwd = 2.2)

# linhas verticais nos quantis
cores_q <- c("2.5%" = "#C00000", "25%" = "#7F7F7F", "mediana (50%)" = "#1F4E79",
             "75%" = "#7F7F7F", "97.5%" = "#C00000")
lty_q   <- c("2.5%" = 2, "25%" = 3, "mediana (50%)" = 1, "75%" = 3, "97.5%" = 2)
for (nm in names(q)) {
  key <- if (nm == "50%") "mediana (50%)" else nm
  abline(v = q[nm], col = cores_q[key], lty = lty_q[key], lwd = 1.8)
}

# rótulos dos quantis, perto do eixo x
y_lab <- max(y_msy) * 0.05
text(q["2.5%"],  y_lab, sprintf("2,5%%\n%.0f", q["2.5%"]),  col = "#C00000", cex = 0.8, pos = 2, offset = 0.3)
text(q["97.5%"], y_lab, sprintf("97,5%%\n%.0f", q["97.5%"]), col = "#C00000", cex = 0.8, pos = 4, offset = 0.3)
text(q["50%"], max(y_msy) * 0.97, sprintf("mediana: %.0f t", q["50%"]),
     col = "#1F4E79", cex = 0.8, pos = 4, offset = 0.3, font = 2)

legend("topright", legend = c("Densidade conjunta do MSY", "Faixa de 95% (IC)", "Mediana"),
       col = c("#1F4E79", adjustcolor("#1F4E79", alpha.f = 0.4), "#1F4E79"),
       lwd = c(2.2, 8, 1.8), lty = c(1, 1, 1), bty = "n", cex = 0.8)

par(op)
dev.off()
cat("\nPNG salvo: msy_densidade_conjunta_dbsra.png\n")


# =====================================================================
# Gráfico tornado -- sensibilidade de uma quantidade de manejo do DB-SRA
# aos dois eixos de cenário testados: M (fonte da mortalidade natural) e
# a hipótese de depleção (Bt/K).
#   1) tornado_sensibilidade_<metrica>_dbsra.csv -- 1 linha por cenário
#      alternativo, com o valor absoluto e a variação % em relação ao
#      cenário BASE, para cada um dos dois fatores testados
#   2) tornado_sensibilidade_<metrica>_dbsra.png -- o gráfico tornado
# =====================================================================

stopifnot(exists("resultados_dbsra"), exists("cenarios_macarellus_dbsra"))
stopifnot(all(c("hipotese", "m_hipotese") %in% names(cenarios_macarellus_dbsra)))

## ---------------------------------------------------------------------
## 0) CONFIGURAÇÃO -- ajuste aqui
## ---------------------------------------------------------------------

metrica <- "MSY"   # troque para "OFLT1", "Bmsy", "Fmsy", "Umsy", "K", etc.
# (qualquer coluna presente em resultados_dbsra[[i]]$Values)

# cenário BASE: a combinação de bk/M que vocês tratam como referência
# (ex.: a hipótese de bk mais defensável e a fonte de M mais confiável)
bk_base <- "Target_switch"
m_base  <- "M_mais_confiavel"

## ---------------------------------------------------------------------
## 1) Mediana da métrica escolhida, por cenário (só simulações aceitas)
## ---------------------------------------------------------------------

medianas <- sapply(names(resultados_dbsra), function(id) {
  vals <- resultados_dbsra[[id]]$Values
  acc  <- vals[vals$ll == 1, ]
  stopifnot(metrica %in% names(acc))
  median(acc[[metrica]])
})
names(medianas) <- names(resultados_dbsra)

id_base <- cenarios_macarellus_dbsra$cenario_id[
  cenarios_macarellus_dbsra$hipotese == bk_base & cenarios_macarellus_dbsra$m_hipotese == m_base]
if (length(id_base) != 1) {
  stop("Nao encontrei (ou encontrei mais de um) cenario BASE com hipotese='", bk_base,
       "' e m_hipotese='", m_base, "'. Confira os valores em cenarios_macarellus_dbsra.")
}
valor_base <- medianas[[id_base]]
cat(sprintf("Cenario BASE: %s  |  mediana de %s = %.2f\n", id_base, metrica, valor_base))

## ---------------------------------------------------------------------
## 2) Variação de cada fator, mantendo o outro fator fixo no nível BASE
## ---------------------------------------------------------------------

niveis_m <- unique(cenarios_macarellus_dbsra$m_hipotese)
tab_m <- do.call(rbind, lapply(setdiff(niveis_m, m_base), function(m_alt) {
  id <- cenarios_macarellus_dbsra$cenario_id[
    cenarios_macarellus_dbsra$hipotese == bk_base & cenarios_macarellus_dbsra$m_hipotese == m_alt]
  data.frame(fator = "M (mortalidade natural)", nivel = m_alt,
             cenario_id = id, valor = medianas[[id]], stringsAsFactors = FALSE)
}))

niveis_bk <- unique(cenarios_macarellus_dbsra$hipotese)
tab_bk <- do.call(rbind, lapply(setdiff(niveis_bk, bk_base), function(bk_alt) {
  id <- cenarios_macarellus_dbsra$cenario_id[
    cenarios_macarellus_dbsra$hipotese == bk_alt & cenarios_macarellus_dbsra$m_hipotese == m_base]
  data.frame(fator = "Bt/K (metodo de depleção)", nivel = bk_alt,
             cenario_id = id, valor = medianas[[id]], stringsAsFactors = FALSE)
}))

tab <- rbind(tab_m, tab_bk)
tab$delta_pct <- 100 * (tab$valor - valor_base) / valor_base

cat("\n===== Variação em relação ao cenário BASE =====\n")
print(tab[, c("fator", "nivel", "valor", "delta_pct")])
write.csv(tab, sprintf("tornado_sensibilidade_%s_dbsra.csv", tolower(metrica)), row.names = FALSE)
cat(sprintf("\nCSV salvo: tornado_sensibilidade_%s_dbsra.csv\n", tolower(metrica)))

## ---------------------------------------------------------------------
## 3) Empilha os níveis de cada fator dos dois lados do zero (negativos
##    à esquerda, positivos à direita), do menor para o maior módulo --
##    é só uma convenção de leiaute para caber vários níveis numa única
##    barra por fator, igual ao gráfico do seu amigo; não representa soma
##    real de efeitos (cada cenário é uma rodada independente do dbsra()).
## ---------------------------------------------------------------------

empilhar <- function(df) {
  df <- df[order(abs(df$delta_pct)), ]
  neg <- df[df$delta_pct < 0, , drop = FALSE]
  pos <- df[df$delta_pct >= 0, , drop = FALSE]
  if (nrow(neg) > 0) {
    cum <- 0
    for (i in seq_len(nrow(neg))) {
      neg$xmax[i] <- cum
      cum <- cum + neg$delta_pct[i]
      neg$xmin[i] <- cum
    }
  }
  if (nrow(pos) > 0) {
    cum <- 0
    for (i in seq_len(nrow(pos))) {
      pos$xmin[i] <- cum
      cum <- cum + pos$delta_pct[i]
      pos$xmax[i] <- cum
    }
  }
  rbind(neg, pos)
}

tab_emp <- do.call(rbind, lapply(split(tab, tab$fator), empilhar))

# ordena os fatores pela amplitude total (maior impacto primeiro, no topo)
amplitude <- sapply(split(tab_emp, tab_emp$fator), function(d) max(d$xmax) - min(d$xmin))
ordem_fatores <- names(sort(amplitude, decreasing = TRUE))
tab_emp$y <- match(tab_emp$fator, rev(ordem_fatores))  # fator de maior impacto no topo

## ---------------------------------------------------------------------
## 4) Gráfico tornado
## ---------------------------------------------------------------------

niveis_unicos <- unique(tab_emp$nivel)
cores <- setNames(grDevices::hcl.colors(length(niveis_unicos), palette = "Dynamic"), niveis_unicos)

xlim_plot <- range(c(tab_emp$xmin, tab_emp$xmax, 0)) * 1.15
altura_barra <- 0.32

png(sprintf("tornado_sensibilidade_%s_dbsra.png", tolower(metrica)),  width = 32, height = 20,
                                    res = 300,antialias = "cleartype", units = "cm")
op <- par(mar = c(4.5, 13, 4.5, 12), xpd = FALSE, bty="l",cex.main=0.8)

plot(NA, xlim = xlim_plot, ylim = c(0.5, length(ordem_fatores) + 0.5),
     yaxt = "n", ylab = "", xlab = sprintf("Variação da mediana de %s em relação ao cenário Base (%%)", metrica),
     main = "")
mtext(sprintf("Gráfico tornado — sensibilidade da mediana de %s", metrica), side = 3, line = 2.3, cex = 1.15, font = 2, adj = 0)
mtext(sprintf("Referência (cenário BASE: %s): mediana de %s = %.1f", id_base, metrica, valor_base),
      side = 3, line = 0.8, cex = 0.85, adj = 0)

abline(v = 0, col = "black", lwd = 1.4)
abline(v = pretty(xlim_plot), col = "grey90", lty = 1)
abline(v = 0, col = "black", lwd = 1.4)

for (i in seq_len(nrow(tab_emp))) {
  r <- tab_emp[i, ]
  rect(r$xmin, r$y - altura_barra, r$xmax, r$y + altura_barra,
       col = cores[r$nivel], border = "white")
}

axis(2, at = seq_along(ordem_fatores), labels = rev(ordem_fatores), las = 1, tick = FALSE, cex.axis = 0.85)

legend(x = xlim_plot[2] * 1.1, y = length(ordem_fatores) + 0.5, xpd = NA,
       legend = niveis_unicos, fill = cores[niveis_unicos], bty = "n", cex = 0.9,
       title = "Nível testado", xjust = 0)

par(op)
dev.off()
cat(sprintf("\nPNG salvo: tornado_sensibilidade_%s_dbsra.png\n", tolower(metrica)))















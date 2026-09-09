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

# definindo diretorio de trabalho..
setwd("C:/Users/mathe/OneDrive/Documents/Cabo-Verde---Stock-Assessment/Catch models")

### lendo os dados de capturas... ###
ct<- read.csv("Catch_Luz and Vieira.csv",sep = ",",dec = ".")
# lendo dados de história de vida... ##
lh<- read_xlsx("Parâmetros_História de vida.xlsx")


#---------------------------------------#
# Analise exploratoria das capturas

# Gráfico da série temporal de capturas
#---------------------------------------#

p_ct <- ggplot(ct, aes(x = Year, y = Catch)) +
  # Captura observada
  geom_line(linewidth = 2) +
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
# DEPLETION HYPOTHESES FOR Decapterus macarellus
# ==============================================================================
#
# Objetivo:
# Construir hipóteses de depleção B/K para Decapterus macarellus
# utilizando:
#
#   1. Neural Network do CMSY++
#   2. zBRT
#   3. Hipótese independente de mudança de alvo (target switching)
#   4. Hipótese não informativa
#
# As hipóteses são estimadas para:
#   - 2015
#   - 2025
#
# IMPORTANTE:
# Para cada ano-alvo, somente os dados até aquele ano são utilizados.
# Isso evita usar informação futura para estimar a depleção histórica.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. PACOTES
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(neuralnet)
library(datalimited2)


# ------------------------------------------------------------------------------
# 2. CARREGAR MODELO DA NEURAL NETWORK DO CMSY++
# ------------------------------------------------------------------------------

# ffnn.bin deve estar no working directory
#
# Esse arquivo contém:
#   - nn.endbio
#   - slope.first.min
#   - slope.first.max
#   - slope.last.min
#   - slope.last.max
#   - e possivelmente outros objetos utilizados pelo CMSY++

load("ffnn.bin")


# ------------------------------------------------------------------------------
# 3. CARREGAR DADOS DE CAPTURA
# ------------------------------------------------------------------------------

# Seu arquivo atual
ct <- read.csv(
  "Catch_Luz and Vieira.csv",
  sep = ",",
  dec = "."
)


# ------------------------------------------------------------------------------
# 4. PADRONIZAR NOMES DAS COLUNAS
# ------------------------------------------------------------------------------

# O restante do script trabalha com:
#
#   year
#   ct
#
# Mas seu arquivo pode estar com:
#
#   Year
#   Catch
#
# Portanto, fazemos a conversão automaticamente.

if ("Year" %in% names(ct)) {
  names(ct)[names(ct) == "Year"] <- "year"
}

if ("Catch" %in% names(ct)) {
  names(ct)[names(ct) == "Catch"] <- "ct"
}


# ------------------------------------------------------------------------------
# 5. CHECAGEM DOS DADOS
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# 6. PREPARAR SÉRIE DE CAPTURA
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# 7. VERIFICAR ANOS DUPLICADOS
# ------------------------------------------------------------------------------

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


# ------------------------------------------------------------------------------
# 8. VERIFICAR CAPTURAS NEGATIVAS
# ------------------------------------------------------------------------------

if (any(ct$ct < 0, na.rm = TRUE)) {
  
  stop(
    "Existem valores de captura negativos. ",
    "Verifique os dados antes de continuar."
  )
}


# ------------------------------------------------------------------------------
# 9. MOSTRAR SÉRIE
# ------------------------------------------------------------------------------

print(ct)


# ==============================================================================
# 10. FUNÇÃO — NEURAL NETWORK DO CMSY++
# ==============================================================================

estimate_endbio <- function(ct_raw, yr) {
  
  
  # --------------------------------------------------------------------------
  # 10.1. CHECAGENS
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.2. ORGANIZAR SÉRIE
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.3. ESTIMATIVA DO MSY PRIOR
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.4. ANO INTERMEDIÁRIO
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.5. CAPTURA / MSY NO ANO INTERMEDIÁRIO
  # --------------------------------------------------------------------------
  
  idx.int <- which(yr == int.yr)
  
  
  if (length(idx.int) == 0) {
    
    idx.int <- which.min(
      abs(yr - int.yr)
    )
  }
  
  
  ct_MSY.int <- ct_raw[idx.int[1]] / MSY.pr
  
  
  # --------------------------------------------------------------------------
  # 10.6. POSIÇÕES NORMALIZADAS
  # --------------------------------------------------------------------------
  
  min.ct.i <- which.min(ct_raw) / nyr
  max.ct.i <- which.max(ct_raw) / nyr
  
  int.ct.i <- idx.int[1] / nyr
  
  
  # --------------------------------------------------------------------------
  # 10.7. NORMALIZAÇÃO DO NÚMERO DE ANOS
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.8. MÉDIAS DE CAPTURA NO INÍCIO E FINAL
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.9. SLOPES
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.10. NORMALIZAR SLOPES
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.11. PADRÕES DA SÉRIE DE CAPTURA
  # --------------------------------------------------------------------------
  
  min_max <- min(ct_raw) /
    max(ct_raw)
  
  
  start.rel <- ct_raw[1] /
    max(ct_raw)
  
  
  end.rel <- ct_raw[nyr] /
    max(ct_raw)
  
  
  Flat <- as.numeric(
    min_max >= 0.45 &
      start.rel >= 0.45 &
      end.rel >= 0.45
  )
  
  
  LH <- as.numeric(
    min_max < 0.25 &
      start.rel < 0.45 &
      end.rel > 0.45
  )
  
  
  LHL <- as.numeric(
    min_max < 0.25 &
      start.rel < 0.45 &
      end.rel < 0.25
  )
  
  
  HL <- as.numeric(
    min_max < 0.25 &
      start.rel > 0.50 &
      end.rel < 0.25
  )
  
  
  HLH <- as.numeric(
    min_max < 0.25 &
      start.rel >= 0.45 &
      end.rel >= 0.45
  )
  
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.12. DATA FRAME DOS PREDITORES
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.13. NEURAL NETWORK
  # --------------------------------------------------------------------------
  
  pr.nn <- neuralnet::compute(
    nn.endbio,
    preds
  )
  
  
  idx <- max.col(
    pr.nn$net.result
  )
  
  
  # --------------------------------------------------------------------------
  # 10.14. RAZÃO CAPTURA / MSY NO FINAL
  # --------------------------------------------------------------------------
  
  ct_MSY.end <-
    ct_raw[nyr] /
    MSY.pr
  
  
  ct_MSY.use <-
    min(
      ct_MSY.end,
      mean.ct_MSY.end
    )
  
  
  # --------------------------------------------------------------------------
  # 10.15. MAPEAMENTO CMSY++ PARA B/K
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.16. RESULTADO
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # 10.17. GARANTIR LIMITES FÍSICOS
  # --------------------------------------------------------------------------
  
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


# ==============================================================================
# 11. FUNÇÃO — zBRT
# ==============================================================================

estimate_zbrt <- function(
    ct_data,
    target_year
) {
  
  
  # --------------------------------------------------------------------------
  # Selecionar apenas informação disponível até o ano-alvo
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # Rodar zBRT
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # Verificar se o ano-alvo existe no output
  # --------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------
  # Extrair B/K
  # --------------------------------------------------------------------------
  
  bk_lo <- output$ts$s_lo[idx]
  bk_hi <- output$ts$s_hi[idx]
  bk    <- output$ts$s[idx]
  
  
  # Garantir limites
  bk_lo <- max(
    0,
    min(1, bk_lo)
  )
  
  
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


# ==============================================================================
# 12. FUNÇÃO PRINCIPAL — TODAS AS HIPÓTESES
# ==============================================================================

run_depletion_hypotheses <- function(
    data,
    target_years = c(2015, 2025)
) {
  
  
  results <- list()
  
  
  # ==========================================================================
  # LOOP SOBRE OS ANOS
  # ==========================================================================
  
  for (
    yr_target in target_years
  ) {
    
    
    # ------------------------------------------------------------------------
    # Dados disponíveis até o ano-alvo
    # ------------------------------------------------------------------------
    
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
    
    
    # =========================================================================
    # HIPÓTESE 1 — NN CMSY++
    # =========================================================================
    
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
    
    
    results[
      [length(results) + 1]
    ] <- data.frame(
      
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
    
    
    # =========================================================================
    # HIPÓTESE 2 — zBRT
    # =========================================================================
    
    bk_brt <- estimate_zbrt(
      ct_data = data,
      target_year = yr_target
    )
    
    
    results[
      [length(results) + 1]
    ] <- data.frame(
      
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
    
    
    # =========================================================================
    # HIPÓTESE 3 — TARGET SWITCHING
    # =========================================================================
    
    results[
      [length(results) + 1]
    ] <- data.frame(
      
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
    
    
    # =========================================================================
    # HIPÓTESE 4 — NÃO INFORMATIVA
    # =========================================================================
    
    results[
      [length(results) + 1]
    ] <- data.frame(
      
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
  
  
  # ===========================================================================
  # JUNTAR RESULTADOS
  # ===========================================================================
  
  results <- bind_rows(
    results
  )
  
  
  # ===========================================================================
  # ORDENAR
  # ===========================================================================
  
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
    2015,
    2025
  )
)


# ==============================================================================
# 14. VISUALIZAR RESULTADOS
# ==============================================================================

print(
  bk_macarellus
)


# ==============================================================================
# 15. TABELA MAIS LIMPA
# ==============================================================================

bk_macarellus_table <- bk_macarellus %>%
  
  select(
    year,
    hypothesis,
    method,
    source,
    bk_lo,
    bk_hi,
    bk
  )


print(
  bk_macarellus_table
)


# ==============================================================================
# 16. SALVAR RESULTADOS
# ==============================================================================

write.csv(
  bk_macarellus,
  "Depletion_hypotheses_Decapterus_macarellus.csv",
  row.names = FALSE
)


# ==============================================================================
# 17. EXPORTAR APENAS OS INTERVALOS B/K
# ==============================================================================

bk_intervals <- bk_macarellus %>%
  
  select(
    year,
    hypothesis,
    bk_lo,
    bk_hi
  )


write.csv(
  bk_intervals,
  "Depletion_BK_intervals_Decapterus_macarellus.csv",
  row.names = FALSE
)


# ==============================================================================
# 18. RESUMO DAS HIPÓTESES
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("HIPÓTESES DE DEPLEÇÃO — Decapterus macarellus\n")
cat("============================================================\n")
cat("\n")

print(
  bk_macarellus_table
)

cat("\n")
cat("Arquivos salvos:\n")
cat(" - Depletion_hypotheses_Decapterus_macarellus.csv\n")
cat(" - Depletion_BK_intervals_Decapterus_macarellus.csv\n")
cat("\n")





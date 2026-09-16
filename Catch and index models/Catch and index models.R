#--------------------------------------------------------------------------------------------------------#
#     Este script contem toda analise baseada em captura e abundância do Decapterus macarellus           #
#   Analises possuem foco para o Decapterus macarellus no arquipelago de Cabo verde                      #
#              utilizando o JABBA- Just Another Bayesian Biomass Assessment (Winker et al., 2018)        #
#                Bayesian state-space surplus production model — BSSPM                                   #
#       Hipóteses de estoques consideradas h1: estoque unitário para o arquipélago                       #
#         h2: estoques unitarios para cada ilha ou ilhas barlavento e sotavento                          #
#                  PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta)                           #  
#                  Frota de cerco industrial / semi-industrial de Cabo Verde                             #
#   tentativa de padronização de CPUE, inclusão de mudança de direcionamento/mudança de alvo no tempo    #
#    ESFORÇO NÃO É ESPECÍFICO DA ESPÉCIE. O denominador são dias de pesca de uma frota multiespecífica;  #
#    parte deles não foi gasta procurando cavala. Efeito: hiperdepleção aparente depois de 2014.         #
#              Desmembramento de CPUE entre as espécies do D. macarellus e Auxis spp.                    #  
# utilização de modelos regressivos (Generalized Linear models - GLM) e Generalized additive models (GAM)#
#   Distribuições candidatas (resposta = toneladas, contínua com zeros):                                 #
#   D1  Tweedie (compound Poisson-Gamma, 1 < p < 2), ligação log                                         #  
#   D2  Hurdle-Gamma (ziGamma + ziformula), ligação log                                                  #
#   D3  Delta-lognormal clássico (binomial + lognormal nas positivas)                                    #  
#   D4  Gamma só nas positivas      — comparação: ignora os zeros                                        #  
#   D5  Gaussiana em log(CPUE + c)  — comparação: a má prática                                           #
#                                                                                                        #  
#         Seleção: backward a partir do modelo cheio, LRT para modelos aninhados                         #  
#         (glmmTMB ajusta por ML, então o LRT dos efeitos fixos é válido),                               #
#         com AIC e BIC como critérios complementares. Distribuição escolhida                            #  
#         por AIC + resíduos simulados (DHARMa), não por AIC sozinho.                                    #      
#         Índice: médias marginais do fator `ano` via emmeans, com offset = 0                            #  
#         (taxa por dia), retrotransformadas com correção de viés, série                                 #  
#         normalizada para média 1 (convenção de entrada do JABBA) e CV anual                            #
#         derivado do erro-padrão na escala do preditor linear.                                          #  
#                                                                                                         #
#      # Codificacao criada por Silva, MLS ; 09/09/2026, Instituto do Mar- IMar, Mindelo                 #
#--------------------------------------------------------------------------------------------------------#

# limpando ambiente de trabalho...
rm(list = ls())

#@pacotes..
#install.packages("readxl")
library(readxl)
#installed.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("writexl")
library(writexl)
#install.packages("mgcv")
library(mgcv)
#-----------------------------
#instalando o JABBA (Just Another Bayesian Biomass Assessment (Winker et al., 2018) )
#install.packages("pak")
library(pak)
#pak::pak("jabbamodel/JABBA")
library(JABBA)
#------------------------------------

# definindo diretorio de trabalho..
setwd("C:/Users/mathe/OneDrive/Documents/Cabo-Verde---Stock-Assessment/Catch and index models")

### lendo os dados de capturas... ###
ct<- read.csv("Catch_Luz and Vieira.csv",sep = ",",dec = ".")
# lendo dados de história de vida... ##
lh<- read_xlsx("Parametros_Historia_de_vida.xlsx")




## =====================================================================
## A) ESTRUTURA DE DADOS QUE ESTE SCRIPT PRODUZ
## =====================================================================
## `viagens` — uma linha por VIAGEM de cada barco:
##   id_viagem, barco, ilha, ano, mes, data_saida, data_volta, dias,
##   cap_macarellus, cap_auxis, cap_katsuwonus, cap_euthynnus, cap_outras
##
## `esforco_estrato` — dias de pesca agregados por ilha x mes x ano,
##   que é o formato em que o IMar disse ter o esforço. Serve para o
##   cenário em que NÃO se consegue o esforço por viagem (ver parte 03).
##
## `verdade` — a trajetória de biomassa usada para gerar os dados.
##   EXISTE SÓ NA SIMULAÇÃO. É o gabarito: permite checar se a correção
##   recupera o sinal, que é exatamente o teste que Hinton & Maunder
##   (2003) recomendam antes de confiar num método de padronização.
## =====================================================================

set.seed(20260916)

anos  <- 1995:2023
ilhas <- c("Sao Vicente", "Santiago", "Sal", "Sao Nicolau", "Maio")
sp    <- c("macarellus", "auxis", "katsuwonus", "euthynnus", "outras")

## ---------------------------------------------------------------------
## 1) A VERDADE OCULTA — biomassa relativa por espécie
## ---------------------------------------------------------------------
# Cavala: declínio moderado e contínuo (cai ~45% em 29 anos). NÃO tem
# degrau em 2014. Todo degrau que aparecer na CPUE nominal é artefato.
B_mac <- seq(1.25, 0.68, length.out = length(anos)) *
  exp(cumsum(rnorm(length(anos), 0, 0.055)))
# Auxis: aumento real moderado (parte do "boom" das capturas é esforço
# redirecionado, parte é disponibilidade que de fato subiu)
B_aux <- seq(0.50, 1.60, length.out = length(anos)) *
  exp(cumsum(rnorm(length(anos), 0, 0.040)))
B_kat <- rep(1, length(anos)) * exp(cumsum(rnorm(length(anos), 0, 0.05)))
B_eut <- seq(1.15, 0.80, length.out = length(anos)) *
  exp(cumsum(rnorm(length(anos), 0, 0.05)))

B <- data.frame(ano = anos, macarellus = B_mac, auxis = B_aux,
                katsuwonus = B_kat, euthynnus = B_eut, outras = 1)
for (s in c("macarellus", "auxis", "katsuwonus", "euthynnus")) {
  B[[s]] <- B[[s]] / mean(B[[s]])       # média 1, como o índice do JABBA
}
verdade <- B

# "anos bons e anos ruins": variação de DISPONIBILIDADE (mexe no q, não em
# B). É o que produz os picos e quedas da série de captura sem que o
# estoque mude — e é também o ruído que nenhum índice de CPUE consegue
# separar da abundância. Mantido modesto de propósito.
disp <- data.frame(ano = anos)
for (s_ in c("macarellus", "auxis", "katsuwonus", "euthynnus", "outras"))
  disp[[s_]] <- exp(rnorm(length(anos), 0, 0.10))

## ---------------------------------------------------------------------
## 2) FROTA — barcos entram e saem ao longo da série (problema P5)
## ---------------------------------------------------------------------
n_barcos <- 26
frota <- data.frame(
  barco   = sprintf("CV-%03d", seq_len(n_barcos)),
  entrada = sample(anos[1:20], n_barcos, replace = TRUE),
  poder   = exp(rnorm(n_barcos, 0, 0.28)),          # poder de pesca (P4)
  base    = sample(ilhas, n_barcos, replace = TRUE,
                   prob = c(.34, .30, .16, .12, .08)),
  stringsAsFactors = FALSE
)
frota$saida <- pmin(frota$entrada + rpois(n_barcos, 17) + 5, max(anos) + 1)

## ---------------------------------------------------------------------
## 3) ESFORÇO ANUAL — forma da série real (Tabela 3 da FAO, 2026):
##    baixo nos anos 90, crescimento forte 2005-2010, queda depois
## ---------------------------------------------------------------------
esforco_alvo <- approx(
  x = c(1995, 2000, 2004, 2007, 2010, 2013, 2016, 2019, 2023),
  y = c(1500, 2100, 2700, 4900, 7200, 6300, 5400, 3500, 2700),
  xout = anos)$y
esforco_alvo <- esforco_alvo * exp(rnorm(length(anos), 0, 0.07))

## ---------------------------------------------------------------------
## 4) TROCA DE ALVO — transição logística suave centrada em 2014
##    (mais realista que um degrau: a frota migrou ao longo de ~3 anos)
## ---------------------------------------------------------------------
p_cavala <- 0.14 + (0.78 - 0.14) / (1 + exp((anos - 2013.6) / 0.9))
p_auxis  <- 0.62 - (0.62 - 0.12) / (1 + exp((anos - 2013.6) / 0.9))
p_mista  <- pmax(0, 1 - p_cavala - p_auxis)
prob_tatica <- data.frame(ano = anos, cavala = p_cavala,
                          auxis = p_auxis, mista = p_mista)

## ---------------------------------------------------------------------
## 5) CAPTURABILIDADE — quanto cada tática, ilha e mês mexe no q
## ---------------------------------------------------------------------
# q base por espécie (t por dia de pesca, numa viagem "média")
q0 <- c(macarellus = 0.24, auxis = 0.30, katsuwonus = 0.16,
        euthynnus = 0.07, outras = 0.09)

# multiplicador de capturabilidade por tática (linhas) x espécie (colunas)
q_tatica <- rbind(
  cavala = c(macarellus = 1.90, auxis = 0.20, katsuwonus = 0.55, euthynnus = 0.70, outras = 0.80),
  auxis  = c(macarellus = 0.42, auxis = 4.60, katsuwonus = 1.25, euthynnus = 0.80, outras = 0.90),
  mista  = c(macarellus = 1.00, auxis = 1.00, katsuwonus = 1.00, euthynnus = 1.00, outras = 1.00)
)

# densidade relativa por ilha (P6: as ilhas não são equivalentes)
dens_ilha <- rbind(
  `Sao Vicente` = c(macarellus = 1.35, auxis = 0.80, katsuwonus = 1.10, euthynnus = 1.00, outras = 1),
  Santiago      = c(macarellus = 0.95, auxis = 1.30, katsuwonus = 1.05, euthynnus = 1.10, outras = 1),
  Sal           = c(macarellus = 0.80, auxis = 1.25, katsuwonus = 1.20, euthynnus = 0.85, outras = 1),
  `Sao Nicolau` = c(macarellus = 1.20, auxis = 0.90, katsuwonus = 0.85, euthynnus = 0.95, outras = 1),
  Maio          = c(macarellus = 0.70, auxis = 1.10, katsuwonus = 0.75, euthynnus = 1.05, outras = 1)
)

# sazonalidade: pico da cavala no 1º semestre, do Auxis no 2º
saz <- function(mes, pico, amp) 1 + amp * cos(2 * pi * (mes - pico) / 12)
saz_sp <- list(macarellus = function(m) saz(m, 3, 0.35),
               auxis      = function(m) saz(m, 9, 0.40),
               katsuwonus = function(m) saz(m, 7, 0.25),
               euthynnus  = function(m) saz(m, 5, 0.20),
               outras     = function(m) rep(1, length(m)))

## ---------------------------------------------------------------------
## 6) GERAÇÃO DAS VIAGENS
## ---------------------------------------------------------------------
# Captura de cada espécie gerada como Poisson composta com Gamma (cardumes
# encontrados x massa por cardume): produz zeros de forma natural, sem
# precisar de regra ad hoc. Um erro lognormal por viagem é somado por
# cima para que a distribuição geradora NÃO seja exatamente nenhuma das
# candidatas testadas na parte 03 — do contrário a comparação de
# distribuições seria viciada a favor da Tweedie.
gera_captura <- function(mu, peso_cardume = 0.55) {
  n <- length(mu)
  lambda <- mu / peso_cardume
  n_card <- rpois(n, lambda)
  out <- numeric(n)
  pos <- n_card > 0
  if (any(pos)) {
    out[pos] <- vapply(which(pos), function(i)
      sum(rgamma(n_card[i], shape = 2.2, scale = peso_cardume / 2.2)),
      numeric(1))
  }
  out * exp(rnorm(n, -0.5 * 0.22^2, 0.22))
}

linhas <- vector("list", 0)
id <- 0L

for (a in seq_along(anos)) {
  ano_i  <- anos[a]
  ativos <- frota[frota$entrada <= ano_i & frota$saida > ano_i, ]
  if (nrow(ativos) == 0) next
  
  # reparte o esforço do ano entre meses (sazonalidade da atividade) e barcos
  peso_mes <- 1 + 0.30 * cos(2 * pi * (1:12 - 4) / 12)
  dias_mes <- esforco_alvo[a] * peso_mes / sum(peso_mes)
  
  for (m in 1:12) {
    # nº de viagens no mês: dias do mês / duração média, repartido entre barcos
    n_viag <- max(1, rpois(1, dias_mes[m] / 6.8))
    if (n_viag == 0) next
    
    b_idx  <- sample(seq_len(nrow(ativos)), n_viag, replace = TRUE)
    dias   <- pmin(pmax(rpois(n_viag, 4.8) + 2, 2), 13)
    
    # ilha da viagem: normalmente a base do barco, às vezes outra
    ilha <- ifelse(runif(n_viag) < 0.78, ativos$base[b_idx],
                   sample(ilhas, n_viag, replace = TRUE))
    
    # tática da viagem (a variável latente que a parte 02 tenta recuperar)
    pt <- as.numeric(prob_tatica[a, c("cavala", "auxis", "mista")])
    tat <- sample(c("cavala", "auxis", "mista"), n_viag, replace = TRUE, prob = pt)
    
    dia_saida  <- sample(1:26, n_viag, replace = TRUE)
    data_saida <- as.Date(sprintf("%d-%02d-%02d", ano_i, m, dia_saida))
    
    cap <- matrix(0, nrow = n_viag, ncol = length(sp),
                  dimnames = list(NULL, sp))
    for (s in sp) {
      mu <- q0[[s]] * ativos$poder[b_idx] * q_tatica[tat, s] *
        dens_ilha[ilha, s] * saz_sp[[s]](m) * B[[s]][a] * disp[[s]][a] * dias
      cap[, s] <- gera_captura(mu)
    }
    
    id <- id + n_viag
    linhas[[length(linhas) + 1]] <- data.frame(
      id_viagem  = (id - n_viag + 1L):id,
      barco      = ativos$barco[b_idx],
      ilha       = ilha,
      ano        = ano_i,
      mes        = m,
      data_saida = data_saida,
      data_volta = data_saida + dias,
      dias       = dias,
      tatica_real = tat,          # NÃO existe nos dados reais — é gabarito
      cap_macarellus = round(cap[, "macarellus"], 3),
      cap_auxis      = round(cap[, "auxis"], 3),
      cap_katsuwonus = round(cap[, "katsuwonus"], 3),
      cap_euthynnus  = round(cap[, "euthynnus"], 3),
      cap_outras     = round(cap[, "outras"], 3),
      stringsAsFactors = FALSE
    )
  }
}

viagens <- do.call(rbind, linhas)
rownames(viagens) <- NULL

## ---------------------------------------------------------------------
## 7) ESFORÇO AGREGADO POR ESTRATO — o formato que o IMar disse ter
## ---------------------------------------------------------------------
esforco_estrato <- aggregate(dias ~ ilha + ano + mes, data = viagens, FUN = sum)
names(esforco_estrato)[names(esforco_estrato) == "dias"] <- "dias_estrato"
esforco_estrato <- esforco_estrato[order(esforco_estrato$ano,
                                         esforco_estrato$mes,
                                         esforco_estrato$ilha), ]
rownames(esforco_estrato) <- NULL

## ---------------------------------------------------------------------
## 8) CONFERÊNCIA
## ---------------------------------------------------------------------
cat("\n===== DADOS SIMULADOS (formato IMar) =====\n")
cat(sprintf("Viagens: %d | Barcos: %d | Anos: %d-%d | Ilhas: %d\n",
            nrow(viagens), length(unique(viagens$barco)),
            min(viagens$ano), max(viagens$ano), length(unique(viagens$ilha))))
cat(sprintf("Esforco total: %.0f dias | Duracao media da viagem: %.1f dias\n",
            sum(viagens$dias), mean(viagens$dias)))
cat(sprintf("Viagens com ZERO de cavala: %.1f%%  (zeros de direcionamento, P3)\n",
            100 * mean(viagens$cap_macarellus == 0)))

res_ano <- aggregate(cbind(cap_macarellus, cap_auxis, dias) ~ ano,
                     data = viagens, FUN = sum)
cat("\nCaptura (t) e esforco (dias) por ano — primeiras e ultimas linhas:\n")
print(round(head(res_ano, 4), 0))
cat("   ...\n")
print(round(tail(res_ano, 6), 0))

cat("\nEstrutura de `viagens`:\n")
str(viagens)

cat("\nOBS: as colunas `tatica_real` e o objeto `verdade` sao GABARITO da\n")
cat("simulacao e NAO existirao nos dados reais do IMar. Servem so para\n")
cat("checar se a padronizacao recupera o sinal que sabemos ser verdadeiro.\n")




#======================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde
#----------------------------------------------------------------------
# Autor : Silva, MLS
# Data  : 2026-09-16
# Parte : 02 de 03 — análise exploratória e INFERÊNCIA DA TÁTICA DE PESCA
#         (a covariável de alvo, que é o coração da correção)
#
# Pré-requisito: objeto `viagens` no ambiente (rodar a parte 01, ou
# carregar os dados reais do IMar no mesmo formato de colunas).
#
# Este script NÃO depende de nenhum pacote além do R base — de propósito,
# para que a parte exploratória rode em qualquer máquina.
#
# Produz:
#   viagens        — acrescido de: cpue_mac, prop_*, alvo, PC1, PC2
#   comp_bm        — matriz de composição por barco-mês (insumo do cluster)
#   esforco_dirigido — dias alocados por tática (cenário H3)
#   6 figuras exploratórias em PNG
#======================================================================

stopifnot(exists("viagens"))

cols_cap <- c("cap_macarellus", "cap_auxis", "cap_katsuwonus",
              "cap_euthynnus", "cap_outras")
sp_nomes <- c("D. macarellus", "Auxis spp.", "K. pelamis",
              "E. alletteratus", "Outras")
stopifnot(all(cols_cap %in% names(viagens)))

## cores da casa
cor_sp   <- hcl.colors(5, palette = "Dark 3")
COR_MAC  <- "#1F4E79"; COR_AUX <- "#C0501B"; COR_NEU <- "#7F7F7F"

## =====================================================================
## 1) VARIÁVEIS DERIVADAS
## =====================================================================
viagens$cap_total <- rowSums(viagens[, cols_cap])
viagens$cpue_mac  <- viagens$cap_macarellus / viagens$dias   # t/dia
viagens$pos_mac   <- as.integer(viagens$cap_macarellus > 0)

# composição proporcional da viagem (a "impressão digital" da tática)
for (i in seq_along(cols_cap)) {
  p <- ifelse(viagens$cap_total > 0,
              viagens[[cols_cap[i]]] / viagens$cap_total, 0)
  viagens[[sub("cap_", "prop_", cols_cap[i])]] <- p
}
cols_prop <- sub("cap_", "prop_", cols_cap)

# fatores (é assim que as variáveis entram no modelo da parte 03: como
# FATORES, para depois serem marginalizadas e sobrar o efeito de ano)
viagens$fano   <- factor(viagens$ano)
viagens$fmes   <- factor(viagens$mes)
viagens$filha  <- factor(viagens$ilha)
viagens$fbarco <- factor(viagens$barco)
viagens$barco_mes <- paste(viagens$barco, viagens$ano, viagens$mes, sep = "_")

## =====================================================================
## 2) FILTROS — documentar SEMPRE o que cada um descarta
##    (regras na linha de Hoyle et al. 2015/2018 e Sant'Ana et al. 2020)
## =====================================================================
n0 <- nrow(viagens)
filtro_log <- data.frame(regra = character(), removidas = integer(),
                         restantes = integer(), stringsAsFactors = FALSE)
reg <- function(regra, antes) {
  filtro_log[nrow(filtro_log) + 1, ] <<-
    list(regra, antes - nrow(viagens), nrow(viagens)); invisible(NULL)
}

a <- nrow(viagens); viagens <- viagens[viagens$dias > 0, ];        reg("esforço > 0", a)
a <- nrow(viagens); viagens <- viagens[viagens$cap_total > 0, ];   reg("captura total > 0 (viagem com alguma pesca)", a)

# barcos com amostra mínima (senão o efeito aleatório de barco não estima)
tb <- table(viagens$barco)
a <- nrow(viagens); viagens <- viagens[viagens$barco %in% names(tb[tb >= 30]), ]
reg("barcos com >= 30 viagens", a)

# estratos ano x ilha com amostra mínima
te <- table(paste(viagens$ano, viagens$ilha))
a <- nrow(viagens)
viagens <- viagens[paste(viagens$ano, viagens$ilha) %in% names(te[te >= 5]), ]
reg("estratos ano x ilha com >= 5 viagens", a)

viagens$fano <- droplevels(factor(viagens$ano))
viagens$filha <- droplevels(factor(viagens$ilha))
viagens$fbarco <- droplevels(factor(viagens$barco))

cat("\n========= FILTROS APLICADOS =========\n"); print(filtro_log)
cat(sprintf("Retidas %d de %d viagens (%.1f%%)\n",
            nrow(viagens), n0, 100 * nrow(viagens) / n0))

## =====================================================================
## 3) SÉRIES ANUAIS E CPUE NOMINAL
##    CPUE nominal = captura da espécie / esforço TOTAL da frota.
##    É o cenário S1 — e é o que a FAO (2026) usou.
## =====================================================================
ser <- aggregate(cbind(cap_macarellus, cap_auxis, cap_katsuwonus,
                       cap_euthynnus, cap_outras, cap_total, dias) ~ ano,
                 data = viagens, FUN = sum)
ser$cpue_nom_mac <- ser$cap_macarellus / ser$dias
ser$cpue_nom_aux <- ser$cap_auxis / ser$dias
ser$prop_zero    <- tapply(viagens$pos_mac == 0, viagens$ano, mean)[as.character(ser$ano)]
ser$n_viagens    <- as.numeric(table(viagens$ano)[as.character(ser$ano)])

cat("\n===== SÉRIE ANUAL (captura t, esforço dias, CPUE nominal t/dia) =====\n")
print(data.frame(ano = ser$ano,
                 cap_mac = round(ser$cap_macarellus, 0),
                 cap_aux = round(ser$cap_auxis, 0),
                 dias    = round(ser$dias, 0),
                 cpue_mac = round(ser$cpue_nom_mac, 3),
                 zeros_mac = sprintf("%.0f%%", 100 * ser$prop_zero)),
      row.names = FALSE)

## =====================================================================
## 4) FIGURAS EXPLORATÓRIAS
## =====================================================================
eixo_anos <- function(x) axis(1, at = pretty(x), labels = pretty(x))

## --- Fig 1: capturas por espécie e esforço -----------------------------
png("exp1_capturas_esforco.png", width = 26, height = 13, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)

matplot(ser$ano, ser[, paste0("cap_", c("macarellus", "auxis", "katsuwonus",
                                        "euthynnus", "outras"))],
        type = "l", lty = 1, lwd = 2.2, col = cor_sp,
        xlab = "Ano", ylab = "Captura (t)",
        main = "A. Captura anual por espécie")
abline(v = 2014, col = "grey70", lty = 2)
legend("topleft", sp_nomes, col = cor_sp, lwd = 2.2, bty = "n", cex = 0.72)

plot(ser$ano, ser$dias, type = "l", lwd = 2.4, col = COR_NEU,
     xlab = "Ano", ylab = "Esforço (dias de pesca)",
     main = "B. Esforço total da frota", ylim = c(0, max(ser$dias) * 1.05))
abline(v = 2014, col = "grey70", lty = 2)
text(2014, max(ser$dias) * 1.02, " troca de alvo", adj = 0, cex = 0.7, col = "grey40")
par(op); dev.off()
cat("\nPNG salvo: exp1_capturas_esforco.png\n")

## --- Fig 2: composição proporcional da captura -------------------------
comp_ano <- as.matrix(ser[, paste0("cap_", c("macarellus", "auxis", "katsuwonus",
                                             "euthynnus", "outras"))])
comp_ano <- comp_ano / rowSums(comp_ano)

png("exp2_composicao_anual.png", width = 24, height = 12, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mar = c(4.2, 4.6, 3, 8.5), bty = "l", cex.main = 0.95, cex = 0.85, xpd = FALSE)
acum <- t(apply(comp_ano, 1, cumsum))
plot(NA, xlim = range(ser$ano), ylim = c(0, 1), xlab = "Ano",
     ylab = "Proporção da captura total",
     main = "Composição da captura da frota — a troca de alvo fica visível")
for (k in ncol(acum):1) {
  polygon(c(ser$ano, rev(ser$ano)), c(acum[, k], rep(0, nrow(acum))),
          col = cor_sp[k], border = NA)
}
abline(v = 2014, col = "white", lty = 2, lwd = 1.6)
par(xpd = TRUE)
legend(max(ser$ano) + 0.6, 0.85, sp_nomes, fill = cor_sp, border = NA,
       bty = "n", cex = 0.75)
par(op); dev.off()
cat("PNG salvo: exp2_composicao_anual.png\n")

## --- Fig 3: CPUE nominal e zeros ---------------------------------------
png("exp3_cpue_nominal.png", width = 26, height = 13, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)

plot(ser$ano, ser$cpue_nom_mac, type = "l", lwd = 2.4, col = COR_MAC,
     xlab = "Ano", ylab = "CPUE nominal (t/dia)",
     main = "A. CPUE nominal — cavala e Auxis",
     ylim = c(0, max(c(ser$cpue_nom_mac, ser$cpue_nom_aux)) * 1.05))
lines(ser$ano, ser$cpue_nom_aux, lwd = 2.4, col = COR_AUX)
abline(v = 2014, col = "grey70", lty = 2)
legend("topleft", c("D. macarellus", "Auxis spp."),
       col = c(COR_MAC, COR_AUX), lwd = 2.4, bty = "n", cex = 0.78)

plot(ser$ano, 100 * ser$prop_zero, type = "l", lwd = 2.4, col = COR_MAC,
     xlab = "Ano", ylab = "% de viagens sem cavala na captura",
     main = "B. Zeros de direcionamento (problema P3)",
     ylim = c(0, max(100 * ser$prop_zero) * 1.1))
abline(v = 2014, col = "grey70", lty = 2)
par(op); dev.off()
cat("PNG salvo: exp3_cpue_nominal.png\n")

## --- Fig 4: sazonalidade e efeito de ilha ------------------------------
png("exp4_mes_ilha.png", width = 26, height = 12, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(5.5, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
boxplot(cpue_mac ~ mes, data = viagens, outline = FALSE, col = "#8FAADC",
        xlab = "Mês", ylab = "CPUE da cavala (t/dia)", lwd = 1,
        main = "A. Sazonalidade (justifica o termo `mes`)")
boxplot(cpue_mac ~ ilha, data = viagens, outline = FALSE, col = "#74C476",
        xlab = "", ylab = "CPUE da cavala (t/dia)", lwd = 1, xaxt = "n",
        main = "B. Efeito de ilha (justifica o termo `ilha`)")
ilhas_u <- levels(factor(viagens$ilha))
axis(1, at = seq_along(ilhas_u), labels = FALSE)
text(seq_along(ilhas_u), par("usr")[3], labels = ilhas_u, srt = 30,
     adj = 1, xpd = NA, cex = 0.78)
par(op); dev.off()
cat("PNG salvo: exp4_mes_ilha.png\n")

## =====================================================================
## 5) INFERÊNCIA DA TÁTICA DE PESCA (a covariável de alvo)
## ---------------------------------------------------------------------
## Nenhum banco registra "o que o mestre queria pescar". A composição da
## captura é a melhor impressão digital disponível dessa intenção.
## Seguimos Winker et al. (2013) e Sant'Ana et al. (2020):
##   (i)  agregar a composição por BARCO-MÊS, não por viagem — dados de
##        viagem carregam ruído de encontro fortuito com cardume, o que
##        gera má-alocação de táticas (Hoyle et al. 2018);
##   (ii) transformar as proporções por raiz quadrada, para que espécies
##        menos abundantes também contribuam para a similaridade;
##   (iii) PCA sobre essa matriz;
##   (iv) DUAS representações concorrentes da tática (hipótese H2):
##        - discreta  : cluster k-means / Ward  -> fator `alvo`
##        - contínua  : escores PC1, PC2 ("DPC")
## =====================================================================

## 5.1 matriz de composição por barco-mês
agg <- aggregate(viagens[, cols_cap], by = list(barco_mes = viagens$barco_mes),
                 FUN = sum)
tot <- rowSums(agg[, cols_cap])
comp_bm <- as.matrix(agg[, cols_cap] / tot)
rownames(comp_bm) <- agg$barco_mes
colnames(comp_bm) <- sub("cap_", "", cols_cap)
comp_bm <- comp_bm[is.finite(rowSums(comp_bm)), , drop = FALSE]
comp_sqrt <- sqrt(comp_bm)            # transformação de Winker et al. (2013)

cat(sprintf("\nMatriz de composição: %d barcos-mês x %d espécies\n",
            nrow(comp_sqrt), ncol(comp_sqrt)))

## 5.2 PCA
pca <- prcomp(comp_sqrt, center = TRUE, scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cat("Variância explicada pelos eixos (%): ",
    paste(sprintf("PC%d=%.1f", 1:4, var_exp[1:4]), collapse = "  "), "\n")

## 5.3 número de grupos: silhueta média (implementada em R base para não
##     exigir o pacote `cluster`; calculada em subamostra por velocidade)
silhueta_media <- function(X, cl, n_sub = 1200) {
  idx <- if (nrow(X) > n_sub) sample(nrow(X), n_sub) else seq_len(nrow(X))
  Xs <- X[idx, , drop = FALSE]; cs <- cl[idx]
  D <- as.matrix(dist(Xs))
  gr <- unique(cs)
  s <- vapply(seq_along(cs), function(i) {
    mesmo <- cs == cs[i]; mesmo[i] <- FALSE
    if (!any(mesmo)) return(0)
    ai <- mean(D[i, mesmo])
    bi <- min(vapply(setdiff(gr, cs[i]),
                     function(g) mean(D[i, cs == g]), numeric(1)))
    (bi - ai) / max(ai, bi)
  }, numeric(1))
  mean(s)
}

# Se o analista quiser fixar o número de táticas (por conhecimento da
# pescaria, ou para comparar com a literatura — Winker et al. 2013 e
# Sant'Ana et al. 2020 acharam 4), basta preencher `k_forcado`.
k_forcado <- NA          # NA = escolher pela silhueta

esc <- pca$x[, 1:3, drop = FALSE]
set.seed(42)
ks <- 2:6
sil <- vapply(ks, function(k)
  silhueta_media(esc, kmeans(esc, centers = k, nstart = 25, iter.max = 50)$cluster),
  numeric(1))
k_otimo <- if (is.na(k_forcado)) ks[which.max(sil)] else k_forcado
cat("Silhueta média por k: ",
    paste(sprintf("k=%d: %.3f", ks, sil), collapse = "  "), "\n")
cat(sprintf("k escolhido: %d%s\n", k_otimo,
            if (is.na(k_forcado)) " (pela silhueta)" else " (fixado pelo analista)"))
if (is.na(k_forcado) && k_otimo == 2)
  cat("[NOTA] A silhueta tende a favorecer k pequeno. Se o painel B da Fig. 5\n",
      "      mostrar uma nuvem CONTÍNUA cortada ao meio (e não grupos separados),\n",
      "      isso é evidência de que a tática é um GRADIENTE, não uma classe --\n",
      "      exatamente o argumento de Winker et al. (2013) para preferir os\n",
      "      escores contínuos (PC1/PC2) ao cluster. A parte 03 testa as duas.\n")

## 5.4 k-means com o k escolhido  +  Ward como método concorrente
km <- kmeans(esc, centers = k_otimo, nstart = 50, iter.max = 100)
sub  <- sample(nrow(esc), min(2500, nrow(esc)))
ward <- cutree(hclust(dist(esc[sub, ]), method = "ward.D2"), k = k_otimo)
tb_kw <- table(kmeans = km$cluster[sub], ward = ward)
conc_kw <- sum(apply(tb_kw, 1, max)) / sum(tb_kw)
cat(sprintf("Concordância k-means x Ward (subamostra de %d barcos-mês): %.1f%%\n",
            length(sub), 100 * conc_kw))
if (conc_kw < 0.80)
  cat("[NOTA] Concordância baixa entre os dois algoritmos: o agrupamento não é\n",
      "      estável, mais um motivo para tratar a tática como gradiente contínuo.\n")

## 5.5 nomear os clusters pela espécie dominante no centróide
cent <- t(vapply(seq_len(k_otimo), function(g)
  colMeans(comp_bm[km$cluster == g, , drop = FALSE]), numeric(ncol(comp_bm))))
colnames(cent) <- colnames(comp_bm)
dom <- colnames(cent)[apply(cent, 1, which.max)]
nome_cl <- make.unique(dom, sep = "_")
cat("\nComposição média de cada grupo (proporção da captura):\n")
print(round(cbind(cent, n = as.numeric(table(km$cluster))), 3))
cat("Rótulos atribuídos: ", paste(nome_cl, collapse = ", "), "\n")

## 5.6 devolver o rótulo e os escores contínuos a cada VIAGEM
mapa <- data.frame(barco_mes = rownames(comp_bm),
                   alvo = factor(nome_cl[km$cluster], levels = unique(nome_cl)),
                   PC1  = pca$x[, 1], PC2 = pca$x[, 2],
                   stringsAsFactors = FALSE)
viagens <- merge(viagens, mapa, by = "barco_mes", all.x = TRUE, sort = FALSE)
viagens <- viagens[!is.na(viagens$alvo), ]
viagens$alvo <- droplevels(viagens$alvo)

cat(sprintf("\nViagens com tática atribuída: %d (%.1f%%)\n",
            nrow(viagens), 100 * nrow(viagens) / n0))
print(table(viagens$alvo))

## 5.7 VALIDAÇÃO — só possível porque estes dados são simulados.
##     Com os dados reais do IMar este bloco não roda (não há gabarito);
##     a validação passa a ser indireta (composição média dos grupos faz
##     sentido biológico? o gráfico temporal reproduz a mudança conhecida?)
if ("tatica_real" %in% names(viagens)) {
  tb <- table(real = viagens$tatica_real, inferida = viagens$alvo)
  acerto <- sum(apply(tb, 1, max)) / sum(tb)
  cat("\n--- VALIDAÇÃO CONTRA O GABARITO DA SIMULAÇÃO ---\n")
  print(tb)
  cat(sprintf("Concordância (melhor pareamento por linha): %.1f%%\n", 100 * acerto))
}

## --- Fig 5: PCA, silhueta e composição dos grupos ----------------------
cor_cl <- hcl.colors(k_otimo, palette = "Dark 3")
png("exp5_clusters_alvo.png", width = 27, height = 10.5, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 3), mar = c(4.4, 4.4, 3, 1), oma = c(0, 0, 0, 4.6),
          bty = "l", cex.main = 0.95, cex = 0.85)

plot(ks, sil, type = "b", pch = 19, lwd = 2, col = COR_MAC,
     xlab = "Número de grupos (k)", ylab = "Silhueta média",
     main = "A. Escolha de k")
points(k_otimo, max(sil), pch = 21, bg = COR_AUX, cex = 1.8)

plot(pca$x[, 1], pca$x[, 2], col = adjustcolor(cor_cl[km$cluster], 0.55),
     pch = 16, cex = 0.5, xlab = sprintf("PC1 (%.0f%%)", var_exp[1]),
     ylab = sprintf("PC2 (%.0f%%)", var_exp[2]),
     main = "B. Táticas no espaço de composição")
points(km$centers[, 1], km$centers[, 2], pch = 21, bg = cor_cl, cex = 1.9, lwd = 1.5)
legend("topright", nome_cl, col = cor_cl, pch = 16, bty = "n", cex = 0.72)

bp <- barplot(t(cent), beside = FALSE, col = cor_sp, border = NA,
              names.arg = nome_cl, las = 2, cex.names = 0.7,
              ylab = "Proporção média da captura",
              main = "C. Composição média de cada tática")
legend(max(bp) + 0.75, 1, rev(sp_nomes), fill = rev(cor_sp), border = NA,
       bty = "n", cex = 0.7, xpd = NA)
par(op); dev.off()
cat("PNG salvo: exp5_clusters_alvo.png\n")

## =====================================================================
## 6) A TROCA DE ALVO AO LONGO DO TEMPO  +  ESFORÇO DIRIGIDO (H3)
## =====================================================================
tat_ano <- prop.table(table(viagens$ano, viagens$alvo), margin = 1)

# esforço alocado a cada tática dentro do estrato (ano x mês x ilha),
# proporcional aos DIAS das viagens de cada tática — é isto que converte
# "esforço total" em "esforço dirigido" (hipótese H3)
esforco_dirigido <- aggregate(dias ~ ano + mes + ilha + alvo,
                              data = viagens, FUN = sum)
names(esforco_dirigido)[names(esforco_dirigido) == "dias"] <- "dias_dirigidos"

dias_alvo_ano <- aggregate(dias_dirigidos ~ ano + alvo,
                           data = esforco_dirigido, FUN = sum)

png("exp6_troca_de_alvo.png", width = 26, height = 12, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)

matplot(as.numeric(rownames(tat_ano)), tat_ano, type = "l", lty = 1, lwd = 2.4,
        col = cor_cl, ylim = c(0, 1), xlab = "Ano",
        ylab = "Proporção das viagens", main = "A. Tática de pesca ao longo do tempo")
abline(v = 2014, col = "grey70", lty = 2)
legend("topleft", nome_cl, col = cor_cl, lwd = 2.4, bty = "n", cex = 0.72)

alvos <- levels(viagens$alvo)
mat_d <- sapply(alvos, function(g) {
  x <- dias_alvo_ano[dias_alvo_ano$alvo == g, ]
  d <- setNames(rep(0, length(ser$ano)), ser$ano)
  d[as.character(x$ano)] <- x$dias_dirigidos; d
})
matplot(ser$ano, mat_d, type = "l", lty = 1, lwd = 2.4, col = cor_cl,
        ylim = c(0, max(c(mat_d, ser$dias)) * 1.05),
        xlab = "Ano", ylab = "Dias de pesca dirigidos",
        main = "B. Esforço dirigido por tática (cenário H3)")
lines(ser$ano, ser$dias, lwd = 2, lty = 2, col = COR_NEU)
abline(v = 2014, col = "grey70", lty = 2)
legend("topleft", c(alvos, "esforço total"), col = c(cor_cl, COR_NEU),
       lwd = 2.2, lty = c(rep(1, length(alvos)), 2), bty = "n", cex = 0.72)
par(op); dev.off()
cat("PNG salvo: exp6_troca_de_alvo.png\n")

cat("\n===== PARTE 02 CONCLUÍDA =====\n")
cat("Objetos prontos para a parte 03: `viagens` (com alvo, PC1, PC2),\n")
cat("`ser` (séries anuais + CPUE nominal), `esforco_dirigido`.\n")




#======================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde
#----------------------------------------------------------------------
# Autor : Silva, MLS
# Data  : 2026-09-16
# Parte : 03 de 03 — modelagem, seleção, diagnóstico e extração dos
#         índices de abundância relativa para entrar no JABBA
#
# Pré-requisito: rodar as partes 01 e 02 (objetos `viagens`, `ser`,
# `esforco_dirigido` no ambiente).
#
# Pacotes: glmmTMB, emmeans, DHARMa (opcionalmente writexl).
#   install.packages(c("glmmTMB","emmeans","DHARMa","writexl"))
#
# AVISO DE TEMPO: com ~16 mil viagens, cada GLMM Tweedie com efeito
# aleatório leva de 1 a 5 minutos. O script guarda os ajustes em `fits`.
#======================================================================
#
# ---------------------------------------------------------------------
# DECISÕES METODOLÓGICAS DESTE SCRIPT (e a justificativa de cada uma)
# ---------------------------------------------------------------------
# (1) TODAS as covariáveis de estrato entram como FATORES (ano, mês,
#     ilha, alvo). Isso permite marginalizá-las depois — a média marginal
#     do fator `ano`, com as demais mediadas, É o índice. O fator `ano`
#     nunca entra em seleção e nunca entra em interação: interação com o
#     ano invalida o efeito de ano como índice (Hinton & Maunder, 2003).
#
# (2) O esforço entra como OFFSET, log(dias), com coeficiente fixado em 1.
#     Modelar captura com offset é equivalente a modelar a taxa, mas
#     preserva a estrutura de erro da captura observada — dividir antes e
#     modelar a razão trata o denominador como se fosse constante conhecida.
#
# (3) SELEÇÃO: backward a partir do modelo cheio, com teste da razão de
#     verossimilhança (LRT) para modelos aninhados. glmmTMB ajusta por ML
#     por padrão, então o LRT dos efeitos fixos é válido. AIC e BIC são
#     reportados como critérios complementares, não como juiz único.
#     Backward em vez de forward porque termo cujo efeito só aparece
#     depois de ajustar outro é sistematicamente perdido no forward.
#     Lembrar do dilema de Hinton & Maunder: covariável de menos deixa a
#     variação de q vazar para o ano; covariável demais come o sinal de
#     abundância. Os dois enviesam, em direções opostas.
#
# (4) DISTRIBUIÇÕES: a resposta é captura em toneladas — contínua,
#     positiva e com massa em zero. Testamos cinco candidatas (D1-D5) e
#     escolhemos por AIC **e** resíduos simulados (DHARMa), nunca por AIC
#     sozinho: AIC compara ajuste, não compara adequação distribucional.
#
# (5) INCERTEZA: o CV anual sai do erro-padrão do efeito de ano na escala
#     do preditor linear. Nos modelos delta, o CV do componente binomial
#     costuma ser excluído por causa de separação perfeita (Sant'Ana et
#     al., 2020) — aqui reportamos as duas versões e dizemos qual foi
#     usada, porque a versão sem o binomial SUBESTIMA a incerteza.
#
# (6) O índice final é normalizado para MÉDIA 1. É a convenção de entrada
#     do JABBA (e do SS): o modelo estima q livremente, então só a forma
#     relativa da série importa.
#
# (7) RESSALVA QUE ACOMPANHA O RESULTADO: sem variáveis operacionais
#     (potência, porão, tempo de busca), isto NÃO é uma padronização
#     completa de CPUE. É uma correção parcial para troca de alvo. As
#     séries nominal e corrigida entram no JABBA como CENÁRIOS, não como
#     "a" estimativa de abundância.
#======================================================================

stopifnot(exists("viagens"), exists("ser"))
suppressPackageStartupMessages({
  library(glmmTMB); library(emmeans); library(DHARMa)
})
tem_writexl <- requireNamespace("writexl", quietly = TRUE)

COR_MAC <- "#1F4E79"; COR_AUX <- "#C0501B"; COR_NEU <- "#7F7F7F"
COR_S2  <- "#2E8B57"; COR_S3 <- "#7030A0"

## variável-resposta e checagens
viagens$captura <- viagens$cap_macarellus
viagens$ldias   <- log(viagens$dias)
stopifnot(all(is.finite(viagens$captura)), all(viagens$dias > 0))

cat("\n===== DADOS PARA A MODELAGEM =====\n")
cat(sprintf("Viagens: %d | Anos: %d | Ilhas: %d | Barcos: %d | Táticas: %d\n",
            nrow(viagens), nlevels(viagens$fano), nlevels(viagens$filha),
            nlevels(viagens$fbarco), nlevels(viagens$alvo)))
cat(sprintf("Proporção de zeros na resposta: %.1f%%\n",
            100 * mean(viagens$captura == 0)))

## =====================================================================
## 1) ESTRUTURAS CANDIDATAS (E0-E5) — a Tweedie é a distribuição de
##    trabalho nesta etapa; a escolha da distribuição vem depois, na
##    melhor estrutura (assim não se compara tudo contra tudo às cegas)
## =====================================================================
# Construtor explícito de fórmulas. Usar update(f, . ~ . - offset(ldias))
# NÃO remove o offset (ele não é um `term.label`), e reformulate() quebra o
# termo aleatório `(1|fbarco)`. Montar a string é o caminho seguro.
monta_formula <- function(resposta, termos, aleatorio = TRUE, offset = TRUE) {
  rhs <- paste(c("fano", termos,
                 if (aleatorio) "(1 | fbarco)",
                 if (offset) "offset(ldias)"), collapse = " + ")
  stats::as.formula(paste(resposta, "~", rhs))
}

estruturas <- list(
  E0 = captura ~ fano + offset(ldias),
  E1 = captura ~ fano + fmes + offset(ldias),
  E2 = captura ~ fano + fmes + filha + offset(ldias),
  E3 = captura ~ fano + fmes + filha + alvo + offset(ldias),
  E4 = captura ~ fano + fmes + filha + alvo + (1 | fbarco) + offset(ldias),
  E5 = captura ~ fano + fmes + filha + PC1 + PC2 + (1 | fbarco) + offset(ldias)
)

cat("\n===== 1) AJUSTE DAS ESTRUTURAS (família Tweedie) =====\n")
fits <- list()
for (nm in names(estruturas)) {
  t0 <- Sys.time()
  fits[[nm]] <- try(glmmTMB(estruturas[[nm]], family = tweedie(link = "log"),
                            data = viagens), silent = TRUE)
  ok <- !inherits(fits[[nm]], "try-error")
  cat(sprintf("  %s: %s  (%.1f min)\n", nm,
              if (ok) sprintf("AIC = %.1f", AIC(fits[[nm]])) else "FALHOU",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
fits <- fits[!vapply(fits, inherits, logical(1), "try-error")]

tab_est <- data.frame(
  modelo  = names(fits),
  df      = vapply(fits, function(m) attr(logLik(m), "df"), numeric(1)),
  logLik  = vapply(fits, function(m) as.numeric(logLik(m)), numeric(1)),
  AIC     = vapply(fits, AIC, numeric(1)),
  BIC     = vapply(fits, BIC, numeric(1)), row.names = NULL)
tab_est$dAIC <- round(tab_est$AIC - min(tab_est$AIC), 1)
cat("\nComparação das estruturas:\n")
print(transform(tab_est, logLik = round(logLik, 1), AIC = round(AIC, 1),
                BIC = round(BIC, 1)), row.names = FALSE)

## LRT dos termos aninhados (E0 < E1 < E2 < E3 < E4)
cat("\nLRT sequencial (cada termo contra o modelo sem ele):\n")
seq_lrt <- c("E0", "E1", "E2", "E3", "E4")
seq_lrt <- seq_lrt[seq_lrt %in% names(fits)]
for (i in seq_len(length(seq_lrt) - 1)) {
  a <- anova(fits[[seq_lrt[i]]], fits[[seq_lrt[i + 1]]])
  cat(sprintf("  %s vs %s : Chisq = %.1f, df = %d, p = %s\n",
              seq_lrt[i], seq_lrt[i + 1], a$Chisq[2], a$`Chi Df`[2],
              format.pval(a$`Pr(>Chisq)`[2], digits = 3, eps = 1e-16)))
}

## =====================================================================
## 2) SELEÇÃO BACKWARD COM LRT (o fator `ano` é intocável)
## =====================================================================
cat("\n===== 2) SELEÇÃO BACKWARD (LRT, alfa = 0.01) =====\n")
termos_sel <- c("fmes", "filha", "alvo")   # `fano` NUNCA entra na seleção
m_atual <- fits$E4
repeat {
  if (length(termos_sel) == 0) break
  p_vals <- setNames(numeric(length(termos_sel)), termos_sel)
  for (tm in termos_sel) {
    f_red <- monta_formula("captura", setdiff(termos_sel, tm))
    m_red <- try(glmmTMB(f_red, family = tweedie(link = "log"), data = viagens),
                 silent = TRUE)
    p_vals[tm] <- if (inherits(m_red, "try-error")) 0 else
      anova(m_red, m_atual)[["Pr(>Chisq)"]][2]
  }
  cat("  p-valores da remoção: ",
      paste(sprintf("%s=%.3g", names(p_vals), p_vals), collapse = "  "), "\n")
  pior <- names(which.max(p_vals))
  if (max(p_vals) > 0.01) {
    cat(sprintf("  -> removendo `%s` (p = %.3g)\n", pior, max(p_vals)))
    termos_sel <- setdiff(termos_sel, pior)
    m_atual <- glmmTMB(monta_formula("captura", termos_sel),
                       family = tweedie(link = "log"), data = viagens)
  } else {
    cat("  -> todos os termos significativos; estrutura final mantida\n"); break
  }
}
formula_atual <- monta_formula("captura", termos_sel)
cat("Estrutura selecionada: "); print(formula_atual)
m_estrutura <- m_atual

## =====================================================================
## 3) DISTRIBUIÇÕES CANDIDATAS (D1-D5) NA ESTRUTURA SELECIONADA
## =====================================================================
cat("\n===== 3) DISTRIBUIÇÕES CANDIDATAS =====\n")
f_fix <- formula_atual
dist_fits <- list()

# D1 — Tweedie (compound Poisson-Gamma): contínua com massa em zero
dist_fits$D1_tweedie <- m_estrutura

# D2 — Hurdle-Gamma: um modelo binomial para a presença + Gamma para a
#      magnitude, ajustados conjuntamente. Em glmmTMB, `ziGamma` com
#      `ziformula` é literalmente um hurdle (a Gamma não tem massa em 0).
zi_f <- ~ fano + alvo
dist_fits$D2_hurdle_gamma <- try(
  glmmTMB(f_fix, ziformula = zi_f, family = ziGamma(link = "log"), data = viagens),
  silent = TRUE)

# D3 — Delta-lognormal clássico (o padrão nos grupos do ICCAT): um
#      binomial para a PRESENÇA e um lognormal para a MAGNITUDE das
#      positivas; o índice é o produto dos dois efeitos de ano,
#      (1 - w) x (y | y > 0), como em Sant'Ana et al. (2020).
viagens_pos <- viagens[viagens$captura > 0, ]
viagens_pos$lcpue <- log(viagens_pos$captura / viagens_pos$dias)
f_bin <- monta_formula("pos_mac", termos_sel, offset = FALSE)
f_pos <- monta_formula("lcpue",   termos_sel, offset = FALSE)
m_bin <- try(glmmTMB(f_bin, family = binomial(), data = viagens), silent = TRUE)
m_pos <- try(glmmTMB(f_pos, family = gaussian(), data = viagens_pos), silent = TRUE)
if (!inherits(m_pos, "try-error")) dist_fits$D3_lognormal_positivas <- m_pos

# D4 — Gamma só nas positivas: ignora os zeros de propósito, para
#      QUANTIFICAR o viés que essa prática comum introduz
dist_fits$D4_gamma_positivas <- try(
  glmmTMB(f_fix, family = Gamma(link = "log"), data = viagens_pos), silent = TRUE)

# D5 — Gaussiana em log(CPUE + c): a má prática, incluída como referência.
#      A resposta JÁ é uma taxa, então este modelo vai SEM offset — é
#      justamente o erro que se comete ao manter o offset por inércia.
c_add <- 0.5 * min(viagens$captura[viagens$captura > 0])
viagens$lcpue_c <- log((viagens$captura + c_add) / viagens$dias)
f_d5 <- monta_formula("lcpue_c", termos_sel, offset = FALSE)
dist_fits$D5_log_mais_c <- try(
  glmmTMB(f_d5, family = gaussian(), data = viagens), silent = TRUE)

dist_fits <- dist_fits[!vapply(dist_fits, inherits, logical(1), "try-error")]

cat("\nAIC/BIC das distribuições (D4 e D5 NÃO são comparáveis por AIC:\n",
    "usam dados/escalas diferentes — estão aqui para comparação de ÍNDICE):\n")
print(data.frame(
  modelo = names(dist_fits),
  n      = vapply(dist_fits, function(m) nrow(model.frame(m)), numeric(1)),
  AIC    = round(vapply(dist_fits, AIC, numeric(1)), 1),
  BIC    = round(vapply(dist_fits, BIC, numeric(1)), 1),
  row.names = NULL))

## Escolha da distribuição: só D1 e D2 são comparáveis por AIC (mesmos
## dados, mesma escala da resposta). D3/D4/D5 usam subconjunto ou escala
## diferentes e entram só como comparação de ÍNDICE.
comparaveis <- intersect(c("D1_tweedie", "D2_hurdle_gamma"), names(dist_fits))
melhor_dist <- comparaveis[which.min(vapply(dist_fits[comparaveis], AIC, numeric(1)))]
m_final <- dist_fits[[melhor_dist]]
cat(sprintf("\nDistribuição escolhida por AIC (confirmar no diagnóstico): %s\n",
            melhor_dist))
if (melhor_dist == "D2_hurdle_gamma")
  cat("[ATENÇÃO] Em modelo hurdle, o emmeans devolve apenas o componente\n",
      "         condicional. Para o índice use o cenário delta (S2d), que\n",
      "         multiplica presença x magnitude explicitamente.\n")

## =====================================================================
## 4) DIAGNÓSTICO DE RESÍDUOS (DHARMa: resíduos simulados, escalonados)
##    Para GLMM não-gaussiano, resíduos de Pearson enganam — os testes
##    formais devem ser feitos sobre resíduos simulados.
## =====================================================================
cat("\n===== 4) DIAGNÓSTICO DE RESÍDUOS =====\n")
diagnostica <- function(m, nome, n_sim = 250) {
  set.seed(1)
  r <- simulateResiduals(m, n = n_sim, plot = FALSE)
  u  <- testUniformity(r, plot = FALSE)       # KS: a distribuição está certa?
  d  <- testDispersion(r, plot = FALSE)       # sobre/subdispersão
  o  <- testOutliers(r, plot = FALSE)         # excesso de extremos
  q  <- try(testQuantiles(r, plot = FALSE), silent = TRUE)  # homocedasticidade
  z  <- try(testZeroInflation(r, plot = FALSE), silent = TRUE)
  data.frame(modelo = nome,
             KS_p   = signif(u$p.value, 3),
             disp_p = signif(d$p.value, 3),
             disp_ratio = signif(as.numeric(d$statistic), 3),
             outlier_p = signif(o$p.value, 3),
             quantis_p = if (inherits(q, "try-error")) NA else signif(q$p.value, 3),
             zeros_p = if (inherits(z, "try-error")) NA else signif(z$p.value, 3),
             row.names = NULL)
}
diag_tab <- do.call(rbind, lapply(names(dist_fits), function(nm)
  tryCatch(diagnostica(dist_fits[[nm]], nm),
           error = function(e) data.frame(modelo = nm, KS_p = NA, disp_p = NA,
                                          disp_ratio = NA, outlier_p = NA,
                                          quantis_p = NA, zeros_p = NA))))
print(diag_tab, row.names = FALSE)
cat("\nLeitura: p > 0.05 em KS, dispersão e quantis = sem evidência contra o\n")
cat("modelo. Com n grande, testes formais rejeitam por desvios triviais —\n")
cat("olhe TAMBÉM os gráficos abaixo, não só os p-valores.\n")

png("diag_residuos_modelo_final.png", width = 26, height = 13, res = 300,
    antialias = "cleartype", units = "cm")
set.seed(1); r_fin <- simulateResiduals(m_final, n = 250, plot = FALSE)
plot(r_fin)      # o próprio DHARMa monta o painel QQ + resíduo vs predito
dev.off()
cat("PNG salvo: diag_residuos_modelo_final.png\n")

# Testes clássicos, aplicáveis apenas ao componente lognormal (D5/delta):
if ("D5_log_mais_c" %in% names(dist_fits)) {
  res5 <- residuals(dist_fits$D5_log_mais_c)
  amostra <- sample(res5, min(5000, length(res5)))
  cat(sprintf("\nShapiro-Wilk nos resíduos do modelo em log (n=%d): p = %.3g\n",
              length(amostra), shapiro.test(amostra)$p.value))
  cat("(normalidade é premissa do componente lognormal, não da Tweedie)\n")
}

## =====================================================================
## 5) EXTRAÇÃO DO ÍNDICE — médias marginais do fator `ano` (emmeans)
## ---------------------------------------------------------------------
## `offset = 0` é ESSENCIAL: sem isso o emmeans prediz na média do
## log(dias) e a série sai em toneladas por viagem-média, não por dia.
## `weights`: "equal" dá o mesmo peso a cada estrato (padrão, e o que se
## quer num índice); "proportional" pondera pelo n observado, o que
## devolve parte do desbalanceamento que estamos tentando remover.
## =====================================================================
extrai_indice <- function(m, nome, lognormal = FALSE, pesos = "equal",
                          usa_offset = TRUE) {
  args <- list(object = m, specs = ~ fano, weights = pesos)
  if (usa_offset) args$offset <- 0     # sem isto a série sai por viagem-média
  em <- do.call(emmeans, args)         # emmeans devolve na escala do link
  s_link <- as.data.frame(summary(em))
  mu <- s_link$emmean; se <- s_link$SE
  # Retrotransformação. Em GLM com ligação log (Tweedie, Gamma) a média
  # marginal já está na escala do preditor: basta exp(). Em modelo
  # gaussiano ajustado sobre log(y), a média da resposta exige a correção
  # de viés lognormal exp(sigma^2/2) — são coisas diferentes, e trocar
  # uma pela outra é um erro comum.
  idx <- if (lognormal) exp(mu + 0.5 * sigma(m)^2) else exp(mu)
  data.frame(ano = as.numeric(as.character(s_link$fano)),
             cenario = nome, indice_bruto = idx, se_log = se,
             cv = sqrt(exp(se^2) - 1), row.names = NULL)
}

# Índice delta (binomial x lognormal), com a variância combinada pelo
# método delta: em escala log, d log(p)/d logit(p) = 1 - p.
indice_delta <- function(m_bin, m_pos, nome, pesos = "equal") {
  sb <- as.data.frame(summary(emmeans(m_bin, ~ fano, weights = pesos)))
  sp <- as.data.frame(summary(emmeans(m_pos, ~ fano, weights = pesos)))
  p  <- plogis(sb$emmean)
  se_logp <- (1 - p) * sb$SE
  idx <- p * exp(sp$emmean + 0.5 * sigma(m_pos)^2)
  se_tot <- sqrt(se_logp^2 + sp$SE^2)
  data.frame(ano = as.numeric(as.character(sb$fano)), cenario = nome,
             indice_bruto = idx, se_log = se_tot,
             cv = sqrt(exp(se_tot^2) - 1),
             cv_sem_binomial = sqrt(exp(sp$SE^2) - 1),
             row.names = NULL)
}
normaliza <- function(d) { d$indice <- d$indice_bruto / mean(d$indice_bruto); d }

cat("\n===== 5) ÍNDICES POR CENÁRIO =====\n")

## S1 — CPUE NOMINAL (captura total da espécie / esforço total da frota)
S1 <- data.frame(ano = ser$ano, cenario = "S1 nominal",
                 indice_bruto = ser$cpue_nom_mac,
                 se_log = NA_real_, cv = NA_real_)
# CV empírico da nominal: dispersão entre viagens dentro do ano
cv_emp <- tapply(viagens$cpue_mac, viagens$ano,
                 function(x) sd(x) / (mean(x) * sqrt(length(x))))
S1$cv <- as.numeric(cv_emp[as.character(S1$ano)])
S1 <- normaliza(S1)

## S2 — CPUE CORRIGIDA (modelo selecionado, com a covariável de alvo)
S2 <- normaliza(extrai_indice(m_final, "S2 corrigida (alvo discreto)"))

## S2d — mesma correção pelo caminho DELTA-LOGNORMAL (padrão ICCAT):
##       índice = P(captura > 0) x média das positivas. Serve para checar
##       se a escolha da distribuição muda a forma da série — se mudar,
##       isso vira incerteza estrutural a reportar.
S2d <- NULL
if (!inherits(m_bin, "try-error") && !inherits(m_pos, "try-error")) {
  S2d <- normaliza(indice_delta(m_bin, m_pos, "S2d delta-lognormal"))
  cat(sprintf("S2d: CV médio COM o binomial = %.3f | SEM = %.3f",
              mean(S2d$cv), mean(S2d$cv_sem_binomial)))
  cat("  <- a versão sem o binomial subestima a incerteza\n")
  S2d$cv_sem_binomial <- NULL
}

## S2b — alvo contínuo (hipótese H2: DPC de Winker et al. 2013)
S2b <- NULL
if ("E5" %in% names(fits)) {
  S2b <- normaliza(extrai_indice(fits$E5, "S2b corrigida (alvo contínuo PC1/PC2)"))
}

## S0 — sem covariável de alvo (isola QUANTO o alvo desloca o índice;
##      é o "influence plot" de Bentley et al. feito à mão)
S0 <- normaliza(extrai_indice(fits$E2, "S0 sem alvo"))

## S3 — esforço dirigido (hipótese H3): só viagens da tática da cavala,
##      com o esforço dirigido no denominador. SENSIBILIDADE, não cenário
##      principal: subsetting tende a gerar hiperestabilidade.
alvo_cavala <- grep("macarellus", levels(viagens$alvo), value = TRUE)[1]
S3 <- NULL
if (!is.na(alvo_cavala)) {
  v3 <- viagens[viagens$alvo == alvo_cavala, ]
  if (nlevels(droplevels(v3$fano)) >= 0.7 * nlevels(viagens$fano)) {
    v3$fano <- droplevels(v3$fano); v3$filha <- droplevels(v3$filha)
    v3$fbarco <- droplevels(v3$fbarco)
    m3 <- try(glmmTMB(monta_formula("captura", setdiff(termos_sel, "alvo")),
                      family = tweedie(link = "log"), data = v3), silent = TRUE)
    if (!inherits(m3, "try-error")) S3 <- normaliza(extrai_indice(m3, "S3 esforço dirigido"))
  } else {
    cat("[AVISO] A tática da cavala não cobre anos suficientes depois de 2014;\n")
    cat("        o cenário S3 fica sem estimativa em parte da série — é\n")
    cat("        exatamente a limitação prevista para o subsetting.\n")
  }
}

indices <- do.call(rbind, Filter(Negate(is.null), list(S1, S0, S2, S2d, S2b, S3)))

cat("\nÍndices (média 1) e CV:\n")
print(transform(indices[, c("ano", "cenario", "indice", "cv")],
                indice = round(indice, 3), cv = round(cv, 3)), row.names = FALSE)

## =====================================================================
## 6) TESTE DAS HIPÓTESES
## =====================================================================
cat("\n===== 6) TESTE DAS HIPÓTESES =====\n")
razao <- function(d, ini, fim) mean(d$indice[d$ano >= ini & d$ano <= fim])
queda <- function(d) razao(d, 2015, max(d$ano)) / razao(d, 2005, 2013)

qn <- queda(S1); qc <- queda(S2)
cat(sprintf("H1 — queda pós-2014 (média 2015+ / média 2005-2013):\n"))
cat(sprintf("     CPUE nominal   : %.2f  (queda de %.0f%%)\n", qn, 100 * (1 - qn)))
cat(sprintf("     CPUE corrigida : %.2f  (queda de %.0f%%)\n", qc, 100 * (1 - qc)))
cat(sprintf("     Diferença: %.0f pontos percentuais de queda atribuíveis\n",
            100 * (qc - qn)))
cat(sprintf("     à realocação de esforço  ->  H1 %s\n",
            if (qc > qn * 1.10) "SUSTENTADA" else "NÃO sustentada"))

corr_s0s2 <- cor(S0$indice, S2$indice)
cat(sprintf("\nH0 — correlação entre índice com e sem a covariável de alvo: %.3f\n",
            corr_s0s2))
cat(sprintf("     %s\n", if (corr_s0s2 > 0.98)
  "Praticamente idênticos: a padronização NÃO está corrigindo nada (reportar!)" else
    "A covariável de alvo desloca o índice de forma relevante"))

if (!is.null(S2b)) {
  cat(sprintf("\nH2 — alvo discreto vs contínuo: correlação %.3f | dAIC = %.1f\n",
              cor(S2$indice, S2b$indice), AIC(fits$E5) - AIC(m_final)))
  cat("     (dAIC negativo favorece o alvo contínuo, como em Winker et al. 2013)\n")
}

## Comparação com a VERDADE — só existe na simulação
if (exists("verdade")) {
  v <- verdade[verdade$ano %in% S2$ano, c("ano", "macarellus")]
  v$verdadeiro <- v$macarellus / mean(v$macarellus)
  mare <- function(d) {
    m <- merge(d, v, by = "ano")
    median(abs(m$indice - m$verdadeiro) / m$verdadeiro)
  }
  cat("\n--- CHECAGEM CONTRA O GABARITO (só possível na simulação) ---\n")
  cat("Erro absoluto relativo mediano (MARE) de cada cenário:\n")
  for (cen in unique(indices$cenario)) {
    d <- indices[indices$cenario == cen, ]
    cat(sprintf("  %-38s MARE = %.3f | r = %.3f\n", cen, mare(d),
                cor(merge(d, v, by = "ano")$indice, merge(d, v, by = "ano")$verdadeiro)))
  }
}

## =====================================================================
## 7) FIGURAS DOS ÍNDICES
## =====================================================================
cores_cen <- c("S1 nominal" = COR_AUX, "S0 sem alvo" = COR_NEU,
               "S2 corrigida (alvo discreto)" = COR_MAC,
               "S2d delta-lognormal" = "#00A0B0",
               "S2b corrigida (alvo contínuo PC1/PC2)" = COR_S2,
               "S3 esforço dirigido" = COR_S3)

png("indices_cpue_cenarios.png", width = 26, height = 14, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)

cens <- unique(indices$cenario)
ylim <- c(0, max(indices$indice, na.rm = TRUE) * 1.12)
plot(NA, xlim = range(indices$ano), ylim = ylim, xlab = "Ano",
     ylab = "Índice de abundância relativa (média = 1)",
     main = "A. Cenários de índice para o JABBA")
abline(v = 2014, col = "grey80", lty = 2)
for (cen in cens) {
  d <- indices[indices$cenario == cen, ]
  lines(d$ano, d$indice, lwd = 2.4, col = cores_cen[cen])
}
if (exists("verdade")) {
  lines(v$ano, v$verdadeiro, lwd = 2.2, lty = 3, col = "black")
}
legend("topright", c(cens, if (exists("verdade")) "biomassa verdadeira"),
       col = c(cores_cen[cens], if (exists("verdade")) "black"),
       lwd = 2.3, lty = c(rep(1, length(cens)), if (exists("verdade")) 3),
       bty = "n", cex = 0.68)

# painel B: índice principal com faixa de incerteza
d2 <- S2
lo <- d2$indice * exp(-1.96 * d2$se_log); hi <- d2$indice * exp(1.96 * d2$se_log)
plot(d2$ano, d2$indice, type = "n", ylim = c(0, max(hi, na.rm = TRUE) * 1.05),
     xlab = "Ano", ylab = "Índice (média = 1)",
     main = "B. Índice corrigido com IC 95%")
polygon(c(d2$ano, rev(d2$ano)), c(lo, rev(hi)),
        col = adjustcolor(COR_MAC, 0.18), border = NA)
lines(d2$ano, d2$indice, lwd = 2.6, col = COR_MAC)
lines(S1$ano, S1$indice, lwd = 2.0, col = COR_AUX, lty = 2)
abline(v = 2014, col = "grey80", lty = 2)
legend("topright", c("corrigida (IC 95%)", "nominal"),
       col = c(COR_MAC, COR_AUX), lwd = c(2.6, 2), lty = c(1, 2),
       bty = "n", cex = 0.72)
par(op); dev.off()
cat("\nPNG salvo: indices_cpue_cenarios.png\n")

## Influence plot artesanal: quanto cada termo desloca o efeito de ano
png("influencia_covariaveis.png", width = 24, height = 12, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mar = c(4.4, 4.6, 3, 1), bty = "l", cex.main = 0.95, cex = 0.85)
passos <- intersect(c("E0", "E1", "E2", "E3", "E4"), names(fits))
cor_p <- hcl.colors(length(passos), palette = "Dark 3")
plot(NA, xlim = range(S2$ano), ylim = c(0, 2.6), xlab = "Ano",
     ylab = "Índice (média = 1)",
     main = "Deslocamento do índice à medida que cada termo entra")
abline(v = 2014, col = "grey80", lty = 2); abline(h = 1, col = "grey85")
for (i in seq_along(passos)) {
  di <- normaliza(extrai_indice(fits[[passos[i]]], passos[i]))
  lines(di$ano, di$indice, lwd = 2.2, col = cor_p[i])
}
legend("topright", passos, col = cor_p, lwd = 2.2, bty = "n", cex = 0.72)
par(op); dev.off()
cat("PNG salvo: influencia_covariaveis.png\n")

## =====================================================================
## 8) EXPORTAÇÃO PARA O JABBA
##    O JABBA espera uma coluna de ano e uma coluna por índice, mais um
##    objeto de mesmo formato com os CVs (ou SE). Média 1 por convenção.
## =====================================================================
jabba_idx <- data.frame(Yr = S1$ano)
jabba_se  <- data.frame(Yr = S1$ano)
for (cen in c("S1 nominal", "S2 corrigida (alvo discreto)")) {
  d <- indices[indices$cenario == cen, ]
  nm <- if (grepl("nominal", cen)) "cpue_nominal" else "cpue_corrigida"
  jabba_idx[[nm]] <- d$indice[match(jabba_idx$Yr, d$ano)]
  jabba_se[[nm]]  <- d$cv[match(jabba_se$Yr, d$ano)]
}
# piso de CV: CV estatístico do modelo subestima a incerteza real (não
# inclui erro de processo nem o componente binomial); 0.15-0.25 é o piso
# usual em avaliações com índice fishery-dependent
jabba_se[, -1] <- lapply(jabba_se[, -1, drop = FALSE],
                         function(x) pmax(x, 0.20, na.rm = TRUE))

write.csv(jabba_idx, "jabba_indices_macarellus.csv", row.names = FALSE)
write.csv(jabba_se,  "jabba_cv_macarellus.csv", row.names = FALSE)
write.csv(indices,   "indices_todos_cenarios.csv", row.names = FALSE)
if (tem_writexl) {
  writexl::write_xlsx(list(indices = jabba_idx, cv = jabba_se,
                           todos = indices, estruturas = tab_est,
                           diagnostico = diag_tab),
                      path = "padronizacao_cpue_macarellus.xlsx")
  cat("XLSX salvo: padronizacao_cpue_macarellus.xlsx\n")
}
cat("CSV salvos: jabba_indices_macarellus.csv, jabba_cv_macarellus.csv,\n")
cat("            indices_todos_cenarios.csv\n")

cat("\n===== COMO REPORTAR =====\n")
cat("1. Apresentar S1 e S2 como CENÁRIOS ALTERNATIVOS de entrada no JABBA,\n")
cat("   nunca uma como 'a' série correta (Hoyle et al. 2024, boas práticas\n")
cat("   17 e 18: não emendar índices conflitantes num só).\n")
cat("2. Reportar a tabela de seleção, os diagnósticos e o gráfico de\n")
cat("   influência — inclusive se os pressupostos foram violados.\n")
cat("3. Declarar que, sem variáveis operacionais, isto é uma CORREÇÃO\n")
cat("   PARCIAL PARA TROCA DE ALVO, não uma padronização completa.\n")
cat("4. O CV exportado tem piso de 0,20: o CV do modelo é de processo\n")
cat("   estatístico e subestima a incerteza da abundância.\n")















































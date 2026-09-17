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
#install.packages("DHARMa")
library(DHARMa)
#install.packages("emmeans")
library(emmeans)
#install.packages("glmmTMB")
library(glmmTMB)

#------------------------------------

# definindo diretorio de trabalho..
setwd("C:/Users/mathe/OneDrive/Documents/Cabo-Verde---Stock-Assessment/Catch and index models")

### lendo os dados de capturas... ###
ct<- read.csv("Catch_Luz and Vieira.csv",sep = ",",dec = ".")
# lendo dados de história de vida... ##
lh<- read_xlsx("Parametros_Historia_de_vida.xlsx")



#======================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde
# Frota industrial / semi-industrial — dados do IMar
#----------------------------------------------------------------------
# Autor : Silva, MLS
# Data  : 2026-09-16
# Parte : 01 de 03 — LEITURA, LIMPEZA E MONTAGEM DA TABELA DE VIAGENS
#         a partir do arquivo INDUSTRIAL_<ano>.csv entregue pelo IMar
#======================================================================
#
# ---------------------------------------------------------------------
# O QUE ESTE SCRIPT FAZ, EM UMA FRASE
# ---------------------------------------------------------------------
# O arquivo do IMar vem em formato LONGO: uma linha por ESPÉCIE capturada
# em cada evento de amostragem. Este script converte isso para o formato
# LARGO que a análise precisa: uma linha por VIAGEM, com uma coluna de
# captura por espécie, mais as variáveis operacionais daquela viagem.
#
# ---------------------------------------------------------------------
# AVISO ESTRUTURAL SOBRE O ANO DE 2024 SOZINHO  (leia antes de tudo)
# ---------------------------------------------------------------------
# O índice de abundância de uma padronização de CPUE É o efeito de ANO.
# Com um único ano de dados esse efeito não existe — não há o que estimar.
# Portanto, com o arquivo de 2024 apenas, este pipeline NÃO produz índice
# de abundância. O que ele produz, e que já é muito útil, é:
#   (a) a verificação de que a estrutura do arquivo do IMar suporta toda
#       a análise planejada (e quais campos servem para quê);
#   (b) a caracterização da frota e das táticas de pesca em 2024;
#   (c) todos os modelos ajustados e diagnosticados, usando MÊS como
#       fator temporal no lugar de ano — o que valida a máquina;
#   (d) a lista objetiva do que ainda falta pedir ao IMar.
# Quando a série completa chegar, o objeto `fator_tempo` (parte 02) muda
# sozinho de `fmes` para `fano` e o restante do código não muda.
#
# ---------------------------------------------------------------------
# DECISÕES DE LEITURA E LIMPEZA (cada uma comentada no código abaixo)
# ---------------------------------------------------------------------
# L1. Encoding: o CSV vem em latin1 (acentos de "S. ANTÃO" etc.). Lemos
#     com readLines + iconv para UTF-8 em vez de `fileEncoding`, porque
#     `fileEncoding` falha dependendo do locale da máquina.
# L2. `quote = ""`: o arquivo não usa aspas; deixar o padrão faz o R
#     interpretar apóstrofos como aspas e truncar a leitura.
# L3. `read.csv2`: separador ";" e decimal com VÍRGULA ("76,71").
# L4. Uma linha por espécie: `Amostragem` é o identificador da viagem
#     (verificado: cada Amostragem tem 1 embarcação e 1 data de partida).
# L5. FILTRAMOS PARA REDE DE CERCO. Das 120 viagens que capturaram cavala
#     preta em 2024, 118 foram de cerco. Manter linha-de-mão e covos no
#     denominador acrescenta esforço que nunca teve chance de pegar
#     cavala — é o problema P1 (esforço não específico) ampliado de graça.
#     A arte é, além disso, o fator de capturabilidade mais óbvio.
# L6. `Profundidade_pesca_engenho1 == 0` é AUSÊNCIA DE REGISTRO, não
#     "zero metros" (31% das viagens). Vira NA.
# L7. Local de pesca: usamos `Nome_banco_pesca` (onde se pescou), não o
#     porto/ilha de desembarque (onde se descarregou). São coisas
#     diferentes e só a primeira é uma covariável de densidade.
# L8. `Tipo_embarcacao` tem I (industrial), A e R. R aparece em 4 viagens
#     de 1 barco — agrupado em "outro" para não criar nível instável.
#======================================================================

## =====================================================================
## 0) PARÂMETROS QUE VOCÊ PODE QUERER MEXER
## =====================================================================
ARQUIVO      <- "INDUSTRIAL_2024.csv"   # ou vetor de arquivos, um por ano
ARTE_ALVO    <- "REDE DE CERCO"          # ver decisão L5
MIN_VIAG_BANCO <- 10                     # bancos com menos viagens viram "OUTROS"
ESPECIE_FOCO <- "DECAPTERUS MACARELLUS"  # nome científico, como no arquivo

## =====================================================================
## 1) LEITURA ROBUSTA (decisões L1-L3)
## =====================================================================
le_imar <- function(caminho) {
  linhas <- readLines(caminho, warn = FALSE)
  # L1: latin1 -> UTF-8. `sub="?"` evita erro fatal se houver byte inválido.
  linhas <- iconv(linhas, from = "latin1", to = "UTF-8", sub = "?")
  # L2 e L3: sem aspas, separador ";", decimal ","
  d <- read.csv2(text = linhas, quote = "", stringsAsFactors = FALSE,
                 strip.white = TRUE)
  names(d) <- trimws(names(d))
  # sobra espaço à direita em quase todo campo de texto do arquivo
  for (j in seq_along(d)) if (is.character(d[[j]])) d[[j]] <- trimws(d[[j]])
  d
}

bruto <- do.call(rbind, lapply(ARQUIVO, le_imar))

cat("======================================================\n")
cat("LEITURA DO ARQUIVO DO IMar\n")
cat("======================================================\n")
cat(sprintf("Linhas lidas (1 por espécie por amostragem): %d\n", nrow(bruto)))
cat(sprintf("Colunas: %d\n", ncol(bruto)))

## =====================================================================
## 2) RELATÓRIO DE QUALIDADE — antes de filtrar qualquer coisa
##    (o que se descarta tem de ser sempre visível e justificado)
## =====================================================================
cat("\n--- qualidade dos campos-chave ---\n")
qual <- function(x, nome) cat(sprintf("  %-28s %5.1f%% ausente/zero\n", nome,
                                      100 * mean(is.na(x) | x == "" | (is.numeric(x) & !is.na(x) & x == 0))))
qual(bruto$Amostragem, "Amostragem (id da viagem)")
qual(bruto$Quantidade, "Quantidade (kg)")
qual(bruto$Num_dias, "Num_dias")
qual(bruto$Num_horas, "Num_horas")
qual(bruto$Numero_pescadores, "Numero_pescadores")
qual(bruto$Profundidade_pesca_engenho1, "Profundidade (0 = ausente)")
qual(bruto$Nome_banco_pesca, "Nome_banco_pesca")

## Linhas sem `Amostragem`: são registros de esforço sem captura associada
## (barco saiu e voltou, mas nenhuma espécie foi lançada na planilha). Não
## dá para saber se foi viagem sem captura ou falha de digitação, então
## saem da análise — mas o número fica registrado.
sem_id <- sum(is.na(bruto$Amostragem) | bruto$Amostragem == "")
cat(sprintf("\nLinhas sem identificador de amostragem: %d (descartadas)\n", sem_id))
d <- bruto[!is.na(bruto$Amostragem) & bruto$Amostragem != "", ]

## =====================================================================
## 3) FILTRO DE ARTE DE PESCA (decisão L5)
## =====================================================================
cat("\n--- viagens por arte de pesca (antes do filtro) ---\n")
viag_arte <- tapply(d$Amostragem, d$Nome_engenho, function(x) length(unique(x)))
print(sort(viag_arte, decreasing = TRUE))

foco <- d[toupper(d$Nome_cientifico) == ESPECIE_FOCO, ]
cat(sprintf("\n%s: %d viagens, %.1f t no total\n", ESPECIE_FOCO,
            length(unique(foco$Amostragem)), sum(foco$Quantidade, na.rm = TRUE) / 1000))
cat("distribuição dessas viagens por arte:\n")
print(tapply(foco$Amostragem, foco$Nome_engenho, function(x) length(unique(x))))
cat(sprintf("-> mantendo apenas `%s` (ver decisão L5 no cabeçalho)\n", ARTE_ALVO))

d <- d[d$Nome_engenho == ARTE_ALVO, ]

## =====================================================================
## 4) ESPÉCIES: quais viram coluna própria
## ---------------------------------------------------------------------
## Usamos NOME CIENTÍFICO, nunca o nome comum: no arquivo, "CAVALA PRETA"
## é Decapterus macarellus, mas "CAVALA BRANCA" é D. punctatus e "CAVALA
## DE RABO VERMELHA" é D. tabl — três espécies do mesmo gênero que seriam
## fundidas se a chave fosse o nome popular.
## As espécies que entram como coluna própria são as que sustentam a
## composição da captura (a impressão digital da tática); o resto vai
## para "outras".
## =====================================================================
kg_sp <- tapply(d$Quantidade, toupper(d$Nome_cientifico), sum, na.rm = TRUE)
kg_sp <- sort(kg_sp, decreasing = TRUE)
cat("\n--- capturas por espécie no engenho selecionado (t) ---\n")
print(round(head(kg_sp, 12) / 1000, 1))

SP_FOCAIS <- c(
  macarellus = "DECAPTERUS MACARELLUS",
  auxis      = "AUXIS SP",
  katsuwonus = "KATSUWONUS PELAMIS",
  carangideo = "CARANX CRYSOS",
  selar      = "SELAR CRUMENOPHTHALMUS",
  sardinella = "SARDINELLA MADERENSIS",
  punctatus  = "DECAPTERUS PUNCTATUS",
  thunnus    = "THUNNUS ALBACARES",
  trachurus  = "TRACHURUS SP"
)
d$sp_col <- "outras"
for (k in names(SP_FOCAIS)) d$sp_col[toupper(d$Nome_cientifico) == SP_FOCAIS[k]] <- k
cat(sprintf("\nEspécies com coluna própria: %d (+ 'outras' com %.1f%% da captura)\n",
            length(SP_FOCAIS),
            100 * sum(d$Quantidade[d$sp_col == "outras"], na.rm = TRUE) /
              sum(d$Quantidade, na.rm = TRUE)))

## =====================================================================
## 5) LONGO -> LARGO: uma linha por VIAGEM
## ---------------------------------------------------------------------
## Captura convertida de kg para TONELADAS (o arquivo vem em kg; a soma
## de 2024 no cerco bate com a ordem de grandeza dos relatórios da FAO).
## As variáveis operacionais são constantes dentro da viagem, então basta
## pegar a primeira ocorrência — mas conferimos isso antes.
## =====================================================================
chk <- tapply(d$Nome_embarcacao, d$Amostragem, function(x) length(unique(x)))
stopifnot(all(chk == 1))   # cada amostragem = 1 embarcação; se falhar, investigar

cap <- tapply(d$Quantidade / 1000, list(d$Amostragem, d$sp_col), sum)
cap[is.na(cap)] <- 0                      # espécie não registrada = 0 t
cap <- as.data.frame(cap)
todas_col <- c(names(SP_FOCAIS), "outras")
for (k in setdiff(todas_col, names(cap))) cap[[k]] <- 0
cap <- cap[, todas_col]
names(cap) <- paste0("cap_", names(cap))
cap$Amostragem <- rownames(cap)

primeiro <- function(col) tapply(d[[col]], d$Amostragem, function(x) x[1])
viagens <- data.frame(
  Amostragem = names(primeiro("Ano")),
  ano        = as.integer(primeiro("Ano")),
  mes        = as.integer(primeiro("Mes")),
  barco      = primeiro("Nome_embarcacao"),
  barco_id   = primeiro("Embarcacao"),
  tipo_emb   = primeiro("Tipo_embarcacao"),
  npesc      = as.numeric(primeiro("Numero_pescadores")),
  dias       = as.numeric(primeiro("Num_dias")),
  horas      = as.numeric(primeiro("Num_horas")),
  prof       = as.numeric(primeiro("Profundidade_pesca_engenho1")),
  banco      = primeiro("Nome_banco_pesca"),
  ilha_desemb = primeiro("Nome_ilha"),
  porto_desemb = primeiro("Nome_porto_desembarque"),
  data_saida = as.Date(primeiro("Data_partida"), format = "%m/%d/%Y"),
  data_volta = as.Date(primeiro("Data_chegada"), format = "%m/%d/%Y"),
  stringsAsFactors = FALSE
)
viagens <- merge(viagens, cap, by = "Amostragem", sort = FALSE)

## =====================================================================
## 6) LIMPEZA DAS VARIÁVEIS OPERACIONAIS (decisões L6-L8)
## =====================================================================

## L6 — profundidade: 0 é ausência de registro, não profundidade zero.
n_prof0 <- sum(viagens$prof == 0, na.rm = TRUE)
viagens$prof[viagens$prof == 0] <- NA
cat(sprintf("\nProfundidade: %d viagens (%.0f%%) com 0 -> convertidas para NA\n",
            n_prof0, 100 * n_prof0 / nrow(viagens)))

## Esforço — duas medidas concorrentes, e vale testar as duas:
##   `dias`  = Num_dias. Confere com (chegada - partida) em todos os casos
##             de 2+ dias; nas viagens de 1 dia o campo conta o dia de
##             saída, então é 1 mesmo quando saída e chegada são no mesmo
##             dia. É a medida "oficial" e comparável com a série da FAO.
##   `horas` = Num_horas, tempo no mar. Bem mais fino (mediana ~31 h, com
##             viagens de 1 dia variando de ~10 a ~28 h). Tende a ser o
##             melhor offset justamente por discriminar dentro do dia.
## A parte 03 ajusta os dois e compara por AIC — é comparação válida,
## porque a resposta e as linhas são as mesmas, só muda o offset.
viagens$dur_datas <- as.numeric(viagens$data_volta - viagens$data_saida)
cat(sprintf("Num_dias == (chegada - partida): %.1f%% das viagens\n",
            100 * mean(viagens$dias == viagens$dur_datas, na.rm = TRUE)))
cat("  (a diferença está toda nas viagens de 1 dia, em que o campo conta\n")
cat("   o dia de saída — comportamento consistente, não é erro)\n")

n_horas_0 <- sum(viagens$horas <= 0 | is.na(viagens$horas))
if (n_horas_0 > 0) cat(sprintf("Num_horas <= 0 ou ausente: %d viagens\n", n_horas_0))

## L7 — banco de pesca: agrupar níveis raros evita coeficientes instáveis.
tb_banco <- table(viagens$banco)
raros <- names(tb_banco[tb_banco < MIN_VIAG_BANCO])
viagens$banco_gr <- ifelse(viagens$banco %in% raros, "OUTROS", viagens$banco)
cat(sprintf("\nBancos de pesca: %d distintos; %d com < %d viagens agrupados em 'OUTROS'\n",
            length(tb_banco), length(raros), MIN_VIAG_BANCO))
cat(sprintf("  -> %d níveis usados no modelo\n", length(unique(viagens$banco_gr))))

## L8 — tipo de embarcação: I = industrial; A e R são minoria. Mantemos I
## e A como níveis e mandamos o resto para "outro".
cat("\nTipo de embarcação (viagens):\n"); print(table(viagens$tipo_emb))
viagens$tipo_emb[!viagens$tipo_emb %in% c("I", "A")] <- "outro" #Industrial mantida, o resto= Outros

## Tripulação — proxy de poder de pesca. É a única variável operacional
## contínua disponível além da profundidade, então é valiosa: é ela que
## permite separar "barco maior" de "mais peixe no mar".
cat(sprintf("\nTripulação: mediana %.0f, intervalo %.0f-%.0f pescadores\n",
            median(viagens$npesc, na.rm = TRUE), min(viagens$npesc, na.rm = TRUE),
            max(viagens$npesc, na.rm = TRUE)))

## =====================================================================
## 7) TABELA FINAL E CONFERÊNCIA
## =====================================================================
cols_cap <- paste0("cap_", todas_col)
viagens$cap_total <- rowSums(viagens[, cols_cap])

cat("\n======================================================\n")
cat("TABELA `viagens` PRONTA\n")
cat("======================================================\n")
cat(sprintf("Viagens: %d | Embarcações: %d | Meses: %d | Ano(s): %s\n",
            nrow(viagens), length(unique(viagens$barco)),
            length(unique(viagens$mes)),
            paste(sort(unique(viagens$ano)), collapse = ", ")))
cat(sprintf("Esforço total: %.0f dias | %.0f horas no mar\n",
            sum(viagens$dias, na.rm = TRUE), sum(viagens$horas, na.rm = TRUE)))
cat(sprintf("Captura total: %.1f t | cavala preta: %.1f t (%.1f%%)\n",
            sum(viagens$cap_total), sum(viagens$cap_macarellus),
            100 * sum(viagens$cap_macarellus) / sum(viagens$cap_total)))
cat(sprintf("Viagens SEM cavala preta: %d (%.1f%%) <- zeros de direcionamento\n",
            sum(viagens$cap_macarellus == 0),
            100 * mean(viagens$cap_macarellus == 0)))
cat(sprintf("Viagens que registraram UMA só espécie: %.1f%%\n",
            100 * mean(rowSums(viagens[, cols_cap] > 0) == 1)))
cat("  ^ este número importa: composição de UMA espécie não identifica\n")
cat("    tática nenhuma. É por isso que a parte 02 agrega a composição\n")
cat("    por BARCO-MÊS antes de agrupar (Hoyle et al. 2018).\n")

cat("\nEstrutura:\n"); str(viagens, give.attr = FALSE)

if (length(unique(viagens$ano)) == 1) {
  cat("\n**********************************************************\n")
  cat("ATENÇÃO: um único ano na base. O efeito de ANO — que É o índice\n")
  cat("de abundância — não pode ser estimado. A parte 02 vai usar MÊS\n")
  cat("como fator temporal para demonstrar o pipeline. Assim que a série\n")
  cat("histórica chegar, o script troca sozinho para ano.\n")
  cat("**********************************************************\n")
}









#======================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde
#----------------------------------------------------------------------
# Autor : Silva, MLS
# Data  : 2026-09-16
# Parte : 02 de 03 — EXPLORATÓRIA E INFERÊNCIA DA TÁTICA DE PESCA
#
# Pré-requisito: objeto `viagens` (rodar a parte 01).
# Só usa R base — roda em qualquer máquina, sem instalar nada.
#
# O QUE ESTE SCRIPT DECIDE, E POR QUÊ
# -----------------------------------
# A análise inteira gira em torno de uma variável que NÃO existe no banco:
# a espécie que o mestre pretendia pescar. Ninguém registra intenção. O
# que existe é a COMPOSIÇÃO DA CAPTURA, que funciona como impressão
# digital da tática usada. Este script transforma essa composição em uma
# covariável utilizável, de duas formas concorrentes:
#
#   (a) DISCRETA  — agrupa as unidades de esforço por similaridade de
#                   composição (k-means / Ward) e usa o rótulo do grupo
#                   como FATOR `alvo`. É o padrão nos grupos do ICCAT
#                   (Sant'Ana et al. 2020).
#   (b) CONTÍNUA  — usa os escores dos primeiros eixos de uma PCA da
#                   composição (PC1, PC2) como preditores contínuos. É o
#                   "DPC" de Winker et al. (2013), que em dois estudos
#                   independentes ajustou melhor que o cluster.
#
# As duas são construídas aqui e comparadas formalmente na parte 03
# (hipótese H2). Não se escolhe uma por gosto: escolhe-se por ajuste.
#
# POR QUE AGREGAR POR BARCO-MÊS ANTES DE AGRUPAR
# ----------------------------------------------
# Em 2024, ~57% das viagens de cerco registraram UMA única espécie. Uma
# viagem com uma espécie só tem composição degenerada (100% daquela
# espécie) e não distingue "tática dirigida" de "sorte num lance". Hoyle
# et al. (2018) recomendam agregar por barco-mês exatamente por isso: a
# agregação dilui o encontro fortuito com o cardume e deixa o padrão de
# estratégia aparecer. É o que fazemos.
#======================================================================

stopifnot(exists("viagens"))

## =====================================================================
## 0) FATOR TEMPORAL — ano se houver série, mês se só houver um ano
## ---------------------------------------------------------------------
## O índice de abundância É o efeito do fator temporal. Com um ano só,
## não há efeito de ano; usamos mês para demonstrar o pipeline. Quando a
## série completa chegar, esta linha troca sozinha e nada mais muda.
## =====================================================================
UM_ANO_SO <- length(unique(viagens$ano)) == 1
fator_tempo  <- if (UM_ANO_SO) "fmes" else "fano"
rotulo_tempo <- if (UM_ANO_SO) "Mês (2024)" else "Ano"
cat(sprintf("\n>> Fator temporal desta rodada: `%s` (%s)\n", fator_tempo,
            if (UM_ANO_SO) "DEMONSTRAÇÃO — um ano só" else "série histórica"))

cols_cap <- grep("^cap_", names(viagens), value = TRUE)
cols_cap <- setdiff(cols_cap, "cap_total")
sp_nomes <- c(cap_macarellus = "D. macarellus", cap_auxis = "Auxis sp.",
              cap_katsuwonus = "K. pelamis",    cap_carangideo = "C. crysos",
              cap_selar = "S. crumenoph.", cap_sardinella = "S. maderensis",
              cap_punctatus = "D. punctatus",   cap_thunnus = "T. albacares",
              cap_trachurus = "Trachurus sp.",  cap_outras = "Outras")
sp_nomes <- sp_nomes[cols_cap]

cor_sp  <- hcl.colors(length(cols_cap), palette = "Dark 3")
COR_MAC <- "#1F4E79"; COR_AUX <- "#C0501B"; COR_NEU <- "#7F7F7F"

## =====================================================================
## 1) VARIÁVEIS DERIVADAS
## ---------------------------------------------------------------------
## Duas CPUE nominais, uma por medida de esforço. Elas NÃO entram no
## modelo (lá o esforço entra como offset); servem para a exploratória e
## para o cenário S1, que reproduz o que a FAO (2026) fez.
## =====================================================================
viagens$cpue_dia  <- viagens$cap_macarellus / viagens$dias
viagens$cpue_hora <- viagens$cap_macarellus / pmax(viagens$horas, 1)
viagens$pos_mac   <- as.integer(viagens$cap_macarellus > 0)

# composição proporcional DA VIAGEM (usada só para descrever; o cluster
# usa a composição por barco-mês, ver seção 5)
for (cc in cols_cap) {
  viagens[[sub("cap_", "prop_", cc)]] <-
    ifelse(viagens$cap_total > 0, viagens[[cc]] / viagens$cap_total, 0)
}

# Fatores. TUDO que vai ser marginalizado depois precisa ser fator — é
# assim que o emmeans consegue tirar a média sobre os níveis e devolver
# o efeito do tempo "limpo" dos demais.
viagens$fano     <- factor(viagens$ano)
viagens$fmes     <- factor(viagens$mes, levels = 1:12)
viagens$fbanco   <- factor(viagens$banco_gr)
viagens$filha    <- factor(viagens$ilha_desemb)
viagens$fbarco   <- factor(viagens$barco)
viagens$ftipo    <- factor(viagens$tipo_emb)
viagens$barco_mes <- paste(viagens$barco, viagens$ano, viagens$mes, sep = "_")

## =====================================================================
## 2) FILTROS — cada um registrado, com o que custou
## ---------------------------------------------------------------------
## Filtro é decisão analítica, não faxina: muda a população amostrada e
## precisa ser reportado no texto. Por isso o log fica num data.frame.
## Os limiares aqui são folgados de propósito, porque com um ano só a
## amostra é pequena; com a série completa dá para apertar (Sant'Ana et
## al. 2020 usam >= 50 lances por embarcação, por exemplo).
## =====================================================================
n0 <- nrow(viagens)
filtro_log <- data.frame(regra = character(), removidas = integer(),
                         restantes = integer(), stringsAsFactors = FALSE)
reg <- function(regra, antes)
  filtro_log[nrow(filtro_log) + 1, ] <<- list(regra, antes - nrow(viagens),
                                              nrow(viagens))

a <- nrow(viagens); viagens <- viagens[!is.na(viagens$dias) & viagens$dias > 0, ]
reg("esforço (dias) válido e > 0", a)
a <- nrow(viagens); viagens <- viagens[!is.na(viagens$horas) & viagens$horas > 0, ]
reg("esforço (horas) válido e > 0", a)
a <- nrow(viagens); viagens <- viagens[viagens$cap_total > 0, ]
reg("alguma captura registrada na viagem", a)

# Embarcações com pouquíssimas viagens não sustentam efeito aleatório e
# ainda desequilibram o cruzamento com o fator temporal.
MIN_VIAG_BARCO <- 5
tb <- table(viagens$barco)
a <- nrow(viagens)
viagens <- viagens[viagens$barco %in% names(tb[tb >= MIN_VIAG_BARCO]), ]
reg(sprintf("embarcações com >= %d viagens", MIN_VIAG_BARCO), a)

# Estratos tempo x banco muito ralos geram coeficientes instáveis.
te <- table(paste(viagens[[sub("^f", "", fator_tempo)]], viagens$banco_gr))
a <- nrow(viagens)
viagens <- viagens[paste(viagens[[sub("^f", "", fator_tempo)]],
                         viagens$banco_gr) %in% names(te[te >= 3]), ]
reg("estratos tempo x banco com >= 3 viagens", a)

for (f in c("fano", "fmes", "fbanco", "filha", "fbarco", "ftipo"))
  viagens[[f]] <- droplevels(viagens[[f]])

cat("\n===== FILTROS =====\n"); print(filtro_log, row.names = FALSE)
cat(sprintf("Retidas %d de %d viagens (%.1f%%)\n", nrow(viagens), n0,
            100 * nrow(viagens) / n0))

## =====================================================================
## 3) SÉRIE NOMINAL (cenário S1) — captura da espécie / esforço TOTAL
## ---------------------------------------------------------------------
## Este é o índice que a FAO usou e o que queremos submeter à prova: o
## denominador é o esforço de TODA a frota de cerco, inclusive as viagens
## que estavam atrás de Auxis. É justamente por isso que ele é suspeito.
## =====================================================================
tempo_num <- viagens[[sub("^f", "", fator_tempo)]]
viagens$tempo <- tempo_num

ser <- aggregate(viagens[, c(cols_cap, "dias", "horas")],
                 by = list(tempo = viagens$tempo), FUN = sum)
ser$cpue_nom_dia  <- ser$cap_macarellus / ser$dias
ser$cpue_nom_hora <- ser$cap_macarellus / ser$horas
ser$prop_zero <- tapply(viagens$pos_mac == 0, viagens$tempo, mean)[as.character(ser$tempo)]
ser$n_viagens <- as.numeric(table(viagens$tempo)[as.character(ser$tempo)])

cat(sprintf("\n===== SÉRIE POR %s =====\n", toupper(rotulo_tempo)))
print(data.frame(tempo = ser$tempo,
                 cap_macarellus_t = round(ser$cap_macarellus, 1),
                 cap_auxis_t = round(ser$cap_auxis, 1),
                 dias = ser$dias, horas = ser$horas,
                 cpue_t_dia = round(ser$cpue_nom_dia, 3),
                 zeros = sprintf("%.0f%%", 100 * ser$prop_zero),
                 n = ser$n_viagens), row.names = FALSE)

## =====================================================================
## 4) FIGURAS EXPLORATÓRIAS
## =====================================================================

## --- Fig 1: captura por espécie e esforço ------------------------------
png("exp1_capturas_esforco.png", width = 26, height = 13, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
matplot(ser$tempo, ser[, cols_cap], type = "l", lty = 1, lwd = 2.2,
        col = cor_sp, xlab = rotulo_tempo, ylab = "Captura (t)",
        main = "A. Captura por espécie")
legend("topleft", sp_nomes, col = cor_sp, lwd = 2.2, bty = "n", cex = 0.62)
plot(ser$tempo, ser$dias, type = "l", lwd = 2.4, col = COR_NEU,
     xlab = rotulo_tempo, ylab = "Esforço (dias de pesca)",
     main = "B. Esforço da frota de cerco", ylim = c(0, max(ser$dias) * 1.05))
par(op); dev.off()
cat("\nPNG salvo: exp1_capturas_esforco.png\n")

## --- Fig 2: composição da captura --------------------------------------
comp <- as.matrix(ser[, cols_cap]); comp <- comp / rowSums(comp)
png("exp2_composicao.png", width = 24, height = 12, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mar = c(4.2, 4.6, 3, 9), bty = "l", cex.main = 0.95, cex = 0.85)
acum <- t(apply(comp, 1, cumsum))
plot(NA, xlim = range(ser$tempo), ylim = c(0, 1), xlab = rotulo_tempo,
     ylab = "Proporção da captura",
     main = "Composição da captura da frota de cerco")
for (k in ncol(acum):1)
  polygon(c(ser$tempo, rev(ser$tempo)), c(acum[, k], rep(0, nrow(acum))),
          col = cor_sp[k], border = NA)
par(xpd = TRUE)
legend(max(ser$tempo) + diff(range(ser$tempo)) * 0.03, 0.9, sp_nomes,
       fill = cor_sp, border = NA, bty = "n", cex = 0.68)
par(op); dev.off()
cat("PNG salvo: exp2_composicao.png\n")

## --- Fig 3: CPUE nominal e zeros ---------------------------------------
## O painel B é o diagnóstico mais importante desta figura: com ~87% de
## viagens sem cavala, qualquer modelo que não trate os zeros direito
## (delta ou Tweedie) vai dar resultado errado.
png("exp3_cpue_nominal.png", width = 26, height = 13, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
plot(ser$tempo, ser$cpue_nom_dia, type = "l", lwd = 2.4, col = COR_MAC,
     xlab = rotulo_tempo, ylab = "CPUE nominal da cavala (t/dia)",
     main = "A. CPUE nominal", ylim = c(0, max(ser$cpue_nom_dia) * 1.05))
plot(ser$tempo, 100 * ser$prop_zero, type = "l", lwd = 2.4, col = COR_AUX,
     xlab = rotulo_tempo, ylab = "% de viagens sem cavala",
     main = "B. Zeros de direcionamento", ylim = c(0, 100))
par(op); dev.off()
cat("PNG salvo: exp3_cpue_nominal.png\n")

## --- Fig 4: as covariáveis operacionais valem a pena? -------------------
## COMO LER: com ~87% de viagens sem cavala, gráfico de CPUE contra
## covariável vira uma parede de zeros e não informa nada. A leitura
## correta para dado assim é decompor, que é exatamente o que o modelo
## delta faz: (i) a PROBABILIDADE de a viagem pegar cavala e (ii) QUANTO
## ela pega, dado que pegou. Uma covariável pode atuar só na primeira
## (tem a ver com onde/como se procura), só na segunda (tem a ver com
## capacidade de captura) ou nas duas. Estes painéis são a justificativa
## empírica de cada termo do modelo.
classe <- function(x, n = 4) {
  q <- unique(quantile(x, seq(0, 1, length.out = n + 1), na.rm = TRUE))
  cut(x, breaks = q, include.lowest = TRUE, dig.lab = 3)
}
barra_prop <- function(prop, n, titulo, xlab, cex_nome = 0.7, cor = "#8FAADC") {
  bp <- barplot(100 * prop, col = cor, border = NA, las = 2,
                cex.names = cex_nome, ylab = "% de viagens com cavala",
                xlab = xlab, main = titulo,
                ylim = c(0, max(100 * prop, na.rm = TRUE) * 1.25))
  text(bp, 100 * prop, labels = paste0("n=", n), pos = 3, cex = 0.55,
       col = "#52514E", xpd = NA)
  invisible(bp)
}

png("exp4_covariaveis.png", width = 26, height = 20, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(2, 2), mar = c(7.5, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)

## A — onde se pega cavala (componente de presença, por banco de pesca)
pr_b <- tapply(viagens$pos_mac, viagens$fbanco, mean)
n_b  <- table(viagens$fbanco)
ord  <- order(pr_b, decreasing = TRUE)
barra_prop(pr_b[ord], n_b[ord], "A. Presença de cavala por banco de pesca",
           "", cex_nome = 0.55)

## B — quanto se pega, dado que pegou (componente de magnitude)
pos <- viagens[viagens$pos_mac == 1, ]
bancos_ok <- names(which(table(pos$fbanco) >= 5))
pos_b <- pos[pos$fbanco %in% bancos_ok, ]
if (nrow(pos_b) > 0) {
  pos_b$fbanco <- droplevels(pos_b$fbanco)
  bx <- boxplot(cpue_dia ~ fbanco, data = pos_b, outline = FALSE, plot = FALSE)
  boxplot(cpue_dia ~ fbanco, data = pos_b, outline = FALSE, col = "#74C476",
          xaxt = "n", xlab = "", ylab = "CPUE (t/dia) entre as positivas",
          lwd = 1, main = "B. Magnitude, só nas viagens com cavala")
  axis(1, at = seq_along(bx$names), labels = FALSE)
  text(seq_along(bx$names), par("usr")[3], labels = bx$names, srt = 45,
       adj = 1, xpd = NA, cex = 0.6)
}

## C — tripulação: proxy de poder de pesca. Em classes, porque a relação
##     não tem de ser linear e classes mostram a forma sem impor nada.
cl_p <- classe(viagens$npesc)
barra_prop(tapply(viagens$pos_mac, cl_p, mean), table(cl_p),
           "C. Tripulação (proxy de poder de pesca)", "Nº de pescadores",
           cor = "#B497D6")

## D — profundidade, só onde foi registrada (os 0 viraram NA na parte 01)
ok <- !is.na(viagens$prof)
if (sum(ok) > 30) {
  cl_d <- classe(viagens$prof[ok])
  barra_prop(tapply(viagens$pos_mac[ok], cl_d, mean), table(cl_d),
             sprintf("D. Profundidade (n=%d com registro)", sum(ok)),
             "Profundidade (m)", cor = "#E8A33D")
}
par(op); dev.off()
cat("PNG salvo: exp4_covariaveis.png\n")

## Tipo de embarcação tem poucos níveis: tabela é mais honesta que gráfico
cat("\nPresença de cavala por tipo de embarcação:\n")
print(round(cbind(n = table(viagens$ftipo),
                  prop_com_cavala = tapply(viagens$pos_mac, viagens$ftipo, mean),
                  cpue_media = tapply(viagens$cpue_dia, viagens$ftipo, mean)), 3))

## =====================================================================
## 5) INFERÊNCIA DA TÁTICA — composição por barco-mês
## =====================================================================

## 5.1 matriz de composição agregada por BARCO-MÊS (ver justificativa no
##     cabeçalho). Proporções em peso; transformação raiz quadrada para
##     que espécies menos abundantes também pesem na similaridade, como
##     em Winker et al. (2013).
agg <- aggregate(viagens[, cols_cap],
                 by = list(barco_mes = viagens$barco_mes), FUN = sum)
tot <- rowSums(agg[, cols_cap])
comp_bm <- as.matrix(agg[, cols_cap][tot > 0, , drop = FALSE] / tot[tot > 0])
rownames(comp_bm) <- agg$barco_mes[tot > 0]
comp_sqrt <- sqrt(comp_bm)

n_sp_bm <- rowSums(comp_bm > 0)
cat(sprintf("\nComposição: %d barcos-mês x %d espécies\n",
            nrow(comp_sqrt), ncol(comp_sqrt)))
cat(sprintf("Barcos-mês com uma só espécie: %.1f%% (era %.1f%% por viagem)\n",
            100 * mean(n_sp_bm == 1),
            100 * mean(rowSums(viagens[, cols_cap] > 0) == 1)))
cat("  ^ a agregação por barco-mês reduz a degeneração da composição;\n")
cat("    se ainda assim a maioria tiver uma espécie só, o cluster fica\n")
cat("    fraco e a representação contínua (PCA) é a mais defensável.\n")

## 5.2 PCA. Sem `scale.` porque as colunas já estão na mesma unidade
##     (proporções transformadas) — escalonar daria peso igual a espécies
##     raras e dominantes, o que não é o que queremos.
pca <- prcomp(comp_sqrt, center = TRUE, scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cat("Variância explicada: ",
    paste(sprintf("PC%d=%.1f%%", 1:min(4, length(var_exp)),
                  var_exp[1:min(4, length(var_exp))]), collapse = "  "), "\n")
cat("Cargas de PC1 (o que este eixo separa):\n")
print(round(sort(pca$rotation[, 1]), 3))

## Quantos eixos levar para o modelo? Guardamos os primeiros PCs até
## acumular ~70% da variação da composição (teto de 4, para não inflar o
## modelo). Com poucas espécies dominantes, 1-2 eixos bastam; numa
## pescaria com composição realmente multidimensional — como a de 2024,
## em que PC1 explica só ~25% — usar apenas PC1/PC2 jogaria fora a maior
## parte do sinal de tática.
cum_var <- cumsum(var_exp)
n_pc <- min(4, max(2, which(cum_var >= 70)[1]), ncol(pca$x))
if (is.na(n_pc)) n_pc <- min(4, ncol(pca$x))
cat(sprintf("Eixos retidos para o modelo: %d (%.0f%% da variação)\n",
            n_pc, cum_var[n_pc]))

## 5.3 Número de grupos pela silhueta média.
##     Implementada aqui em R base para não exigir o pacote `cluster`.
##     ATENÇÃO à leitura: a silhueta favorece sistematicamente k pequeno.
##     Ela responde "os grupos estão separados?", não "quantas táticas
##     existem?". Se o gráfico de PC1xPC2 mostrar uma nuvem contínua, a
##     resposta honesta é que não há grupos — há um gradiente.
silhueta_media <- function(X, cl, n_sub = 1200) {
  idx <- if (nrow(X) > n_sub) sample(nrow(X), n_sub) else seq_len(nrow(X))
  Xs <- X[idx, , drop = FALSE]; cs <- cl[idx]
  D <- as.matrix(dist(Xs)); gr <- unique(cs)
  if (length(gr) < 2) return(NA_real_)
  mean(vapply(seq_along(cs), function(i) {
    mesmo <- cs == cs[i]; mesmo[i] <- FALSE
    if (!any(mesmo)) return(0)
    ai <- mean(D[i, mesmo])
    bi <- min(vapply(setdiff(gr, cs[i]), function(g) mean(D[i, cs == g]), numeric(1)))
    (bi - ai) / max(ai, bi)
  }, numeric(1)))
}

k_forcado <- NA          # preencha para fixar k por conhecimento da pescaria
esc <- pca$x[, 1:min(3, ncol(pca$x)), drop = FALSE]
set.seed(42)
ks  <- 2:6
sil <- vapply(ks, function(k)
  silhueta_media(esc, kmeans(esc, centers = k, nstart = 25, iter.max = 50)$cluster),
  numeric(1))
k_otimo <- if (is.na(k_forcado)) ks[which.max(sil)] else k_forcado
cat("\nSilhueta média: ",
    paste(sprintf("k=%d: %.3f", ks, sil), collapse = "  "), "\n")
cat(sprintf("k adotado: %d%s\n", k_otimo,
            if (is.na(k_forcado)) " (silhueta)" else " (fixado)"))

## 5.4 k-means (partição) e Ward (hierárquico) — dois algoritmos com
##     lógicas diferentes. Se discordarem, o agrupamento não é estável e
##     isso é informação, não um contratempo.
km  <- kmeans(esc, centers = k_otimo, nstart = 50, iter.max = 100)
sub <- sample(nrow(esc), min(2500, nrow(esc)))
ward <- cutree(hclust(dist(esc[sub, , drop = FALSE]), method = "ward.D2"), k = k_otimo)
tb_kw <- table(kmeans = km$cluster[sub], ward = ward)
conc_kw <- sum(apply(tb_kw, 1, max)) / sum(tb_kw)
cat(sprintf("Concordância k-means x Ward: %.1f%%%s\n", 100 * conc_kw,
            if (conc_kw < 0.80) "  <- baixa: agrupamento instável" else ""))

## 5.5 Nomear cada grupo pela espécie dominante do seu centróide, para
##     que o fator `alvo` seja legível no output do modelo.
cent <- t(vapply(seq_len(k_otimo), function(g)
  colMeans(comp_bm[km$cluster == g, , drop = FALSE]), numeric(ncol(comp_bm))))
colnames(cent) <- sub("cap_", "", colnames(comp_bm))
nome_cl <- make.unique(colnames(cent)[apply(cent, 1, which.max)], sep = "_")
cat("\nComposição média de cada tática inferida:\n")
print(round(cbind(cent, n = as.numeric(table(km$cluster))), 3))

## Qual dessas táticas é a "da cavala"? NÃO se identifica pelo nome do
## grupo: o nome vem da espécie dominante do centróide, e a cavala pode
## não dominar nenhum grupo (em 2024 ela é ~7% da captura do cerco). A
## tática relevante é aquela com a MAIOR proporção de cavala no
## centróide, seja qual for o nome dela. Isso é guardado para o cenário
## H3 da parte 03.
alvo_cavala <- nome_cl[which.max(cent[, "macarellus"])]
cat(sprintf("\nTática com maior fração de cavala no centróide: '%s' (%.1f%% da captura)\n",
            alvo_cavala, 100 * max(cent[, "macarellus"])))
if (max(cent[, "macarellus"]) < 0.30)
  cat("  [NOTA] Nenhuma tática é dominada pela cavala. Ela é capturada\n",
      "        acompanhando outras espécies, não como alvo exclusivo —\n",
      "        o que enfraquece o cenário de 'esforço dirigido' (H3) e\n",
      "        reforça o uso da composição como covariável (H1/H2).\n")

## 5.6 Levar o rótulo e os escores de volta para cada VIAGEM.
mapa <- data.frame(barco_mes = rownames(comp_bm),
                   alvo = factor(nome_cl[km$cluster], levels = unique(nome_cl)),
                   stringsAsFactors = FALSE)
PCs <- paste0("PC", seq_len(n_pc))
for (i in seq_len(n_pc)) mapa[[PCs[i]]] <- pca$x[, i]
viagens <- merge(viagens, mapa, by = "barco_mes", all.x = TRUE, sort = FALSE)
viagens <- viagens[!is.na(viagens$alvo), ]
viagens$alvo <- droplevels(viagens$alvo)
cat(sprintf("\nViagens com tática atribuída: %d\n", nrow(viagens)))
print(table(viagens$alvo))

## 5.7 Validação possível SEM gabarito: a tática inferida tem de separar
##     a captura de cavala. Se não separar, ela não está medindo alvo.
cat("\nCPUE média da cavala por tática inferida (t/dia):\n")
print(round(tapply(viagens$cpue_dia, viagens$alvo, mean), 3))
cat("Proporção de viagens com cavala, por tática:\n")
print(round(tapply(viagens$pos_mac, viagens$alvo, mean), 3))

## --- Fig 5: PCA, silhueta e composição das táticas ---------------------
cor_cl <- hcl.colors(k_otimo, palette = "Dark 3")
png("exp5_taticas.png", width = 27, height = 10.5, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 3), mar = c(4.6, 4.4, 3, 1), oma = c(0, 0, 0, 5),
          bty = "l", cex.main = 0.95, cex = 0.85)
plot(ks, sil, type = "b", pch = 19, lwd = 2, col = COR_MAC,
     xlab = "Número de grupos (k)", ylab = "Silhueta média", main = "A. Escolha de k")
points(k_otimo, sil[ks == k_otimo], pch = 21, bg = COR_AUX, cex = 1.9)
plot(pca$x[, 1], pca$x[, 2], col = adjustcolor(cor_cl[km$cluster], 0.6),
     pch = 16, cex = 0.7, xlab = sprintf("PC1 (%.0f%%)", var_exp[1]),
     ylab = sprintf("PC2 (%.0f%%)", var_exp[2]),
     main = "B. Táticas no espaço de composição")
points(km$centers[, 1], km$centers[, 2], pch = 21, bg = cor_cl, cex = 2, lwd = 1.5)
legend("topright", nome_cl, col = cor_cl, pch = 16, bty = "n", cex = 0.7)
bp <- barplot(t(cent), col = cor_sp, border = NA, names.arg = nome_cl,
              las = 2, cex.names = 0.65, ylab = "Proporção média da captura",
              main = "C. Composição de cada tática")
legend(max(bp) + 0.8, 1, rev(sp_nomes), fill = rev(cor_sp), border = NA,
       bty = "n", cex = 0.62, xpd = NA)
par(op); dev.off()
cat("PNG salvo: exp5_taticas.png\n")

## =====================================================================
## 6) ESFORÇO DIRIGIDO (insumo do cenário H3 na parte 03)
## ---------------------------------------------------------------------
## Reparte os dias de pesca de cada estrato entre as táticas, na
## proporção dos dias das viagens de cada uma. É isto que converte
## "esforço total da frota" em "esforço dirigido à cavala" — o campo que
## o IMar não tem e que a composição da captura permite reconstruir.
## =====================================================================
esforco_dirigido <- aggregate(cbind(dias, horas) ~ tempo + banco_gr + alvo,
                              data = viagens, FUN = sum)
dias_alvo <- aggregate(dias ~ tempo + alvo, data = esforco_dirigido, FUN = sum)

png("exp6_esforco_dirigido.png", width = 26, height = 12, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
tat <- prop.table(table(viagens$tempo, viagens$alvo), margin = 1)
matplot(as.numeric(rownames(tat)), tat, type = "l", lty = 1, lwd = 2.4,
        col = cor_cl, ylim = c(0, 1), xlab = rotulo_tempo,
        ylab = "Proporção das viagens", main = "A. Mistura de táticas")
legend("topleft", nome_cl, col = cor_cl, lwd = 2.4, bty = "n", cex = 0.7)
mat_d <- sapply(levels(viagens$alvo), function(g) {
  x <- dias_alvo[dias_alvo$alvo == g, ]
  v <- setNames(rep(0, nrow(ser)), ser$tempo)
  v[as.character(x$tempo)] <- x$dias; v
})
matplot(ser$tempo, mat_d, type = "l", lty = 1, lwd = 2.4, col = cor_cl,
        ylim = c(0, max(c(mat_d, ser$dias)) * 1.05), xlab = rotulo_tempo,
        ylab = "Dias de pesca", main = "B. Esforço dirigido por tática")
lines(ser$tempo, ser$dias, lwd = 2, lty = 2, col = COR_NEU)
legend("topleft", c(levels(viagens$alvo), "esforço total"),
       col = c(cor_cl, COR_NEU), lwd = 2.2,
       lty = c(rep(1, nlevels(viagens$alvo)), 2), bty = "n", cex = 0.7)
par(op); dev.off()
cat("PNG salvo: exp6_esforco_dirigido.png\n")

cat("\n===== PARTE 02 CONCLUÍDA =====\n")
cat("Objetos para a parte 03: `viagens` (com alvo, PC1, PC2, tempo),\n")
cat("`ser`, `esforco_dirigido`, `fator_tempo`, `rotulo_tempo`.\n")



















#======================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde
#----------------------------------------------------------------------
# Autor : Silva, MLS
# Data  : 2026-09-16
# Parte : 03 de 03 — MODELAGEM, SELEÇÃO, DIAGNÓSTICO E ÍNDICES
#
# Pré-requisito: partes 01 e 02 (objetos `viagens`, `ser`, `fator_tempo`,
# `rotulo_tempo`, `PCs`, `alvo_cavala`, `esforco_dirigido`).
# Pacotes: glmmTMB, emmeans, DHARMa (writexl é opcional).
#======================================================================
#
# ======================= O QUE ESTAMOS TESTANDO =======================
# Produzir índices alternativos de abundância relativa da cavala preta
# para entrar no JABBA como CENÁRIOS CONCORRENTES:
#   S1 nominal   — captura da cavala / esforço total da frota de cerco
#                  (o que a FAO 2026 fez)
#   S2 corrigida — efeito do fator temporal num modelo que controla
#                  banco de pesca, tripulação, tipo de embarcação,
#                  embarcação e TÁTICA DE PESCA
# A distância entre as duas É a medida do viés de direcionamento.
#
# ================= PROBLEMAS DESTA CPUE (resumo) ======================
# P1 esforço não é específico da espécie (denominador multiespecífico)
# P2 troca de alvo ao longo do tempo (Schirripa & Goodyear 2010)
# P3 zeros de direcionamento — 87% das viagens de cerco de 2024 não
#    registraram cavala. Excluí-los ou somar constante enviesa o índice.
# P4 poucas variáveis operacionais — MAS o arquivo do IMar tem mais do
#    que esperávamos: tripulação, profundidade, banco de pesca, tipo e
#    identidade da embarcação. Isso eleva a análise de "correção de
#    alvo" para algo mais perto de uma padronização de fato.
# P5 composição da frota muda ao longo do tempo
# P6 desbalanceamento espacial e sazonal
#
# ========================== HIPÓTESES =================================
# H0 a CPUE nominal é aceitável (o índice não muda ao controlar o resto)
# H1 parte da variação da CPUE nominal é troca de alvo, não abundância
# H2 tática contínua (PCA) ajusta melhor que tática discreta (cluster)
# H3 usar só o esforço dirigido reproduz H1 (sensibilidade; tende a
#    gerar hiperestabilidade, então nunca é o cenário principal)
#
# ===================== MODELOS QUE SERÃO TESTADOS =====================
# ESTRUTURAS (o fator temporal NUNCA entra na seleção — ele é o índice):
#   E0  tempo
#   E1  tempo + banco
#   E2  tempo + banco + tipo_emb + npesc
#   E3  tempo + banco + tipo_emb + npesc + alvo          (tática discreta)
#   E4  E3 + (1 | barco)                                  [GLMM]
#   E5  tempo + banco + tipo_emb + npesc + PC1..PCn + (1 | barco)
#                                                         (tática contínua)
# OFFSET: log(dias) contra log(horas) — comparável por AIC (mesma
#   resposta, mesmas linhas; só muda o termo de esforço).
# DISTRIBUIÇÕES (resposta = toneladas, contínua, com 87% de zeros):
#   D1 Tweedie             D2 Hurdle-Gamma
#   D3 Lognormal (positivas, componente do delta)
#   D4 Gamma (positivas)   D5 Gaussiana em log(CPUE + c)
#
# ---------------------------------------------------------------------
# SOBRE COMPARAR AIC/BIC ENTRE AS DISTRIBUIÇÕES  (ponto importante)
# ---------------------------------------------------------------------
# AIC e BIC só podem ser comparados entre modelos ajustados às MESMAS
# observações e à MESMA variável-resposta. Por isso a tabela de
# distribuições traz uma coluna `grupo`, e o dAIC é calculado DENTRO de
# cada grupo — nunca entre grupos:
#   grupo A (D1, D2): mesma resposta (toneladas), todas as linhas.
#                     -> AIC/BIC comparáveis entre si. É aqui que se
#                        escolhe a distribuição do modelo principal.
#   grupo B (D3, D4): só as viagens com captura positiva (13% das
#                     linhas). Comparáveis entre si, NÃO com o grupo A.
#   grupo C (D5):     a resposta é log(CPUE + c), outra escala; a
#                     verossimilhança nem é da mesma quantidade.
#                     Sozinho no grupo — serve só para comparar o
#                     ÍNDICE que produz, nunca o AIC.
# O AIC dos grupos B e C aparece na tabela porque é informação útil
# dentro do grupo; o que não se pode fazer é olhar a coluna inteira e
# escolher o menor número. A coluna `grupo` existe exatamente para
# impedir esse erro.
#
# E vale o princípio geral: AIC/BIC medem ajuste, não adequação. A
# escolha final considera também os resíduos simulados (DHARMa) e o
# comportamento do índice. Um modelo pode ganhar no AIC e ter resíduo
# ruim — nesse caso ele não é o escolhido.
#======================================================================

stopifnot(exists("viagens"), exists("ser"), exists("fator_tempo"))
suppressPackageStartupMessages({
  library(glmmTMB); library(emmeans); library(DHARMa)
})
tem_writexl <- requireNamespace("writexl", quietly = TRUE)

COR_MAC <- "#1F4E79"; COR_AUX <- "#C0501B"; COR_NEU <- "#7F7F7F"
COR_S2  <- "#2E8B57"; COR_S3 <- "#7030A0"; COR_S2D <- "#00A0B0"

## =====================================================================
## 0) PREPARO DA RESPOSTA E DO ESFORÇO
## =====================================================================
viagens$captura <- viagens$cap_macarellus     # toneladas
viagens$ldias   <- log(viagens$dias)
viagens$lhoras  <- log(viagens$horas)

## Tripulação centrada: com `npesc` bruto, o intercepto do modelo passa a
## ser "a captura de um barco com ZERO pescadores", que não existe e
## deixa a escala do intercepto sem sentido. Centrar não muda o ajuste,
## mas torna os coeficientes legíveis.
viagens$npesc_c <- as.numeric(scale(viagens$npesc, center = TRUE, scale = FALSE))

cat("\n===== DADOS PARA A MODELAGEM =====\n")
cat(sprintf("Viagens: %d | níveis de %s: %d | bancos: %d | barcos: %d | táticas: %d\n",
            nrow(viagens), fator_tempo, nlevels(viagens[[fator_tempo]]),
            nlevels(viagens$fbanco), nlevels(viagens$fbarco),
            nlevels(viagens$alvo)))
cat(sprintf("Zeros na resposta: %.1f%%  (P3 — decisivo para a escolha da distribuição)\n",
            100 * mean(viagens$captura == 0)))

## PROFUNDIDADE — decisão explícita.
## `prof` tem ~31% de ausentes. Se entrasse no modelo principal, o
## glmmTMB descartaria essas linhas e os modelos passariam a ser
## ajustados a conjuntos DIFERENTES de dados — o que invalidaria toda a
## comparação por AIC e por LRT. Então: `prof` fica FORA do conjunto
## principal e é testada à parte, no subconjunto completo, como
## sensibilidade. Imputar não é opção razoável aqui porque o ausente
## quase certamente não é aleatório (depende do amostrador/porto).
cat(sprintf("Profundidade ausente em %.0f%% das viagens -> fora do modelo\n",
            100 * mean(is.na(viagens$prof))))
cat("  principal; testada depois, isolada, no subconjunto completo.\n")

## =====================================================================
## 1) CONSTRUTOR DE FÓRMULAS
## ---------------------------------------------------------------------
## Montar a fórmula como string é o caminho seguro: `update()` NÃO
## remove `offset()` (ele não é um term.label) e `reformulate()` quebra o
## termo aleatório `(1|fbarco)`. Já custou bug antes.
## =====================================================================
monta_formula <- function(resposta, termos, aleatorio = TRUE,
                          offset_var = "ldias") {
  rhs <- paste(c(fator_tempo, termos,
                 if (aleatorio) "(1 | fbarco)",
                 if (!is.null(offset_var)) sprintf("offset(%s)", offset_var)),
               collapse = " + ")
  stats::as.formula(paste(resposta, "~", rhs))
}

## =====================================================================
## 2) ESTRUTURAS CANDIDATAS (E0-E5), todas com família Tweedie
## ---------------------------------------------------------------------
## A Tweedie é a distribuição de TRABALHO nesta etapa: fixamos a
## distribuição para comparar estruturas, e só depois, já na melhor
## estrutura, comparamos distribuições. Comparar tudo contra tudo (6
## estruturas x 5 distribuições) multiplicaria testes sem necessidade e
## tornaria o resultado dependente da ordem em que se olha.
## =====================================================================
termos_E <- list(
  E0 = character(0),
  E1 = c("fbanco"),
  E2 = c("fbanco", "ftipo", "npesc_c"),
  E3 = c("fbanco", "ftipo", "npesc_c", "alvo"),
  E4 = c("fbanco", "ftipo", "npesc_c", "alvo"),
  E5 = c("fbanco", "ftipo", "npesc_c", PCs)
)
aleat_E <- c(E0 = FALSE, E1 = FALSE, E2 = FALSE, E3 = FALSE, E4 = TRUE, E5 = TRUE)

cat("\n===== 1) ESTRUTURAS (Tweedie, offset = log dias) =====\n")
fits <- list()
for (nm in names(termos_E)) {
  t0 <- Sys.time()
  f  <- monta_formula("captura", termos_E[[nm]], aleatorio = aleat_E[[nm]])
  fits[[nm]] <- try(glmmTMB(f, family = tweedie(link = "log"), data = viagens),
                    silent = TRUE)
  ok <- !inherits(fits[[nm]], "try-error")
  cat(sprintf("  %-3s %-62s %s (%.1f min)\n", nm,
              paste(deparse(f), collapse = ""),
              if (ok) sprintf("AIC=%.1f", AIC(fits[[nm]])) else "FALHOU",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
fits <- fits[!vapply(fits, inherits, logical(1), "try-error")]

tab_est <- data.frame(
  modelo = names(fits),
  df     = vapply(fits, function(m) attr(logLik(m), "df"), numeric(1)),
  logLik = round(vapply(fits, function(m) as.numeric(logLik(m)), numeric(1)), 1),
  AIC    = round(vapply(fits, AIC, numeric(1)), 1),
  BIC    = round(vapply(fits, BIC, numeric(1)), 1), row.names = NULL)
tab_est$dAIC <- round(tab_est$AIC - min(tab_est$AIC), 1)
cat("\nTodas ajustadas às mesmas linhas e à mesma resposta -> AIC comparável:\n")
print(tab_est, row.names = FALSE)

## LRT sequencial: o ganho de cada termo, um a um.
cat("\nLRT sequencial (cada modelo contra o anterior, mais simples):\n")
seqm <- intersect(c("E0", "E1", "E2", "E3", "E4"), names(fits))
for (i in seq_len(length(seqm) - 1)) {
  a <- anova(fits[[seqm[i]]], fits[[seqm[i + 1]]])
  cat(sprintf("  %s -> %s : Chisq=%.1f, df=%d, p=%s\n", seqm[i], seqm[i + 1],
              a$Chisq[2], a$`Chi Df`[2],
              format.pval(a$`Pr(>Chisq)`[2], digits = 3, eps = 1e-16)))
}

## =====================================================================
## 3) MEDIDA DE ESFORÇO: dias contra horas no mar
## ---------------------------------------------------------------------
## Comparação legítima por AIC: mesma resposta, mesmas linhas, só muda o
## offset. `horas` costuma ganhar porque discrimina dentro do dia (uma
## viagem de 1 dia pode ter 10 h ou 28 h no mar), mas isso tem de ser
## verificado, não assumido.
## =====================================================================
cat("\n===== 2) QUAL MEDIDA DE ESFORÇO USAR =====\n")
melhor_E <- tab_est$modelo[which.min(tab_est$AIC)]
f_dias  <- monta_formula("captura", termos_E[[melhor_E]], aleat_E[[melhor_E]], "ldias")
f_horas <- monta_formula("captura", termos_E[[melhor_E]], aleat_E[[melhor_E]], "lhoras")
m_dias  <- fits[[melhor_E]]
m_horas <- try(glmmTMB(f_horas, family = tweedie(link = "log"), data = viagens),
               silent = TRUE)
if (!inherits(m_horas, "try-error")) {
  cat(sprintf("  offset log(dias)  : AIC = %.1f\n", AIC(m_dias)))
  cat(sprintf("  offset log(horas) : AIC = %.1f\n", AIC(m_horas)))
  usar_horas <- AIC(m_horas) < AIC(m_dias)
} else usar_horas <- FALSE
OFFSET <- if (usar_horas) "lhoras" else "ldias"
UNID   <- if (usar_horas) "t/hora no mar" else "t/dia de pesca"
cat(sprintf("  -> adotado: offset(%s)  [%s]\n", OFFSET, UNID))

## =====================================================================
## 4) SELEÇÃO BACKWARD COM LRT
## ---------------------------------------------------------------------
## Backward, não forward: um termo cujo efeito só aparece depois de
## ajustar outro é sistematicamente perdido no forward. Partimos do
## modelo cheio e removemos o termo com maior p, enquanto houver algum
## acima de alfa. O fator temporal não entra na seleção — ele É o índice,
## e removê-lo seria remover o objeto da análise.
## Alfa = 0,01 (e não 0,05) porque com muitas viagens o LRT detecta
## efeitos irrelevantes; o custo de manter covariável inútil é comer
## sinal de abundância (o dilema de Hinton & Maunder, 2003).
## =====================================================================
cat("\n===== 3) SELEÇÃO BACKWARD (LRT, alfa = 0,01) =====\n")
termos_sel <- termos_E$E4
m_atual <- glmmTMB(monta_formula("captura", termos_sel, TRUE, OFFSET),
                   family = tweedie(link = "log"), data = viagens)
repeat {
  if (length(termos_sel) == 0) break
  p_vals <- setNames(numeric(length(termos_sel)), termos_sel)
  for (tm in termos_sel) {
    m_red <- try(glmmTMB(monta_formula("captura", setdiff(termos_sel, tm), TRUE, OFFSET),
                         family = tweedie(link = "log"), data = viagens), silent = TRUE)
    p_vals[tm] <- if (inherits(m_red, "try-error")) 0 else
      anova(m_red, m_atual)[["Pr(>Chisq)"]][2]
  }
  cat("  p de remoção: ",
      paste(sprintf("%s=%.3g", names(p_vals), p_vals), collapse = "  "), "\n")
  if (max(p_vals) > 0.01) {
    pior <- names(which.max(p_vals))
    cat(sprintf("  -> remove `%s` (p=%.3g)\n", pior, max(p_vals)))
    termos_sel <- setdiff(termos_sel, pior)
    m_atual <- glmmTMB(monta_formula("captura", termos_sel, TRUE, OFFSET),
                       family = tweedie(link = "log"), data = viagens)
  } else { cat("  -> todos os termos retidos\n"); break }
}
f_fix <- monta_formula("captura", termos_sel, TRUE, OFFSET)
cat("Estrutura final: "); print(f_fix)

## Sensibilidade da profundidade, no subconjunto com registro (ver seção 0)
cat("\n--- profundidade (sensibilidade, subconjunto completo) ---\n")
v_prof <- viagens[!is.na(viagens$prof), ]
m_sp <- try(glmmTMB(monta_formula("captura", termos_sel, TRUE, OFFSET),
                    family = tweedie(link = "log"), data = v_prof), silent = TRUE)
m_cp <- try(glmmTMB(monta_formula("captura", c(termos_sel, "prof"), TRUE, OFFSET),
                    family = tweedie(link = "log"), data = v_prof), silent = TRUE)
if (!inherits(m_sp, "try-error") && !inherits(m_cp, "try-error")) {
  a <- anova(m_sp, m_cp)
  cat(sprintf("  n = %d viagens | LRT da profundidade: p = %s | dAIC = %.1f\n",
              nrow(v_prof), format.pval(a$`Pr(>Chisq)`[2], digits = 3),
              AIC(m_cp) - AIC(m_sp)))
  cat("  (mesmo se significativa, NÃO entra no modelo principal: usá-la\n")
  cat("   custaria 31% das viagens e impediria comparar com os demais)\n")
}

## ---------------------------------------------------------------------
## Guarda contra modelo mais complexo que os dados suportam.
## Os modelos ajustados só às capturas POSITIVAS trabalham com uma fração
## pequena das viagens (em 2024, ~13%). Manter ali um fator com 20 níveis
## de banco significaria estimar 20 coeficientes com pouquíssimas
## observações por nível: o ajuste não converge ou converge para lixo.
## Regra: um termo categórico só entra se houver pelo menos
## MIN_POR_NIVEL observações por nível.
## ---------------------------------------------------------------------
MIN_POR_NIVEL <- 10
termos_viaveis <- function(dados, termos, min_por_nivel = MIN_POR_NIVEL) {
  manter <- vapply(termos, function(tm) {
    x <- dados[[tm]]
    if (is.factor(x) || is.character(x)) {
      k <- length(unique(as.character(x[!is.na(x)])))
      nrow(dados) >= min_por_nivel * k
    } else TRUE
  }, logical(1))
  if (any(!manter))
    cat(sprintf("  [ajuste] termos retirados por amostra insuficiente: %s\n",
                paste(termos[!manter], collapse = ", ")))
  termos[manter]
}

## =====================================================================
## 5) DISTRIBUIÇÕES CANDIDATAS
## ---------------------------------------------------------------------
## Ver a discussão sobre grupos de comparabilidade no cabeçalho.
## =====================================================================
cat("\n===== 4) DISTRIBUIÇÕES =====\n")
dist_fits <- list(); dist_grupo <- c()

## D1 — Tweedie: Poisson composta com Gamma. Tem massa em zero e cauda
##      contínua positiva num modelo só. É a escolha natural aqui.
dist_fits$D1_tweedie <- m_atual;                       dist_grupo["D1_tweedie"] <- "A"

## D2 — Hurdle-Gamma: um binomial para a presença e um Gamma para a
##      magnitude. Em glmmTMB, `ziGamma` + `ziformula` é literalmente um
##      hurdle (a Gamma não tem massa em zero). Deixa a probabilidade de
##      zero depender das covariáveis de direcionamento, o que é
##      conceitualmente o que acontece na pescaria.
zi_f <- as.formula(paste("~", paste(c(fator_tempo,
                                      intersect("alvo", termos_sel)), collapse = " + ")))
dist_fits$D2_hurdle_gamma <- try(
  glmmTMB(f_fix, ziformula = zi_f, family = ziGamma(link = "log"), data = viagens),
  silent = TRUE);                                      dist_grupo["D2_hurdle_gamma"] <- "A"

## D3/D4 — só as viagens com captura positiva. São o componente de
##      magnitude do delta clássico (D3, lognormal) e sua alternativa
##      Gamma (D4). Ajustados a 13% das linhas -> grupo B.
viagens_pos <- viagens[viagens$captura > 0, ]
viagens_pos$lcpue <- log(viagens_pos$captura /
                           (if (OFFSET == "lhoras") viagens_pos$horas else viagens_pos$dias))
viagens_pos[[fator_tempo]] <- droplevels(viagens_pos[[fator_tempo]])
for (f in c("fbanco", "fbarco", "ftipo", "alvo"))
  if (f %in% names(viagens_pos)) viagens_pos[[f]] <- droplevels(viagens_pos[[f]])
cat(sprintf("  positivas: %d viagens (%.1f%% do total)\n", nrow(viagens_pos),
            100 * nrow(viagens_pos) / nrow(viagens)))
termos_pos <- termos_viaveis(viagens_pos, termos_sel)
f_bin <- monta_formula("pos_mac", termos_sel, TRUE, NULL)
f_pos <- monta_formula("lcpue",   termos_pos, TRUE, NULL)
m_bin <- try(glmmTMB(f_bin, family = binomial(), data = viagens), silent = TRUE)
m_pos <- try(glmmTMB(f_pos, family = gaussian(), data = viagens_pos), silent = TRUE)
if (!inherits(m_pos, "try-error")) {
  dist_fits$D3_lognormal_pos <- m_pos;                 dist_grupo["D3_lognormal_pos"] <- "B"
}
dist_fits$D4_gamma_pos <- try(
  glmmTMB(monta_formula("captura", termos_pos, TRUE, OFFSET),
          family = Gamma(link = "log"), data = viagens_pos),
  silent = TRUE);                                      dist_grupo["D4_gamma_pos"] <- "B"

## D5 — gaussiana em log(CPUE + c): a prática que queremos justamente
##      evitar, incluída para MEDIR o estrago. Responde numa escala
##      diferente -> grupo C, sozinha. Sem offset: a resposta já é taxa.
c_add <- 0.5 * min(viagens$captura[viagens$captura > 0])
viagens$lcpue_c <- log((viagens$captura +
                          c_add) / (if (OFFSET == "lhoras") viagens$horas else viagens$dias))
dist_fits$D5_log_mais_c <- try(
  glmmTMB(monta_formula("lcpue_c", termos_sel, TRUE, NULL),
          family = gaussian(), data = viagens),
  silent = TRUE);                                      dist_grupo["D5_log_mais_c"] <- "C"

dist_fits <- dist_fits[!vapply(dist_fits, inherits, logical(1), "try-error")]
dist_grupo <- dist_grupo[names(dist_fits)]

tab_dist <- data.frame(
  modelo = names(dist_fits),
  grupo  = dist_grupo,
  resposta = c(D1_tweedie = "captura (t)", D2_hurdle_gamma = "captura (t)",
               D3_lognormal_pos = "log CPUE | >0", D4_gamma_pos = "captura (t) | >0",
               D5_log_mais_c = "log(CPUE + c)")[names(dist_fits)],
  n      = vapply(dist_fits, function(m) nrow(model.frame(m)), numeric(1)),
  AIC    = round(vapply(dist_fits, AIC, numeric(1)), 1),
  BIC    = round(vapply(dist_fits, BIC, numeric(1)), 1), row.names = NULL)
# dAIC calculado DENTRO de cada grupo — comparar entre grupos não faz sentido
tab_dist$dAIC_no_grupo <- ave(tab_dist$AIC, tab_dist$grupo,
                              FUN = function(x) round(x - min(x), 1))
cat("\nA / B / C são os grupos de comparabilidade (ver cabeçalho).\n")
cat("Compare AIC apenas DENTRO do mesmo grupo:\n")
print(tab_dist, row.names = FALSE)

## Escolha do modelo principal: dentro do grupo A, que é o único que
## responde à pergunta "qual distribuição descreve melhor a captura,
## zeros incluídos".
gA <- tab_dist$modelo[tab_dist$grupo == "A"]
melhor_dist <- gA[which.min(tab_dist$AIC[tab_dist$grupo == "A"])]
m_final <- dist_fits[[melhor_dist]]
cat(sprintf("\nModelo principal (menor AIC no grupo A): %s\n", melhor_dist))
cat("Confirmar no diagnóstico antes de aceitar — AIC mede ajuste, não\n")
cat("adequação da distribuição.\n")
if (melhor_dist == "D2_hurdle_gamma")
  cat("[ATENÇÃO] Em hurdle, o emmeans devolve só o componente condicional.\n",
      "         Use o cenário S2d (delta) para o índice.\n")

## =====================================================================
## 6) DIAGNÓSTICO DE RESÍDUOS (DHARMa)
## ---------------------------------------------------------------------
## Em GLMM não-gaussiano o resíduo de Pearson engana: a relação
## média-variância não é constante e o gráfico "parece" ruim mesmo
## quando o modelo está certo. O DHARMa simula do modelo ajustado e
## transforma os resíduos para a escala uniforme, onde a leitura é a
## mesma para qualquer distribuição.
## O que cada teste responde:
##   KS         a distribuição assumida está certa?
##   dispersão  há sobre/subdispersão?
##   outliers   há mais extremos do que o modelo consegue gerar?
##   quantis    a variância é homogênea ao longo do predito?
##              (é o teste de homocedasticidade aqui)
##   zeros      o modelo gera a quantidade certa de zeros?
## Com n grande, p pequeno aparece por desvio trivial — olhe os gráficos.
## =====================================================================
cat("\n===== 5) DIAGNÓSTICO =====\n")
diagnostica <- function(m, nome, n_sim = 250) {
  set.seed(1)
  r <- simulateResiduals(m, n = n_sim, plot = FALSE)
  u <- testUniformity(r, plot = FALSE); d <- testDispersion(r, plot = FALSE)
  o <- testOutliers(r, plot = FALSE)
  q <- try(testQuantiles(r, plot = FALSE), silent = TRUE)
  z <- try(testZeroInflation(r, plot = FALSE), silent = TRUE)
  data.frame(modelo = nome, KS_p = signif(u$p.value, 3),
             disp_p = signif(d$p.value, 3),
             disp_ratio = signif(as.numeric(d$statistic), 3),
             outlier_p = signif(o$p.value, 3),
             quantis_p = if (inherits(q, "try-error")) NA else signif(q$p.value, 3),
             zeros_p = if (inherits(z, "try-error")) NA else signif(z$p.value, 3),
             row.names = NULL)
}
diag_tab <- do.call(rbind, lapply(names(dist_fits), function(nm)
  tryCatch(diagnostica(dist_fits[[nm]], nm), error = function(e)
    data.frame(modelo = nm, KS_p = NA, disp_p = NA, disp_ratio = NA,
               outlier_p = NA, quantis_p = NA, zeros_p = NA))))
print(diag_tab, row.names = FALSE)

png("diag_residuos.png", width = 26, height = 13, res = 300,
    antialias = "cleartype", units = "cm")
set.seed(1); r_fin <- simulateResiduals(m_final, n = 250, plot = FALSE)
plot(r_fin)
dev.off()
cat("PNG salvo: diag_residuos.png\n")

## Normalidade só faz sentido onde a premissa existe: no componente
## lognormal. Não se testa normalidade de resíduo de Tweedie.
if ("D3_lognormal_pos" %in% names(dist_fits)) {
  res <- residuals(dist_fits$D3_lognormal_pos)
  am <- sample(res, min(5000, length(res)))
  cat(sprintf("\nShapiro-Wilk no componente lognormal (n=%d): p = %.3g\n",
              length(am), shapiro.test(am)$p.value))
}

## =====================================================================
## 7) EXTRAÇÃO DO ÍNDICE (emmeans)
## ---------------------------------------------------------------------
## É aqui que "entrar como fator" se paga: o emmeans calcula a média
## marginal do fator temporal MEDIANDO os demais fatores — isto é, a
## taxa de captura esperada num banco médio, com tripulação média, numa
## mistura média de táticas. O que sobra é o efeito do tempo.
##
## Dois detalhes que mudam o resultado e passam despercebidos:
##  (1) `offset = 0`. Sem isso o emmeans usa a MÉDIA do offset e a série
##      sai por "viagem média", não por dia/hora. Com isso, sai em
##      unidade de esforço, que é o que o JABBA espera.
##  (2) `weights`. "equal" dá o mesmo peso a cada nível (padrão, e o
##      certo para um índice: queremos a média sobre estratos, não sobre
##      as viagens que por acaso foram amostradas). "proportional"
##      pondera pelo n observado e devolve parte do desbalanceamento que
##      estamos justamente tentando tirar.
## =====================================================================
extrai_indice <- function(m, nome, lognormal = FALSE, pesos = "equal",
                          usa_offset = TRUE) {
  esp <- stats::as.formula(paste("~", fator_tempo))
  args <- list(object = m, specs = esp, weights = pesos)
  if (usa_offset) args$offset <- 0
  s <- as.data.frame(summary(do.call(emmeans, args)))
  mu <- s$emmean; se <- s$SE
  # retrotransformação: ligação log de GLM -> exp(); gaussiana ajustada
  # em log(y) -> exp(mu + sigma^2/2). Trocar uma pela outra é erro comum.
  idx <- if (lognormal) exp(mu + 0.5 * sigma(m)^2) else exp(mu)
  data.frame(tempo = as.numeric(as.character(s[[fator_tempo]])),
             cenario = nome, indice_bruto = idx, se_log = se,
             cv = sqrt(exp(se^2) - 1), row.names = NULL)
}
## Índice delta: P(captura > 0) x média das positivas. Variância
## combinada pelo método delta — em escala log, d log(p)/d logit(p) = 1-p.
indice_delta <- function(m_bin, m_pos, nome, pesos = "equal") {
  esp <- stats::as.formula(paste("~", fator_tempo))
  sb <- as.data.frame(summary(emmeans(m_bin, esp, weights = pesos)))
  sp <- as.data.frame(summary(emmeans(m_pos, esp, weights = pesos)))
  p <- plogis(sb$emmean); se_logp <- (1 - p) * sb$SE
  idx <- p * exp(sp$emmean + 0.5 * sigma(m_pos)^2)
  se_tot <- sqrt(se_logp^2 + sp$SE^2)
  data.frame(tempo = as.numeric(as.character(sb[[fator_tempo]])), cenario = nome,
             indice_bruto = idx, se_log = se_tot, cv = sqrt(exp(se_tot^2) - 1),
             cv_sem_binomial = sqrt(exp(sp$SE^2) - 1), row.names = NULL)
}
normaliza <- function(d) { d$indice <- d$indice_bruto / mean(d$indice_bruto); d }

cat("\n===== 6) ÍNDICES POR CENÁRIO =====\n")

## S1 — nominal: captura agregada / esforço agregado, sem modelo nenhum.
esf_col <- if (OFFSET == "lhoras") "horas" else "dias"
S1 <- data.frame(tempo = ser$tempo, cenario = "S1 nominal",
                 indice_bruto = ser$cap_macarellus / ser[[esf_col]],
                 se_log = NA_real_, cv = NA_real_)
# CV empírico: erro-padrão relativo da CPUE entre viagens do mesmo período
cv_emp <- tapply(viagens$captura / viagens[[esf_col]], viagens$tempo,
                 function(x) sd(x) / (mean(x) * sqrt(length(x))))
S1$cv <- as.numeric(cv_emp[as.character(S1$tempo)])
S1 <- normaliza(S1)

## S0 — modelo SEM a covariável de tática. A diferença S0 - S2 isola o
##      efeito do direcionamento: é o "influence plot" de Bentley et al.
##      reduzido ao termo que interessa.
E_sem_alvo <- setdiff(termos_sel, c("alvo", PCs))
m_S0 <- glmmTMB(monta_formula("captura", E_sem_alvo, TRUE, OFFSET),
                family = tweedie(link = "log"), data = viagens)
S0 <- normaliza(extrai_indice(m_S0, "S0 sem tática"))

## S2 — cenário principal: modelo selecionado, com tática discreta.
S2 <- normaliza(extrai_indice(m_final, "S2 corrigida (tática discreta)"))

## S2d — mesma correção pelo caminho delta-lognormal (padrão ICCAT).
##       Se a forma da série mudar muito entre S2 e S2d, isso é
##       incerteza estrutural e tem de ser reportada como tal.
S2d <- NULL
if (!inherits(m_bin, "try-error") && !inherits(m_pos, "try-error")) {
  S2d <- normaliza(indice_delta(m_bin, m_pos, "S2d delta-lognormal"))
  cat(sprintf("S2d: CV médio com binomial = %.3f | sem = %.3f (a versão sem subestima)\n",
              mean(S2d$cv), mean(S2d$cv_sem_binomial)))
  S2d$cv_sem_binomial <- NULL
}

## S2b — tática contínua (H2).
S2b <- if ("E5" %in% names(fits))
  normaliza(extrai_indice(fits$E5, "S2b tática contínua (PCs)")) else NULL

## S3 — esforço dirigido (H3): só as viagens da tática com maior fração
##      de cavala. `alvo_cavala` vem da parte 02 e é escolhido pelo
##      CENTRÓIDE, não pelo nome do grupo.
S3 <- NULL
if (exists("alvo_cavala") && alvo_cavala %in% levels(viagens$alvo)) {
  v3 <- viagens[viagens$alvo == alvo_cavala, ]
  v3[[fator_tempo]] <- droplevels(v3[[fator_tempo]])
  cobre <- nlevels(v3[[fator_tempo]]) >= 0.7 * nlevels(viagens[[fator_tempo]])
  if (cobre) {
    for (f in c("fbanco", "fbarco", "ftipo")) v3[[f]] <- droplevels(v3[[f]])
    cat(sprintf("  S3: %d viagens da tática '%s'\n", nrow(v3), alvo_cavala))
    m3 <- try(glmmTMB(monta_formula("captura", termos_viaveis(v3, E_sem_alvo),
                                    TRUE, OFFSET),
                      family = tweedie(link = "log"), data = v3), silent = TRUE)
    if (!inherits(m3, "try-error"))
      S3 <- normaliza(extrai_indice(m3, "S3 esforço dirigido"))
  } else {
    cat("[AVISO] a tática da cavala não cobre períodos suficientes;\n")
    cat("        S3 fica sem estimativa — limitação esperada do subsetting.\n")
  }
}

indices <- do.call(rbind, Filter(Negate(is.null), list(S1, S0, S2, S2d, S2b, S3)))
cat("\nÍndices (média 1) e CV:\n")
print(transform(indices[, c("tempo", "cenario", "indice", "cv")],
                indice = round(indice, 3), cv = round(cv, 3)), row.names = FALSE)

## =====================================================================
## 8) TESTE DAS HIPÓTESES
## =====================================================================
cat("\n===== 7) HIPÓTESES =====\n")
r_S0S2 <- cor(S0$indice, S2$indice)
cat(sprintf("H0 — correlação entre índice com e sem tática: %.3f\n", r_S0S2))
cat(sprintf("     %s\n", if (r_S0S2 > 0.98)
  "praticamente idênticos: a correção não está mudando nada (reportar!)" else
    "a tática desloca o índice de forma relevante"))
amp <- function(d) max(d$indice) / min(d$indice)
cat(sprintf("H1 — amplitude (máx/mín) do índice: nominal %.2f | corrigida %.2f\n",
            amp(S1), amp(S2)))
cat(sprintf("     correlação nominal x corrigida: %.3f\n",
            cor(S1$indice[match(S2$tempo, S1$tempo)], S2$indice)))
if (!is.null(S2b))
  cat(sprintf("H2 — discreta x contínua: r = %.3f | dAIC (contínua - discreta) = %.1f\n",
              cor(S2$indice, S2b$indice), AIC(fits$E5) - AIC(m_final)))

## =====================================================================
## 9) FIGURA DOS ÍNDICES
## =====================================================================
cores_cen <- c("S1 nominal" = COR_AUX, "S0 sem tática" = COR_NEU,
               "S2 corrigida (tática discreta)" = COR_MAC,
               "S2d delta-lognormal" = COR_S2D,
               "S2b tática contínua (PCs)" = COR_S2,
               "S3 esforço dirigido" = COR_S3)
png("indices_cenarios.png", width = 26, height = 14, res = 300,
    antialias = "cleartype", units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
cens <- unique(indices$cenario)
plot(NA, xlim = range(indices$tempo), ylim = c(0, max(indices$indice, na.rm = TRUE) * 1.12),
     xlab = rotulo_tempo, ylab = "Índice relativo (média = 1)",
     main = "A. Cenários de índice")
for (cen in cens) {
  d <- indices[indices$cenario == cen, ]
  lines(d$tempo, d$indice, lwd = 2.4, col = cores_cen[cen])
}
legend("topright", cens, col = cores_cen[cens], lwd = 2.3, bty = "n", cex = 0.66)

lo <- S2$indice * exp(-1.96 * S2$se_log); hi <- S2$indice * exp(1.96 * S2$se_log)
plot(S2$tempo, S2$indice, type = "n", ylim = c(0, max(hi, na.rm = TRUE) * 1.05),
     xlab = rotulo_tempo, ylab = "Índice (média = 1)",
     main = "B. Índice corrigido com IC 95%")
polygon(c(S2$tempo, rev(S2$tempo)), c(lo, rev(hi)),
        col = adjustcolor(COR_MAC, 0.18), border = NA)
lines(S2$tempo, S2$indice, lwd = 2.6, col = COR_MAC)
lines(S1$tempo, S1$indice, lwd = 2, col = COR_AUX, lty = 2)
legend("topright", c("corrigida (IC 95%)", "nominal"), col = c(COR_MAC, COR_AUX),
       lwd = c(2.6, 2), lty = c(1, 2), bty = "n", cex = 0.72)
par(op); dev.off()
cat("\nPNG salvo: indices_cenarios.png\n")

## =====================================================================
## 10) EXPORTAÇÃO
## ---------------------------------------------------------------------
## Formato do JABBA: uma coluna de tempo e uma coluna por índice, mais
## uma tabela equivalente de CV. Série com média 1 (o JABBA estima q).
## Piso de CV em 0,20: o CV que sai do modelo é de processo estatístico
## e ignora erro de processo, erro de reporte e a incerteza da própria
## inferência de tática. Entregar CV de 0,04 ao JABBA faria o modelo
## confiar no índice muito mais do que ele merece.
## =====================================================================
saida_tempo <- sort(unique(indices$tempo))
jabba_idx <- data.frame(tempo = saida_tempo)
jabba_cv  <- data.frame(tempo = saida_tempo)
for (cen in c("S1 nominal", "S2 corrigida (tática discreta)")) {
  d <- indices[indices$cenario == cen, ]
  nm <- if (grepl("nominal", cen)) "cpue_nominal" else "cpue_corrigida"
  jabba_idx[[nm]] <- d$indice[match(saida_tempo, d$tempo)]
  jabba_cv[[nm]]  <- d$cv[match(saida_tempo, d$tempo)]
}
names(jabba_idx)[1] <- names(jabba_cv)[1] <- if (fator_tempo == "fano") "Yr" else "Mes"
jabba_cv[, -1] <- lapply(jabba_cv[, -1, drop = FALSE],
                         function(x) pmax(x, 0.20, na.rm = TRUE))

write.csv(jabba_idx, "jabba_indices_macarellus.csv", row.names = FALSE)
write.csv(jabba_cv,  "jabba_cv_macarellus.csv", row.names = FALSE)
write.csv(indices,   "indices_todos_cenarios.csv", row.names = FALSE)
if (tem_writexl)
  writexl::write_xlsx(list(indices = jabba_idx, cv = jabba_cv, todos = indices,
                           estruturas = tab_est, distribuicoes = tab_dist,
                           diagnostico = diag_tab),
                      path = "padronizacao_cpue_macarellus.xlsx")
cat("Arquivos salvos: jabba_indices_macarellus.csv, jabba_cv_macarellus.csv,\n")
cat("                 indices_todos_cenarios.csv",
    if (tem_writexl) ", padronizacao_cpue_macarellus.xlsx\n" else "\n")

## =====================================================================
## 11) COMO REPORTAR
## =====================================================================
cat("\n===== COMO REPORTAR =====\n")
cat("1. S1 e S2 são CENÁRIOS ALTERNATIVOS de entrada no JABBA, não uma\n")
cat("   série certa e outra errada (Hoyle et al. 2024, boas práticas 17-18).\n")
cat("2. Reportar as tabelas de seleção e de diagnóstico, inclusive quando\n")
cat("   os pressupostos forem violados.\n")
cat("3. Declarar o que a tática é: uma variável INFERIDA da composição da\n")
cat("   captura, não observada. A incerteza dessa inferência não está no CV.\n")
cat("4. A profundidade ficou fora do modelo principal por ter 31% de\n")
cat("   ausentes; o efeito dela está reportado como sensibilidade.\n")
if (exists("UM_ANO_SO") && UM_ANO_SO) {
  cat("\n5. ATENÇÃO: com um único ano, o que saiu aqui é um índice MENSAL\n")
  cat("   de 2024. Ele NÃO é um índice de abundância interanual e não\n")
  cat("   deve entrar no JABBA como tal. Serve para validar o pipeline e\n")
  cat("   para caracterizar a pescaria. O índice de verdade depende da\n")
  cat("   série histórica.\n")
}

















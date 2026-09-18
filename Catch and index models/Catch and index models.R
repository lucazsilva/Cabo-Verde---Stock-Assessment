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
#                                                                                                        #  
#         Seleção: backward a partir do modelo cheio, LRT para modelos aninhados                         #  
#         (glmmTMB ajusta por ML, então o LRT dos efeitos fixos é válido),                               #
#         com AIC e BIC como critérios complementares. Distribuição escolhida                            #  
#         por AIC + resíduos simulados (DHARMa), não por AIC sozinho.                                    #      
#         Índice: médias marginais do fator `ano` via emmeans, com offset = 0                            #  
#         (taxa por dia), retrotransformadas com correção de viés, série                                 #  
#         normalizada para média 1 (convenção de entrada do JABBA) e CV anual                            #
#         derivado do erro-padrão na escala do preditor linear.                                          #  
#                                                                                                        #
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


#===============================================================================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde Frota industrial — dados do IMar, série 2019-2025
# (1) A série agora é 2019-2025 (7 anos completos). O efeito de ANO — que É o índice de abundância — passa a ser estimável. A parte 02 troca
#     sozinha o fator temporal de mês para ano.
# (2) PROFUNDIDADE foi REMOVIDA do pipeline. Ela é  lida e o diagnóstico é impresso, mas não vai para a tabela final:
#     28% das viagens têm 0 (= ausência de registro, não zero metro) e o ausente não é aleatório. Ver decisão L6.
# (3) MÊS vira TRIMESTRE como fator sazonal (decisão L9) e TRIPULAÇÃO vira FATOR em classes (decisão L10). Motivo em cada decisão.
# (4) A identidade da embarcação passa a ser o CÓDIGO e não o nome  (decisão L11) — Os nomes se repetem entre barcos
#     diferentes que operam ao mesmo tempo.
# (5) `Tipo_embarcacao` foi INVESTIGADO e DESCARTADO como covariável (decisão L12).
#================================================================================================================================
# --------------------------------------------------------------------------------------------------------
# DECISÕES DE LEITURA E LIMPEZA (cada uma comentada no código abaixo)
# --------------------------------------------------------------------------------------------------------
# L1. Encoding: o CSV vem em latin1 (acentos de "S. ANTÃO"). Lemos com readLines + iconv para UTF-8 
#em vez de `fileEncoding`, porque `fileEncoding` falha dependendo do locale da máquina.
# L2. `quote = ""`: o arquivo não usa aspas; deixar o padrão faz o R interpretar apóstrofos de nomes próprios como aspas e truncar.
# L3. `read.csv2`: separador ";" e decimal com VÍRGULA ("32,92").
# L4. Todo campo de texto vem preenchido com espaços à direita (largura fixa). Sem `trimws` 
# em TUDO, "DECAPTERUS MACARELLUS      " nunca casa com "DECAPTERUS MACARELLUS" e a espécie some da análise.
# L5. FILTRAMOS PARA REDE DE CERCO. É a arte que captura a cavala preta. Manter linha-de-mão e
#     covos no denominador acrescentaria esforço que nunca teve chance de pegar cavala 
# — é o problema P1 (esforço não específico)  ampliado de graça. A arte é, além disso, o fator de capturabilidade
#     mais óbvio que existe.
# L6. PROFUNDIDADE SAI. `Profundidade_pesca_engenho1 == 0` é ausência de registro, não "zero metros" 
#— 28% das viagens. Como o ausente depende do amostrador/porto (não é aleatório), imputar enviesaria;
#     e usar só o subconjunto completo custaria 28% das viagens e  impediria comparar modelos 
#  por AIC (conjuntos diferentes de linhas). Fica de fora, e o diagnóstico é impresso para o texto.
# L7. Local de pesca: usamos `Nome_banco_pesca` (onde se pescou), não o porto/ilha de 
# desembarque (onde se descarregou). São coisas diferentes e só a primeira é covariável de densidade. 
# A ilha de desembarque é guardada só para descrição — NÃO entra no modelo, porque a cobertura 
# dela muda ao longo da série (S. Nicolau só aparece a partir de 2021), o que a deixa confundida com o ano.
# L8. Bancos com menos de MIN_VIAG_BANCO viagens viram "OUTROS": nível com 3 viagens em 7 anos não 
# estima coeficiente, só instabiliza.
# L9. MÊS -> TRIMESTRE. A sazonalidade da cavala é forte (presença cai de ~23% em janeiro para ~7% em agosto) 
# e precisa estar no modelo. Mas 12 níveis de mês, cruzados com 7 anos e dezenas de bancos, deixam
#     muitas células quase vazias e o coeficiente vira ruído. Trimestre mantém a forma sazonal 
# com 4 níveis bem povoados.
# L10. TRIPULAÇÃO em CLASSES, não contínua. `Numero_pescadores` é proxy de poder de pesca. Entrar 
# como contínua impõe relação monotônica e linear na escala do log, o que não há motivo para supor; 
# e a cauda (1-9 e 44-50 pescadores) puxa a reta. Classes deixam a forma aparecer e são robustas 
# a esses extremos. Cortes em <=12 / 13-15 /16-17 / >=18 (quartis aproximados, ~19/29/30/23% das viagens).
# L11. EMBARCAÇÃO = CÓDIGO (`Embarcacao`), não nome. Na varredura, 21  nomes aparecem com 2-3 códigos diferentes 
# e as séries se SOBREPÕEM no tempo (ex.: "BARBARA" código 18 de 2019-2025 e código
#     3413 de 2020-2025, simultâneos). São barcos distintos com o mesmo nome. Usar o nome fundiria 
# históricos de barcos diferentes; o código é 1:1 com o nome (nenhum código tem dois nomes).
# L12. `Tipo_embarcacao` (I / S) NÃO É USADO. Parecia um estrato de frota, mas a varredura 
# mostrou três coisas: (a) 30 dos 34 barcos que aparecem como "S" também aparecem como "I" — é o MESMO barco
#   com dois códigos em viagens diferentes; (b) "S" praticamente só existe em 2019 (452 viagens) e 
# some depois (12, 7, 1, 0, 0, 0), ou seja, está quase perfeitamente confundido com o ano — justo o ano
# de maior captura de cavala; (c) os perfis de I e S em 2019 são iguais (mesmos dias no mar, 
#mesma tripulação, mesma captura média). É um código administrativo de registro, não um tipo de embarcação.
#Usá-lo como covariável roubaria sinal do efeito de ano.
# L13. Viagens com MAIS DE UMA ARTE (43 viagens) são descartadas: os dias de mar delas cobrem também 
#o que foi pescado com outra arte, então o esforço atribuído ao cerco ficaria inflado.
#======================================================================

## =====================================================================
## 0) PARÂMETROS 
## =====================================================================
ARQUIVO        <- "INDUSTRIAL_2019_2025_atualizado_17.09.2026.csv"
ARTE_ALVO      <- "REDE DE CERCO"        # decisão L5
MIN_VIAG_BANCO <- 100                    # decisão L8 (45 bancos + "OUTROS")
ESPECIE_FOCO   <- "DECAPTERUS MACARELLUS"
CORTES_NPESC   <- c(0, 12, 15, 17, Inf)  # decisão L10

## =====================================================================
## 1) LEITURA ROBUSTA (decisões L1-L4)
## =====================================================================
le_imar <- function(caminho) {
  linhas <- readLines(caminho, warn = FALSE)
  # L1: latin1 -> UTF-8. `sub="?"` evita erro fatal se houver byte inválido.
  linhas <- iconv(linhas, from = "latin1", to = "UTF-8", sub = "?")
  # L2 e L3: sem aspas, separador ";", decimal ","
  d <- read.csv2(text = linhas, quote = "", stringsAsFactors = FALSE,
                 strip.white = TRUE)
  names(d) <- trimws(names(d))
  # L4: o arquivo é de largura fixa disfarçada de CSV — sobra espaço à
  # direita em TODO campo de texto. Sem este laço, nenhuma comparação de
  # string funciona.
  for (j in seq_along(d)) if (is.character(d[[j]])) d[[j]] <- trimws(d[[j]])
  d
}

bruto <- do.call(rbind, lapply(ARQUIVO, le_imar))

cat("======================================================\n")
cat("LEITURA DO ARQUIVO DO IMar\n")
cat("======================================================\n")
cat(sprintf("Arquivo(s): %s\n", paste(ARQUIVO, collapse = ", ")))
cat(sprintf("Linhas lidas (1 por espécie por amostragem): %d\n", nrow(bruto)))
cat(sprintf("Colunas: %d\n", ncol(bruto)))

## =====================================================================
## 2) RELATÓRIO DE QUALIDADE — ANTES de filtrar qualquer coisa
## ---------------------------------------------------------------------
## Princípio: o que se descarta tem de ser sempre visível e justificado.
## Esta seção é o que vai virar o parágrafo de "dados" do artigo.
## =====================================================================
cat("\n--- qualidade dos campos-chave (sobre TODAS as linhas) ---\n")
qual <- function(x, nome) {
  ausente <- is.na(x) | (is.character(x) & !is.na(x) & x == "")
  zero    <- is.numeric(x) & !is.na(x) & x == 0
  cat(sprintf("  %-30s %5.1f%% ausente  %5.1f%% zero\n", nome,
              100 * mean(ausente), 100 * mean(zero)))
}
qual(bruto$Amostragem,                  "Amostragem (id da viagem)")
qual(bruto$Quantidade,                  "Quantidade (kg)")
qual(bruto$Num_dias,                    "Num_dias")
qual(bruto$Num_horas,                   "Num_horas")
qual(bruto$Numero_pescadores,           "Numero_pescadores")
qual(bruto$Profundidade_pesca_engenho1, "Profundidade (0 = ausente)")
qual(bruto$Nome_banco_pesca,            "Nome_banco_pesca")

## Cobertura temporal: confirma que a série é contínua e quais anos estão
## completos. Um ano com poucos meses amostrados não é comparável com os
## demais e apareceria como "queda de abundância" que é só amostragem.
cat("\n--- cobertura temporal (linhas por ano) ---\n")
print(table(bruto$Ano, useNA = "ifany"))
cat("\n--- meses amostrados em cada ano ---\n")
print(tapply(bruto$Mes, bruto$Ano, function(x) length(unique(x[!is.na(x)]))))

## Linhas sem `Amostragem`: são registros de esforço SEM captura
## associada — o barco saiu e voltou, mas nenhuma espécie foi lançada na
## planilha. Não dá para saber se foi viagem sem captura (informação
## legítima e valiosa) ou falha de digitação, e elas não têm identificador
## para entrar na tabela de viagens. Saem da análise, com o número
## registrado. São 0,1% das linhas — irrelevante para o resultado.
sem_id <- sum(is.na(bruto$Amostragem) | bruto$Amostragem == "")
cat(sprintf("\nLinhas sem identificador de amostragem: %d (%.2f%%) -> descartadas\n",
            sem_id, 100 * sem_id / nrow(bruto)))
d <- bruto[!is.na(bruto$Amostragem) & bruto$Amostragem != "", ]

## Conferência de que `Amostragem` é mesmo um identificador de viagem:
## cada uma tem de ter UM ano, UMA embarcação e UMA data. Se falhar, todo
## o resto do script está construído sobre chão falso.
chk_ano  <- tapply(d$Ano,         d$Amostragem, function(x) length(unique(x)))
chk_emb  <- tapply(d$Embarcacao,  d$Amostragem, function(x) length(unique(x)))
chk_data <- tapply(d$Data_partida, d$Amostragem, function(x) length(unique(x)))
cat(sprintf("Amostragem é id único de viagem? ano:%s embarcação:%s data:%s\n",
            all(chk_ano == 1), all(chk_emb == 1), all(chk_data == 1)))
stopifnot(all(chk_ano == 1), all(chk_emb == 1), all(chk_data == 1))
cat(sprintf("Viagens (Amostragem distintas) no arquivo: %d\n",
            length(unique(d$Amostragem))))

## =====================================================================
## 3) FILTRO DE ARTE DE PESCA (decisões L5 e L13)
## =====================================================================
cat("\n--- viagens por arte de pesca (antes do filtro) ---\n")
viag_arte <- tapply(d$Amostragem, d$Nome_engenho, function(x) length(unique(x)))
print(sort(viag_arte, decreasing = TRUE))

## Onde a cavala preta é capturada? Esta tabela é a JUSTIFICATIVA
## empírica do filtro de arte — não filtramos por hábito, filtramos
## porque a captura está concentrada numa arte só.
foco <- d[toupper(d$Nome_cientifico) == ESPECIE_FOCO, ]
cat(sprintf("\n%s: %d viagens, %.1f t no total (todas as artes)\n", ESPECIE_FOCO,
            length(unique(foco$Amostragem)),
            sum(as.numeric(foco$Quantidade), na.rm = TRUE) / 1000))
cat("distribuição dessas viagens por arte:\n")
print(sort(tapply(foco$Amostragem, foco$Nome_engenho,
                  function(x) length(unique(x))), decreasing = TRUE))
cat(sprintf("-> mantendo apenas `%s` (decisão L5)\n", ARTE_ALVO))

## L13 — viagens que usaram mais de uma arte. O campo de esforço
## (`Num_dias`) é da VIAGEM inteira; se parte dela foi pescada com outra
## arte, o esforço que atribuiríamos ao cerco estaria inflado e a CPUE,
## deprimida. São poucas — saem.
n_artes <- tapply(d$Engenho, d$Amostragem, function(x) length(unique(x)))
viag_multi <- names(n_artes[n_artes > 1])
cat(sprintf("Viagens com mais de uma arte: %d -> descartadas (decisão L13)\n",
            length(viag_multi)))
d <- d[!d$Amostragem %in% viag_multi, ]

d <- d[d$Nome_engenho == ARTE_ALVO, ]
cat(sprintf("Linhas após o filtro de arte: %d (%d viagens)\n",
            nrow(d), length(unique(d$Amostragem))))

## =====================================================================
## 4) ESPÉCIES: quais viram coluna própria
## ---------------------------------------------------------------------
## Usamos NOME CIENTÍFICO, nunca o nome comum. No arquivo, "CAVALA PRETA"
## é Decapterus macarellus, mas "CAVALA BRANCA" é D. punctatus e "CAVALA
## DE RABO VERMELHA" é D. tabl — três espécies do mesmo gênero que seriam
## fundidas se a chave fosse o nome popular. Esse erro seria invisível no
## resultado final e catastrófico para a interpretação.
##
## As espécies que ganham coluna própria são as que sustentam a
## COMPOSIÇÃO da captura — a impressão digital da tática de pesca. O
## resto vai para "outras". O critério é cobrir a grande maioria do peso
## desembarcado com poucas colunas: espécie rara não informa tática, só
## acrescenta ruído e zeros à matriz de composição.
## =====================================================================
d$Quantidade <- as.numeric(d$Quantidade)
kg_sp <- tapply(d$Quantidade, toupper(d$Nome_cientifico), sum, na.rm = TRUE)
kg_sp <- sort(kg_sp, decreasing = TRUE)
cat("\n--- 15 maiores capturas por espécie no cerco (t) ---\n")
print(round(head(kg_sp, 15) / 1000, 1))
cat(sprintf("Espécies distintas registradas no cerco: %d\n", length(kg_sp)))
cat(sprintf("As 11 maiores cobrem %.1f%% do peso desembarcado\n",
            100 * sum(head(kg_sp, 11)) / sum(kg_sp)))

## O gênero Auxis aparece com três grafias ("Auxis sp", "Auxis rochei",
## "Auxis thazard") que são o mesmo complexo comercial (cachorrinha/
## judeu/melva) e não são separadas de forma consistente no desembarque.
## Tratá-las como espécies distintas fragmentaria o eixo de composição
## mais importante da pescaria. Por isso o gênero inteiro é uma coluna só
## — casado por PREFIXO, não por igualdade exata.
SP_FOCAIS <- c(
  macarellus = "DECAPTERUS MACARELLUS",
  katsuwonus = "KATSUWONUS PELAMIS",
  selar      = "SELAR CRUMENOPHTHALMUS",
  carangideo = "CARANX CRYSOS",
  punctatus  = "DECAPTERUS PUNCTATUS",
  sardinella = "SARDINELLA MADERENSIS",
  thunnus    = "THUNNUS ALBACARES",
  spicara    = "SPICARA MELANURUS",
  elagatis   = "ELAGATIS BIPINNULATA",
  apsilus    = "APSILUS FUSCUS"
)
SP_PREFIXO <- c(auxis = "AUXIS")   # casa Auxis sp / rochei / thazard

d$sp_col <- "outras"
for (k in names(SP_FOCAIS))
  d$sp_col[toupper(d$Nome_cientifico) == SP_FOCAIS[k]] <- k
for (k in names(SP_PREFIXO))
  d$sp_col[startsWith(toupper(d$Nome_cientifico), SP_PREFIXO[k])] <- k

cat(sprintf("\nEspécies com coluna própria: %d (+ 'outras', que ficou com %.1f%% do peso)\n",
            length(SP_FOCAIS) + length(SP_PREFIXO),
            100 * sum(d$Quantidade[d$sp_col == "outras"], na.rm = TRUE) /
              sum(d$Quantidade, na.rm = TRUE)))

## =====================================================================
## 5) LONGO -> LARGO: uma linha por VIAGEM
## ---------------------------------------------------------------------
## Captura convertida de kg para TONELADAS (o arquivo vem em kg).
## `tapply` cruzando viagem x coluna-de-espécie produz a matriz de
## capturas; o NA que sobra significa "espécie não registrada nesta
## viagem", que é captura ZERO — e esses zeros são o dado mais importante
## da análise (são eles que carregam a informação de direcionamento).
## =====================================================================
cap <- tapply(d$Quantidade / 1000, list(d$Amostragem, d$sp_col), sum)
cap[is.na(cap)] <- 0
cap <- as.data.frame(cap)
todas_col <- c(names(SP_FOCAIS), names(SP_PREFIXO), "outras")
for (k in setdiff(todas_col, names(cap))) cap[[k]] <- 0   # espécie ausente no período
cap <- cap[, todas_col]
names(cap) <- paste0("cap_", names(cap))
cap$Amostragem <- rownames(cap)

## As variáveis operacionais são constantes dentro da viagem (conferido
## na seção 2), então basta pegar a primeira ocorrência de cada uma.
primeiro <- function(col) tapply(d[[col]], d$Amostragem, function(x) x[1])
viagens <- data.frame(
  Amostragem   = names(primeiro("Ano")),
  ano          = as.integer(primeiro("Ano")),
  mes          = as.integer(primeiro("Mes")),
  barco_id     = primeiro("Embarcacao"),          # decisão L11
  barco_nome   = primeiro("Nome_embarcacao"),     # só para leitura humana
  tipo_emb     = primeiro("Tipo_embarcacao"),     # só para o diagnóstico L12
  npesc        = as.numeric(primeiro("Numero_pescadores")),
  dias         = as.numeric(primeiro("Num_dias")),
  horas        = as.numeric(primeiro("Num_horas")),
  prof         = as.numeric(primeiro("Profundidade_pesca_engenho1")),
  banco        = primeiro("Nome_banco_pesca"),
  ilha_desemb  = primeiro("Nome_ilha"),
  porto_desemb = primeiro("Nome_porto_desembarque"),
  data_saida   = as.Date(primeiro("Data_partida"), format = "%m/%d/%Y"),
  data_volta   = as.Date(primeiro("Data_chegada"), format = "%m/%d/%Y"),
  stringsAsFactors = FALSE
)
viagens <- merge(viagens, cap, by = "Amostragem", sort = FALSE)

## =====================================================================
## 6) DIAGNÓSTICOS QUE JUSTIFICAM AS DECISÕES L6, L11 e L12
## ---------------------------------------------------------------------
## Estes blocos não alteram os dados: eles IMPRIMEM a evidência que
## sustenta cada decisão do cabeçalho, para o texto do artigo e para
## quem for revisar o script.
## =====================================================================

## --- L6: por que a profundidade sai ----------------------------------
n_prof0 <- sum(viagens$prof == 0, na.rm = TRUE)
cat(sprintf("\n[L6] Profundidade: %d viagens (%.1f%%) com valor 0 = sem registro.\n",
            n_prof0, 100 * n_prof0 / nrow(viagens)))
cat("     Ausência não aleatória (depende do amostrador/porto). Usá-la\n")
cat("     custaria essas viagens e impediria comparar modelos por AIC\n")
cat("     (conjuntos de linhas diferentes). REMOVIDA do pipeline.\n")
viagens$prof <- NULL

## --- L11: nome de embarcação não serve como identidade ---------------
nomes_multi <- names(which(tapply(viagens$barco_id, viagens$barco_nome,
                                  function(x) length(unique(x))) > 1))
cat(sprintf("\n[L11] Nomes de embarcação usados por mais de um código: %d\n",
            length(nomes_multi)))
if (length(nomes_multi) > 0) {
  ex <- nomes_multi[1]
  sub_ex <- viagens[viagens$barco_nome == ex, ]
  cat(sprintf("      Exemplo '%s': códigos %s, anos %s\n", ex,
              paste(unique(sub_ex$barco_id), collapse = "/"),
              paste(range(sub_ex$ano), collapse = "-")))
  cat("      Como as séries se sobrepõem no tempo, são barcos distintos\n")
  cat("      com o mesmo nome. A identidade usada é o CÓDIGO.\n")
}

## --- L12: por que `Tipo_embarcacao` não pode ser covariável ----------
cat("\n[L12] Tipo_embarcacao x ano (viagens):\n")
print(table(viagens$ano, viagens$tipo_emb))
barcos_I <- unique(viagens$barco_nome[viagens$tipo_emb == "I"])
barcos_S <- unique(viagens$barco_nome[viagens$tipo_emb == "S"])
cat(sprintf("      Barcos que aparecem como 'S': %d | também como 'I': %d\n",
            length(barcos_S), length(intersect(barcos_I, barcos_S))))
cat("      => mesmo barco com dois códigos + código quase restrito a 2019\n")
cat("      => confundido com o ano. NÃO entra como covariável (decisão L12).\n")

## =====================================================================
## 7) VARIÁVEIS DERIVADAS E AGRUPAMENTOS (decisões L8, L9, L10)
## =====================================================================

## --- L9: trimestre ---------------------------------------------------
## A sazonalidade entra como trimestre, não como mês. Guardamos o mês
## também, porque a exploratória da parte 02 usa o perfil mensal para
## MOSTRAR que a sazonalidade existe e que 4 níveis a descrevem bem.
viagens$trimestre <- (viagens$mes - 1) %/% 3 + 1

## --- L10: tripulação em classes --------------------------------------
viagens$npesc_cat <- cut(viagens$npesc, breaks = CORTES_NPESC,
                         labels = c("<=12", "13-15", "16-17", ">=18"),
                         include.lowest = TRUE)
cat("\n[L10] Tripulação em classes (proxy de poder de pesca):\n")
print(table(viagens$npesc_cat, useNA = "ifany"))
cat(sprintf("      contínua: mediana %.0f, intervalo %.0f-%.0f pescadores\n",
            median(viagens$npesc, na.rm = TRUE), min(viagens$npesc, na.rm = TRUE),
            max(viagens$npesc, na.rm = TRUE)))

## --- L8: bancos raros agrupados --------------------------------------
tb_banco <- table(viagens$banco)
raros <- names(tb_banco[tb_banco < MIN_VIAG_BANCO])
viagens$banco_gr <- ifelse(viagens$banco %in% raros, "OUTROS", viagens$banco)
cat(sprintf("\n[L8] Bancos de pesca: %d distintos; %d com < %d viagens agrupados\n",
            length(tb_banco), length(raros), MIN_VIAG_BANCO))
cat(sprintf("     em 'OUTROS' (%.1f%% das viagens) -> %d níveis no modelo\n",
            100 * mean(viagens$banco %in% raros),
            length(unique(viagens$banco_gr))))

## --- esforço: duas medidas concorrentes ------------------------------
## `dias`  = Num_dias. Íntegro (sem NA, sem zero), inteiro, e é a unidade
##           das séries oficiais (comparável com o relatório da FAO).
## `horas` = Num_horas. Mais fino, mas com problemas: 60 viagens com
##           valor <= 0 e ~21% com horas > 24 x dias, o que mostra que os
##           dois campos não foram derivados de forma consistente.
## Conferência com as datas (partida/chegada) abaixo: `dias` bate em 100%
## das viagens, contando o dia de saída nas viagens de 1 dia.
viagens$dur_datas <- as.numeric(viagens$data_volta - viagens$data_saida)
cat(sprintf("\n[esforço] Num_dias == (chegada - partida): %.1f%% das viagens;\n",
            100 * mean(viagens$dias == viagens$dur_datas, na.rm = TRUE)))
cat(sprintf("          == (chegada - partida) + 1: %.1f%%; incoerentes: %.1f%%\n",
            100 * mean(viagens$dias == viagens$dur_datas + 1, na.rm = TRUE),
            100 * mean(viagens$dias != viagens$dur_datas &
                         viagens$dias != viagens$dur_datas + 1, na.rm = TRUE)))
cat(sprintf("          Num_horas <= 0 ou ausente: %d viagens\n",
            sum(viagens$horas <= 0 | is.na(viagens$horas))))
cat(sprintf("          Num_horas > 24 x Num_dias: %.1f%% das viagens\n",
            100 * mean(viagens$horas > 24 * viagens$dias, na.rm = TRUE)))
cat("          -> `dias` é a medida padrão; `horas` é comparada por AIC\n")
cat("             na parte 03 (mesma resposta, mesmas linhas, só muda o\n")
cat("             offset), mas com esta ressalva registrada.\n")

## =====================================================================
## 8) TABELA FINAL E CONFERÊNCIA
## =====================================================================
cols_cap <- paste0("cap_", todas_col)
viagens$cap_total <- rowSums(viagens[, cols_cap])

cat("\n======================================================\n")
cat("TABELA `viagens` PRONTA\n")
cat("======================================================\n")
cat(sprintf("Viagens: %d | Embarcações: %d | Anos: %s\n",
            nrow(viagens), length(unique(viagens$barco_id)),
            paste(range(viagens$ano), collapse = "-")))
cat(sprintf("Esforço total: %.0f dias de pesca\n", sum(viagens$dias, na.rm = TRUE)))
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

## Resumo por ano — é a primeira olhada no que a série tem a dizer.
cat("\n--- resumo por ano (frota de cerco) ---\n")
res_ano <- data.frame(
  ano       = sort(unique(viagens$ano)),
  viagens   = as.numeric(table(viagens$ano)),
  dias      = round(tapply(viagens$dias, viagens$ano, sum)),
  cap_total = round(tapply(viagens$cap_total, viagens$ano, sum)),
  cavala_t  = round(tapply(viagens$cap_macarellus, viagens$ano, sum), 1),
  pct_cavala = round(100 * tapply(viagens$cap_macarellus, viagens$ano, sum) /
                       tapply(viagens$cap_total, viagens$ano, sum), 1),
  pct_viag_com_cavala = round(100 * tapply(viagens$cap_macarellus > 0,
                                           viagens$ano, mean), 1),
  row.names = NULL)
print(res_ano, row.names = FALSE)
cat("\nLeitura: a cavala cai de ~18% para ~2-7% da captura do cerco ao\n")
cat("longo da série, e a fração de viagens que a encontram cai junto.\n")
cat("Saber quanto disso é abundância e quanto é direcionamento é\n")
cat("exatamente o que as partes 02 e 03 vão tentar separar.\n")

cat("\nEstrutura da tabela:\n"); str(viagens, give.attr = FALSE)



#===================================================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde
# Parte 2
# A análise inteira gira em torno de uma variável que NÃO existe no banco: a espécie que o mestre
#pretendia pescar. Ninguém registra intenção. O que existe é a COMPOSIÇÃO DA CAPTURA, que funciona 
# como impressão digital da tática usada. Este script transforma essa composição em covariável 
# utilizável, de duas formas CONCORRENTES:
#
#   (a) DISCRETA  — agrupa as unidades de esforço por similaridade de composição 
#                   (k-means, conferido com Ward) e usa o rótulo do grupo como FATOR `alvo`.
#                   É o padrão nos grupos do ICCAT (Sant'Ana et al. 2020).
#   (b) CONTÍNUA  — usa os escores dos primeiros eixos de uma PCA dacomposição (PC1..PCn) como 
#                   preditores contínuos. É o "DPC" de Winker et al. (2013), que em estudos
#                   independentes ajustou melhor que o cluster.
#
#   * a PCA não agrupa. Ela reescreve a matriz de composição em eixos novos, ortogonais,
#     ordenados por variância explicada. Cada unidade recebe uma NOTA CONTÍNUA em cada eixo. 
#     É um gradiente.
#   * o k-means agrupa. Ele parte os pontos em k grupos discretos. Aquiele opera SOBRE OS 
#     ESCORES DA PCA (não sobre a composição bruta), porque a PCA já filtrou ruído e dimensões 
#     redundantes.
#   * logo: PCA -> k-means é sequencial NA CONSTRUÇÃO; mas os dois PRODUTOS (o fator `alvo` 
#     e os escores PC1..PCn) entram na parte 03 como representações ALTERNATIVAS da mesma 
#     coisa, comparadas por AIC/BIC. É a hipótese H2.
#
# POR QUE AGREGAR POR BARCO-MÊS ANTES DE AGRUPAR
# ----------------------------------------------
# ~59% das viagens de cerco registram UMA única espécie. Uma viagem com uma espécie só tem composição
# degenerada (100% daquela espécie) e não distingue "tática dirigida" de "sorte num lance". 
# Hoyle et al. (2018) recomendam agregar por barco-mês exatamente por isso: a agregação
# dilui o encontro fortuito com o cardume e deixa o padrão de estratégia aparecer.
#==================================================================================================

stopifnot(exists("viagens"))

## =====================================================================
## 0) FATOR TEMPORAL
## ---------------------------------------------------------------------
## O índice de abundância É o efeito do fator temporal. Com a série
## 2019-2025 completa, esse fator é o ANO. O código continua detectando
## sozinho o caso de um ano só (que cairia para mês) para não quebrar se
## alguém rodar um recorte.
## =====================================================================
UM_ANO_SO    <- length(unique(viagens$ano)) == 1
fator_tempo  <- if (UM_ANO_SO) "fmes" else "fano"
rotulo_tempo <- if (UM_ANO_SO) "Mês" else "Ano"
cat(sprintf("\n>> Fator temporal desta rodada: `%s` (%s)\n", fator_tempo,
            if (UM_ANO_SO) "DEMONSTRAÇÃO — um ano só" else
              sprintf("série %s", paste(range(viagens$ano), collapse = "-"))))

cols_cap <- setdiff(grep("^cap_", names(viagens), value = TRUE), "cap_total")
sp_nomes <- c(cap_macarellus = "D. macarellus", cap_katsuwonus = "K. pelamis",
              cap_selar      = "S. crumenoph.", cap_carangideo = "C. crysos",
              cap_punctatus  = "D. punctatus",  cap_sardinella = "S. maderensis",
              cap_thunnus    = "T. albacares",  cap_spicara    = "S. melanurus",
              cap_elagatis   = "E. bipinnulata", cap_apsilus   = "A. fuscus",
              cap_auxis      = "Auxis spp.",    cap_outras     = "Outras")
sp_nomes <- sp_nomes[cols_cap]

## Antialiasing dos PNG: "cleartype" só existe no Windows; em Linux/Mac
## o R para com erro. Esta linha deixa o script rodar nos dois.
AA <- if (.Platform$OS.type == "windows") "cleartype" else "default"

cor_sp  <- hcl.colors(length(cols_cap), palette = "Dark 3")
COR_MAC <- "#1F4E79"; COR_AUX <- "#C0501B"; COR_NEU <- "#7F7F7F"

## =====================================================================
## 1) VARIÁVEIS DERIVADAS
## ---------------------------------------------------------------------
## As CPUE nominais NÃO entram no modelo (lá o esforço entra como
## offset); servem para a exploratória e para o cenário S1, que
## reproduz o que a FAO (2026) fez.
## =====================================================================
viagens$cpue_dia  <- viagens$cap_macarellus / viagens$dias
viagens$cpue_hora <- viagens$cap_macarellus / pmax(viagens$horas, 1)
viagens$pos_mac   <- as.integer(viagens$cap_macarellus > 0)

# composição proporcional DA VIAGEM (só para descrever; o cluster usa a
# composição por barco-mês, ver seção 5)
for (cc in cols_cap)
  viagens[[sub("cap_", "prop_", cc)]] <-
  ifelse(viagens$cap_total > 0, viagens[[cc]] / viagens$cap_total, 0)

## Fatores. TUDO que vai ser marginalizado depois precisa ser fator — é
## assim que o emmeans consegue tirar a média sobre os níveis e devolver
## o efeito do tempo "limpo" dos demais.
## Note o que NÃO está aqui: `tipo_emb` (decisão L12 da parte 01) e
## `ilha_desemb` (decisão L7) — ambos confundidos com o ano.
viagens$fano   <- factor(viagens$ano)
viagens$fmes   <- factor(viagens$mes, levels = 1:12)
viagens$ftri   <- factor(viagens$trimestre, levels = 1:4,
                         labels = c("T1", "T2", "T3", "T4"))
viagens$fbanco <- factor(viagens$banco_gr)
viagens$fbarco <- factor(viagens$barco_id)          # código, não nome (L11)
viagens$barco_mes <- paste(viagens$barco_id, viagens$ano, viagens$mes, sep = "_")

## =====================================================================
## 2) FILTROS — cada um registrado, com o que custou
## ---------------------------------------------------------------------
## Filtro é decisão analítica, não faxina: muda a população amostrada e
## precisa ser reportado no texto. Por isso o log fica num data.frame.
## =====================================================================
MIN_VIAG_BARCO  <- 10   # com 7 anos dá para exigir mais que antes
MIN_POR_ESTRATO <- 3

n0 <- nrow(viagens)
filtro_log <- data.frame(regra = character(), removidas = integer(),
                         restantes = integer(), stringsAsFactors = FALSE)
reg <- function(regra, antes)
  filtro_log[nrow(filtro_log) + 1, ] <<- list(regra, antes - nrow(viagens),
                                              nrow(viagens))

## (i) esforço válido. `dias` é a medida principal; `horas` só é exigida
## porque a parte 03 compara os dois offsets por AIC e a comparação só
## vale se as LINHAS forem as mesmas nos dois modelos.
a <- nrow(viagens); viagens <- viagens[!is.na(viagens$dias) & viagens$dias > 0, ]
reg("esforço (dias) válido e > 0", a)
a <- nrow(viagens); viagens <- viagens[!is.na(viagens$horas) & viagens$horas > 0, ]
reg("esforço (horas) válido e > 0", a)

## (ii) viagem sem captura nenhuma não tem composição e, portanto, não
## tem tática atribuível. (Não confundir com viagem sem CAVALA: essas
## são a informação central e ficam.)
a <- nrow(viagens); viagens <- viagens[viagens$cap_total > 0, ]
reg("alguma captura registrada na viagem", a)

## (iii) embarcações com pouquíssimas viagens não sustentam efeito
## aleatório e desequilibram o cruzamento com o fator temporal.
tb <- table(viagens$barco_id)
a <- nrow(viagens)
viagens <- viagens[viagens$barco_id %in% names(tb[tb >= MIN_VIAG_BARCO]), ]
reg(sprintf("embarcações com >= %d viagens", MIN_VIAG_BARCO), a)

## (iv) estratos tempo x banco muito ralos geram coeficientes instáveis
## (e, no limite, níveis que só existem em um ano — que o modelo
## confundiria com efeito de ano).
te <- table(paste(viagens$ano, viagens$banco_gr))
a <- nrow(viagens)
viagens <- viagens[paste(viagens$ano, viagens$banco_gr) %in%
                     names(te[te >= MIN_POR_ESTRATO]), ]
reg(sprintf("estratos ano x banco com >= %d viagens", MIN_POR_ESTRATO), a)

for (f in c("fano", "fmes", "ftri", "fbanco", "fbarco"))
  viagens[[f]] <- droplevels(viagens[[f]])

cat("\n================== FILTROS ================\n"); print(filtro_log, row.names = FALSE)
cat(sprintf("Retidas %d de %d viagens (%.1f%%)\n", nrow(viagens), n0,
            100 * nrow(viagens) / n0))
cat(sprintf("Após filtros: %d anos, %d bancos, %d embarcações\n",
            nlevels(viagens$fano), nlevels(viagens$fbanco),
            nlevels(viagens$fbarco)))
cat("Viagens por ano após filtros:\n"); print(table(viagens$ano))

## =====================================================================
## 3) SÉRIE NOMINAL (cenário S1) — captura da espécie / esforço TOTAL
## ---------------------------------------------------------------------
## Este é o índice que a FAO usou e o que queremos submeter à prova: o
## denominador é o esforço de TODA a frota de cerco, inclusive as viagens
## que estavam atrás de Auxis. É justamente por isso que ele é suspeito.
## =====================================================================
viagens$tempo <- viagens[[sub("^f", "", fator_tempo)]]

ser <- aggregate(viagens[, c(cols_cap, "dias", "horas")],
                 by = list(tempo = viagens$tempo), FUN = sum)
ser$cpue_nom_dia  <- ser$cap_macarellus / ser$dias
ser$cpue_nom_hora <- ser$cap_macarellus / ser$horas
ser$prop_zero <- tapply(viagens$pos_mac == 0, viagens$tempo, mean)[as.character(ser$tempo)]
ser$n_viagens <- as.numeric(table(viagens$tempo)[as.character(ser$tempo)])

cat(sprintf("\n===== SÉRIE POR %s =====\n", toupper(rotulo_tempo)))
print(data.frame(tempo = ser$tempo,
                 cavala_t = round(ser$cap_macarellus, 1),
                 auxis_t  = round(ser$cap_auxis, 1),
                 dias = ser$dias,
                 cpue_t_dia = round(ser$cpue_nom_dia, 3),
                 zeros = sprintf("%.0f%%", 100 * ser$prop_zero),
                 n = ser$n_viagens), row.names = FALSE)

## =====================================================================
## 4) FIGURAS EXPLORATÓRIAS
## =====================================================================

## --- Fig 1: captura por espécie e esforço ------------------------------
png("exp1_capturas_esforco.png", width = 26, height = 13, res = 300,
    antialias = AA, units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
matplot(ser$tempo, ser[, cols_cap], type = "l", lty = 1, lwd = 2.2,
        col = cor_sp, xlab = rotulo_tempo, ylab = "Captura (t)",
        main = "A. Captura por especie")
legend("topleft", sp_nomes, col = cor_sp, lwd = 2.2, bty = "n", cex = 0.6)
plot(ser$tempo, ser$dias, type = "b", pch = 19, lwd = 2.4, col = COR_NEU,
     xlab = rotulo_tempo, ylab = "Esforco (dias de pesca)",
     main = "B. Esforco amostrado da frota de cerco",
     ylim = c(0, max(ser$dias) * 1.05))
par(op); dev.off()
cat("\nPNG salvo: exp1_capturas_esforco.png\n")

## --- Fig 2: composição da captura --------------------------------------
## É a figura que conta a história central: a fração da cavala encolhe e
## a do Auxis cresce. Isso é compatível COM TROCA DE ALVO e TAMBÉM com
## queda real de abundância — a figura levanta a questão, não a resolve.
comp <- as.matrix(ser[, cols_cap]); comp <- comp / rowSums(comp)
png("exp2_composicao.png", width = 24, height = 12, res = 300,
    antialias = AA, units = "cm")
op <- par(mar = c(4.2, 4.6, 3, 9), bty = "l", cex.main = 0.95, cex = 0.85)
acum <- t(apply(comp, 1, cumsum))
plot(NA, xlim = range(ser$tempo), ylim = c(0, 1), xlab = rotulo_tempo,
     ylab = "Proporcao da captura",
     main = "Composicao da captura da frota de cerco")
for (k in ncol(acum):1)
  polygon(c(ser$tempo, rev(ser$tempo)), c(acum[, k], rep(0, nrow(acum))),
          col = cor_sp[k], border = NA)
par(xpd = TRUE)
legend(max(ser$tempo) + diff(range(ser$tempo)) * 0.03, 0.95, sp_nomes,
       fill = cor_sp, border = NA, bty = "n", cex = 0.66)
par(op); dev.off()
cat("PNG salvo: exp2_composicao.png\n")

## --- Fig 3: CPUE nominal e zeros ---------------------------------------
## O painel B é o diagnóstico mais importante desta figura: com ~88% de
## viagens sem cavala, qualquer modelo que não trate os zeros direito
## (Tweedie ou hurdle) vai dar resultado errado.
png("exp3_cpue_nominal.png", width = 26, height = 13, res = 300,
    antialias = AA, units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
plot(ser$tempo, ser$cpue_nom_dia, type = "b", pch = 19, lwd = 2.4, col = COR_MAC,
     xlab = rotulo_tempo, ylab = "CPUE nominal da cavala (t/dia)",
     main = "A. CPUE nominal", ylim = c(0, max(ser$cpue_nom_dia) * 1.05))
plot(ser$tempo, 100 * ser$prop_zero, type = "b", pch = 19, lwd = 2.4, col = COR_AUX,
     xlab = rotulo_tempo, ylab = "% de viagens sem cavala",
     main = "B. Zeros de direcionamento", ylim = c(0, 100))
par(op); dev.off()
cat("PNG salvo: exp3_cpue_nominal.png\n")

## --- Fig 4: as covariáveis operacionais valem a pena? -------------------
## COMO LER: com ~88% de viagens sem cavala, gráfico de CPUE contra
## covariável vira uma parede de zeros e não informa nada. A leitura
## correta para dado assim é DECOMPOR, que é exatamente o que o modelo
## hurdle faz: (i) a PROBABILIDADE de a viagem pegar cavala e (ii) QUANTO
## ela pega, dado que pegou. Uma covariável pode atuar só na primeira
## (tem a ver com onde/como se procura), só na segunda (tem a ver com
## capacidade de captura) ou nas duas. Estes painéis são a justificativa
## empírica de cada termo do modelo.
## `n_vertical = TRUE` escreve o n de pé: com muitas barras de altura
## parecida os rótulos deitados se sobrepõem e viram borrão.
barra_prop <- function(prop, n, titulo, xlab, cex_nome = 0.7,
                       cor = "#8FAADC", n_vertical = FALSE) {
  bp <- barplot(100 * prop, col = cor, border = NA, names.arg = NA,
                ylab = "% de viagens com cavala",
                xlab = xlab, main = titulo,
                ylim = c(0, max(100 * prop, na.rm = TRUE) * 1.3))
  text(bp, par("usr")[3], labels = names(prop), srt = 45, adj = 1,
       xpd = NA, cex = cex_nome)
  text(bp, 100 * prop, labels = paste0("n=", n), pos = 3, cex = 0.55,
       col = "#52514E", xpd = NA, srt = if (n_vertical) 90 else 0,
       offset = if (n_vertical) 0.9 else 0.5)
  invisible(bp)
}
png("exp4_covariaveis.png", width = 26, height = 20, res = 300,
    antialias = AA, units = "cm")
op <- par(mfrow = c(2, 2), mar = c(7.5, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)

## A — ONDE se pega cavala (componente de presença, por banco de pesca).
##     Só os 20 bancos com mais viagens, senão o eixo fica ilegível.
pr_b <- tapply(viagens$pos_mac, viagens$fbanco, mean)
n_b  <- table(viagens$fbanco)
top_b <- names(sort(n_b, decreasing = TRUE))[1:min(20, length(n_b))]
ord  <- top_b[order(pr_b[top_b], decreasing = TRUE)]
barra_prop(pr_b[ord], n_b[ord], "A. Presenca de cavala (20 maiores bancos)",
           "", cex_nome = 0.55, n_vertical = TRUE)

## B — QUANTO se pega, dado que pegou (componente de magnitude).
pos <- viagens[viagens$pos_mac == 1, ]
## Só os 15 bancos com mais viagens POSITIVAS: abaixo disso a caixa é
## desenhada sobre 3-4 pontos e não descreve distribuição nenhuma.
n_pos_b <- sort(table(pos$fbanco), decreasing = TRUE)
bancos_ok <- names(n_pos_b[n_pos_b >= 10])[1:min(15, sum(n_pos_b >= 10))]
pos_b <- pos[pos$fbanco %in% bancos_ok, ]
if (nrow(pos_b) > 0) {
  pos_b$fbanco <- droplevels(pos_b$fbanco)
  bx <- boxplot(cpue_dia ~ fbanco, data = pos_b, outline = FALSE, plot = FALSE)
  boxplot(cpue_dia ~ fbanco, data = pos_b, outline = FALSE, col = "#74C476",
          xaxt = "n", xlab = "", ylab = "CPUE (t/dia) entre as positivas",
          lwd = 1, main = "B. Magnitude, so nas viagens com cavala")
  axis(1, at = seq_along(bx$names), labels = FALSE)
  text(seq_along(bx$names), par("usr")[3], labels = bx$names, srt = 45,
       adj = 1, xpd = NA, cex = 0.6)
}

## C — tripulação em classes (decisão L10). Se as barras forem
##     praticamente iguais, a covariável não está medindo poder de pesca
##     e isso também é resultado.
barra_prop(tapply(viagens$pos_mac, viagens$npesc_cat, mean),
           table(viagens$npesc_cat),
           "C. Tripulacao (proxy de poder de pesca)", "Numero de pescadores",
           cor = "#B497D6")

## D — trimestre (decisão L9). É a variável sazonal que vai para o modelo.
barra_prop(tapply(viagens$pos_mac, viagens$ftri, mean), table(viagens$ftri),
           "D. Trimestre (sazonalidade)", "Trimestre", cor = "#E8A33D")
par(op); dev.off()
cat("PNG salvo: exp4_covariaveis.png\n")

## --- Fig 5: CPUE no tempo, na frota e na sazonalidade fina -------------
## Esta figura existe para responder quatro perguntas que a exploratória
## precisa responder ANTES de modelar:
##   A. o padrão sazonal é o mesmo em todos os anos? (se não for, um
##      efeito aditivo de trimestre não basta e seria preciso interação)
##   B. a sazonalidade mensal tem forma que 4 níveis descrevem bem?
##   C. a captura de cavala está concentrada em poucos barcos? (se
##      estiver, o efeito aleatório de embarcação é indispensável)
##   D. a frota mudou ao longo da série? (é o problema P5: composição de
##      frota variando no tempo contamina o efeito de ano)
png("exp5_cpue_frota_tempo.png", width = 26, height = 20, res = 300,
    antialias = AA, units = "cm")
op <- par(mfrow = c(2, 2), mar = c(4.6, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)

## A — CPUE nominal por ano, separada por trimestre
cpue_at <- tapply(viagens$cap_macarellus, list(viagens$ano, viagens$ftri), sum) /
  tapply(viagens$dias,           list(viagens$ano, viagens$ftri), sum)
cor_tri <- hcl.colors(4, palette = "Zissou 1")
matplot(as.numeric(rownames(cpue_at)), cpue_at, type = "b", pch = 19, lty = 1,
        lwd = 2, col = cor_tri, xlab = "Ano", ylab = "CPUE (t/dia)",
        main = "A. CPUE nominal por ano e trimestre")
legend("topright", colnames(cpue_at), col = cor_tri, lwd = 2, pch = 19,
       bty = "n", cex = 0.72)

## B — sazonalidade mensal: presença e magnitude na mesma figura
pr_m <- tapply(viagens$pos_mac, viagens$mes, mean)
cp_m <- tapply(viagens$cap_macarellus, viagens$mes, sum) /
  tapply(viagens$dias, viagens$mes, sum)
plot(as.numeric(names(pr_m)), 100 * pr_m, type = "b", pch = 19, lwd = 2.2,
     col = COR_MAC, xlab = "Mes", ylab = "% de viagens com cavala",
     main = "B. Sazonalidade mensal", xaxt = "n",
     ylim = c(0, max(100 * pr_m) * 1.15))
axis(1, at = 1:12)
abline(v = c(3.5, 6.5, 9.5), lty = 3, col = COR_NEU)   # limites de trimestre
par(new = TRUE)
plot(as.numeric(names(cp_m)), cp_m, type = "b", pch = 17, lty = 2, lwd = 2,
     col = COR_AUX, axes = FALSE, xlab = "", ylab = "",
     ylim = c(0, max(cp_m) * 1.15))
axis(4, col.axis = COR_AUX)
mtext("CPUE (t/dia)", side = 4, line = 2.2, cex = 0.75, col = COR_AUX)
legend("topright", c("% com cavala", "CPUE (eixo dir.)"),
       col = c(COR_MAC, COR_AUX), lwd = 2, pch = c(19, 17), lty = c(1, 2),
       bty = "n", cex = 0.7)

## C — CPUE por embarcação (só as com >= 30 viagens, ordenadas).
##     A dispersão entre barcos é a justificativa do termo (1 | fbarco).
nb <- table(viagens$fbarco)
bok <- names(nb[nb >= 30])
cp_b <- tapply(viagens$cap_macarellus[viagens$fbarco %in% bok],
               droplevels(viagens$fbarco[viagens$fbarco %in% bok]), sum) /
  tapply(viagens$dias[viagens$fbarco %in% bok],
         droplevels(viagens$fbarco[viagens$fbarco %in% bok]), sum)
cp_b <- sort(cp_b, decreasing = TRUE)
## Sem rótulo por barra: o código do barco individual não informa nada a
## quem lê; o que importa é a FORMA do perfil (quão desigual é a frota).
barplot(cp_b, col = COR_MAC, border = NA, names.arg = rep("", length(cp_b)),
        ylim = c(0, max(cp_b) * 1.12),
        ylab = "CPUE de cavala (t/dia)",
        xlab = sprintf("Embarcacoes ordenadas (n = %d)", length(cp_b)),
        main = "C. CPUE por embarcacao (>=30 viagens)")
media_frota <- sum(viagens$cap_macarellus) / sum(viagens$dias)
abline(h = media_frota, lty = 2, col = COR_AUX, lwd = 2)
text(par("usr")[2], media_frota+0.15*media_frota, "media da frota", 
     pos = 2,adj = c(1, -5), offset = 0.5, cex = 0.9, col = COR_AUX)

## D — a frota muda? barcos ativos por ano e concentração do esforço
nb_ano <- tapply(viagens$barco_id, viagens$ano, function(x) length(unique(x)))
plot(as.numeric(names(nb_ano)), nb_ano, type = "b", pch = 19, lwd = 2.2,
     col = COR_MAC, xlab = "Ano", ylab = "Embarcacoes ativas",
     main = "D. Composicao da frota ao longo da serie",
     ylim = c(0, max(nb_ano) * 1.15))
par(new = TRUE)
vpb <- as.numeric(table(viagens$ano)) / nb_ano
plot(as.numeric(names(nb_ano)), vpb, type = "b", pch = 17, lty = 2, lwd = 2,
     col = COR_AUX, axes = FALSE, xlab = "", ylab = "",
     ylim = c(0, max(vpb) * 1.15))
axis(4, col.axis = COR_AUX)
mtext("Viagens por embarcacao", side = 4, line = 2.2, cex = 0.75, col = COR_AUX)
legend("bottomright", c("barcos ativos", "viagens/barco (dir.)"),
       col = c(COR_MAC, COR_AUX), lwd = 2, pch = c(19, 17), lty = c(1, 2),
       bty = "n", cex = 0.7)
par(op); dev.off()
cat("PNG salvo: exp5_cpue_frota_tempo.png\n")

## Tabelas que acompanham a figura 5 (números exatos para o texto)
cat("\n--- CPUE nominal (t/dia) por ano e trimestre ---\n")
print(round(cpue_at, 3))
cat("\n--- presença de cavala por classe de tripulação ---\n")
print(round(cbind(n = as.numeric(table(viagens$npesc_cat)),
                  prop_com_cavala = tapply(viagens$pos_mac, viagens$npesc_cat, mean),
                  cpue_media = tapply(viagens$cpue_dia, viagens$npesc_cat, mean)), 3))


## =====================================================================
## 5) INFERÊNCIA DA TÁTICA — composição por barco-mês
## =====================================================================

## 5.1 matriz de composição agregada por BARCO-MÊS (ver justificativa no
##     cabeçalho). Proporções em peso; transformação raiz quadrada para
##     que espécies menos abundantes também pesem na similaridade, como
##     em Winker et al. (2013). Sem a raiz, a distância entre unidades
##     seria decidida quase só pelo Auxis, que domina o peso.
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
cat("    fraco pois depende de agrupamento de viagens semelhantes\n")
cat(" Então a representação contínua (PCA) pode ser a mais defensável.\n")

## 5.2 PCA. Sem `scale.` porque as colunas já estão na mesma unidade
##     (proporções transformadas) — escalonar daria peso igual a
##     espécies raras e dominantes, que não é o que queremos.
pca <- prcomp(comp_sqrt, center = TRUE, scale. = FALSE)
var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cat("\nVariância explicada pelos eixos da PCA: ",
    paste(sprintf("PC%d=%.1f%%", 1:min(4, length(var_exp)),
                  var_exp[1:min(4, length(var_exp))]), collapse = "  "), "\n")
cat("Cargas de PC1 (o que este eixo separa — negativo x positivo):\n")
print(round(sort(pca$rotation[, 1]), 3))
cat("Cargas de PC2:\n")
print(round(sort(pca$rotation[, 2]), 3))

## Quantos eixos levar para o modelo? Guardamos os primeiros PCs até
## acumular ~70% da variação da composição (teto de 4, para não inflar o
## modelo). Numa pescaria com composição realmente multidimensional,
## usar apenas PC1/PC2 jogaria fora a maior parte do sinal de tática.
cum_var <- cumsum(var_exp)
n_pc <- min(4, max(2, which(cum_var >= 70)[1]), ncol(pca$x))
if (is.na(n_pc)) n_pc <- min(4, ncol(pca$x))
cat(sprintf("Eixos retidos para o modelo: %d (%.0f%% da variação)\n",
            n_pc, cum_var[n_pc]))

## 5.3 Número de grupos pela silhueta média.
##     Implementada em R base para não exigir o pacote `cluster`.
##     ATENÇÃO à leitura: a silhueta favorece sistematicamente k pequeno.
##     Ela responde "os grupos estão separados?", não "quantas táticas
##     existem?". Se o gráfico de PC1xPC2 mostrar uma nuvem contínua, a
##     resposta honesta é que não há grupos — há um gradiente, e a
##     representação contínua (PCs) é a mais fiel.
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
##     isso é informação, não contratempo: é mais um argumento a favor
##     da representação contínua.
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
## não dominar nenhum grupo. A tática relevante é aquela com a MAIOR
## proporção de cavala no centróide, seja qual for o nome dela. Isso é
## guardado para o cenário S3 (esforço dirigido) da parte 03.
alvo_cavala <- nome_cl[which.max(cent[, "macarellus"])]
frac_cavala <- max(cent[, "macarellus"])
cat(sprintf("\nTática com maior fração de cavala no centróide: '%s' (%.1f%% da captura)\n",
            alvo_cavala, 100 * frac_cavala))
if (frac_cavala < 0.30) {
  cat("  [NOTA] Nenhuma tática é dominada pela cavala. Ela é capturada\n")
  cat("         acompanhando outras espécies, não como alvo exclusivo —\n")
  cat("         o que enfraquece o cenário de 'esforço dirigido' (S3) e\n")
  cat("         reforça o uso da composição como covariável (H1/H2).\n")
} else {
  cat("  [OK] Existe uma tática em que a cavala domina a composição: o\n")
  cat("       cenário de esforço dirigido (S3) tem base empírica.\n")
}

## 5.6 Levar o rótulo e os escores de volta para cada VIAGEM.
##     Repare que as DUAS coisas vão juntas: `alvo` (discreto, do
##     k-means) e PC1..PCn (contínuos, da PCA). Elas NÃO são usadas ao
##     mesmo tempo no mesmo modelo — são alternativas comparadas na
##     parte 03 (estruturas E3/E4 contra E5).
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
##     a captura de cavala. Se não separar, ela não está medindo alvo —
##     e aí o termo `alvo` não vai sobreviver à seleção da parte 03.
cat("\nCPUE média da cavala por tática inferida (t/dia):\n")
print(round(tapply(viagens$cpue_dia, viagens$alvo, mean), 3))
cat("Proporção de viagens com cavala, por tática:\n")
print(round(tapply(viagens$pos_mac, viagens$alvo, mean), 3))

## A tática mudou ao longo da série? Esta tabela é o coração da hipótese
## H1: se a mistura de táticas mudou, a CPUE nominal mistura mudança de
## comportamento com mudança de abundância.
cat("\nDistribuição das táticas por ano (proporção das viagens):\n")
print(round(prop.table(table(viagens$ano, viagens$alvo), margin = 1), 3))

## --- Fig 6: PCA, silhueta e composição das táticas ---------------------
cor_cl <- hcl.colors(k_otimo, palette = "Dark 3")
png("exp6_taticas.png", width = 27, height = 10.5, res = 300,
    antialias = AA, units = "cm")
op <- par(mfrow = c(1, 3), mar = c(4.6, 4.4, 3, 1), oma = c(0, 0, 0, 5),
          bty = "l", cex.main = 0.95, cex = 0.85)
plot(ks, sil, type = "b", pch = 19, lwd = 2, col = COR_MAC,
     xlab = "Numero de grupos (k)", ylab = "Silhueta media", main = "A. Escolha de k")
points(k_otimo, sil[ks == k_otimo], pch = 21, bg = COR_AUX, cex = 1.9)
plot(pca$x[, 1], pca$x[, 2], col = adjustcolor(cor_cl[km$cluster], 0.6),
     pch = 16, cex = 0.7, xlab = sprintf("PC1 (%.0f%%)", var_exp[1]),
     ylab = sprintf("PC2 (%.0f%%)", var_exp[2]),
     main = "B. Taticas no espaco de composicao")
points(km$centers[, 1], km$centers[, 2], pch = 21, bg = cor_cl, cex = 2, lwd = 1.5)
legend("topright", nome_cl, col = cor_cl, pch = 16, bty = "n", cex = 0.7)
bp <- barplot(t(cent), col = cor_sp, border = NA, names.arg = nome_cl,
              las = 2, cex.names = 0.65, ylab = "Proporcao media da captura",
              main = "C. Composicao de cada tatica")
legend(max(bp) + 0.8, 1, rev(sp_nomes), fill = rev(cor_sp), border = NA,
       bty = "n", cex = 0.6, xpd = NA)
par(op); dev.off()
cat("PNG salvo: exp6_taticas.png\n")

## =====================================================================
## 6) ESFORÇO DIRIGIDO (insumo do cenário S3 na parte 03)
## ---------------------------------------------------------------------
## Reparte os dias de pesca de cada estrato entre as táticas, na
## proporção dos dias das viagens de cada uma. É isto que converte
## "esforço total da frota" em "esforço dirigido à cavala" — o campo que
## o IMar não tem e que a composição da captura permite reconstruir.
## RESSALVA: se nenhuma tática for dominada pela cavala (ver 5.5), este
## cenário mede "esforço da tática que mais encontra cavala", que é
## menos do que "esforço dirigido à cavala". A diferença tem de estar no
## texto.
## =====================================================================
esforco_dirigido <- aggregate(cbind(dias, horas) ~ tempo + banco_gr + alvo,
                              data = viagens, FUN = sum)
dias_alvo <- aggregate(dias ~ tempo + alvo, data = esforco_dirigido, FUN = sum)

png("exp7_esforco_dirigido.png", width = 26, height = 12, res = 300,
    antialias = AA, units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
tat <- prop.table(table(viagens$tempo, viagens$alvo), margin = 1)
matplot(as.numeric(rownames(tat)), tat, type = "b", pch = 19, lty = 1, lwd = 2.2,
        col = cor_cl, ylim = c(0, 1), xlab = rotulo_tempo,
        ylab = "Proporcao das viagens", main = "A. Mistura de taticas")
legend("topleft", nome_cl, col = cor_cl, lwd = 2.2, pch = 19, bty = "n", cex = 0.68)
mat_d <- sapply(levels(viagens$alvo), function(g) {
  x <- dias_alvo[dias_alvo$alvo == g, ]
  v <- setNames(rep(0, nrow(ser)), ser$tempo)
  v[as.character(x$tempo)] <- x$dias; v
})
matplot(ser$tempo, mat_d, type = "l", lty = 1, lwd = 2.4, col = cor_cl,
        ylim = c(0, max(c(mat_d, ser$dias)) * 1.05), xlab = rotulo_tempo,
        ylab = "Dias de pesca", main = "B. Esforco dirigido por tatica")
lines(ser$tempo, ser$dias, lwd = 2, lty = 2, col = COR_NEU)
legend("topleft", c(levels(viagens$alvo), "esforco total"),
       col = c(cor_cl, COR_NEU), lwd = 2.2,
       lty = c(rep(1, nlevels(viagens$alvo)), 2), bty = "n", cex = 0.68)
par(op); dev.off()
cat("PNG salvo: exp7_esforco_dirigido.png\n")

cat("\n===== PARTE 02 CONCLUÍDA =====\n")
cat("Objetos para a parte 03: `viagens` (com alvo, PC1..PCn, tempo),\n")
cat("`ser`, `esforco_dirigido`, `fator_tempo`, `rotulo_tempo`, `PCs`,\n")
cat("`alvo_cavala`, `frac_cavala`.\n")



#====================================================================================================
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde Produzir índices 
# alternativos de abundância relativa da cavala preta para entrar no JABBA como CENÁRIOS CONCORRENTES:
#   S1 nominal   — captura da cavala / esforço total da frota de cerco
#                  (o que a FAO 2026 fez)
#   S2 corrigida — efeito do fator temporal num modelo que controla
#                  banco de pesca, trimestre, tripulação, embarcação e
#                  TÁTICA DE PESCA
# A distância entre as duas É a medida do viés de direcionamento.
#
# ================= PROBLEMAS DESTA CPUE (resumo) ===================================================
# P1 esforço não é específico da espécie (denominador multiespecífico)
# P2 troca de alvo ao longo do tempo (Schirripa & Goodyear 2010). A parte 02 mostrou que isso É REAL 
#    nesta série: a proporção de viagens da tática "selar" (a que mais encontra cavala) vai de ~42%
#    em 2019 a ~18% em 2021, volta a ~36% em 2023-24 e cai a ~17% em2025 — e a CPUE nominal da 
#    cavala acompanha esse movimento.
# P3 zeros de direcionamento — ~88% das viagens de cerco não registram cavala. Excluí-los ou somar 
#    constante enviesa o índice.
# P4 poucas variáveis operacionais, mas o arquivo do IMar tem mais do que era esperado: tripulação,
#    banco de pesca e identidade da embarcação. Profundidade e tipo de embarcação foram descartados
#    na parte 01 (decisões L6 e L12) — ver lá o porquê.
# P5 composição da frota muda ao longo do tempo
# P6 desbalanceamento espacial e sazonal
#
# ========================== HIPÓTESES ========================================================
# H0 a CPUE nominal é aceitável (o índice não muda ao controlar o resto)
# H1 parte da variação da CPUE nominal é troca de alvo, não abundância
# H2 tática contínua (PCA) ajusta melhor que tática discreta (cluster)
# H3 usar só o esforço dirigido reproduz H1 (sensibilidade; tende a gerar hiperestabilidade, 
#    então nunca é o cenário principal)
#
# ===================== MODELOS QUE SERÃO TESTADOS ============================================
# ESTRUTURAS (o fator temporal NUNCA entra na seleção — ele É o índice):
#   E0  tempo
#   E1  tempo + banco
#   E2  tempo + banco + trimestre + tripulação(classes)
#   E3  E2 + alvo                                   (tática DISCRETA)
#   E4  E3 + (1 | barco)                            [GLMM]
#   E5  E2 + PC1..PCn + (1 | barco)                 (tática CONTÍNUA)
# OFFSET: log(dias) contra log(horas).
# DISTRIBUIÇÕES (resposta = toneladas, contínua, com ~88% de zeros):
#   D1 Tweedie            D2 Hurdle-Gamma
#
# Duas famílias que modelam a captura inteira, zeros incluídos, nas MESMAS linhas:
#   D1 Tweedie      — Poisson composta com Gamma; massa em zero e cauda contínua positiva 
#                     num modelo só.
#   D2 Hurdle-Gamma — dois processos explícitos: um binomial para a presença e um Gamma para 
#                     a magnitude. Em glmmTMB,`ziGamma` + `ziformula` é literalmente um hurdle
#                     (a Gamma não tem massa em zero).
# Como as duas usam a MESMA resposta e as MESMAS linhas, AIC e BIC são diretamente comparáveis 
# entre elas. 
#
# ---------------------------------------------------------------------------------
# Sem LRT ; Diagnóstico e seleção de modelos
# ---------------------------------------------------------------------------------
# A seleção é por AIC/BIC + diagnóstico de resíduos (DHARMa). O LRT foi retirado 
# por três motivos: (i) com ~6 mil viagens ele acusa significância em efeito 
# irrelevante, o que empurra o modelo a comer sinal de abundância (o dilema de 
# Hinton & Maunder 2003); (ii) ele só vale entre modelos aninhados, enquanto 
# o AIC compara também os não aninhados que aparecem aqui (tática discreta 
# contra contínua); (iii) AIC e BIC já dão a ordenação, e quem decide adequação 
# é o resíduo, não o p-valor.
# Princípio que continua valendo: AIC/BIC medem AJUSTE, não ADEQUAÇÃO Um modelo 
# pode ganhar no AIC e ter resíduo ruim — nesse caso ele não é o escolhido.
# A tabela de diagnóstico é parte da decisão, não um anexo.
#==================================================================================

stopifnot(exists("viagens"), exists("ser"), exists("fator_tempo"))
suppressPackageStartupMessages({
  library(glmmTMB); library(emmeans); library(DHARMa)
})
tem_writexl <- requireNamespace("writexl", quietly = TRUE)

AA <- if (.Platform$OS.type == "windows") "cleartype" else "default"
COR_MAC <- "#1F4E79"; COR_AUX <- "#C0501B"; COR_NEU <- "#7F7F7F"
COR_S2  <- "#2E8B57"; COR_S3 <- "#7030A0"; COR_S2D <- "#00A0B0"

## =====================================================================
## 0) PREPARO DA RESPOSTA E DO ESFORÇO
## ---------------------------------------------------------------------
## A resposta é a CAPTURA em toneladas, não a CPUE. O esforço entra como
## OFFSET (coeficiente fixo em 1 na escala log), o que é equivalente a
## modelar a taxa mas preserva a estrutura de erro da captura — inclusive
## os zeros, que desapareceriam se dividíssemos.
## =====================================================================
viagens$captura <- viagens$cap_macarellus     # toneladas
viagens$ldias   <- log(viagens$dias)
viagens$lhoras  <- log(viagens$horas)

cat("\n===== DADOS PARA A MODELAGEM =====\n")
cat(sprintf("Viagens: %d | níveis de %s: %d | bancos: %d | barcos: %d | táticas: %d\n",
            nrow(viagens), fator_tempo, nlevels(viagens[[fator_tempo]]),
            nlevels(viagens$fbanco), nlevels(viagens$fbarco),
            nlevels(viagens$alvo)))
cat(sprintf("Zeros na resposta: %.1f%%  (P3 — decisivo para a escolha da distribuição)\n",
            100 * mean(viagens$captura == 0)))
cat(sprintf("Captura da cavala: %.1f t em %.0f dias de pesca\n",
            sum(viagens$captura), sum(viagens$dias)))

## Covariáveis que NÃO entram, e por quê (fica registrado no output para
## quem for ler o log do script sem ler a parte 01):
cat("\nCovariáveis descartadas na parte 01:\n")
cat("  profundidade   — 28% de ausentes codificados como 0 (decisão L6)\n")
cat("  tipo_embarcacao— mesmo barco com dois códigos e código quase\n")
cat("                   restrito a 2019 => confundido com ano (decisão L12)\n")
cat("  ilha/porto de desembarque — cobertura muda ao longo da série (L7)\n")

## =====================================================================
## 1) CONSTRUTOR DE FÓRMULAS
## ---------------------------------------------------------------------
## Montar a fórmula como string para evitar perder termos nos testes
## =====================================================================
monta_formula <- function(resposta, termos, aleatorio = TRUE,
                          offset_var = "ldias") {
  rhs <- paste(c(fator_tempo, termos,
                 if (aleatorio) "(1 | fbarco)",
                 if (!is.null(offset_var)) sprintf("offset(%s)", offset_var)),
               collapse = " + ")
  stats::as.formula(paste(resposta, "~", rhs))
}

## Guarda contra modelo mais complexo do que os dados sustentam: um termo
## categórico só entra se houver pelo menos MIN_POR_NIVEL observações por
## nível. Importa nos subconjuntos (cenário S3), onde `fbanco` tem
## dezenas de níveis e sobram poucas centenas de viagens.
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

## ==================================================================================
## 2) ESTRUTURAS CANDIDATAS (E0-E5), todas com família Tweedie
## ----------------------------------------------------------------------------------
## A Tweedie é a distribuição de TRABALHO nesta etapa: fixamos a distribuição para 
# comparar estruturas e só depois, já na melhor estrutura, comparamos as duas 
# distribuições. Comparar tudo contra tudo multiplicaria ajustes sem necessidade 
# e tornaria o resultado dependente da ordem em que se olha.
##
## As seis estruturas são uma escada: cada degrau acrescenta um tipo deexplicação 
#alternativa à abundância.
##   E0 só o tempo         -> o índice "cru" do modelo
##   E1 + banco            -> onde se pescou (P6, desbalanceamento espacial)
##   E2 + trimestre + trip.-> quando se pescou e com que poder de pesca
##   E3 + alvo             -> o que se estava tentando pescar (P2)
##   E4 + (1|barco)        -> quem pescou (P5, composição de frota)
##   E5 tática contínua    -> a alternativa da H2
## Todas ajustadas às MESMAS linhas e à MESMA resposta => AIC comparável.
## ===============================================================================
termos_E <- list(
  E0 = character(0),
  E1 = c("fbanco"),
  E2 = c("fbanco", "ftri", "npesc_cat"),
  E3 = c("fbanco", "ftri", "npesc_cat", "alvo"),
  E4 = c("fbanco", "ftri", "npesc_cat", "alvo"),
  E5 = c("fbanco", "ftri", "npesc_cat", PCs)
)
aleat_E <- c(E0 = FALSE, E1 = FALSE, E2 = FALSE, E3 = FALSE, E4 = TRUE, E5 = TRUE)

cat("\n===== 1) ESTRUTURAS (Tweedie, offset = log dias) =====\n")
fits <- list()
for (nm in names(termos_E)) {
  t0 <- Sys.time()
  f  <- monta_formula("captura", termos_E[[nm]], aleatorio = aleat_E[[nm]])
  fits[[nm]] <- try(glmmTMB(f, family = tweedie(link = "log"), data = viagens),
                    silent = FALSE)
  ok <- !inherits(fits[[nm]], "try-error")
  cat(sprintf("  %-3s %-64s %s (%.1f min)\n", nm,
              paste(deparse(f), collapse = ""),
              if (ok) sprintf("AIC=%.1f", AIC(fits[[nm]])) else "FALHOU",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
fits <- fits[!vapply(fits, inherits, logical(1), "try-error")]
stopifnot(length(fits) > 0)

# ---------------------------------------------------------------------
# Teste rápido de convergência dos modelos Tweedie em `fits`.
# Verifica pdHess (Hessiana definida positiva no ótimo) e o maior
# |gradiente| em valor absoluto — os dois sinais que realmente importam
# quando o otimizador solta avisos tipo "singular convergence".
# ---------------------------------------------------------------------
LIMIAR_GRADIENTE <- 1e-2  # acima disso, vale desconfiar do ajuste

cat("\n=== Checagem de convergencia (Tweedie) ===\n")
for (nome in names(fits)) {
  m      <- fits[[nome]]
  sdr    <- m$sdr
  pdhess <- isTRUE(sdr$pdHess)
  grad   <- sdr$gradient.fixed
  grad_max <- if (is.null(grad)) NA_real_ else max(abs(grad))
  msg_otim <- m$fit$message
  ok <- pdhess && !is.na(grad_max) && grad_max < LIMIAR_GRADIENTE
  
  cat(sprintf(
    "%-4s | pdHess=%-5s | max|grad|=%.2e | %-28s | %s\n",
    nome, pdhess, grad_max, msg_otim,
    if (ok) "OK" else "VERIFICAR"
  ))
}
cat("===========================================\n")

#*** O modelo E3 tem direção plana longa no otimizador .Isso sugere que a superfície 
#de verossimilhança do E3 tem alguma direção quase plana — bem possivelmente colinearidade 
#parcial entre alvo (que vem da composição de espécies) e fbanco/fano 
#(já que a composição de espécies varia sistematicamente por banco e por ano) — 
#e que adicionar o efeito aleatório de embarcação (1 | fbarco) ajuda a "resolver" 
#essa direção achatada, absorvendo parte dessa variação compartilhada
#--------------------------------------------------------------------

#tabela comparativa dos modelos com índices de ajuste
tab_est <- data.frame(
  modelo = names(fits),
  df     = vapply(fits, function(m) attr(logLik(m), "df"), numeric(1)),
  logLik = round(vapply(fits, function(m) as.numeric(logLik(m)), numeric(1)), 1),
  AIC    = round(vapply(fits, AIC, numeric(1)), 1),
  BIC    = round(vapply(fits, BIC, numeric(1)), 1), row.names = NULL)
tab_est$dAIC <- round(tab_est$AIC - min(tab_est$AIC), 1)
tab_est$dBIC <- round(tab_est$BIC - min(tab_est$BIC), 1)
cat("\nTodas ajustadas às mesmas linhas e à mesma resposta -> AIC/BIC comparáveis.\n")
cat("dAIC/dBIC são a distância para o melhor da coluna; <2 = empate técnico.\n")
print(tab_est, row.names = FALSE)
cat("\nAIC e BIC podem discordar: o BIC pune mais a complexidade, então\n")
cat("tende a preferir estrutura menor. Se discordarem, vale reportar os\n")
cat("dois e checar se o ÍNDICE muda — se não muda, a discordância é\n")
cat("acadêmica.\n")

## =====================================================================
## 3) MEDIDA DE ESFORÇO: dias contra horas no mar
## ---------------------------------------------------------------------
## Comparação legítima por AIC: mesma resposta, mesmas linhas, só muda o
## offset. A comparação é feita na estrutura E4 (a discreta completa),
## que é a estrutura do modelo principal — e não na "melhor da escada",
## para que a escolha do offset não dependa de qual degrau ganhou.
## RESSALVA REGISTRADA NA PARTE 01: `Num_horas` é inconsistente com
## `Num_dias` em ~21% das viagens (horas > 24 x dias), o que sugere que
## os dois campos não foram derivados do mesmo jeito. Por isso, se a
## diferença de AIC for pequena, a escolha fica com `dias`, que é
## íntegro e é a unidade das séries oficiais.
## =====================================================================
cat("\n===== 2) QUAL MEDIDA DE ESFORÇO USAR =====\n")
MARGEM_AIC <- 2      # abaixo disso é empate técnico
m_dias  <- if ("E4" %in% names(fits)) fits[["E4"]] else
  try(glmmTMB(monta_formula("captura", termos_E$E4, TRUE, "ldias"),
              family = tweedie(link = "log"), data = viagens), silent = TRUE)
m_horas <- try(glmmTMB(monta_formula("captura", termos_E$E4, TRUE, "lhoras"),
                       family = tweedie(link = "log"), data = viagens),
               silent = FALSE)
if (!inherits(m_dias, "try-error") && !inherits(m_horas, "try-error")) {
  cat(sprintf("  offset log(dias)  : AIC = %.1f\n", AIC(m_dias)))
  cat(sprintf("  offset log(horas) : AIC = %.1f\n", AIC(m_horas)))
  ## `isTRUE` protege contra AIC não finito (modelo que convergiu mal):
  ## nesse caso a comparação devolveria NA e o script pararia aqui.
  usar_horas <- isTRUE((AIC(m_dias) - AIC(m_horas)) > MARGEM_AIC)
} else usar_horas <- FALSE
OFFSET <- if (usar_horas) "lhoras" else "ldias"
UNID   <- if (usar_horas) "t/hora no mar" else "t/dia de pesca"
cat(sprintf("  -> adotado: offset(%s)  [%s]\n", OFFSET, UNID))
if (!usar_horas)
  cat("     (dias por integridade do campo; ver ressalva acima; Melhor ajuste)\n")

## =====================================================================
## 4) REFINAMENTO DA ESTRUTURA POR AIC (sem LRT)
## ---------------------------------------------------------------------
## QUAL ESTRUTURA É REFINADA, E POR QUÊ NÃO É "A MELHOR DA ESCADA":
## o modelo principal é deliberadamente o da TÁTICA DISCRETA (E4). Não
## porque ele ganhe sempre, mas porque as duas representações de tática
## precisam continuar sendo DUAS COISAS DIFERENTES até o fim do script:
## S2 (discreta) e S2b (contínua) são os dois lados da hipótese H2 e só 
# sãocomparados formalmente na seção 8. Se deixássemos a escada escolher,
## e ela escolhesse E5, os dois cenários (S2 e S2b) virariam o MESMO modelo com
## dois nomes — e a comparação H2 (hipotese discreta vs contínua) se 
# tornaria vazia. Precisamos manter a comparação discreta vs contínua viva
## A estrutura contínua é refinada junto, com a mesma base de termos,
## para que a comparação seja de TÁTICA e não de conjunto de covariáveis.
##
## O refinamento em si: partimos da estrutura cheia e testamos a REMOÇÃO
## de cada termo, um por vez. Um termo só sai se removê-lo MELHORAR o AIC
## em mais de MARGEM_AIC — a carga da prova é para retirar, não para
## manter. É a versão por AIC do backward clássico.
## Backward, e não forward: um termo cujo efeito só aparece depois de
## ajustar outro é sistematicamente perdido no forward.
## O fator temporal nunca entra na seleção: ele É o índice, e removê-lo
## seria remover o objeto da análise.
## =====================================================================
cat("\n===== 3) REFINAMENTO BACKWARD POR AIC =====\n")
usa_aleat <- TRUE                       # E4/E5 têm (1 | fbarco)
ajusta <- function(termos, dados = viagens, aleat = usa_aleat)
  try(glmmTMB(monta_formula("captura", termos, aleat, OFFSET),
              family = tweedie(link = "log"), data = dados), silent = FALSE)

## Os PCs formam UM bloco: são eixos da mesma PCA e remover PC2 mantendo
## PC3 não tem interpretação. Esta função devolve a lista de blocos
## removíveis de um conjunto de termos.
blocos_de <- function(termos) {
  pcs <- intersect(PCs, termos)
  outros <- setdiff(termos, PCs)
  bl <- as.list(outros); names(bl) <- outros
  if (length(pcs) > 0) bl[["PCs"]] <- pcs
  bl
}

backward_aic <- function(termos0) {
  termos <- termos0
  m <- ajusta(termos)
  if (inherits(m, "try-error")) return(list(termos = termos, fit = m))
  cat(sprintf("  partida: %s | AIC = %.1f\n",
              paste(c(fator_tempo, termos), collapse = " + "), AIC(m)))
  repeat {
    bl <- blocos_de(termos)
    if (length(bl) == 0) break
    aic_sem <- vapply(bl, function(b) {
      mr <- ajusta(setdiff(termos, b))
      if (inherits(mr, "try-error")) Inf else AIC(mr)
    }, numeric(1))
    ganho <- AIC(m) - aic_sem            # positivo = remover melhora
    cat("  dAIC ao remover: ",
        paste(sprintf("%s=%+.1f", names(ganho), ganho), collapse = "  "), "\n")
    ## guarda: se nenhum ajuste reduzido convergiu, `ganho` é todo NA e a
    ## comparação devolveria NA — o modelo atual fica como está.
    if (all(is.na(ganho)) || !isTRUE(max(ganho, na.rm = TRUE) > MARGEM_AIC)) {
      cat("  -> nenhum termo melhora o AIC ao sair: todos retidos\n"); break
    }
    pior <- names(which.max(ganho))
    cat(sprintf("  -> remove `%s` (AIC melhora %.1f)\n",
                paste(bl[[pior]], collapse = "+"), max(ganho, na.rm = TRUE)))
    termos <- setdiff(termos, bl[[pior]])
    m <- ajusta(termos)
  }
  list(termos = termos, fit = m)
}

cat("\n-- estrutura DISCRETA (E4: fator `alvo`) --\n")
sel_disc   <- backward_aic(termos_E$E4)
termos_sel <- sel_disc$termos
m_atual    <- sel_disc$fit
f_fix      <- monta_formula("captura", termos_sel, usa_aleat, OFFSET)
cat("Estrutura discreta final: ")
cat("Melhor modelo dentro da família discreta: ")
print(f_fix, showEnv = FALSE)

## A estrutura contínua reaproveita EXATAMENTE a mesma base de termos,
## trocando o fator `alvo` pelos eixos da PCA. É o que torna a
## comparação H2 uma comparação de REPRESENTAÇÃO DE TÁTICA.
#Se cada estrutura passasse pelo seu próprio backward, poderia 
#acontecer de o backward do contínuo decidir manter npesc_cat enquanto
#o backward do discreto decidiu tirar (hipoteticamente) — e aí, 
#quando você comparasse o AIC final de um contra o outro na seção 8, 
#a diferença estaria misturando dois efeitos: (a) discreto vs. contínuo E (b) 
##conjuntos de covariáveis diferentes. Você não conseguiria separar 
#qual dos dois motivos está gerando a diferença de AIC.
termos_cont <- c(setdiff(termos_sel, "alvo"), PCs)
cat("\n-- estrutura CONTÍNUA (E5: escores PC1..PCn) --\n")
m_cont <- ajusta(termos_cont)
if (!inherits(m_cont, "try-error"))
  cat(sprintf("  %s | AIC = %.1f\n",
              paste(c(fator_tempo, termos_cont), collapse = " + "), AIC(m_cont)))

## =====================================================================
## 5) AS DUAS DISTRIBUIÇÕES
## ---------------------------------------------------------------------
## Mesma estrutura, mesma resposta, mesmas linhas -> AIC/BIC comparáveis
## diretamente. Não há mais grupos de comparabilidade para administrar.
## =====================================================================
cat("\n===== 4) DISTRIBUIÇÕES =====\n")
dist_fits <- list()

## D1 — Tweedie.
dist_fits$D1_tweedie <- m_atual

## D2 — Hurdle-Gamma. O `ziformula` modela a probabilidade de a viagem
##      NÃO registrar cavala. Ele recebe o fator temporal, o trimestre e
##      a representação de tática que estiver no modelo: são essas que
##      descrevem a decisão de procurar (ou não) a espécie. Pôr o banco
##      aqui também seria defensável, mas multiplicaria parâmetros no
##      componente que tem menos informação.
##      A função abaixo monta o `ziformula` A PARTIR dos termos do modelo,
##      para que o componente de zeros acompanhe a estrutura escolhida
##      (com `alvo` no modelo discreto, com os PCs no contínuo).
monta_zi <- function(termos) {
  z <- c(fator_tempo, intersect("ftri", termos),
         intersect(c("alvo", PCs), termos))
  stats::as.formula(paste("~", paste(z, collapse = " + ")))
}
zi_f <- monta_zi(termos_sel)
cat(sprintf("  ziformula (probabilidade de zero): %s\n",
            paste(deparse(zi_f), collapse = "")))
dist_fits$D2_hurdle_gamma <- try(
  glmmTMB(f_fix, ziformula = zi_f, family = ziGamma(link = "log"),
          data = viagens), silent = FALSE)

## Por que a Tweedie pode ganhar aqui: o hurdle-Gamma reporta AIC/BIC = NA
## porque o componente `zi` (probabilidade de zero) não convergiu de forma
## confiável — Hessiana não-positiva-definida, provavelmente por
## quase-separação no termo `alvo`: a tática "macarellus" foi construída a
## partir da própria composição de captura, então dentro desse nível quase
## não existem viagens com captura zero, e o coeficiente logístico
## correspondente tenta ir para o infinito. Isso invalida a verossimilhança
## reportada, então o AIC não pode ser calculado (nem deveria ser usado).
## A Tweedie não sofre disso porque modela a massa de zeros e os valores
## positivos numa densidade só (via seu parâmetro de potência), sem precisar
## de um sub-modelo logístico separado que dependa de `alvo` — por isso ela
## tende a "vencer" por default quando o hurdle-Gamma quebra dessa forma,
## e isso é reportado como argumento a favor da Tweedie, não só um acaso.

dist_fits <- dist_fits[!vapply(dist_fits, inherits, logical(1), "try-error")]

tab_dist <- data.frame(
  modelo   = names(dist_fits),
  resposta = "captura (t), zeros incluídos",
  n        = vapply(dist_fits, function(m) nrow(model.frame(m)), numeric(1)),
  df       = vapply(dist_fits, function(m) attr(logLik(m), "df"), numeric(1)),
  AIC      = round(vapply(dist_fits, AIC, numeric(1)), 1),
  BIC      = round(vapply(dist_fits, BIC, numeric(1)), 1), row.names = NULL)
tab_dist$dAIC <- round(tab_dist$AIC - min(tab_dist$AIC), 1)
tab_dist$dBIC <- round(tab_dist$BIC - min(tab_dist$BIC), 1)
print(tab_dist, row.names = FALSE)

melhor_dist <- tab_dist$modelo[which.min(tab_dist$AIC)]
m_final <- dist_fits[[melhor_dist]]
cat(sprintf("\nMenor AIC: %s\n", melhor_dist))
cat("Confirmar no diagnóstico antes de aceitar — AIC mede ajuste, não\n")
cat("adequação da distribuição. Se o resíduo do vencedor for ruim e o do\n")
cat("outro for bom, o outro é o escolhido (e isso se reporta).\n")

## A partir daqui, TODO ajuste usa a distribuição vencedora — inclusive
## a estrutura contínua e o subconjunto do cenário S3. É o que garante
## que as diferenças entre cenários venham da ESTRUTURA e não de estarmos
## comparando famílias diferentes sem querer.
ajusta_final <- function(termos, dados = viagens, aleat = usa_aleat) {
  f <- monta_formula("captura", termos, aleat, OFFSET)
  if (melhor_dist == "D2_hurdle_gamma")
    try(glmmTMB(f, ziformula = monta_zi(termos), family = ziGamma(link = "log"),
                data = dados), silent = TRUE)
  else
    try(glmmTMB(f, family = tweedie(link = "log"), data = dados), silent = TRUE)
}

## Estrutura contínua reajustada na distribuição vencedora: é ESTE
## modelo que responde à hipótese H2, porque só difere do principal na
## representação da tática.
m_cont_final <- ajusta_final(termos_cont)
if (!inherits(m_cont_final, "try-error"))
  cat(sprintf("Estrutura contínua na mesma distribuição: AIC = %.1f (discreta: %.1f)\n",
              AIC(m_cont_final), AIC(m_final)))

## =====================================================================
## 6) DIAGNÓSTICO DE RESÍDUOS (DHARMa)
## ---------------------------------------------------------------------
## Em GLMM não-gaussiano o resíduo de Pearson engana: a relação
## média-variância não é constante e o gráfico "parece" ruim mesmo quando
## o modelo está certo. O DHARMa simula do modelo ajustado e transforma
## os resíduos para a escala uniforme, onde a leitura é a mesma para
## qualquer distribuição.
## O que cada teste responde:
##   KS         a distribuição assumida está certa?
##   dispersão  há sobre/subdispersão?
##   outliers   há mais extremos do que o modelo consegue gerar?
##   quantis    a variância é homogênea ao longo do predito?
##              (é o teste de homocedasticidade aqui)
##   zeros      o modelo gera a quantidade certa de zeros?
## Com n grande, p pequeno aparece por desvio trivial — por isso os
## gráficos são salvos: é neles que se vê se o desvio é grande ou só
## detectável.
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
cat("Leitura rápida: disp_ratio perto de 1 é bom; p pequeno em KS ou\n")
cat("quantis com n grande pede olhar o gráfico antes de condenar.\n")

png("diag_residuos.png", width = 26, height = 13, res = 300,
    antialias = AA, units = "cm")
set.seed(1); r_fin <- simulateResiduals(m_final, n = 250, plot = FALSE)
plot(r_fin)
dev.off()
cat("PNG salvo: diag_residuos.png (modelo escolhido)\n")

if (length(dist_fits) > 1) {
  outro <- setdiff(names(dist_fits), melhor_dist)[1]
  png("diag_residuos_alternativo.png", width = 26, height = 13, res = 300,
      antialias = AA, units = "cm")
  set.seed(1); plot(simulateResiduals(dist_fits[[outro]], n = 250, plot = FALSE))
  dev.off()
  cat(sprintf("PNG salvo: diag_residuos_alternativo.png (%s)\n", outro))
}

## =====================================================================
## 7) EXTRAÇÃO DO ÍNDICE (emmeans)
## ---------------------------------------------------------------------
## É aqui que "entrar como fator" se paga: o emmeans calcula a média
## marginal do fator temporal MEDIANDO os demais fatores — isto é, a taxa
## de captura esperada num banco médio, num trimestre médio, com
## tripulação média, numa mistura média de táticas. O que sobra é o
## efeito do tempo.
##
## Dois detalhes que mudam o resultado e passam despercebidos:
##  (1) `offset = 0`. Sem isso o emmeans usa a MÉDIA do offset e a série
##      sai por "viagem média", não por dia. Com isso, sai em unidade de
##      esforço, que é o que o JABBA espera.
##  (2) `weights`. "equal" dá o mesmo peso a cada nível — o certo para um
##      índice: queremos a média sobre ESTRATOS, não sobre as viagens que
##      por acaso foram amostradas. "proportional" pondera pelo n
##      observado e devolve parte do desbalanceamento (P6) que estamos
##      justamente tentando remover.
## =====================================================================
especificacao <- stats::as.formula(paste("~", fator_tempo))

extrai_indice <- function(m, nome, pesos = "equal", usa_offset = TRUE) {
  args <- list(object = m, specs = especificacao, weights = pesos)
  if (usa_offset) args$offset <- 0
  s <- as.data.frame(summary(do.call(emmeans, args)))
  data.frame(tempo = as.numeric(as.character(s[[fator_tempo]])),
             cenario = nome, indice_bruto = exp(s$emmean), se_log = s$SE,
             cv = sqrt(exp(s$SE^2) - 1), row.names = NULL)
}

## Índice de um modelo HURDLE. Este é o ponto técnico mais delicado do
## script: num hurdle, o `emmeans` devolve por padrão só o componente
## CONDICIONAL — a captura esperada DADO que houve captura. Usar isso
## como índice ignoraria completamente a mudança na probabilidade de
## encontrar a espécie, que é exatamente onde está o sinal aqui (a
## proporção de viagens com cavala cai de 24% para 7%).
## O valor esperado correto é:
##      E[Y] = (1 - p_zero) x E[Y | Y > 0]
## Na escala log:  log E[Y] = log(1 - p) + eta_cond
## e, pelo método delta, com p = plogis(z):
##      d log(1-p)/dz = -p   =>   contribuição ao erro-padrão = p * SE(z)
## Os dois componentes são tratados como independentes (premissa usual).
indice_hurdle <- function(m, nome, pesos = "equal") {
  s_c <- as.data.frame(summary(emmeans(m, especificacao, component = "cond",
                                       offset = 0, weights = pesos)))
  s_z <- as.data.frame(summary(emmeans(m, especificacao, component = "zi",
                                       weights = pesos)))
  p   <- stats::plogis(s_z$emmean)            # P(zero estrutural)
  idx <- (1 - p) * exp(s_c$emmean)
  se  <- sqrt(s_c$SE^2 + (p * s_z$SE)^2)
  data.frame(tempo = as.numeric(as.character(s_c[[fator_tempo]])),
             cenario = nome, indice_bruto = idx, se_log = se,
             cv = sqrt(exp(se^2) - 1), row.names = NULL)
}

## Despachante: usa o caminho certo conforme a família do modelo.
## Num glmmTMB sem `ziformula`, o campo guardado é `~0`; com hurdle, é a
## fórmula que passamos. O teste abaixo cobre também o caso `~1` (hurdle
## com probabilidade constante), que `all.vars()` sozinho deixaria passar
## como se não fosse hurdle.
indice_de <- function(m, nome) {
  zf <- m$modelInfo$allForm$ziformula
  eh_hurdle <- !is.null(zf) &&
    !identical(gsub("\\s", "", paste(deparse(zf), collapse = "")), "~0")
  if (eh_hurdle) indice_hurdle(m, nome) else extrai_indice(m, nome)
}
normaliza <- function(d) { d$indice <- d$indice_bruto / mean(d$indice_bruto); d }

cat("\n===== 6) ÍNDICES POR CENÁRIO =====\n")

## S1 — nominal: captura agregada / esforço agregado, sem modelo nenhum.
##      O CV é empírico (erro-padrão relativo da CPUE entre as viagens do
##      mesmo período), porque não há modelo de onde tirar variância.
esf_col <- if (OFFSET == "lhoras") "horas" else "dias"
S1 <- data.frame(tempo = ser$tempo, cenario = "S1 nominal",
                 indice_bruto = ser$cap_macarellus / ser[[esf_col]],
                 se_log = NA_real_, cv = NA_real_)
cv_emp <- tapply(viagens$captura / viagens[[esf_col]], viagens$tempo,
                 function(x) sd(x) / (mean(x) * sqrt(length(x))))
S1$cv <- as.numeric(cv_emp[as.character(S1$tempo)])
S1 <- normaliza(S1)

## S0 — modelo SEM a covariável de tática, com o resto igual. A diferença
##      S0 - S2 isola o efeito do direcionamento: é o "influence plot" de
##      Bentley et al. reduzido ao termo que interessa. Sem este cenário
##      não dá para afirmar que a correção veio da tática e não de outro
##      termo qualquer.
E_sem_alvo <- setdiff(termos_sel, c("alvo", PCs))
m_S0 <- ajusta_final(E_sem_alvo)
S0 <- if (!inherits(m_S0, "try-error"))
  normaliza(indice_de(m_S0, "S0 sem tática")) else NULL

## S2 — cenário principal: modelo selecionado (estrutura + distribuição
##      vencedoras), com tática discreta.
S2 <- normaliza(indice_de(m_final, "S2 corrigida (tática discreta)"))

## S2t — a MESMA estrutura na outra distribuição. Não é redundância: se o
##      índice muda de forma ao trocar Tweedie por hurdle, isso é
##      incerteza ESTRUTURAL e tem de ir para o texto; se não muda, é um
##      argumento forte de robustez.
outro_nome <- setdiff(names(dist_fits), melhor_dist)
S2t <- if (length(outro_nome) > 0)
  normaliza(indice_de(dist_fits[[outro_nome[1]]],
                      sprintf("S2t %s", outro_nome[1]))) else NULL

## S2b — tática CONTÍNUA (hipótese H2): os escores da PCA no lugar do
##      fator de cluster, mesma base de termos e MESMA distribuição.
S2b <- if (!inherits(m_cont_final, "try-error"))
  normaliza(indice_de(m_cont_final, "S2b tática contínua (PCs)")) else NULL

## S3 — esforço dirigido (H3): ajusta o modelo SÓ nas viagens da tática
##      com maior fração de cavala. `alvo_cavala` vem da parte 02 e é
##      escolhido pelo CENTRÓIDE, não pelo nome do grupo.
##      RESSALVA: se nenhuma tática for dominada pela cavala (o caso
##      desta série — a melhor tem ~15%), este cenário mede "esforço da
##      tática que mais encontra cavala", que é MENOS do que "esforço
##      dirigido à cavala". Continua sendo sensibilidade útil, mas não é
##      o cenário principal e tende a hiperestabilidade (ao restringir às
##      viagens que encontram a espécie, a queda fica achatada).
S3 <- NULL
if (exists("alvo_cavala") && alvo_cavala %in% levels(viagens$alvo)) {
  v3 <- viagens[viagens$alvo == alvo_cavala, ]
  v3[[fator_tempo]] <- droplevels(v3[[fator_tempo]])
  cobre <- nlevels(v3[[fator_tempo]]) >= 0.7 * nlevels(viagens[[fator_tempo]])
  if (cobre) {
    for (f in c("fbanco", "fbarco", "ftri", "npesc_cat"))
      if (f %in% names(v3)) v3[[f]] <- droplevels(v3[[f]])
    cat(sprintf("  S3: %d viagens da tática '%s' (%.0f%% de cavala no centróide)\n",
                nrow(v3), alvo_cavala, 100 * frac_cavala))
    m3 <- ajusta_final(termos_viaveis(v3, E_sem_alvo), dados = v3)
    if (!inherits(m3, "try-error"))
      S3 <- normaliza(indice_de(m3, "S3 esforço dirigido"))
  } else {
    cat("[AVISO] a tática da cavala não cobre períodos suficientes;\n")
    cat("        S3 fica sem estimativa — limitação esperada do subsetting.\n")
  }
}

indices <- do.call(rbind, Filter(Negate(is.null), list(S1, S0, S2, S2t, S2b, S3)))
cat("\nÍndices (média 1) e CV:\n")
print(transform(indices[, c("tempo", "cenario", "indice", "cv")],
                indice = round(indice, 3), cv = round(cv, 3)), row.names = FALSE)

## =====================================================================
## 8) TESTE DAS HIPÓTESES
## =====================================================================
cat("\n===== 7) HIPÓTESES =====\n")
casa <- function(a, b) b$indice[match(a$tempo, b$tempo)]
if (!is.null(S0)) {
  r_S0S2 <- cor(S0$indice, casa(S0, S2), use = "complete.obs")
  cat(sprintf("H0 — correlação entre índice com e sem tática: %.3f\n", r_S0S2))
  cat(sprintf("     %s\n", if (r_S0S2 > 0.98)
    "praticamente idênticos: a correção não muda nada (reportar!)" else
      "a tática desloca o índice de forma relevante"))
}
amp <- function(d) max(d$indice, na.rm = TRUE) / min(d$indice, na.rm = TRUE)
cat(sprintf("H1 — amplitude (máx/mín): nominal %.2f | corrigida %.2f\n",
            amp(S1), amp(S2)))
cat(sprintf("     correlação nominal x corrigida: %.3f\n",
            cor(S2$indice, casa(S2, S1), use = "complete.obs")))
cat("     amplitude menor na corrigida = parte da variação nominal era\n")
cat("     comportamento de frota, não abundância (o efeito esperado).\n")
if (!is.null(S2b)) {
  d_h2 <- AIC(m_cont_final) - AIC(m_final)   # negativo = contínua ganha
  cat(sprintf("H2 — discreta x contínua: r = %.3f | dAIC (contínua - discreta) = %+.1f\n",
              cor(S2$indice, casa(S2, S2b), use = "complete.obs"), d_h2))
  cat(sprintf("     %s\n", if (d_h2 < -MARGEM_AIC)
    "a representação CONTÍNUA (PCs) ajusta melhor — como em Winker et al. (2013)"
    else if (d_h2 > MARGEM_AIC)
      "a representação DISCRETA (cluster) ajusta melhor"
    else "empate técnico: as duas descrevem a tática igualmente bem"))
}
if (!is.null(S3))
  cat(sprintf("H3 — esforço dirigido x corrigida: r = %.3f | amplitude %.2f\n",
              cor(S2$indice, casa(S2, S3), use = "complete.obs"), amp(S3)))

## =====================================================================
## 9) FIGURA DOS ÍNDICES
## =====================================================================
cores_cen <- setNames(
  c(COR_AUX, COR_NEU, COR_MAC, COR_S2D, COR_S2, COR_S3),
  c("S1 nominal", "S0 sem tática", "S2 corrigida (tática discreta)",
    if (length(outro_nome) > 0) sprintf("S2t %s", outro_nome[1]) else "S2t",
    "S2b tática contínua (PCs)", "S3 esforço dirigido"))

png("indices_cenarios.png", width = 26, height = 14, res = 300,
    antialias = AA, units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3, 1), bty = "l",
          cex.main = 0.95, cex = 0.85)
cens <- unique(indices$cenario)
plot(NA, xlim = range(indices$tempo),
     ylim = c(0, max(indices$indice, na.rm = TRUE) * 1.12),
     xlab = rotulo_tempo, ylab = "Indice relativo (media = 1)",
     main = "A. Cenarios de indice")
for (cen in cens) {
  d <- indices[indices$cenario == cen, ]
  lines(d$tempo, d$indice, lwd = 2.4, col = cores_cen[cen])
  points(d$tempo, d$indice, pch = 19, cex = 0.8, col = cores_cen[cen])
}
legend("topright", cens, col = cores_cen[cens], lwd = 2.3, bty = "n", cex = 0.62)

lo <- S2$indice * exp(-1.96 * S2$se_log); hi <- S2$indice * exp(1.96 * S2$se_log)
plot(S2$tempo, S2$indice, type = "n", ylim = c(0, max(hi, na.rm = TRUE) * 1.05),
     xlab = rotulo_tempo, ylab = "Indice (media = 1)",
     main = "B. Indice corrigido com IC 95%")
polygon(c(S2$tempo, rev(S2$tempo)), c(lo, rev(hi)),
        col = adjustcolor(COR_MAC, 0.18), border = NA)
lines(S2$tempo, S2$indice, lwd = 2.6, col = COR_MAC)
points(S2$tempo, S2$indice, pch = 19, col = COR_MAC)
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
## Piso de CV em 0,20: o CV que sai do modelo é de processo estatístico e
## ignora erro de processo, erro de reporte e a incerteza da própria
## INFERÊNCIA DE TÁTICA (que não tem variância nenhuma no cálculo).
## Entregar CV de 0,04 ao JABBA faria o modelo confiar no índice muito
## mais do que ele merece.
## =====================================================================
PISO_CV <- 0.20
saida_tempo <- sort(unique(indices$tempo))
cen_export <- c("S1 nominal" = "cpue_nominal",
                "S2 corrigida (tática discreta)" = "cpue_corrigida",
                "S2b tática contínua (PCs)" = "cpue_corrigida_pcs")
jabba_idx <- data.frame(tempo = saida_tempo)
jabba_cv  <- data.frame(tempo = saida_tempo)
for (cen in names(cen_export)) {
  if (!cen %in% indices$cenario) next
  d <- indices[indices$cenario == cen, ]
  jabba_idx[[cen_export[cen]]] <- d$indice[match(saida_tempo, d$tempo)]
  jabba_cv[[cen_export[cen]]]  <- d$cv[match(saida_tempo, d$tempo)]
}
names(jabba_idx)[1] <- names(jabba_cv)[1] <- if (fator_tempo == "fano") "Yr" else "Mes"
jabba_cv[, -1] <- lapply(jabba_cv[, -1, drop = FALSE],
                         function(x) pmax(x, PISO_CV, na.rm = TRUE))

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
cat("2. Reportar as tabelas de seleção (estruturas e distribuições) e de\n")
cat("   diagnóstico, inclusive quando os pressupostos forem violados.\n")
cat("3. Declarar o que a tática é: variável INFERIDA da composição da\n")
cat("   captura, não observada. A incerteza dessa inferência não está no\n")
cat("   CV — daí o piso de 0,20.\n")
cat("4. Declarar as covariáveis descartadas e por quê (profundidade,\n")
cat("   tipo de embarcação, ilha de desembarque) — todas por problema de\n")
cat("   dado, não por não serem significativas.\n")
cat(sprintf("5. Esta série cobre %s. Ela é POSTERIOR à mudança de alvo de\n",
            paste(range(viagens$ano), collapse = "-")))
cat("   2014, então descreve a dinâmica da cavala já na condição de\n")
cat("   captura acompanhante. Não é a mesma pergunta que a série\n")
cat("   histórica longa responde, e as duas não devem ser soldadas numa\n")
cat("   única série 'corrigida' contínua.\n")









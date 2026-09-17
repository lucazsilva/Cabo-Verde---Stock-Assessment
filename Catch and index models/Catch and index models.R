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
# L10. TRIPULAÇÃO em CLASSES, não contínua. `Numero_pescadores` é proxy
#     de poder de pesca. Entrar como contínua impõe relação monotônica e
#     linear na escala do log, o que não há motivo para supor; e a cauda
#     (1-9 e 44-50 pescadores) puxa a reta. Classes deixam a forma
#     aparecer e são robustas a esses extremos. Cortes em <=12 / 13-15 /
#     16-17 / >=18 (quartis aproximados, ~19/29/30/23% das viagens).
# L11. EMBARCAÇÃO = CÓDIGO (`Embarcacao`), não nome. Na varredura, 21
#     nomes aparecem com 2-3 códigos diferentes e as séries se
#     SOBREPÕEM no tempo (ex.: "BARBARA" código 18 de 2019-2025 e código
#     3413 de 2020-2025, simultâneos). São barcos distintos com o mesmo
#     nome. Usar o nome fundiria históricos de barcos diferentes; o
#     código é 1:1 com o nome (nenhum código tem dois nomes).
# L12. `Tipo_embarcacao` (I / S) NÃO É USADO. Parecia um estrato de
#     frota, mas a varredura mostrou três coisas: (a) 30 dos 34 barcos
#     que aparecem como "S" também aparecem como "I" — é o MESMO barco
#     com dois códigos em viagens diferentes; (b) "S" praticamente só
#     existe em 2019 (452 viagens) e some depois (12, 7, 1, 0, 0, 0), ou
#     seja, está quase perfeitamente confundido com o ano — justo o ano
#     de maior captura de cavala; (c) os perfis de I e S em 2019 são
#     iguais (mesmos dias no mar, mesma tripulação, mesma captura
#     média). É um código administrativo de registro, não um tipo de
#     embarcação. Usá-lo como covariável roubaria sinal do efeito de ano.
# L13. Viagens com MAIS DE UMA ARTE (43 viagens) são descartadas: os
#     dias de mar delas cobrem também o que foi pescado com outra arte,
#     então o esforço atribuído ao cerco ficaria inflado.
#======================================================================

## =====================================================================
## 0) PARÂMETROS QUE VOCÊ PODE QUERER MEXER
## =====================================================================
ARQUIVO        <- "INDUSTRIAL_2019_2025_atualizado_17.09.2026.csv"
ARTE_ALVO      <- "REDE DE CERCO"        # decisão L5
MIN_VIAG_BANCO <- 30                     # decisão L8 (45 bancos + "OUTROS")
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
















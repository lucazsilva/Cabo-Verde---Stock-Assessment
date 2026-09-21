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
#install.packages("rjags")
library(rjags)
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
# PADRONIZAÇÃO DE CPUE — Decapterus macarellus (cavala preta), Cabo Verde Frota industrial — dados do IMar, série 2015-2025
# (1) A série agora é 2015-2025 (10 anos; 2018 NÃO EXISTE em nenhum dos dois extratos recebidos — decisão L14).
#     O efeito de ANO — que É o índice de abundância — é estimável. A parte 02 troca sozinha o fator temporal de mês para ano.
# (2) PROFUNDIDADE foi REMOVIDA do pipeline. Ela é  lida, mas não vai para a tabela final:
#     28% das viagens têm 0 (= ausência de registro, não zero metro) e o ausente não é aleatório. Ver decisão L6.
# (3) MÊS vira TRIMESTRE como fator sazonal (decisão L9) e TRIPULAÇÃO vira FATOR em classes (decisão L10). Motivo em cada decisão.
# (4) A identidade da embarcação passa a ser o CÓDIGO e não o nome  (decisão L11) — Os nomes se repetem entre barcos
#     diferentes que operam ao mesmo tempo.
# (5) `Tipo_embarcacao` foi INVESTIGADO e DESCARTADO como covariável (decisão L12).
# (6) NOVO: o espaço entra como ILHA DO BANCO (decisão L15) e não mais como banco individual agrupado por limiar de viagens.
#     São 244 bancos distintos nas viagens de cerco — qualquer limiar útil jogava ~40% das viagens num nível "OUTROS"
#     que não significa lugar nenhum. A ilha do banco tem 5 níveis, todos presentes em 9-10 dos 10 anos.
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
# L7. Local de pesca: a origem é `Nome_banco_pesca` (onde se pescou), nunca o porto/ilha de
# desembarque DA VIAGEM (onde se descarregou). São coisas diferentes e só a primeira é covariável
# de densidade. A ilha de desembarque da viagem continua FORA do modelo como variável própria,
# porque a cobertura dela muda ao longo da série e isso a deixa confundida com o ano. Ela é usada
# apenas como INSUMO para construir a ilha do banco (L15), o que é outra coisa — ver lá.
# L8. [SUBSTITUÍDA POR L15] O agrupamento de bancos raros em "OUTROS" por limiar de viagens foi
# abandonado. Motivo empírico: são 244 bancos distintos nas viagens de cerco, e a distribuição é tão
# assimétrica que não existe limiar bom. Com MIN_VIAG_BANCO = 30 sobram 51 níveis (grade do emmeans
# grande demais e células vazias); com 100 sobram 14 níveis, mas 40,5% das viagens caem em "OUTROS"
# — um nível que mistura bancos de ilhas diferentes e não significa lugar nenhum. O parâmetro
# continua no script porque `banco_gr` ainda é usado no efeito ALEATÓRIO opcional (E4b), onde o
# encolhimento resolve a esparsidade sem precisar de limiar agressivo.
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
# L13. Viagens com MAIS DE UMA ARTE (49 viagens) são descartadas: os dias de mar delas cobrem também
#o que foi pescado com outra arte, então o esforço atribuído ao cerco ficaria inflado.
#
# L14. QUATRO FONTES, COLADAS POR ANO, CADA UMA NO SEU PEDAÇO DA SÉRIE. Nenhuma sozinha cobre
#      1989-2025, e elas não têm a mesma granularidade nem o mesmo formato:
#      - "Cavala_desembarques-esforço_pesca industrial_1989-2014.csv" : AGREGADA. Uma linha por ANO,
#        com desembarque e esforço já somados. Não tem viagem, embarcação, tripulação nem composição
#        de captura. Nível = "anual".
#      - "INDUSTRIAL_2015_2017_Sem_validação.csv"  : nível de viagem, formato de referência, 2015-2017.
#      - "INDUSTRIAL_2018_formato_diferente.csv"   : nível de viagem, 14 colunas em vez de 39 (ver L17).
#      - "INDUSTRIAL_2019_2025_atualizado_17.09.2026.csv" : nível de viagem, formato de referência, é
#        ESTE que define o formato-alvo para onde todo o resto é traduzido.
#      REGRA: cada arquivo entra só nos anos em que é a melhor fonte. O script confere sozinho que não
#      há colisão de `Amostragem` entre blocos e para (`stopifnot`) se houver — id repetido fundiria
#      duas viagens diferentes numa só lá na frente.
#      O produto é `esforco_completo` (39 colunas do padrão + 4 de procedência: fonte, nivel,
#      engenhos_agregados, usar_na_serie), escrito em CSV e XLSX. A coluna `nivel` é o que impede o
#      histórico agregado de entrar por engano nas etapas que exigem viagem.
#      COBERTURA RESULTANTE: 1989-2025 sem nenhum ano faltando — 2018 agora existe (era buraco na
#      versão anterior), mas fica fora da série por L17.
#
# L16. FRONTEIRA DA TROCA DE ALVO: pré ≤ 2014, pós ≥ 2015 (`ANO_CORTE_ALVO`). Coincide com o que já
#      estava no cabeçalho ("hiperdepleção aparente depois de 2014") e com a troca de regime dos
#      dados: até 2014 a fonte é agregada, de 2015 em diante é viagem a viagem. É o corte que separa
#      os cenários nominais C2 (pré) e C3 (pós).
#
# L17. 2018 É TRADUZIDO PARA O FORMATO PADRÃO, MAS FICA FORA DA SÉRIE (`INCLUIR_2018 = FALSE`).
#      (mudar 2018 em função dos demais, nunca o contrário) — o que existe vai
#      para a coluna equivalente; o que não existe fica em branco; o que só existe em 2018 é
#      descartado por não ser informação perene do banco:
#        ILHA->Nome_ilha | (coluna sem nome)->Nome_embarcacao | DATA PARTIDA/CHEGADA->Data_partida/
#        Data_chegada | MÊS->Mes | ANO->Ano | PORTO EMBARQUE->Nome_porto_armamento | PORTO
#        DESEMBARQUE->Nome_porto_desembarque | ZONA DE PESCA->Nome_banco_pesca | ENGENHO->
#        Nome_engenho | GRUPO->Grupo | ESPÉCIE CAPTURADA->Especie | NOME CIENTÍFICO->Nome_cientifico
#        | QUANTIDADE CAPTURADA (Kg)->Quantidade.   DESCARTADA: DIA (redundante com Data_chegada).
#      DERIVADOS: `Amostragem` não existe e é construído como barco+partida+chegada (prefixo "2018_",
#      porque não é comparável com os ids do IMar); `Num_dias` não existe e sai de max(chegada-
#      partida, 1) — regra calibrada no arquivo de referência, onde reproduz 100% dos registros.
#      EM BRANCO POR NÃO EXISTIR EM 2018: Embarcacao (o CÓDIGO; só há o nome), Numero_pescadores,
#      Num_horas, Profundidade, Preco, Valor, Familia, Genero.
#      POR QUE FICA FORA DA SÉRIE: 2018 destoa dos vizinhos de forma que não dá para atribuir a
#      abundância com o que se sabe hoje. CPUE nominal de 1,03 t/dia contra 0,28 (2017) e 0,38 (2019);
#      a cavala é 27,7% da captura contra 12,3% e 18,2%; a captura por registro tem mediana de 4,0 t
#      contra 1,5 t em 2019-2025. O esforço e o número de viagens são normais (1.326 viagens, 2.153
#      dias) — o que muda é só a captura. Como o levantamento é de outra origem e outro formato, a
#      hipótese de artefato de amostragem não pode ser descartada, e um artefato aqui viraria um pico
#      de abundância espúrio bem no meio da série. Fica na tabela (auditável, e a tradução está
#      pronta), mas `usar_na_serie = FALSE` e nenhum cenário o usa. Basta `INCLUIR_2018 <- TRUE`
#      para reincorporá-lo se o IMar confirmar os números.
#      FALTA TAMBÉM, mesmo se voltar: sem tripulação e sem código de embarcação, 2018 não pode entrar
#      na CPUE PADRONIZADA — só na nominal.
#
# L18. QUAL ESFORÇO USAR DE 1989-2014. A planilha traz duas colunas: "Rede cerco" e "Total"
#      (cerco + linha de mão). Usamos a do CERCO, para bater com o filtro de arte do resto da série.
#      Nos anos em que a própria planilha anota "engenhos agregados" (1994-1999 e 2014) o número já
#      mistura artes; ali cai para o total e a linha fica marcada com `engenhos_agregados = TRUE`,
#      para a ressalva aparecer no texto em vez de sumir na média. 2013 não tem esforço nenhum
#      ("Sem esforço" na planilha): fica com `usar_na_serie = FALSE` e sem CPUE — a captura de
#      2013 (2.210 t) continua registrada, só não vira índice.
#
# L15. ESPAÇO = ILHA DO BANCO DE PESCA, não o banco individual. Cada BANCO recebe, de uma vez por
#      todas, o nome da ilha de desembarque MAIS FREQUENTE entre as viagens que pescaram nele; essa
#      etiqueta vira a variavel `ilha_banco` (5 níveis).
#      POR QUE ISSO NÃO É A MESMA COISA QUE USAR A ILHA DE DESEMBARQUE DA VIAGEM (o que L7 proíbe):
#      a etiqueta é uma propriedade FIXA DO BANCO, calculada uma vez sobre a série inteira. Duas
#      viagens ao mesmo banco recebem a mesma ilha mesmo que tenham desembarcado em portos
#      diferentes. O porto onde o barco escolheu descarregar — que é o que muda com a cobertura de
#      amostragem ao longo dos anos — deixa de entrar na covariável. O que entra é "em que zona do
#      arquipélago fica este pesqueiro", que é geografia e não muda com o ano.
#      POR QUE A ILHA E NÃO O BANCO: ver L8. E porque 5 níveis bem povoados, presentes em quase todos
#      os anos, estimam coeficiente; 244 níveis (ou 14 + um "OUTROS" com 40% das viagens) não.
#      QUALIDADE DO MAPEAMENTO (impressa pelo script): 197 dos 244 bancos são "puros" — todas as
#      viagens deles desembarcaram na mesma ilha. Ponderando por viagem, 91,2% das viagens estão na
#      ilha modal do seu banco.
#      RESSALVA QUE VAI PARA O TEXTO: os ~9% restantes são topônimos genéricos que se repetem em
#      várias ilhas de Cabo Verde — TARRAFAL (existe em Santiago, S. Nicolau e S. Antão), SANTA MARIA,
#      CALHETA, BAIA, PONTA. Nesses casos a atribuição modal força uma ilha só e erra em parte das
#      viagens. O script imprime a lista dos bancos com pureza < 80% para que ela possa ser conferida
#      com quem conhece a pescaria — é o tipo de coisa que uma tabela de coordenadas dos bancos
#      resolveria de vez, e que vale pedir ao IMar.
#======================================================================

## =====================================================================
## 0) PARÂMETROS E FONTES
## =====================================================================
## As QUATRO fontes do IMar. Cada uma entra só nos anos em que é a melhor
## (ou a única) fonte — ver decisão L14. Quando o IMar mandar atualização,
## é aqui que se mexe.
ARQ_1517 <- "INDUSTRIAL_2015_2017_Sem_validação.csv"
ARQ_1925 <- "INDUSTRIAL_2019_2025_atualizado_17.09.2026.csv"
ARQ_2018 <- "INDUSTRIAL_2018_formato_diferente.csv"
ARQ_HIST <- "Cavala_desembarques-esforço_pesca industrial_1989-2014.csv"

ARTE_ALVO      <- "REDE DE CERCO"        # decisão L5
ESPECIE_FOCO   <- "DECAPTERUS MACARELLUS"
CORTES_NPESC   <- c(0, 12, 15, 17, Inf)  # decisão L10
ANO_CORTE_ALVO <- 2014                   # decisão L16 — fronteira da troca de alvo
INCLUIR_2018   <- FALSE                  # decisão L17 — 2018 fica FORA da série

## MIN_VIAG_BANCO NÃO define mais a covariável espacial (ver L8/L15). Ele
## só agrupa a cauda de bancos raros para o efeito ALEATÓRIO opcional
## `(1 | fbanco)` da estrutura E4b. Como ali o encolhimento já cuida dos
## níveis pequenos, o limiar pode ser brando.
MIN_VIAG_BANCO <- 10
PUREZA_ALERTA  <- 0.80

## =====================================================================
## 0.1) FERRAMENTAS DE LEITURA E ESCRITA (decisões L1-L4)
## =====================================================================
## L1-L4: latin1 -> UTF-8, sem aspas, separador ";", decimal ",", e
## `trimws` em TUDO (o arquivo é de largura fixa disfarçada de CSV).
le_csv_imar <- function(caminho, ...) {
  linhas <- readLines(caminho, warn = FALSE)
  linhas <- iconv(linhas, from = "latin1", to = "UTF-8", sub = "?")
  d <- read.csv2(text = linhas, quote = "", stringsAsFactors = FALSE,
                 strip.white = TRUE, ...)
  names(d) <- trimws(names(d))
  for (j in seq_along(d)) if (is.character(d[[j]])) d[[j]] <- trimws(d[[j]])
  d
}

## Remoção de acentos por substituição de BYTES. Não usa `iconv TRANSLIT`
## (que em algumas máquinas transforma "ç" em "?") nem `chartr` (que
## quebra com multibyte). Como as strings já vieram convertidas para
## UTF-8, cada acentuado é uma sequência fixa de bytes e a troca por byte
## dá o mesmo resultado no Windows, no Linux e no Mac.
sem_acento <- function(x) {
  x <- as.character(x); Encoding(x) <- "UTF-8"
  de <- c("á","à","â","ã","ä","å","é","è","ê",
          "ë","í","ì","î","ï","ó","ò","ô","õ",
          "ö","ú","ù","û","ü","ç","ñ",
          "Á","À","Â","Ã","Ä","Å","É","È","Ê",
          "Ë","Í","Ì","Î","Ï","Ó","Ò","Ô","Õ",
          "Ö","Ú","Ù","Û","Ü","Ç","Ñ")
  para <- c("a","a","a","a","a","a","e","e","e","e","i","i","i","i","o","o","o","o","o",
            "u","u","u","u","c","n","A","A","A","A","A","A","E","E","E","E","I","I","I","I",
            "O","O","O","O","O","U","U","U","U","C","N")
  Encoding(de) <- "UTF-8"
  for (i in seq_along(de)) x <- gsub(de[i], para[i], x, fixed = TRUE, useBytes = TRUE)
  Encoding(x) <- "UTF-8"; x
}
norm_txt <- function(x) toupper(trimws(gsub("\\s+", " ", sem_acento(x), useBytes = TRUE)))

## Escritor de CSV byte-exato. `write.csv(fileEncoding = "UTF-8")`
## reconverte as strings a partir do encoding nativo e, numa máquina cujo
## locale não seja UTF-8, TRUNCA o campo no primeiro acentuado: "S. ANTÃO"
## saía como `"S. ANT` sem fechar aspas e o arquivo ficava irrecuperável
## (relia 8.353 das 26.932 linhas, em silêncio). Aqui as linhas são
## montadas à mão — aspas internas dobradas, como manda o padrão CSV — e
## gravadas como BYTES, sem o R reinterpretar nada.
escreve_csv_utf8 <- function(d, caminho) {
  campo <- function(x) {
    x <- enc2utf8(as.character(x))
    out <- paste0('"', gsub('"', '""', x, fixed = TRUE, useBytes = TRUE), '"')
    out[is.na(x)] <- "NA"
    out
  }
  linhas <- c(paste(campo(names(d)), collapse = ","),
              do.call(paste, c(lapply(d, campo), sep = ",")))
  con <- file(caminho, open = "wb"); on.exit(close(con))
  writeLines(enc2utf8(linhas), con, useBytes = TRUE)
  invisible(nrow(d))
}

## As 39 colunas do formato de referência (o da planilha 2019-2025).
COLS_PADRAO <- c(
  "Amostragem","Ano","Mes","Grupo","Familia","Genero","Nome_cientifico","Especie",
  "Engenho","Nome_engenho","Quantidade","Preco","Valor","Data_amostragem","Tipo_pesca",
  "Regiao","Ilha","Nome_ilha","Concelho","Nome_concelho","Localidade","Nome_localidade",
  "Tipo_embarcacao","Embarcacao","Nome_embarcacao","Porto_armamento","Nome_porto_armamento",
  "Porto_desembarque","Nome_porto_desembarque","Numero_pescadores","Data_partida",
  "Hora_partida","Data_chegada","Num_dias","Hora_chegada","Num_horas","Banco_pesca",
  "Nome_banco_pesca","Profundidade_pesca_engenho1")

## Encaixa qualquer fonte no formato de referência: o que existe é
## copiado, o que não existe fica NA, e o que só existe naquela fonte é
## descartado (não é informação perene do banco). As 4 colunas de
## procedência ficam DEPOIS das 39, para o formato original continuar
## reconhecível.
molda_padrao <- function(d, fonte, nivel, agregados = FALSE, usar = TRUE) {
  out <- as.data.frame(matrix(NA_character_, nrow = nrow(d), ncol = length(COLS_PADRAO)),
                       stringsAsFactors = FALSE)
  names(out) <- COLS_PADRAO
  for (cc in intersect(names(d), COLS_PADRAO)) out[[cc]] <- as.character(d[[cc]])
  out$fonte <- fonte; out$nivel <- nivel
  out$engenhos_agregados <- agregados; out$usar_na_serie <- usar
  out
}

## =====================================================================
## 0.2) MONTAGEM DA SÉRIE COMPLETA 1989-2025 (decisões L14, L17, L18)
## =====================================================================
cat("======================================================\n")
cat("PARTE 0 — SÉRIE DE ESFORÇO COMPLETA 1989-2025\n")
cat("======================================================\n")

## --- (a) 2015-2017 e 2019-2025: já estão no formato de referência -----
b1517 <- le_csv_imar(ARQ_1517); b1925 <- le_csv_imar(ARQ_1925)
b1517 <- b1517[b1517$Ano %in% as.character(2015:2017), , drop = FALSE]
b1925 <- b1925[b1925$Ano %in% as.character(2019:2025), , drop = FALSE]
col_1517 <- unique(b1517$Amostragem[b1517$Amostragem != ""])
col_1925 <- unique(b1925$Amostragem[b1925$Amostragem != ""])
cat(sprintf("\n[a] 2015-2017: %d linhas | 2019-2025: %d linhas\n", nrow(b1517), nrow(b1925)))
cat(sprintf("    ids `Amostragem` colididos entre os dois blocos: %d %s\n",
            length(intersect(col_1517, col_1925)),
            if (length(intersect(col_1517, col_1925)) == 0) "(OK)" else "<<< PARE"))
stopifnot(length(intersect(col_1517, col_1925)) == 0)
P_1517 <- molda_padrao(b1517, "IMar 2015-2017 (sem validacao)", "viagem")
P_1925 <- molda_padrao(b1925, "IMar 2019-2025 (atualizado)",    "viagem")

## --- (b) 2018: formato diferente, traduzido para o de referência ------
## O arquivo de 2018 tem 14 colunas contra 39. A regra (decisão L17) é
## traduzir 2018 PARA o formato maior, nunca o contrário: o que existe
## entra na coluna equivalente, o que não existe fica em branco, e o que
## só existe em 2018 é descartado.
b18 <- le_csv_imar(ARQ_2018, check.names = FALSE)
nm <- names(b18); nm[nm == "" | is.na(nm) | nm %in% c("NA", "X")] <- "BARCO"; names(b18) <- nm
pega <- function(d, alvo) {
  i <- which(norm_txt(names(d)) == norm_txt(alvo))
  if (!length(i)) rep("", nrow(d)) else as.character(d[[i[1]]])
}
d18 <- data.frame(
  ilha = pega(b18,"ILHA"),               barco = pega(b18,"BARCO"),
  dpart= pega(b18,"DATA PARTIDA"),       dcheg = pega(b18,"DATA CHEGADA"),
  mes  = pega(b18,"MES"),                ano   = pega(b18,"ANO"),
  pemb = pega(b18,"PORTO EMBARQUE"),     pdes  = pega(b18,"PORTO DESEMBARQUE"),
  zona = pega(b18,"ZONA DE PESCA"),      eng   = pega(b18,"ENGENHO"),
  grupo= pega(b18,"GRUPO"),              esp   = pega(b18,"ESPECIE CAPTURADA"),
  sci  = pega(b18,"NOME CIENTIFICO"),
  qtd  = pega(b18,"QUANTIDADE CAPTURADA (Kg)"), stringsAsFactors = FALSE)

## `Num_dias` não existe em 2018 e é DERIVADO das datas. A regra foi
## calibrada no arquivo de referência: lá, `Num_dias = max(chegada -
## partida, 1)` reproduz 100% dos registros (viagem que sai e volta no
## mesmo dia conta como 1 dia). Data invertida vira NA, não 0.
dp <- as.Date(d18$dpart, format = "%m/%d/%Y")
dc <- as.Date(d18$dcheg, format = "%m/%d/%Y")
dur <- as.numeric(dc - dp)
d18$num_dias <- ifelse(is.na(dur) | dur < 0, NA_real_, pmax(dur, 1))

## Engenho em 2018 vem com caixa inconsistente ("Rede de cerco", "rede de
## cerco", "cerco") — normalizado para o vocabulário do formato padrão.
eng <- norm_txt(d18$eng)
eng[grepl("CERCO", eng) & grepl("LINHA", eng)] <- "CERCO E LINHA"
eng[grepl("CERCO", eng) & eng != "CERCO E LINHA"] <- ARTE_ALVO
eng[grepl("^LINHA", eng) & !grepl("VARA", eng)] <- "LINHA A MAO"
eng[eng == ""] <- NA_character_

## ILHA: 2018 escreve por extenso ("São Vicente"); o padrão usa a forma
## abreviada e acentuada ("S. VICENTE", "S. ANTÃO"). A grafia canônica
## NÃO é escrita aqui à mão: ela é LIDA dos arquivos padrão e localizada
## pela forma sem acento. Assim o nível do 2018 recebe exatamente os
## mesmos bytes do resto da tabela, sem depender de escapes `\u` — que,
## num locale não-UTF-8, o R guarda como o texto literal "<c3><83>" e
## fazem a correspondência falhar em silêncio.
ilhas_padrao <- unique(c(P_1517$Nome_ilha, P_1925$Nome_ilha))
ilhas_padrao <- ilhas_padrao[!is.na(ilhas_padrao) & ilhas_padrao != ""]
canon_ilha <- setNames(ilhas_padrao, norm_txt(ilhas_padrao))
ALIAS_ILHA <- c("SAO VICENTE"="S. VICENTE", "SANTO ANTAO"="S. ANTAO",
                "SAO NICOLAU"="S. NICOLAU", "SAL"="SAL", "SANTIAGO"="SANTIAGO",
                "BOA VISTA"="BOA VISTA", "MAIO"="MAIO", "FOGO"="FOGO", "BRAVA"="BRAVA")
ilha18  <- norm_txt(d18$ilha)
chave18 <- ifelse(ilha18 %in% names(ALIAS_ILHA), ALIAS_ILHA[ilha18], ilha18)
ilha18  <- ifelse(chave18 %in% names(canon_ilha), canon_ilha[chave18], ilha18)

## Identificador de viagem: 2018 não tem `Amostragem`. Construído a
## partir de barco + data de partida + data de chegada, que é o que
## identifica uma maré. NÃO é comparável com os ids do resto da série —
## por isso recebe o prefixo "2018_".
d18$viagem <- paste(norm_txt(d18$barco), d18$dpart, d18$dcheg, sep = "|")
id18 <- paste0("2018_", match(d18$viagem, unique(d18$viagem)))

cat(sprintf("\n[b] 2018: %d linhas -> %d viagens (barco + partida + chegada)\n",
            nrow(d18), length(unique(id18))))
cat(sprintf("    ilhas casadas com a grafia do padrão: %d de %d\n",
            sum(ilha18 %in% ilhas_padrao), length(ilha18)))
cat(sprintf("    datas invertidas (chegada < partida): %d -> Num_dias = NA\n",
            sum(!is.na(dur) & dur < 0)))
cat(sprintf("    engenho normalizado: %s\n",
            paste(sprintf("%s=%d", names(table(eng)), table(eng)), collapse = " | ")))
cat("    DESCARTADA (só existe em 2018, não é perene): DIA\n")
cat("    SEM EQUIVALENTE (ficam em branco): Embarcacao (código), Numero_pescadores,\n")
cat("                      Num_horas, Profundidade, Preco, Valor, Familia, Genero\n")
cat(sprintf("    usar_na_serie = %s (decisão L17)\n", INCLUIR_2018))

P_2018 <- molda_padrao(data.frame(
  Amostragem = id18, Ano = d18$ano, Mes = d18$mes, Grupo = norm_txt(d18$grupo),
  Nome_cientifico = norm_txt(d18$sci), Especie = norm_txt(d18$esp),
  Nome_engenho = eng, Quantidade = d18$qtd, Tipo_pesca = "INDUSTRIAL",
  Nome_ilha = ilha18, Nome_embarcacao = norm_txt(d18$barco),
  Nome_porto_armamento = norm_txt(d18$pemb),
  Nome_porto_desembarque = norm_txt(d18$pdes),
  Data_partida = d18$dpart, Data_chegada = d18$dcheg, Num_dias = d18$num_dias,
  Nome_banco_pesca = norm_txt(d18$zona), stringsAsFactors = FALSE),
  "IMar 2018 (formato diferente)", "viagem", FALSE, INCLUIR_2018)

## --- (c) 1989-2014: série histórica AGREGADA --------------------------
## Esta planilha não é de viagem: é uma linha por ANO, com desembarque e
## esforço já somados. Ela não tem — e não pode ter — composição de
## captura, embarcação, tripulação ou banco. Entra no formato de
## referência com as colunas que existem preenchidas e o resto em branco,
## marcada com nivel = "anual" para que nenhuma etapa de viagem a use por
## engano.
## DECISÃO L18 (qual coluna de esforço): a planilha traz "Rede cerco" e
## "Total" (cerco + linha de mão). Usamos o cerco, para bater com o
## filtro de arte do resto da série; nos anos em que a própria planilha
## anota "engenhos agregados" o número já mistura artes e aí cai para o
## total, com a marca `engenhos_agregados = TRUE`. 2013 não tem esforço
## nenhum e fica com usar_na_serie = FALSE.
lin <- iconv(readLines(ARQ_HIST, warn = FALSE), "latin1", "UTF-8", sub = "?")
h <- read.csv2(text = lin[-(1:3)], header = FALSE, quote = "",
               stringsAsFactors = FALSE, strip.white = TRUE)
names(h)[1:8] <- c("ano","lm_des","lm_esf","rc_des","rc_esf","tot_des","tot_esf","obs")
h <- h[!is.na(suppressWarnings(as.numeric(h$ano))), , drop = FALSE]
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)))))
for (cc in c("ano","lm_des","lm_esf","rc_des","rc_esf","tot_des","tot_esf")) h[[cc]] <- num(h[[cc]])
h$obs       <- trimws(as.character(h$obs))
h$agregado  <- grepl("agregad", h$obs, ignore.case = TRUE)
h$desem_t   <- ifelse(!is.na(h$rc_des) & !h$agregado, h$rc_des, h$tot_des)
h$esforco_d <- ifelse(!is.na(h$rc_esf) & !h$agregado, h$rc_esf, h$tot_esf)
h$usar      <- !is.na(h$desem_t) & !is.na(h$esforco_d)
cat(sprintf("\n[c] histórico 1989-2014: %d anos | engenhos agregados: %d | sem esforço: %d\n",
            nrow(h), sum(h$agregado), sum(!h$usar)))
if (any(!h$usar))
  cat(sprintf("    ano(s) sem esforço utilizável: %s\n",
              paste(h$ano[!h$usar], collapse = ", ")))

P_HIST <- molda_padrao(data.frame(
  Amostragem = paste0("HIST_", h$ano), Ano = as.character(h$ano),
  Nome_cientifico = ESPECIE_FOCO, Especie = "CAVALA PRETA",
  Nome_engenho = ifelse(h$agregado, "AGREGADO (CERCO + LINHA)", ARTE_ALVO),
  Quantidade = ifelse(is.na(h$desem_t), NA_character_,
                      as.character(round(h$desem_t * 1000, 3))),
  Tipo_pesca = "INDUSTRIAL", Num_dias = as.character(h$esforco_d),
  stringsAsFactors = FALSE),
  "Historico 1989-2014 (agregado)", "anual", h$agregado, h$usar)

## --- (d) colagem e higiene -------------------------------------------
esforco_completo <- rbind(P_HIST, P_2018, P_1517, P_1925)
esforco_completo <- esforco_completo[order(as.numeric(esforco_completo$Ano)), ]
rownames(esforco_completo) <- NULL

## Aspas, ponto-e-vírgula e quebras de linha DENTRO de um campo estragam
## qualquer releitura. `useBytes = TRUE` é obrigatório: sem ele, uma
## classe de caracteres aplicada a uma string UTF-8 num locale C corta a
## string no primeiro acentuado. Como os caracteres removidos são ASCII
## de 1 byte (e um byte ASCII nunca aparece como continuação de um
## caractere UTF-8), a troca por byte é segura.
sujos <- 0
for (cc in COLS_PADRAO) {
  x <- esforco_completo[[cc]]
  ruim <- !is.na(x) & grepl("[\"\r\n;]", x, useBytes = TRUE)
  if (any(ruim)) {
    sujos <- sujos + sum(ruim)
    y <- trimws(gsub("[ \t]+", " ", gsub("[\"\r\n;]", " ", x, useBytes = TRUE), useBytes = TRUE))
    y[!is.na(y) & y == ""] <- NA_character_
    Encoding(y) <- "UTF-8"
    esforco_completo[[cc]] <- y
  }
}
cat(sprintf("\n[d] TABELA UNIFICADA: %d linhas x %d colunas (%d do padrão + 4 de procedência)\n",
            nrow(esforco_completo), ncol(esforco_completo), length(COLS_PADRAO)))
cat(sprintf("    campos higienizados (aspas/;/quebra de linha): %d\n", sujos))
print(table(esforco_completo$fonte, esforco_completo$nivel))
anos_ec <- sort(unique(as.numeric(esforco_completo$Ano)))
buracos <- setdiff(min(anos_ec):max(anos_ec), anos_ec)
cat(sprintf("    cobertura: %d-%d | anos ausentes: %s\n", min(anos_ec), max(anos_ec),
            if (length(buracos)) paste(buracos, collapse = ", ") else "nenhum"))

escreve_csv_utf8(esforco_completo, "esforco_completo_1989_2025.csv")
volta <- read.csv("esforco_completo_1989_2025.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
cat(sprintf("    CSV escrito e relido: %d x %d | íntegro: %s\n", nrow(volta), ncol(volta),
            identical(dim(volta), dim(esforco_completo))))
stopifnot(identical(dim(volta), dim(esforco_completo)))
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(esforco_completo = esforco_completo),
                      "esforco_completo_1989_2025.xlsx")
  cat("    XLSX escrito: esforco_completo_1989_2025.xlsx\n")
}

## --- (e) série anual nominal da cavala --------------------------------
## Esta é a série que alimenta os cenários nominais longos (C1/C2/C3).
## Para 1989-2014 vem pronta da planilha histórica; para 2015-2025 é
## agregada das viagens, com o MESMO filtro de arte e a MESMA exclusão de
## viagem multi-arte (L13) usados no resto do script.
ec <- esforco_completo
ec$Ano_n <- as.numeric(ec$Ano)
ec$Q <- suppressWarnings(as.numeric(ec$Quantidade))
ec$D <- suppressWarnings(as.numeric(ec$Num_dias))

hist_an <- data.frame(ano = h$ano, fonte = "Historico 1989-2014 (agregado)",
                      cavala_t = h$desem_t, dias = h$esforco_d, n_viagens = NA_real_,
                      cap_total_t = NA_real_, engenhos_agregados = h$agregado,
                      usar_na_serie = h$usar, stringsAsFactors = FALSE)

vg <- ec[ec$nivel == "viagem" & !is.na(ec$Amostragem) & ec$Amostragem != "", ]
n_artes <- tapply(vg$Nome_engenho, vg$Amostragem, function(z) length(unique(z[!is.na(z)])))
multi <- names(n_artes)[!is.na(n_artes) & n_artes > 1]
cat(sprintf("\n[e] série anual: %d viagens multi-arte descartadas (decisão L13)\n", length(multi)))
vg <- vg[!(vg$Amostragem %in% multi) & vg$Nome_engenho %in% ARTE_ALVO, ]

por <- function(x, by, f) tapply(x, by, f)
ids <- unique(vg$Amostragem)
V <- data.frame(
  id   = ids,
  ano  = as.numeric(por(vg$Ano_n, vg$Amostragem, function(z) z[1])[ids]),
  dias = as.numeric(por(vg$D,     vg$Amostragem, function(z) z[1])[ids]),
  usar = as.logical(por(vg$usar_na_serie, vg$Amostragem, function(z) z[1])[ids]),
  tot  = as.numeric(por(vg$Q, vg$Amostragem, function(z) sum(z, na.rm = TRUE))[ids]),
  stringsAsFactors = FALSE)
mv <- vg[toupper(vg$Nome_cientifico) == ESPECIE_FOCO, ]
mm <- por(mv$Q, mv$Amostragem, function(z) sum(z, na.rm = TRUE))
V$mac <- as.numeric(mm[V$id]); V$mac[is.na(V$mac)] <- 0
V <- V[!is.na(V$dias) & V$dias > 0 & !is.na(V$ano), ]

viag_an <- do.call(rbind, lapply(split(V, V$ano), function(s) data.frame(
  ano = s$ano[1],
  fonte = if (s$ano[1] == 2018) "IMar 2018 (formato diferente)" else "IMar viagem",
  cavala_t = sum(s$mac) / 1000, dias = sum(s$dias), n_viagens = nrow(s),
  cap_total_t = sum(s$tot) / 1000, engenhos_agregados = FALSE,
  usar_na_serie = all(s$usar), stringsAsFactors = FALSE)))

serie_anual <- rbind(hist_an, viag_an)
serie_anual <- serie_anual[order(serie_anual$ano), ]
serie_anual$cpue_nominal <- serie_anual$cavala_t / serie_anual$dias
serie_anual$periodo <- ifelse(serie_anual$ano <= ANO_CORTE_ALVO,
                              sprintf("pre-alvo (<=%d)", ANO_CORTE_ALVO),
                              sprintf("pos-alvo (>=%d)", ANO_CORTE_ALVO + 1))
rownames(serie_anual) <- NULL
cat("\n--- SÉRIE ANUAL NOMINAL DA CAVALA (captura / dias de mar) ---\n")
print(transform(serie_anual, cavala_t = round(cavala_t, 1),
                cap_total_t = round(cap_total_t, 1),
                cpue_nominal = round(cpue_nominal, 4)), row.names = FALSE)
escreve_csv_utf8(serie_anual, "serie_anual_cavala_1989_2025.csv")
cat("\nCSV escrito: serie_anual_cavala_1989_2025.csv\n")

## =====================================================================
## 0.3) ENTRADA DA PARTE 01
## ---------------------------------------------------------------------
## Daqui em diante o script trabalha com o nível de VIAGEM. O histórico
## agregado (nivel = "anual") fica de fora — ele não tem composição de
## captura, então não pode passar pela inferência de tática nem pela
## padronização. Ele volta no fim, nos cenários nominais longos.
## =====================================================================
bruto <- esforco_completo[esforco_completo$nivel == "viagem" &
                            esforco_completo$usar_na_serie, COLS_PADRAO, drop = FALSE]
rownames(bruto) <- NULL
for (j in seq_along(bruto)) bruto[[j]][is.na(bruto[[j]])] <- ""
cat(sprintf("\n>> Entram na parte 01 (nível viagem): %d linhas, anos %s\n",
            nrow(bruto), paste(range(as.numeric(bruto$Ano)), collapse = "-")))
if (!INCLUIR_2018)
  cat("   (2018 fora por decisão L17 — ver o bloco de decisões)\n")
rm(b1517, b1925, b18, d18, volta)

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

## --- L15: ILHA DO BANCO — a covariável espacial do modelo ------------
## Cada BANCO recebe a ilha de desembarque mais frequente entre as
## viagens que pescaram nele. A etiqueta é calculada sobre a série
## inteira e vale para todas as viagens daquele banco: é propriedade do
## PESQUEIRO, não da viagem. Por isso ela não carrega a variação de
## cobertura de amostragem que proíbe usar a ilha de desembarque crua
## como covariável (decisão L7).
banco_valido <- !is.na(viagens$banco) & viagens$banco != ""
if (any(!banco_valido))
  cat(sprintf("\n[L15] AVISO: %d viagens sem banco registrado -> ilha 'INDEFINIDA'\n",
              sum(!banco_valido)))

## moda da ilha de desembarque dentro de cada banco
mapa_ilha <- tapply(viagens$ilha_desemb[banco_valido],
                    viagens$banco[banco_valido],
                    function(x) { tb <- table(x); names(tb)[which.max(tb)] })
viagens$ilha_banco <- unname(mapa_ilha[viagens$banco])
viagens$ilha_banco[is.na(viagens$ilha_banco)] <- "INDEFINIDA"

## Diagnóstico do mapeamento. "Pureza" de um banco = fração das viagens
## dele que de fato desembarcaram na ilha que ele recebeu. Pureza 1
## significa que o banco só aparece ligado a uma ilha — atribuição sem
## ambiguidade. Pureza baixa significa topônimo repetido entre ilhas.
pureza <- tapply(seq_len(nrow(viagens))[banco_valido],
                 viagens$banco[banco_valido],
                 function(ix) mean(viagens$ilha_desemb[ix] == viagens$ilha_banco[ix]))
n_por_banco <- table(viagens$banco[banco_valido])
pureza      <- pureza[names(n_por_banco)]

cat(sprintf("\n[L15] Bancos distintos: %d -> %d ilhas de pesca\n",
            length(n_por_banco), length(unique(viagens$ilha_banco))))
cat(sprintf("      bancos com pureza 100%% (uma ilha só): %d de %d\n",
            sum(pureza == 1), length(pureza)))
cat(sprintf("      pureza média ponderada por viagem: %.1f%%\n",
            100 * sum(pureza * as.numeric(n_por_banco)) / sum(as.numeric(n_por_banco))))
cat(sprintf("      -> %.1f%% das viagens estão na ilha modal do seu banco\n",
            100 * mean(viagens$ilha_desemb == viagens$ilha_banco)))

## Lista de conferência: topônimos que se repetem entre ilhas. Esta
## tabela é para ser lida por quem conhece a pescaria — a atribuição
## modal força uma ilha só e nesses bancos ela erra em parte das
## viagens. É a ressalva de L15 que vai para o texto.
amb <- which(pureza < PUREZA_ALERTA & as.numeric(n_por_banco) >= 30)
if (length(amb) > 0) {
  cat(sprintf("\n      bancos AMBÍGUOS (pureza < %.0f%% e >= 30 viagens) — conferir:\n",
              100 * PUREZA_ALERTA))
  for (k in amb[order(-as.numeric(n_por_banco)[amb])]) {
    b  <- names(n_por_banco)[k]
    ds <- sort(table(viagens$ilha_desemb[viagens$banco == b]), decreasing = TRUE)
    cat(sprintf("        %-26s n=%4d  pureza=%4.1f%%  ->  %s\n", b,
                as.integer(n_por_banco[k]), 100 * pureza[k],
                paste(sprintf("%s:%d", names(ds), as.integer(ds)), collapse = "  ")))
  }
  cat("      ^ topônimos genéricos que existem em mais de uma ilha de Cabo Verde.\n")
  cat("        Uma tabela de coordenadas dos bancos resolveria isto de vez.\n")
}

cat("\n      distribuição das viagens por ilha do banco:\n")
print(table(viagens$ilha_banco))

## --- L8 (resíduo): banco agrupado, só para o efeito ALEATÓRIO --------
## `banco_gr` não é mais a covariável espacial do modelo. Ele sobrevive
## porque a estrutura E4b testa `(1 | fbanco)` como refinamento DENTRO
## da ilha — e ali o limiar pode ser brando, porque o efeito aleatório
## encolhe sozinho os níveis com pouca informação.
tb_banco <- table(viagens$banco)
raros <- names(tb_banco[tb_banco < MIN_VIAG_BANCO])
viagens$banco_gr <- ifelse(viagens$banco %in% raros, "OUTROS", viagens$banco)
cat(sprintf("\n[L8] `banco_gr` (só para o efeito aleatório E4b): %d bancos, %d com < %d\n",
            length(tb_banco), length(raros), MIN_VIAG_BANCO))
cat(sprintf("     viagens agrupados em 'OUTROS' (%.1f%% das viagens) -> %d níveis\n",
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
## 2015-2025, esse fator é o ANO. O código continua detectando sozinho o
## caso de um ano só (que cairia para mês) para não quebrar se alguém
## rodar um recorte. O ano ausente (2018, decisão L14) simplesmente não
## vira nível do fator — o índice sai com um buraco, que é o correto:
## inventar o valor de 2018 seria pior do que declarar que faltou.
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
## `ilha_desemb` (decisão L7) — ambos confundidos com o ano. `filha` NÃO
## é a ilha de desembarque da viagem: é a ilha ATRIBUÍDA AO BANCO em L15,
## que é geografia do pesqueiro e não escolha de porto.
viagens$fano   <- factor(viagens$ano)
viagens$fmes   <- factor(viagens$mes, levels = 1:12)
viagens$ftri   <- factor(viagens$trimestre, levels = 1:4,
                         labels = c("T1", "T2", "T3", "T4"))
viagens$filha  <- factor(viagens$ilha_banco)        # decisão L15 — covariável espacial
viagens$fbanco <- factor(viagens$banco_gr)          # só para o (1|fbanco) opcional (E4b)
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

## (iv) estratos tempo x ESPAÇO muito ralos geram coeficientes instáveis
## (e, no limite, níveis que só existem em um ano — que o modelo
## confundiria com efeito de ano). O estrato agora é ano x ILHA do banco
## (decisão L15) e não mais ano x banco: como a ilha tem 5 níveis em vez
## de 244, o mesmo critério corta MUITO menos viagem. É exatamente o
## ganho que motivou a troca — o filtro deixa de ser uma peneira e volta
## a ser o que devia ser, uma guarda contra célula vazia.
te <- table(paste(viagens$ano, viagens$ilha_banco))
a <- nrow(viagens)
viagens <- viagens[paste(viagens$ano, viagens$ilha_banco) %in%
                     names(te[te >= MIN_POR_ESTRATO]), ]
reg(sprintf("estratos ano x ilha do banco com >= %d viagens", MIN_POR_ESTRATO), a)

for (f in c("fano", "fmes", "ftri", "filha", "fbanco", "fbarco"))
  viagens[[f]] <- droplevels(viagens[[f]])

cat("\n================== FILTROS ================\n"); print(filtro_log, row.names = FALSE)
cat(sprintf("Retidas %d de %d viagens (%.1f%%)\n", nrow(viagens), n0,
            100 * nrow(viagens) / n0))
cat(sprintf("Após filtros: %d anos, %d ilhas de pesca, %d bancos, %d embarcações\n",
            nlevels(viagens$fano), nlevels(viagens$filha),
            nlevels(viagens$fbanco), nlevels(viagens$fbarco)))
cat("Viagens por ano após filtros:\n"); print(table(viagens$ano))

## Cruzamento ano x ilha do banco: é aqui que se vê se o efeito de ano e
## o efeito de espaço são separáveis. Uma coluna (ano) concentrada numa
## linha só (ilha) significa que, NAQUELE ano, o modelo não consegue
## distinguir "foi um ano ruim" de "só se pescou naquela ilha".
cat("\nViagens por ano x ilha do banco (base da separabilidade ano/espaço):\n")
print(table(viagens$ilha_banco, viagens$ano))
ilhas_por_ano <- colSums(table(viagens$ilha_banco, viagens$ano) > 0)
if (any(ilhas_por_ano <= 1))
  cat(sprintf("AVISO: ano(s) com uma ilha só: %s — o efeito de ano desses anos\n       absorve o efeito de espaço. Reportar como limitação.\n",
              paste(names(ilhas_por_ano)[ilhas_por_ano <= 1], collapse = ", ")))
n_por_ano <- table(viagens$ano)
if (any(n_por_ano < 150))
  cat(sprintf("AVISO: ano(s) com < 150 viagens: %s — índice com IC largo.\n",
              paste(names(n_por_ano)[n_por_ano < 150], collapse = ", ")))

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

## A — ONDE se pega cavala (componente de presença, por ILHA do banco).
##     Agora cabem todos os níveis no eixo: são 5, não 244 (decisão L15).
pr_b <- tapply(viagens$pos_mac, viagens$filha, mean)
n_b  <- table(viagens$filha)
ord  <- names(sort(pr_b, decreasing = TRUE))
barra_prop(pr_b[ord], n_b[ord], "A. Presenca de cavala (ilha do banco)",
           "", cex_nome = 0.70, n_vertical = FALSE)

## B — QUANTO se pega, dado que pegou (componente de magnitude).
pos <- viagens[viagens$pos_mac == 1, ]
## Só as ilhas com ao menos 10 viagens POSITIVAS: abaixo disso a caixa é
## desenhada sobre 3-4 pontos e não descreve distribuição nenhuma.
n_pos_b <- sort(table(pos$filha), decreasing = TRUE)
bancos_ok <- names(n_pos_b[n_pos_b >= 10])
pos_b <- pos[pos$filha %in% bancos_ok, ]
if (nrow(pos_b) > 0) {
  pos_b$filha <- droplevels(pos_b$filha)
  bx <- boxplot(cpue_dia ~ filha, data = pos_b, outline = FALSE, plot = FALSE)
  boxplot(cpue_dia ~ filha, data = pos_b, outline = FALSE, col = "#74C476",
          xaxt = "n", xlab = "", ylab = "CPUE (t/dia) entre as positivas",
          lwd = 1, main = "B. Magnitude, so nas viagens com cavala")
  ## Com 5 ilhas os rótulos cabem na horizontal; a versão rotada a 45
  ## graus era necessária quando aqui havia 15 nomes de banco.
  axis(1, at = seq_along(bx$names), labels = bx$names, cex.axis = 0.72,
       tick = TRUE, mgp = c(3, 0.6, 0))
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
png("exp6_taticas.png", width = 30, height = 10.5, res = 300,
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
legend("topleft", nome_cl, col = cor_cl, pch = 16, bty = "n", cex = 0.8, pt.cex = 1.6)
bp <- barplot(t(cent), col = cor_sp, border = NA, names.arg = nome_cl,
              las = 2, cex.names = 0.65, ylab = "Proporcao media da captura",
              main = "C. Composicao de cada tatica")
legend(max(bp) + 0.8, 1, rev(sp_nomes), fill = rev(cor_sp), border = NA,
       bty = "n", cex = 0.8, xpd = NA,pt.cex = 1.2)
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
## O estrato espacial aqui acompanha o do modelo: ilha do banco (L15).
esforco_dirigido <- aggregate(cbind(dias, horas) ~ tempo + ilha_banco + alvo,
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
#   E1  tempo + ilha do banco                       (espaço; decisão L15)
#   E2  E1 + trimestre + tripulação(classes)
#   E3  E2 + alvo                                   (tática DISCRETA)
#   E4  E3 + (1 | barco)                            [GLMM]
#   E4b E4 + (1 | banco)                            espaço FINO dentro da ilha
#   E5  E2 + PC1..PCn + (1 | barco)                 (tática CONTÍNUA)
#
# SOBRE O E4b — POR QUE ELE EXISTE:
# a ilha (5 níveis) é o gradiente espacial GROSSO. Dentro de cada ilha há bancos com
# produtividades diferentes, e essa variação fina não some só porque não cabe como fator fixo.
# O E4b a coloca como efeito ALEATÓRIO: cada banco ganha um desvio em torno da média da sua
# ilha, estimado com encolhimento (bancos com poucas viagens são puxados para a média da ilha
# em vez de receberem um coeficiente instável). É a estrutura hierarquicamente correta para
# "muitos níveis, alguns com pouca informação" — e, ao contrário do fator fixo, NÃO entra na
# grade de referência do emmeans, então não reaparece o problema de grade que derrubou a
# extração do índice na versão anterior. O AIC decide se essa camada extra se paga.
#
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
cat(sprintf("Viagens: %d | níveis de %s: %d | ilhas: %d | bancos: %d | barcos: %d | táticas: %d\n",
            nrow(viagens), fator_tempo, nlevels(viagens[[fator_tempo]]),
            nlevels(viagens$filha), nlevels(viagens$fbanco),
            nlevels(viagens$fbarco), nlevels(viagens$alvo)))
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
cat("  ilha/porto de desembarque DA VIAGEM — cobertura muda ao longo da\n")
cat("                   série => confundido com ano (decisão L7).\n")
cat("                   NÃO confundir com `filha`, que É usada: aquela é a\n")
cat("                   ilha ATRIBUÍDA AO BANCO (L15), propriedade fixa do\n")
cat("                   pesqueiro, e não o porto escolhido pela viagem.\n")
cat("  banco individual — 244 níveis; entra como efeito ALEATÓRIO se o AIC\n")
cat("                   mandar (seção 1.5), nunca como fator fixo (L8/L15)\n")

## =====================================================================
## 1) CONSTRUTOR DE FÓRMULAS
## ---------------------------------------------------------------------
## Montar a fórmula como string para evitar perder termos nos testes
## =====================================================================
## `aleatorio` aceita três formas, para que a mesma função sirva tanto
## aos modelos sem efeito aleatório quanto ao E4b, que tem dois:
##   TRUE                    -> (1 | fbarco)            [o padrão]
##   FALSE                   -> nenhum efeito aleatório
##   vetor de caracteres     -> exatamente esses termos, ex.:
##                              c("(1 | fbarco)", "(1 | fbanco)")
monta_formula <- function(resposta, termos, aleatorio = TRUE,
                          offset_var = "ldias") {
  ale <- if (isTRUE(aleatorio))  "(1 | fbarco)"
  else if (isFALSE(aleatorio)) character(0)
  else as.character(aleatorio)
  rhs <- paste(c(fator_tempo, termos, ale,
                 if (!is.null(offset_var)) sprintf("offset(%s)", offset_var)),
               collapse = " + ")
  stats::as.formula(paste(resposta, "~", rhs))
}

## Guarda contra modelo mais complexo do que os dados sustentam: um termo
## categórico só entra se houver pelo menos MIN_POR_NIVEL observações por
## nível. Importa nos subconjuntos (cenário S3), onde os fatores perdem
## povoamento e sobram poucas centenas de viagens.
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
## As estruturas são uma escada: cada degrau acrescenta um tipo de explicação
## alternativa à abundância. A leitura da tabela de AIC é literalmente "quanto
## desta série ainda parece abundância depois de descontar isto".
##   E0  só o tempo          -> o índice "cru" do modelo
##   E1  + ilha do banco     -> ONDE se pescou (P6, desbalanceamento espacial)
##   E2  + trimestre + trip. -> QUANDO se pescou e com que poder de pesca
##   E3  + alvo              -> O QUE se estava tentando pescar (P2)
##   E4  + (1|barco)         -> QUEM pescou (P5, composição de frota)
##   E4b + (1|banco)         -> onde DENTRO da ilha (espaço fino, com encolhimento)
##   E5  tática contínua     -> a alternativa da H2
## Todas ajustadas às MESMAS linhas e à MESMA resposta => AIC comparável.
##
## POR QUE O EFEITO DE BARCO É ALEATÓRIO E O DE ILHA É FIXO:
## não é preferência de estilo, são perguntas diferentes. Da ilha queremos o
## CONTRASTE entre níveis, e são 5 níveis bem povoados que precisam aparecer na
## grade do emmeans para a média marginal ser sobre estratos. Do barco não
## queremos contraste nenhum — queremos só descontar que barcos diferentes pescam
## diferente. São ~200 barcos, muitos com poucas viagens; como fator fixo cada um
## gastaria um parâmetro e os raros teriam coeficientes puro ruído. Como efeito
## aleatório eles são descritos por UM parâmetro (a variância entre barcos), e os
## barcos com pouca informação são encolhidos para a média em vez de receberem
## estimativa própria. De quebra, efeito aleatório não entra na grade de
## referência do emmeans — que é o que impede a explosão de grade que travou a
## extração do índice na versão anterior.
## ===============================================================================
termos_E <- list(
  E0  = character(0),
  E1  = c("filha"),
  E2  = c("filha", "ftri", "npesc_cat"),
  E3  = c("filha", "ftri", "npesc_cat", "alvo"),
  E4  = c("filha", "ftri", "npesc_cat", "alvo"),
  E4b = c("filha", "ftri", "npesc_cat", "alvo"),
  E5  = c("filha", "ftri", "npesc_cat", PCs)
)
## `list` e não `c`, porque E4b precisa de um VETOR de termos aleatórios.
aleat_E <- list(E0 = FALSE, E1 = FALSE, E2 = FALSE, E3 = FALSE,
                E4  = TRUE,                                  # (1|fbarco)
                E4b = c("(1 | fbarco)", "(1 | fbanco)"),     # + espaço fino
                E5  = TRUE)

cat("\n===== 1) ESTRUTURAS (Tweedie, offset = log dias) =====\n")
fits <- list()
for (nm in names(termos_E)) {
  t0 <- Sys.time()
  f  <- monta_formula("captura", termos_E[[nm]], aleatorio = aleat_E[[nm]])
  fits[[nm]] <- try(glmmTMB(f, family = tweedie(link = "log"), data = viagens),
                    silent = FALSE)
  ok <- !inherits(fits[[nm]], "try-error")
  cat(sprintf("  %-3s %-70s %s (%.1f min)\n", nm,
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

#*** COMO LER OS AVISOS DE CONVERGÊNCIA:
# "singular convergence (7)" é um código do otimizador (nlminb) dizendo que ele
# parou de conseguir melhorar. NÃO é, por si só, sinal de modelo inválido: se
# pdHess = TRUE e o maior gradiente está na casa de 1e-4, o ótimo encontrado é
# legítimo e o aviso é falso alarme. O que condena o ajuste é pdHess = FALSE
# (Hessiana não-positiva-definida), porque aí não há mínimo bem definido e os
# erros-padrão — logo os CV do índice — não valem nada.
# Na versão anterior o E3 dava esse aviso: a superfície de verossimilhança tinha
# direção quase plana, provavelmente por colinearidade parcial entre `alvo` (que
# vem da composição de espécies) e o espaço/ano (já que a composição varia
# sistematicamente por lugar e por ano); acrescentar (1 | fbarco) no E4 absorvia
# parte dessa variação compartilhada e "resolvia" a direção achatada.
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
## 1.5) A ESTRUTURA ALEATÓRIA: o banco dentro da ilha se paga?
## ---------------------------------------------------------------------
## Esta é a comparação que a decisão L15 deixou em aberto. A ilha entrou
## como fator fixo porque é o gradiente espacial que tem níveis
## suficientes para ser estimado em todos os anos. A pergunta que sobra
## é se a variação ENTRE BANCOS DA MESMA ILHA ainda carrega sinal depois
## disso.
##
## E4  = ... + (1 | fbarco)                  -> só o barco é aleatório
## E4b = ... + (1 | fbarco) + (1 | fbanco)   -> banco também
##
## Os dois têm os MESMOS efeitos fixos e as MESMAS linhas, e o glmmTMB
## ajusta por máxima verossimilhança (não REML), então o AIC é
## diretamente comparável entre eles — inclusive diferindo em efeito
## aleatório, que é o caso em que o REML invalidaria a comparação.
##
## Por que o teste é por AIC e não por LRT: testar variância zero põe o
## parâmetro na BORDA do espaço permitido (uma variância não pode ser
## negativa), e nessa situação o LRT não segue a qui-quadrado usual — o
## p-valor sai conservador. O AIC não depende dessa distribuição nula e
## responde a pergunta que interessa aqui, que é de compromisso entre
## ajuste e parcimônia, não de significância.
##
## A decisão é carregada para TODO o resto do script (backward, as duas
## distribuições, os cenários), para que a estrutura aleatória seja uma
## escolha só, feita uma vez, e não varie de modelo para modelo.
## Margem de AIC abaixo da qual dois modelos são empate técnico. Definida
## aqui porque é usada a partir desta seção (Burnham & Anderson: dAIC < 2
## não distingue modelos).
MARGEM_AIC <- 2

cat("\n===== 1.5) ESTRUTURA ALEATÓRIA: (1|fbarco) x (1|fbarco)+(1|fbanco) =====\n")
usa_aleat <- TRUE     # padrão: só o barco
if (all(c("E4", "E4b") %in% names(fits))) {
  aic_E4  <- AIC(fits[["E4"]]);  aic_E4b <- AIC(fits[["E4b"]])
  cat(sprintf("  E4  (1|fbarco)             : AIC = %.1f\n", aic_E4))
  cat(sprintf("  E4b (1|fbarco) + (1|fbanco): AIC = %.1f\n", aic_E4b))
  cat(sprintf("  dAIC (E4b - E4) = %+.1f\n", aic_E4b - aic_E4))
  if (isTRUE((aic_E4 - aic_E4b) > MARGEM_AIC)) {
    usa_aleat <- c("(1 | fbarco)", "(1 | fbanco)")
    cat("  -> o banco dentro da ilha SE PAGA: a variação fina de pesqueiro\n")
    cat("     ainda carrega sinal depois de controlar a ilha. Adotado E4b.\n")
  } else {
    cat("  -> o banco dentro da ilha NÃO se paga: a ilha já absorve o que\n")
    cat("     havia de espaço. Fica só (1|fbarco) — mais simples e mais\n")
    cat("     estável. Isto é resultado, e vai para o texto.\n")
  }
} else {
  cat("  E4 ou E4b não convergiu — mantido o padrão (1|fbarco).\n")
}

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
## A comparação usa a estrutura aleatória JÁ ESCOLHIDA na seção 1.5, para
## que offset e efeito aleatório não sejam decididos um contra o outro.
cat("\n===== 2) QUAL MEDIDA DE ESFORÇO USAR =====\n")
m_dias  <- try(glmmTMB(monta_formula("captura", termos_E$E4, usa_aleat, "ldias"),
                       family = tweedie(link = "log"), data = viagens), silent = FALSE)
m_horas <- try(glmmTMB(monta_formula("captura", termos_E$E4, usa_aleat, "lhoras"),
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
##
## O QUE O AIC ESTÁ FAZENDO AQUI, EM UMA FRASE:
##   AIC = -2 x log-verossimilhança + 2 x (nº de parâmetros)
## O primeiro termo premia ajuste; o segundo cobra por cada parâmetro
## gasto. Retirar um bloco de termos sempre PIORA a verossimilhança (o
## modelo menor não pode ajustar melhor), mas devolve parâmetros. Se o
## que se devolve vale mais que o que se perde, o AIC cai e o termo sai.
## `MARGEM_AIC = 2` é a convenção de Burnham & Anderson: diferenças
## menores que isso não distinguem modelos, então exigimos folga de 2
## para mexer na estrutura.
##
## POR QUE NÃO EXISTE p-VALOR NESTA SEÇÃO:
## com ~8 mil viagens, qualquer termo com efeito minúsculo sai
## "significativo" num LRT. Num contexto normal isso seria só um
## incômodo; aqui é perigoso, porque cada covariável que entra no modelo
## COME sinal do efeito de ano — e o efeito de ano é o produto final. É o
## dilema de Hinton & Maunder (2003): a covariável que "explica" a queda
## da CPUE pode estar explicando justamente a queda de abundância que
## queremos medir. Por isso a régua é de compromisso (AIC), não de
## significância, e por isso a carga da prova é para RETIRAR.
## =====================================================================
cat("\n===== 3) REFINAMENTO BACKWARD POR AIC =====\n")
## `usa_aleat` vem da seção 1.5 e NÃO é redefinido aqui: a estrutura
## aleatória já foi decidida por AIC e vale para todo o resto do script.
cat(sprintf("  estrutura aleatória em uso: %s\n",
            if (isTRUE(usa_aleat)) "(1 | fbarco)" else paste(usa_aleat, collapse = " + ")))
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
##
## POR QUE O AIC PODE COMPARAR DUAS DISTRIBUIÇÕES DIFERENTES:
## o AIC é comparável sempre que os modelos descrevem A MESMA VARIÁVEL
## RESPOSTA, nas MESMAS OBSERVAÇÕES, e a verossimilhança reportada é a
## completa (com todas as constantes de normalização). As duas condições
## valem aqui: a resposta é `captura` em toneladas nos dois casos, sem
## transformação, e o glmmTMB devolve a log-verossimilhança exata, não
## uma quasi-verossimilhança. O que NÃO seria comparável é modelar
## log(captura) num e captura no outro — aí as respostas são variáveis
## diferentes e o AIC perde o sentido. Mudar só a função de LIGAÇÃO,
## mantendo a distribuição e a resposta, é comparável pelo mesmo motivo.
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

#Comparação entre os modelos tweedie e Hurdle-Gamma
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
## COMO O DHARMa CONSTRÓI O RESÍDUO (vale entender para ler a tabela):
## para cada observação, ele simula centenas de valores a partir do
## modelo ajustado e pergunta em que QUANTIL da distribuição simulada cai
## o valor observado. Se o modelo estiver correto, esses quantis são
## Uniforme(0,1) por construção — qualquer que seja a distribuição
## (Tweedie, Gamma, Poisson). É isso que torna a leitura a mesma para
## todos os modelos e resolve o problema do resíduo de Pearson em GLMM.
##
## O que cada teste responde:
##   KS         os resíduos são mesmo Uniforme(0,1)? É o teste global de
##              adequação da distribuição assumida.
##   dispersão  a variância dos dados bate com a que o modelo prevê?
##              `disp_ratio` = observada/esperada. ~1 é bom; >1 é
##              sobredispersão (o modelo subestima a variabilidade e os
##              IC do índice saem estreitos demais); <1 é subdispersão.
##   outliers   há mais valores fora do envelope simulado do que o
##              modelo consegue gerar?
##   quantis    a variância é homogênea ao longo do predito? Ajusta
##              regressões quantílicas (0,25/0,50/0,75) dos resíduos
##              contra o valor predito: se o modelo está bem
##              especificado, as três curvas saem horizontais. Curva
##              inclinada = o erro depende do nível do predito, que é a
##              heterocedasticidade aqui.
##   zeros      o modelo gera a quantidade certa de zeros? Com ~86% de
##              zeros na resposta, este é o teste que mais importa para
##              escolher entre Tweedie e hurdle.
##
## COMO LER O p-VALOR AQUI (importante): com ~8 mil observações, um
## desvio irrelevante já produz p < 0,001. O p-valor responde "o desvio é
## detectável?", e a resposta é quase sempre sim. A pergunta que
## interessa é "o desvio é GRANDE?", e essa só o tamanho do efeito
## responde — `disp_ratio` perto de 1 e os gráficos. Por isso os PNG são
## salvos e são eles, não a tabela, que decidem se a distribuição serve.
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

## GUARDA DA GRADE DE REFERÊNCIA — o erro que travou a versão anterior.
## Com `weights = "equal"` o emmeans monta o produto cartesiano de TODOS
## os níveis de TODOS os fatores fixos e recusa passar de `rg.limit`
## (10.000 por padrão). Fatores com muitos níveis estouram isso: com o
## banco como fator fixo em 46 níveis a grade pedida era 25.760 e a
## extração do índice morria ali. Com a ilha (5 níveis) a grade cai para
## alguns milhares. Efeito ALEATÓRIO não entra na grade — é a segunda
## razão para `(1 | fbarco)` e `(1 | fbanco)` serem aleatórios.
## Este bloco calcula o tamanho ANTES de chamar o emmeans e avisa, em vez
## de deixar o script morrer 200 linhas adiante com uma mensagem críptica.
tamanho_grade <- function(m) {
  mf <- model.frame(m)
  fx <- names(mf)[vapply(mf, is.factor, logical(1))]
  ## termos aleatórios não entram na grade de referência
  aleat <- unlist(lapply(m$modelInfo$reTrms, function(z) names(z$cnms)))
  fx <- setdiff(fx, aleat)
  if (length(fx) == 0) return(1L)
  prod(vapply(mf[fx], nlevels, numeric(1)))
}
cat("\n--- tamanho da grade de referência do emmeans (limite: 10.000) ---\n")
for (nm in names(dist_fits)) {
  g <- tryCatch(tamanho_grade(dist_fits[[nm]]), error = function(e) NA_real_)
  cat(sprintf("  %-18s %6.0f %s\n", nm, g,
              if (!is.na(g) && g > 10000)
                "<<< vai estourar: reduza níveis ou passe o fator a aleatório"
              else "OK"))
}

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

#----------------------------------------------------------------------
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
#-----------------------------------------------------------------------
## S0 — modelo SEM a covariável de tática, com o resto igual. A diferença
##      S0 - S2 isola o efeito do direcionamento: é o "influence plot" de
##      Bentley et al. reduzido ao termo que interessa. Sem este cenário
##      não dá para afirmar que a correção veio da tática e não de outro
##      termo qualquer.
E_sem_alvo <- setdiff(termos_sel, c("alvo", PCs))
m_S0 <- ajusta_final(E_sem_alvo)
S0 <- if (!inherits(m_S0, "try-error"))
  normaliza(indice_de(m_S0, "S0 sem tática")) else NULL
#-----------------------------------------------------------------------
## S2 — cenário principal: modelo selecionado (estrutura + distribuição
##      vencedoras), com tática discreta.
S2 <- normaliza(indice_de(m_final, "S2 corrigida (tática discreta)"))
#-----------------------------------------------------------------------
## S2t — a MESMA estrutura na outra distribuição. Não é redundância: se o
##      índice muda de forma ao trocar Tweedie por hurdle, isso é
##      incerteza ESTRUTURAL e tem de ir para o texto; se não muda, é um
##      argumento forte de robustez.
outro_nome <- setdiff(names(dist_fits), melhor_dist)
S2t <- if (length(outro_nome) > 0)
  normaliza(indice_de(dist_fits[[outro_nome[1]]],
                      sprintf("S2t %s", outro_nome[1]))) else NULL
#------------------------------------------------------------------------
## S2b — tática CONTÍNUA (hipótese H2): os escores da PCA no lugar do
##      fator de cluster, mesma base de termos e MESMA distribuição.
##      ATENÇÃO: a seção 7.6 mostra que esta série é CIRCULAR — os PCs
##      nascem de uma composição que inclui a própria cavala, então a
##      covariável carrega parte da resposta e o efeito de ano sai
##      achatado. Ela é mantida por ser o RESULTADO do teste (e para a
##      figura mostrar o contraste), mas NÃO é candidata a entrada do
##      JABBA. O contraponto honesto à discreta é o `S2b_limpo`,
##      construído na seção 7.6 a partir da composição sem a cavala.
S2b <- if (!inherits(m_cont_final, "try-error"))
  normaliza(indice_de(m_cont_final, "S2b PCs COM cavala (circular)")) else NULL
#-------------------------------------------------------------------------
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
    for (f in c("filha", "fbanco", "fbarco", "ftri", "npesc_cat"))
      if (f %in% names(v3)) v3[[f]] <- droplevels(v3[[f]])
    cat(sprintf("  S3: %d viagens da tática '%s' (%.0f%% de cavala no centróide)\n",
                nrow(v3), alvo_cavala, 100 * frac_cavala))
    ## O subconjunto é pequeno e pode não sustentar a mesma estrutura
    ## aleatória do modelo principal: um `(1 | fbanco)` com dezenas de
    ## níveis e poucas viagens por nível não é estimável. Se for o caso,
    ## S3 cai para `(1 | fbarco)` e isso é registrado — a alternativa
    ## (o ajuste falhar em silêncio) perderia o cenário inteiro.
    aleat_S3 <- usa_aleat
    if (!isTRUE(usa_aleat) && "(1 | fbanco)" %in% usa_aleat &&
        nrow(v3) < 10 * nlevels(v3$fbanco)) {
      aleat_S3 <- "(1 | fbarco)"
      cat(sprintf("      [S3] (1|fbanco) retirado: %d viagens para %d bancos\n",
                  nrow(v3), nlevels(v3$fbanco)))
    }
    m3 <- ajusta_final(termos_viaveis(v3, E_sem_alvo), dados = v3, aleat = aleat_S3)
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
## 7.6) TESTE DE CIRCULARIDADE DA COVARIÁVEL DE TÁTICA
## ---------------------------------------------------------------------
## POR QUE ISTO É RESULTADO E NÃO RASCUNHO:
## a matriz de composição que gera as DUAS representações de tática
## (`alvo`, do k-means, e PC1..PCn, da PCA) sai de `cols_cap`, que INCLUI
## `cap_macarellus` — a própria espécie-resposta. Isso abre a porta para
## uma circularidade: a covariável passa a conter parte da resposta, e o
## modelo "explica" a captura de cavala com uma variável que é, em parte,
## a captura de cavala. O sintoma é sempre o mesmo par: ajuste
## espetacular e efeito de ANO achatado. É o dilema de Hinton & Maunder
## (2003) descrito na seção 4, aqui na sua forma mais aguda — a
## covariável não apenas compete com o sinal de abundância, ela É o sinal
## de abundância entrando pela porta dos fundos.
##
## O TESTE: reconstruir a composição SEM a cavala e refazer a PCA. O
## direcionamento pode ser descrito inteiramente pelo que sobra na rede
## ("este barco-mês foi 80% Auxis / 20% Katsuwonus"), e como as
## proporções são RENORMALIZADAS, o mix relativo das outras espécies fica
## preservado intacto. Perde-se apenas o número que era a resposta
## disfarçada. Se a vantagem da representação contínua sobreviver à
## limpeza, ela é real; se desaparecer, ela era vazamento.
##
## POR QUE ISSO NÃO DESTRÓI O CENÁRIO DE ESFORÇO DIRIGIDO (S3):
## o teste NÃO altera o pipeline. O k-means, o `alvo`, o `alvo_cavala` e
## o S3 continuam exatamente como estavam; o que se mede aqui é se
## deveriam mudar. O item (4) responde isso com um número.
##
## OS TRÊS MODELOS SÃO AJUSTADOS NAS MESMAS LINHAS. Barcos-mês que só
## pegaram cavala ficam sem composição depois da remoção e saem do teste.
## AIC comparado entre conjuntos de linhas diferentes não significa nada.
## =====================================================================
cat("\n===== 7.6) CIRCULARIDADE DA TÁTICA =====\n")

## (1) quanto a agregação por barco-mês dilui a viagem individual
vpbm <- as.numeric(table(viagens$barco_mes))
cat(sprintf("(1) viagens por barco-mês: mediana %.0f | só uma viagem: %.1f%%\n",
            stats::median(vpbm), 100 * mean(vpbm == 1)))
cat("    quanto mais barcos-mês de uma viagem só, mais a composição do\n")
cat("    barco-mês É a composição da própria viagem — vazamento direto.\n")

## (2) os PCs enxergam a resposta?
esf_v    <- viagens[[if (OFFSET == "lhoras") "horas" else "dias"]]
cor_cpue <- sapply(PCs, function(p) cor(viagens[[p]], viagens$captura / esf_v,
                                        use = "complete.obs"))
cor_pres <- sapply(PCs, function(p) cor(viagens[[p]], as.numeric(viagens$captura > 0),
                                        use = "complete.obs"))
cat("\n(2) correlação dos PCs com a PRÓPRIA resposta:\n")
print(round(rbind(`com a CPUE` = cor_cpue, `com a presença` = cor_pres), 3))
cat("    correlação alta aqui já é o vazamento visível sem ajustar modelo.\n")

## (3) PCA reconstruída SEM a cavala, e os três modelos nas mesmas linhas
cols_sm <- setdiff(cols_cap, "cap_macarellus")
ag2 <- aggregate(viagens[, cols_sm], by = list(barco_mes = viagens$barco_mes), FUN = sum)
tt2 <- rowSums(ag2[, cols_sm]); ok2 <- tt2 > 0
cp2 <- sqrt(as.matrix(ag2[ok2, cols_sm] / tt2[ok2]))
rownames(cp2) <- ag2$barco_mes[ok2]
pca2 <- prcomp(cp2, center = TRUE, scale. = FALSE)
ve2  <- 100 * pca2$sdev^2 / sum(pca2$sdev^2)
np2  <- min(4, max(2, which(cumsum(ve2) >= 70)[1]))
if (is.na(np2)) np2 <- min(4, ncol(pca2$x))
PCs2 <- paste0("PCsm", seq_len(np2))
cat(sprintf("\n(3) PCA sem a cavala: %d eixos retidos (%.0f%% da variação)\n",
            np2, cumsum(ve2)[np2]))

## Os escores novos vão para uma CÓPIA de `viagens`: o objeto que o resto
## do script usa não pode ser poluído por variáveis de teste.
v_circ <- viagens
for (j in seq_len(np2))
  v_circ[[PCs2[j]]] <- as.numeric(pca2$x[match(viagens$barco_mes, rownames(pca2$x)), j])
v_circ <- v_circ[stats::complete.cases(v_circ[, PCs2]), ]
for (f in c("filha", "fbanco", "fbarco", "ftri", "npesc_cat", "alvo", fator_tempo))
  if (f %in% names(v_circ)) v_circ[[f]] <- droplevels(v_circ[[f]])
cat(sprintf("    linhas do teste: %d de %d\n", nrow(v_circ), nrow(viagens)))

m_disc_c <- ajusta_final(termos_sel, dados = v_circ)
m_pc_com <- ajusta_final(c(setdiff(termos_sel, "alvo"), PCs), dados = v_circ)
## `monta_zi()` lê `PCs` do ambiente global. Trocamos temporariamente para
## que os eixos NOVOS entrem também no componente de zeros — senão o
## modelo limpo competiria em desvantagem justamente no lado do hurdle
## onde o vazamento mais pesa, e o teste ficaria viciado a favor dele.
PCs_bkp <- PCs; PCs <- PCs2
m_pc_sem <- ajusta_final(c(setdiff(termos_sel, "alvo"), PCs2), dados = v_circ)
PCs <- PCs_bkp

aic_ok <- function(m) if (inherits(m, "try-error")) NA_real_ else AIC(m)
a_disc <- aic_ok(m_disc_c); a_com <- aic_ok(m_pc_com); a_sem <- aic_ok(m_pc_sem)
d_com  <- a_com - a_disc;   d_sem <- a_sem - a_disc
cat(sprintf("\n    discreto (`alvo`)   AIC = %9.1f\n", a_disc))
cat(sprintf("    PCs COM a cavala    AIC = %9.1f   (dAIC %+.1f)\n", a_com, d_com))
cat(sprintf("    PCs SEM a cavala    AIC = %9.1f   (dAIC %+.1f)\n", a_sem, d_sem))

## (4) o `alvo` DISCRETO depende da espécie-resposta?
## Ele nasce da mesma matriz, então tem a mesma exposição em princípio.
## Mas o k-means joga tudo dentro de um grupo na mesma categoria e
## descarta a variação interna, então vaza muito menos. Aqui isso deixa
## de ser argumento e vira medida: refazemos o agrupamento sobre a PCA
## limpa e comparamos a atribuição.
conc_km <- NA_real_
if (exists("km") && exists("esc") && exists("k_otimo")) {
  set.seed(1)
  esc2  <- pca2$x[, 1:min(3, ncol(pca2$x)), drop = FALSE]
  km2   <- stats::kmeans(esc2, centers = k_otimo, nstart = 50, iter.max = 100)
  bm_c  <- intersect(rownames(esc), rownames(esc2))
  tb_km <- table(com_cavala = km$cluster[match(bm_c, rownames(esc))],
                 sem_cavala = km2$cluster[match(bm_c, rownames(esc2))])
  conc_km <- 100 * sum(apply(tb_km, 2, max)) / sum(tb_km)
  cat(sprintf("\n(4) atribuição de tática com x sem a cavala: %.1f%% de concordância\n",
              conc_km))
  print(tb_km)
  cat("    concordância alta = o `alvo` discreto NÃO depende da espécie-\n")
  cat("    resposta, e o pipeline (k-means, `alvo_cavala`, S3) fica de pé.\n")
}

## (5) a série ainda achata?
## Este é o teste que importa para a avaliação: AIC mede ajuste, mas o
## produto final é a FORMA da série. Se o efeito de ano volta a cair ao
## limpar a covariável, o achatamento era artefato.
S2b_limpo <- if (!inherits(m_pc_sem, "try-error"))
  normaliza(indice_de(m_pc_sem, "S2b PCs SEM cavala")) else NULL
r_sem <- NA_real_; r_com <- NA_real_
if (!is.null(S2b_limpo) && !is.null(S2b)) {
  cmp <- data.frame(ano      = S2b_limpo$tempo,
                    discreta = S2$indice[match(S2b_limpo$tempo, S2$tempo)],
                    pc_sem   = S2b_limpo$indice,
                    pc_com   = S2b$indice[match(S2b_limpo$tempo, S2b$tempo)])
  cat("\n(5) as três séries lado a lado (média = 1):\n")
  print(transform(cmp, discreta = round(discreta, 3), pc_sem = round(pc_sem, 3),
                  pc_com = round(pc_com, 3)), row.names = FALSE)
  r_sem <- cor(cmp$discreta, cmp$pc_sem, use = "complete.obs")
  r_com <- cor(cmp$discreta, cmp$pc_com, use = "complete.obs")
  amp_c <- function(x) max(x, na.rm = TRUE) / min(x, na.rm = TRUE)
  cat(sprintf("    discreta            : amplitude %.2f\n", amp_c(cmp$discreta)))
  cat(sprintf("    PCs SEM a cavala    : amplitude %.2f | r com a discreta = %.3f\n",
              amp_c(cmp$pc_sem), r_sem))
  cat(sprintf("    PCs COM a cavala    : amplitude %.2f | r com a discreta = %.3f\n",
              amp_c(cmp$pc_com), r_com))
}

## VEREDITO — tudo derivado dos números acima, nada fixo no texto.
circularidade <- isTRUE(d_com < -MARGEM_AIC) && isTRUE(d_sem > -MARGEM_AIC)
cat("\nVEREDITO: ")
if (circularidade) {
  cat("CIRCULARIDADE CONFIRMADA.\n")
  cat(sprintf("  A representação contínua parecia superior por %.0f pontos de AIC.\n", -d_com))
  cat(sprintf("  Removida a espécie-resposta da composição, a diferença vira %+.1f\n", d_sem))
  cat(sprintf("  — uma virada de %.0f pontos causada por UMA coluna da matriz.\n", d_sem - d_com))
  cat("  O S2b (PCs COM cavala) NÃO pode ser cenário de entrada do JABBA:\n")
  cat("  o efeito de ano dele está achatado porque a covariável absorveu o\n")
  cat("  sinal de abundância que o índice deveria medir.\n")
  if (isTRUE(d_sem > MARGEM_AIC))
    cat("  Entre as representações SEM vazamento, a DISCRETA ajusta melhor.\n")
  else if (isTRUE(abs(d_sem) <= MARGEM_AIC))
    cat("  Entre as representações SEM vazamento, há empate técnico.\n")
} else if (isTRUE(d_sem < -MARGEM_AIC)) {
  cat("SEM circularidade detectável — a vantagem da representação contínua\n")
  cat("  sobrevive à remoção da espécie-resposta da composição. Ela é real.\n")
} else {
  cat("a representação DISCRETA ajusta melhor entre as opções limpas.\n")
}

tab_circ <- data.frame(
  representacao    = c("discreta (alvo)", "continua PCs COM cavala",
                       "continua PCs SEM cavala"),
  AIC              = round(c(a_disc, a_com, a_sem), 1),
  dAIC_vs_discreta = round(c(0, d_com, d_sem), 1),
  r_com_discreta   = round(c(1, r_com, r_sem), 3),
  n_linhas         = nrow(v_circ),
  concordancia_kmeans_pct = round(conc_km, 1),
  row.names = NULL, stringsAsFactors = FALSE)

## A série limpa entra no painel de cenários: é ela, e não o S2b, que
## serve de contraponto honesto à discreta.
if (!is.null(S2b_limpo))
  indices <- rbind(indices, S2b_limpo[, names(indices), drop = FALSE])

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
  ## Leitura cruzada com a seção 7.6: um `alvo` que quase não desloca o
  ## índice é também um `alvo` que quase não está vazando a resposta —
  ## covariável contaminada desloca MUITO, como o S2b desloca.
}
amp <- function(d) max(d$indice, na.rm = TRUE) / min(d$indice, na.rm = TRUE)
amp_nom <- amp(S1); amp_cor <- amp(S2)
cat(sprintf("H1 — amplitude (máx/mín): nominal %.2f | corrigida %.2f\n",
            amp_nom, amp_cor))
cat(sprintf("     correlação nominal x corrigida: %.3f\n",
            cor(S2$indice, casa(S2, S1), use = "complete.obs")))
## O texto abaixo é CONDICIONAL: a versão anterior afirmava sempre que a
## amplitude corrigida era menor, o que não é um resultado garantido —
## e nesta série não é o que acontece.
if (amp_cor < amp_nom) {
  cat("     amplitude MENOR na corrigida: parte da variação nominal era\n")
  cat("     comportamento de frota e não abundância (o efeito usual).\n")
} else {
  cat("     amplitude MAIOR na corrigida: a padronização NÃO achatou a\n")
  cat("     série. Padronizar não garante amplitude menor — ao remover\n")
  cat("     desbalanceamento a correção pode SEPARAR anos que o cálculo\n")
  cat("     nominal confundia. Não é anomalia, mas o texto NÃO pode usar\n")
  cat("     a frase usual de que a correção removeu variação de frota.\n")
}
if (!is.null(S2b)) {
  d_h2 <- AIC(m_cont_final) - AIC(m_final)   # negativo = contínua ganha
  cat(sprintf("H2 — discreta x contínua: r = %.3f | dAIC bruto (contínua - discreta) = %+.1f\n",
              cor(S2$indice, casa(S2, S2b), use = "complete.obs"), d_h2))
  if (exists("circularidade") && isTRUE(circularidade)) {
    cat("     ESTE dAIC BRUTO É ARTEFATO — não reportar sozinho.\n")
    cat(sprintf("     A seção 7.6 mostra que ele vira %+.1f quando a espécie-\n", d_sem))
    cat("     resposta sai da matriz de composição. A vantagem aparente da\n")
    cat("     representação contínua era a cavala se explicando a si mesma.\n")
    cat(sprintf("     CONCLUSÃO VÁLIDA: %s\n",
                if (isTRUE(d_sem > MARGEM_AIC))
                  "entre representações limpas, a DISCRETA ajusta melhor."
                else if (isTRUE(d_sem < -MARGEM_AIC))
                  "entre representações limpas, a CONTÍNUA ainda ajusta melhor."
                else "entre representações limpas, há empate técnico."))
    cat("     A leitura de Winker et al. (2013) a favor da representação\n")
    cat("     contínua NÃO se sustenta nesta série — e o motivo (a espécie-\n")
    cat("     resposta dentro da composição) é reportável por si só.\n")
  } else if (d_h2 < -MARGEM_AIC) {
    cat("     a representação CONTÍNUA (PCs) ajusta melhor — como em Winker et al. (2013)\n")
  } else if (d_h2 > MARGEM_AIC) {
    cat("     a representação DISCRETA (cluster) ajusta melhor\n")
  } else {
    cat("     empate técnico: as duas descrevem a tática igualmente bem\n")
  }
}
if (!is.null(S3))
  cat(sprintf("H3 — esforço dirigido x corrigida: r = %.3f | amplitude %.2f\n",
              cor(S2$indice, casa(S2, S3), use = "complete.obs"), amp(S3)))
## =====================================================================
## 9) FIGURA DOS ÍNDICES
## =====================================================================
## O cenário circular entra em CINZA: ele está na figura como resultado
## do teste da seção 7.6, não como candidato. A série contínua LIMPA
## (`S2b PCs SEM cavala`) é que serve de contraponto à discreta.
cores_cen <- setNames(
  c(COR_AUX, COR_NEU, COR_MAC, COR_S2D, "grey60", COR_S2, COR_S3),
  c("S1 nominal", "S0 sem tática", "S2 corrigida (tática discreta)",
    if (length(outro_nome) > 0) sprintf("S2t %s", outro_nome[1]) else "S2t",
    "S2b PCs COM cavala (circular)", "S2b PCs SEM cavala",
    "S3 esforço dirigido"))

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
## 10) OS CINCO CENÁRIOS FINAIS (entrada do JABBA)
## ---------------------------------------------------------------------
## Cada cenário é uma SÉRIE ALTERNATIVA de abundância relativa, não uma
## versão "mais certa" da mesma coisa. Eles respondem a perguntas
## diferentes e cobrem períodos diferentes:
##
##   C1  nominal 1989-2025   — a série inteira, como a FAO a construiu.
##                             Junta duas fontes de natureza diferente
##                             (agregada até 2014, viagem a viagem
##                             depois), o que é a sua principal fraqueza.
##   C2  nominal pré-alvo    — só até `ANO_CORTE_ALVO` (L16). É o período
##                             em que a cavala ainda era ALVO da frota,
##                             então o esforço do cerco é um denominador
##                             razoável para ela.
##   C3  nominal pós-alvo    — de `ANO_CORTE_ALVO`+1 em diante. Mesmo
##                             cálculo, período em que a cavala virou
##                             captura acompanhante — é justamente aqui
##                             que o denominador deixa de ser específico
##                             e a CPUE nominal fica suspeita (P1/P2).
##   C4  padronizada SEM tática — modelo com espaço, sazonalidade,
##                             tripulação e embarcação, mas SEM a
##                             covariável de direcionamento. É o cenário
##                             S0 da seção 6.
##   C5  padronizada COM tática — o mesmo modelo MAIS a tática inferida
##                             da composição da captura. É o cenário S2.
##   C6  nominal esforço dirigido — mesmo cálculo simples de C3 (captura/
##                             dias, sem GLM), mas só com as viagens da
##                             tática mais associada à cavala (`alvo_cavala`,
##                             seção 5). Ataca o mesmo problema de C5 por
##                             outro mecanismo: C5 AJUSTA o efeito de
##                             tática com toda viagem dentro do modelo; C6
##                             DESCARTA as viagens não-dirigidas. Só existe
##                             2015-2025 (mesma razão de C4/C5) e carrega a
##                             ressalva da seção 5.5: nenhuma tática é de
##                             fato dominada pela cavala nesta frota.
##
## A distância C4 -> C5 é a medida do viés de direcionamento pela via da
## PADRONIZAÇÃO; a distância C3 -> C6 é a mesma medida pela via do
## SUBSETTING — duas formas independentes de estimar o mesmo viés, e a
## distância C3 -> C5 é o efeito total da padronização.
##
## POR QUE CADA SÉRIE É NORMALIZADA PELA PRÓPRIA MÉDIA: o JABBA estima um
## coeficiente de capturabilidade (q) por índice, então o que ele lê é a
## FORMA da série, não o nível. Normalizar cada cenário pela sua própria
## média deixa isso explícito e evita que um recorte (C2, C3) herde a
## escala da série inteira.
##
## SOBRE O CV DA PARTE HISTÓRICA: a planilha 1989-2014 é agregada — um
## número por ano, sem viagens por trás. Não existe variância amostral
## para extrair dali. Esses anos recebem o PISO de CV, e a coluna
## `cv_origem` registra que é piso, não estimativa. Tratar um CV inventado
## como se fosse medido seria o erro mais caro desta etapa.
## =====================================================================
PISO_CV <- 0.20

norm1 <- function(x) x / mean(x, na.rm = TRUE)

## --- C1/C2/C3: nominais, a partir da série anual da Parte 0 ----------
sa <- serie_anual[serie_anual$usar_na_serie & !is.na(serie_anual$cpue_nominal), ]

## CV empírico do período de viagem (erro-padrão relativo da CPUE entre as
## viagens do ano); nos anos agregados não há como calculá-lo.
cv_emp_ano <- tapply(viagens$captura / viagens[[if (OFFSET == "lhoras") "horas" else "dias"]],
                     viagens$ano, function(z) sd(z) / (mean(z) * sqrt(length(z))))
monta_nominal <- function(d, nome) {
  if (nrow(d) == 0) return(NULL)
  cv <- as.numeric(cv_emp_ano[as.character(d$ano)])
  data.frame(tempo = d$ano, cenario = nome,
             indice = norm1(d$cpue_nominal),
             cv = ifelse(is.na(cv), PISO_CV, cv),
             cv_origem = ifelse(is.na(cv), "piso (fonte agregada)", "empirico"),
             fonte = d$fonte, engenhos_agregados = d$engenhos_agregados,
             row.names = NULL, stringsAsFactors = FALSE)
}
C1 <- monta_nominal(sa, "C1 nominal 1989-2025")
C2 <- monta_nominal(sa[sa$ano <= ANO_CORTE_ALVO, ],
                    sprintf("C2 nominal pre-alvo (<=%d)", ANO_CORTE_ALVO))
C3 <- monta_nominal(sa[sa$ano >  ANO_CORTE_ALVO, ],
                    sprintf("C3 nominal pos-alvo (>=%d)", ANO_CORTE_ALVO + 1))

## --- C4/C5: padronizadas, vindas dos modelos da seção 6 --------------
## S0 = sem a covariável de tática; S2 = com ela. Já vêm normalizadas.
de_modelo <- function(d, nome) {
  if (is.null(d)) return(NULL)
  data.frame(tempo = d$tempo, cenario = nome, indice = d$indice, cv = d$cv,
             cv_origem = "modelo", fonte = "IMar viagem (padronizado)",
             engenhos_agregados = FALSE, row.names = NULL, stringsAsFactors = FALSE)
}
C4 <- de_modelo(S0, "C4 padronizada SEM tatica")
C5 <- de_modelo(S2, "C5 padronizada COM tatica")

## --- C6: nominal, mas com o esforço filtrado só para a tática-cavala --
## Pedido explícito (set/2026): um cenário NOMINAL — sem GLM nenhum,
## captura/dias puro, do mesmo jeito que C1-C3 — só que com o
## DENOMINADOR restrito às viagens que a inferência de tática (seção 5)
## classificou como `alvo_cavala`, em vez do esforço de TODA a frota de
## cerco. A pergunta que C6 testa, que é diferente da de C5: se o
## problema da CPUE pós-2015 é o DENOMINADOR (esforço "genérico",
## incluindo viagens que nunca iriam encontrar cavala), filtrar por
## tática deveria destravar sinal de depleção que a nominal-total (C3)
## não mostra. C5 já ataca isso por outro caminho — CORRIGINDO o efeito
## de tática como covariável num GLM, com toda viagem dentro; C6
## SUBTRAI as viagens não-dirigidas em vez de corrigi-las. São mecanismos
## diferentes (ajustar x descartar) e por isso vale rodar os dois.
##
## SÓ EXISTE 2015-2025 — mesma razão de C4/C5: a tática só é inferível de
## viagem com composição registrada, e só há viagem a viagem a partir de
## 2015. A planilha 1989-2014 é agregada; não tem o que classificar.
##
## RESSALVA DA PRÓPRIA SEÇÃO 5.5, e ela é séria: nenhuma tática desta
## frota é de fato DOMINADA pela cavala (a melhor gira em torno de
## `frac_cavala`*100%, tipicamente ~15%). `alvo_cavala` é "a tática que
## MAIS ENCONTRA cavala", não "a tática exclusiva da cavala" — então C6
## está mais para uma sensibilidade do que para um índice "dirigido" no
## sentido em que a palavra é usada em pescarias com arte realmente
## seletiva. Declarar isso é obrigatório se C6 for reportado.
C6 <- NULL
if (exists("alvo_cavala") && alvo_cavala %in% levels(viagens$alvo)) {
  vd      <- viagens[viagens$alvo == alvo_cavala, ]
  anos_vd <- sort(unique(vd$ano))
  ## mesmo critério de cobertura que o S3 já usa (seção 6): sem isso, um
  ## cenário que só existe em 2-3 anos entraria como se fosse série.
  cobre_c6 <- length(anos_vd) >= 0.7 * length(unique(viagens$ano))
  if (cobre_c6) {
    cavala_ano_vd <- tapply(vd$captura, vd$ano, sum)
    dias_ano_vd   <- tapply(vd$dias,    vd$ano, sum)
    cpue_vd       <- cavala_ano_vd / dias_ano_vd
    ## CV empírico DENTRO do subconjunto dirigido — não reaproveita
    ## cv_emp_ano (esse é calculado com TODAS as viagens do ano).
    cv_emp_vd <- tapply(vd$captura / vd$dias, vd$ano, function(z)
      if (length(z) > 1) sd(z) / (mean(z) * sqrt(length(z))) else NA_real_)
    C6 <- data.frame(tempo = as.numeric(names(cpue_vd)),
                     cenario = "C6 nominal esforco dirigido",
                     indice = norm1(as.numeric(cpue_vd)),
                     cv = pmax(ifelse(is.na(as.numeric(cv_emp_vd)), PISO_CV,
                                      as.numeric(cv_emp_vd)), PISO_CV),
                     cv_origem = "empirico (subconjunto dirigido)",
                     fonte = "IMar viagem (so tatica-cavala)",
                     engenhos_agregados = FALSE, row.names = NULL,
                     stringsAsFactors = FALSE)
    cat(sprintf("\nC6: %d viagens em %d anos na tatica '%s' (%.0f%% de cavala no centroide -- ver ressalva 5.5)\n",
                nrow(vd), length(anos_vd), alvo_cavala, 100 * frac_cavala))
  } else {
    cat("[AVISO] tatica-cavala nao cobre anos suficientes para C6; cenario nao construido.\n")
  }
} else {
  cat("[AVISO] alvo_cavala inexistente/nao encontrado; C6 nao construido.\n")
}

cenarios <- do.call(rbind, Filter(Negate(is.null), list(C1, C2, C3, C4, C5, C6)))

cat("\n===== OS CENARIOS DE INDICE =====\n")
resumo_cen <- do.call(rbind, lapply(split(cenarios, cenarios$cenario), function(s) data.frame(
  cenario = s$cenario[1], anos = sprintf("%d-%d", min(s$tempo), max(s$tempo)),
  n = nrow(s), amplitude = round(max(s$indice, na.rm = TRUE) /
                                   min(s$indice, na.rm = TRUE), 2),
  cv_medio = round(mean(pmax(s$cv, PISO_CV), na.rm = TRUE), 3),
  row.names = NULL, stringsAsFactors = FALSE)))
print(resumo_cen[order(resumo_cen$cenario), ], row.names = FALSE)
## Condicional, e não afirmação fixa: padronizar NÃO garante amplitude
## menor. A versão anterior imprimia a frase do caso usual mesmo quando
## os números diziam o contrário, o que iria direto para o texto errado.
cat("\nAmplitude = máx/mín do índice.\n")
if (!is.null(C3) && !is.null(C5)) {
  amp_C3 <- max(C3$indice, na.rm = TRUE) / min(C3$indice, na.rm = TRUE)
  amp_C5 <- max(C5$indice, na.rm = TRUE) / min(C5$indice, na.rm = TRUE)
  if (amp_C5 < amp_C3) {
    cat("Amplitude MENOR na padronizada (C5) que na nominal do mesmo período\n")
    cat("(C3): parte da variação nominal era comportamento de frota e não\n")
    cat("abundância — é o efeito usual da correção.\n")
  } else {
    cat("Amplitude MAIOR na padronizada (C5) que na nominal do mesmo período\n")
    cat("(C3): a correção não achatou a série. Ao remover desbalanceamento\n")
    cat("ela pode SEPARAR anos que o cálculo nominal confundia — reportar\n")
    cat("assim, sem a frase usual de 'a correção removeu variação de frota'.\n")
  }
}

## Correlação entre os cenários que compartilham período
if (!is.null(C3) && !is.null(C5)) {
  anos_com <- intersect(C3$tempo, C5$tempo)
  if (length(anos_com) > 2) {
    r35 <- cor(C3$indice[match(anos_com, C3$tempo)],
               C5$indice[match(anos_com, C5$tempo)], use = "complete.obs")
    cat(sprintf("\nC3 (nominal pos-alvo) x C5 (padronizada com tatica): r = %.3f em %d anos\n",
                r35, length(anos_com)))
  }
}
if (!is.null(C4) && !is.null(C5)) {
  r45 <- cor(C4$indice, C5$indice[match(C4$tempo, C5$tempo)], use = "complete.obs")
  cat(sprintf("C4 (sem tatica) x C5 (com tatica): r = %.3f", r45))
  cat(sprintf("  -> %s\n", if (r45 > 0.98)
    "a tatica quase nao desloca o indice (reportar!)" else
      "a tatica desloca o indice de forma relevante"))
}
## C3 x C6: o teste direto do que C6 foi criado para responder — será que
## restringir o denominador à tática-cavala muda a FORMA da série nominal
## do período pós-alvo, ou ela reproduz C3 e o denominador não era o
## problema?
if (!is.null(C3) && !is.null(C6)) {
  anos_com36 <- intersect(C3$tempo, C6$tempo)
  if (length(anos_com36) > 2) {
    r36 <- cor(C3$indice[match(anos_com36, C3$tempo)],
               C6$indice[match(anos_com36, C6$tempo)], use = "complete.obs")
    cat(sprintf("C3 (nominal pos-alvo, esforco total) x C6 (nominal, esforco dirigido): r = %.3f em %d anos",
                r36, length(anos_com36)))
    cat(sprintf("  -> %s\n", if (r36 > 0.9)
      "quase identicas: o denominador (esforco total x dirigido) NAO era o problema" else
        "divergem: filtrar o esforco muda a forma da serie, vale investigar por que"))
  }
}

## --- planilhas por cenário, no formato do JABBA ----------------------
## Uma coluna de ano e uma coluna por cenário. Anos que um cenário não
## cobre ficam NA — o JABBA aceita e simplesmente não usa aquele ponto.
anos_saida <- sort(unique(cenarios$tempo))
jabba_idx <- data.frame(Yr = anos_saida)
jabba_cv  <- data.frame(Yr = anos_saida)
rotulo <- c("C1 nominal 1989-2025" = "C1_nominal_total",
            setNames("C2_nominal_pre_alvo",  sprintf("C2 nominal pre-alvo (<=%d)", ANO_CORTE_ALVO)),
            setNames("C3_nominal_pos_alvo",  sprintf("C3 nominal pos-alvo (>=%d)", ANO_CORTE_ALVO + 1)),
            "C4 padronizada SEM tatica" = "C4_padronizada_sem_tatica",
            "C5 padronizada COM tatica" = "C5_padronizada_com_tatica",
            "C6 nominal esforco dirigido" = "C6_nominal_dirigida")
for (cen in names(rotulo)) {
  if (!cen %in% cenarios$cenario) next
  d <- cenarios[cenarios$cenario == cen, ]
  jabba_idx[[rotulo[[cen]]]] <- d$indice[match(anos_saida, d$tempo)]
  jabba_cv[[rotulo[[cen]]]]  <- pmax(d$cv[match(anos_saida, d$tempo)], PISO_CV)
}

cat("\n--- índices (média 1) por cenário ---\n")
## O ano é inteiro; só as colunas de índice levam casas decimais.
mostra <- jabba_idx; mostra[[1]] <- as.integer(mostra[[1]])
for (j in 2:ncol(mostra)) mostra[[j]] <- round(mostra[[j]], 3)
print(mostra, row.names = FALSE)
cat("NA = o cenário não cobre aquele ano (o JABBA apenas não usa o ponto).\n")

escreve_csv_utf8(jabba_idx,  "jabba_indices_macarellus.csv")
escreve_csv_utf8(jabba_cv,   "jabba_cv_macarellus.csv")
escreve_csv_utf8(cenarios,   "cenarios_cpue_macarellus.csv")
escreve_csv_utf8(indices,    "indices_todos_cenarios.csv")
escreve_csv_utf8(tab_circ,   "teste_circularidade_tatica.csv")
if (tem_writexl)
  writexl::write_xlsx(list(indices_jabba = jabba_idx, cv_jabba = jabba_cv,
                           cenarios = cenarios, resumo = resumo_cen,
                           serie_anual = serie_anual, todos_indices = indices,
                           estruturas = tab_est, distribuicoes = tab_dist,
                           diagnostico = diag_tab, circularidade = tab_circ),
                      path = "cpue_macarellus_cenarios.xlsx")
cat("\nArquivos salvos:\n")
cat("  jabba_indices_macarellus.csv   (índices, média 1, uma coluna por cenário)\n")
cat("  jabba_cv_macarellus.csv        (CV correspondente, com piso de", PISO_CV, ")\n")
cat("  cenarios_cpue_macarellus.csv   (formato longo, com fonte e origem do CV)\n")
cat("  indices_todos_cenarios.csv     (os cenários intermediários S0-S3)\n")
cat("  teste_circularidade_tatica.csv (a tabela da seção 7.6)\n")
if (tem_writexl) cat("  cpue_macarellus_cenarios.xlsx  (tudo acima em abas)\n")

## --- figura dos cinco cenários ---------------------------------------
png("cenarios_finais.png", width = 26, height = 13, res = 300, antialias = AA, units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3, 1), bty = "l", cex.main = 0.95, cex = 0.85)
cor_c <- c("C1 nominal 1989-2025" = COR_NEU)
cor_c[sprintf("C2 nominal pre-alvo (<=%d)", ANO_CORTE_ALVO)]    <- COR_AUX
cor_c[sprintf("C3 nominal pos-alvo (>=%d)", ANO_CORTE_ALVO + 1)] <- COR_S2D
cor_c["C4 padronizada SEM tatica"] <- COR_S3
cor_c["C5 padronizada COM tatica"] <- COR_MAC
cor_c["C6 nominal esforco dirigido"] <- COR_S2

plot(NA, xlim = range(cenarios$tempo), ylim = c(0, max(cenarios$indice, na.rm = TRUE) * 1.1),
     xlab = "Ano", ylab = "Indice relativo (media = 1)",
     main = "A. Serie completa 1989-2025")
abline(v = ANO_CORTE_ALVO + 0.5, lty = 3, col = "grey40")
## rótulo embaixo, para não brigar com a legenda no canto superior
text(ANO_CORTE_ALVO + 0.5, max(cenarios$indice, na.rm = TRUE) * 0.04,
     "mudanca de alvo", pos = 4, cex = 0.62, col = "grey30")
for (cen in unique(cenarios$cenario)) {
  d <- cenarios[cenarios$cenario == cen, ]
  lines(d$tempo, d$indice, lwd = 2.2, col = cor_c[cen])
  points(d$tempo, d$indice, pch = 19, cex = 0.6, col = cor_c[cen])
}
legend("topright", names(cor_c), col = cor_c, lwd = 2.2, bty = "n", cex = 0.58)

## painel B: só o período padronizado, onde os cenários são comparáveis
sub <- cenarios[cenarios$tempo > ANO_CORTE_ALVO, ]
plot(NA, xlim = range(sub$tempo), ylim = c(0, max(sub$indice, na.rm = TRUE) * 1.1),
     xlab = "Ano", ylab = "Indice relativo (media = 1)",
     main = sprintf("B. Periodo pos-alvo (>=%d)", ANO_CORTE_ALVO + 1))
for (cen in unique(sub$cenario)) {
  d <- sub[sub$cenario == cen, ]
  lines(d$tempo, d$indice, lwd = 2.4, col = cor_c[cen])
  points(d$tempo, d$indice, pch = 19, cex = 0.8, col = cor_c[cen])
}
legend("topright", unique(sub$cenario), col = cor_c[unique(sub$cenario)],
       lwd = 2.4, bty = "n", cex = 0.62)
par(op); dev.off()
cat("PNG salvo: cenarios_finais.png\n")

## =====================================================================
## 12) SÉRIE DE DESEMBARQUES (CAPTURA) — entrada de CATCH do JABBA
## ---------------------------------------------------------------------
## Até aqui C1-C5 são tudo CPUE (abundância relativa). O JABBA também
## precisa da série de REMOÇÕES (captura em toneladas), que é um insumo
## INDEPENDENTE do índice — e o IMar não forneceu uma série de
## desembarques pronta (a que alimentou o CMSY/DBSRA antes veio da FAO).
##
## A captura de cavala preta já está dentro das MESMAS planilhas de
## esforço usadas para a CPUE — `serie_anual`, montada na seção 0.2 —
## então a série de desembarques sai do mesmo lugar, sem reler nada.
##
## POR QUE A CAPTURA NÃO HERDA OS MESMOS "BURACOS" DA CPUE:
## o JABBA é um modelo de estado-espaço — ele precisa de uma remoção
## (mesmo que pequena) em CADA ano do modelo para atualizar a biomassa;
## o ÍNDICE não, ele só entra nos anos em que existe e o JABBA
## simplesmente pula os anos sem índice. Por isso 2013 (sem esforço) e
## 2018 (formato e levantamento diferentes — decisão L17) ENTRAM aqui
## mesmo estando marcados `usar_na_serie = FALSE` em `serie_anual`: não
## dá para pular um ano de captura sem quebrar o balanço de massa do
## modelo. A marcação de confiabilidade viaja junto na coluna
## `confiavel`, para o texto — não para remover a linha.
## =====================================================================
cat("\n===== 12) SÉRIE DE DESEMBARQUES (CAPTURA) =====\n")

desembarques <- data.frame(
  Yr = serie_anual$ano,
  Catch_macarellus_t = round(serie_anual$cavala_t, 3),
  ## NA em 1989-2014: a planilha histórica é agregada e não abre a
  ## captura total da frota por espécie — só a de cavala.
  Catch_total_frota_t = round(serie_anual$cap_total_t, 3),
  fonte = serie_anual$fonte,
  confiavel = serie_anual$usar_na_serie,
  engenhos_agregados = serie_anual$engenhos_agregados,
  periodo = serie_anual$periodo,
  row.names = NULL, stringsAsFactors = FALSE)

cat(sprintf("Cobertura: %d-%d (%d anos)\n",
            min(desembarques$Yr), max(desembarques$Yr), nrow(desembarques)))
falta <- is.na(desembarques$Catch_macarellus_t)
if (any(falta)) {
  cat(sprintf("Anos SEM captura registrada (%d): %s\n", sum(falta),
              paste(desembarques$Yr[falta], collapse = ", ")))
  cat("  ^ ficam como estão (NA) — o pedido foi para não travar nisso; o\n")
  cat("    JABBA aceita, mas o manual do pacote recomenda decidir entre\n")
  cat("    interpolar ou usar um valor mínimo antes de rodar, porque um\n")
  cat("    NA no meio da série de captura quebra o balanço daquele ano.\n")
} else {
  cat("Nenhum ano sem captura registrada — série completa para o JABBA.\n")
}
nao_conf <- !desembarques$confiavel
if (any(nao_conf))
  cat(sprintf("Anos com captura presente mas marcados NÃO confiáveis para o índice: %s\n",
              paste(desembarques$Yr[nao_conf], collapse = ", ")))
cat("  ^ continuam na série de captura — a ressalva de 2013 (sem esforço)\n")
cat("    e 2018 (levantamento e formato diferentes) vale para este número\n")
cat("    também e deve constar no texto, não justifica remover a linha.\n")

print(desembarques, row.names = FALSE)

escreve_csv_utf8(desembarques, "desembarques_macarellus_1989_2025.csv")
cat("\nCSV escrito: desembarques_macarellus_1989_2025.csv\n")
if (tem_writexl) {
  writexl::write_xlsx(list(desembarques = desembarques),
                      "desembarques_macarellus_1989_2025.xlsx")
  cat("XLSX escrito: desembarques_macarellus_1989_2025.xlsx\n")
}

## --- gráfico -----------------------------------------------------------
png("desembarques_macarellus.png", width = 24, height = 12, res = 300,
    antialias = AA, units = "cm")
op <- par(mar = c(4.4, 4.6, 3, 1), bty = "l", cex.main = 0.95, cex = 0.85)
cor_fonte <- c("Historico 1989-2014 (agregado)" = COR_NEU,
               "IMar viagem" = COR_MAC,
               "IMar 2018 (formato diferente)" = COR_AUX)
plot(desembarques$Yr, desembarques$Catch_macarellus_t, type = "n",
     xlab = "Ano", ylab = "Desembarque de cavala preta (t)",
     ylim = c(0, max(desembarques$Catch_macarellus_t, na.rm = TRUE) * 1.15),
     main = "Serie de desembarques - D. macarellus, frota de cerco (Cabo Verde)", lwd=2)
abline(v = ANO_CORTE_ALVO + 0.5, lty = 3, col = "grey40")
text(ANO_CORTE_ALVO + 3, max(desembarques$Catch_macarellus_t, na.rm = TRUE) * 1,
     "mudanca de alvo", pos = 2, cex = 0.7, col = "grey30")
lines(desembarques$Yr, desembarques$Catch_macarellus_t, lwd =2, col = "grey75")
for (f in names(cor_fonte)) {
  d <- desembarques[desembarques$fonte == f, ]
  if (nrow(d) == 0) next
  points(d$Yr[d$confiavel],  d$Catch_macarellus_t[d$confiavel],
         pch = 19, cex = 1, col = cor_fonte[f])
  points(d$Yr[!d$confiavel], d$Catch_macarellus_t[!d$confiavel],
         pch = 4,  cex = 1.3, lwd = 2, col = cor_fonte[f])
}
legend("topleft", names(cor_fonte), col = cor_fonte, pch = 19, bty = "n", cex = 0.9, pt.cex = 1)
legend("topright", c("confiavel p/ indice", "nao confiavel p/ indice (entra na captura mesmo assim)"),
       pch = c(19, 4), pt.lwd = c(1, 2), col = "grey30", bty = "n", cex = 0.9, pt.cex = 1)
par(op); dev.off()
cat("PNG salvo: desembarques_macarellus.png\n")

## =====================================================================
## 13) COMO REPORTAR
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
cat("   tipo de embarcação, ilha de desembarque da viagem) — todas por\n")
cat("   problema de dado, não por não serem significativas.\n")
cat("5. Declarar como o espaço foi tratado (decisão L15): a covariável é a\n")
cat("   ILHA ATRIBUÍDA AO BANCO pela moda da ilha de desembarque, não o\n")
cat("   banco individual nem a ilha de desembarque da viagem. Reportar a\n")
cat("   pureza do mapeamento e a lista de bancos ambíguos — são topônimos\n")
cat("   que se repetem entre ilhas e é a principal fonte de erro de\n")
cat("   classificação espacial desta análise.\n")
cat("6. Declarar a origem dos dados (decisão L14): QUATRO fontes coladas\n")
cat("   por ano, com granularidades diferentes — agregada por ano até\n")
cat("   2014, viagem a viagem de 2015 em diante. A coluna `fonte` da\n")
cat("   planilha unificada identifica cada trecho.\n")
cat(sprintf("7. Declarar que 2018 EXISTE mas ficou FORA da série (L17): CPUE\n"))
cat("   nominal de 1,03 t/dia contra 0,28 em 2017 e 0,38 em 2019, com\n")
cat("   esforço normal. Como vem de levantamento e formato diferentes,\n")
cat("   não dá para separar abundância de artefato de amostragem. A\n")
cat("   tradução dele está pronta na planilha; basta INCLUIR_2018 <- TRUE\n")
cat("   se o IMar confirmar os números. Sem tripulação e sem código de\n")
cat("   embarcação, 2018 nunca poderá entrar na série PADRONIZADA.\n")
cat("8. Declarar que 2013 tem captura mas NÃO tem esforço (L18): entra na\n")
cat("   tabela de capturas do JABBA e fica fora do índice.\n")
cat("9. NÃO soldar C2 e C3 numa única série contínua. São o mesmo cálculo\n")
cat("   em dois regimes diferentes de pescaria: antes de 2015 a cavala era\n")
cat("   alvo e o esforço do cerco é um denominador razoável para ela;\n")
cat("   depois, ela é captura acompanhante e o mesmo denominador passa a\n")
cat("   medir outra coisa. C1 existe para mostrar a série inteira, mas a\n")
cat("   descontinuidade de fonte e de regime tem de estar no texto.\n")
cat("10. Os CV dos anos 1989-2014 são PISO, não estimativa: a fonte é\n")
cat("   agregada e não tem variância amostral. A coluna `cv_origem` do\n")
cat("   arquivo cenarios_cpue_macarellus.csv marca cada caso.\n")
cat("11. REPORTAR O TESTE DE CIRCULARIDADE (seção 7.6) como resultado de\n")
cat("   método, não como nota de rodapé. A representação CONTÍNUA da\n")
cat("   tática (escores da PCA) parecia superior à discreta por uma\n")
cat("   margem enorme de AIC, e essa vantagem se inverte quando a\n")
cat("   espécie-resposta é removida da matriz de composição que gera os\n")
cat("   eixos. Era a cavala se explicando a si mesma: a covariável\n")
cat("   continha parte da resposta, o ajuste disparava e o efeito de ANO\n")
cat("   saía achatado — que é justamente o índice que se quer medir.\n")
cat("   Consequências a declarar: (a) o cenário contínuo NÃO entra no\n")
cat("   JABBA; (b) o `alvo` discreto foi verificado e a sua atribuição\n")
cat("   praticamente não muda ao limpar a composição, por isso o\n")
cat("   pipeline (k-means, tática da cavala, S3) foi mantido; (c) o\n")
cat("   `alvo` ainda nasce de uma matriz que inclui a cavala, o que é\n")
cat("   limitação a declarar, sustentada pelos números do teste e não\n")
cat("   por argumento. É a demonstração empírica do dilema de Hinton &\n")
cat("   Maunder (2003) na própria série.\n")
cat("12b. C6 (esforco dirigido) so cobre 2015-2025 e carrega a ressalva da\n")
cat("   secao 5.5: nenhuma tatica desta frota e de fato DOMINADA pela\n")
cat("   cavala (a melhor gira em torno da fracao reportada em `frac_cavala`).\n")
cat("   `alvo_cavala` e 'a tatica que MAIS encontra cavala', nao 'a tatica\n")
cat("   exclusiva da cavala' — declarar isso sempre que C6 for reportado.\n")
cat("   Comparar C3 x C6 (o script imprime a correlacao acima): se forem\n")
cat("   quase identicas, o denominador de esforco NAO era o problema da\n")
cat("   CPUE pos-2015 — o mesmo viés de composicao que afeta C3 provavelmente\n")
cat("   afeta C6 tambem, porque a tatica nasce da mesma matriz de composicao\n")
cat("   que inclui a cavala (secao 7.6).\n")
cat("12. A série de CAPTURA (seção 12, `desembarques_macarellus_1989_2025`)\n")
cat("   substitui a série da FAO usada antes no CMSY/DBSRA: agora vem da\n")
cat("   mesma fonte primária (IMar) que a CPUE, extraída das planilhas de\n")
cat("   esforço. Declarar que ela NÃO segue o mesmo filtro de qualidade do\n")
cat("   índice — 2013 e 2018 entram na captura mesmo estando marcados como\n")
cat("   não confiáveis para a CPUE — porque o JABBA precisa de uma remoção\n")
cat("   em todo ano do modelo, e o índice não. As mesmas ressalvas de 2013\n")
cat("   (sem esforço) e 2018 (levantamento e formato diferentes, decisão\n")
cat("   L17) se aplicam ao número de captura desses anos.\n")


#----------------------------------------------------------------------------------------------------------#
#   AVALIACAO DE ESTOQUE — Decapterus macarellus (cavala preta), Cabo Verde                                #
#   JABBA — Just Another Bayesian Biomass Assessment (Winker, Carvalho & Kapur 2018, Fish. Res. 204:275)   #
#                                                                                                          #
#   Este script é a PARTE 04 do pipeline. Ele NAO depende do `Catch_and_index_models.R` estar na           #
#   memória: lê os CSV que aquele script já escreveu no diretório. Pode ser colado no fim dele ou          #
#   rodado sozinho — as poucas funções compartilhadas são redefinidas aqui só se não existirem.            #
#                                                                                                          #
#   ENTRADAS (escritas pelas seções 10 e 12 da parte 03):                                                  #
#     jabba_indices_macarellus.csv        índices de CPUE, média 1, uma coluna por cenário C1..C6          #
#     jabba_cv_macarellus.csv             CV correspondente a cada ponto de índice                         #
#     desembarques_macarellus_1989_2025.csv  captura anual (t) de cavala preta                             #
#                                                                                                          #
#   SAIDAS: uma pasta por rodada + planilhas consolidadas + figuras comparativas.                          #
#----------------------------------------------------------------------------------------------------------#

## =========================================================================
## 0) ANTES DE RODAR — LEIA ISTO
## -------------------------------------------------------------------------
## O JABBA é um modelo BAYESIANO ajustado por MCMC, e quem faz o MCMC é o
## JAGS — um programa EXTERNO ao R, que precisa ser instalado à parte. Esse
## é o erro nº 1 de quem roda JABBA pela primeira vez: `library(JABBA)`
## funciona, `fit_jabba()` morre com uma mensagem críptica sobre não achar
## o JAGS.
##
##   Windows: baixar e instalar o JAGS 4.3.x em
##            https://sourceforge.net/projects/mcmc-jags/files/
##            DEPOIS disso, no R: install.packages(c("rjags","R2jags"))
##            (a ordem importa — o rjags compila contra o JAGS instalado)
##
## O QUE O MODELO FAZ, EM UMA FRASE: ele assume que a biomassa cresce por
## uma curva de produção excedente e encolhe pela captura,
##
##      B[t+1] = B[t] + produção(B[t]) - Captura[t]
##
## e usa os índices de CPUE como observações ruidosas de B[t] (via
## CPUE[t] = q * B[t]). Tudo que ele estima — K, r, MSY, B/Bmsy — sai dessa
## tensão entre o que a captura removeu e o que o índice diz que sobrou.
##
## AS TRES COISAS QUE MAIS DECIDEM O RESULTADO (nesta ordem):
##   1. A SERIE DE CAPTURA. O modelo acredita nela quase cegamente. Se a
##      captura estiver subestimada, K e MSY saem subestimados junto.
##      Ver a seção 3 — tem um problema real a resolver aqui.
##   2. O PRIOR DE r. Stobberup & Erzini (2006), que avaliaram ESTE MESMO
##      estoque em Cabo Verde, reportaram que "a posterior marginal de r
##      era quase idêntica à prior" — ou seja, os dados não informam r, o
##      prior informa. Por isso a seção 5 tem análise de sensibilidade de
##      prior e ela NAO é opcional.
##   3. QUANTOS q O MODELO ESTIMA. Cada índice ganha o seu, e é assim que
##      se acomoda a mudança de alvo de 2015. Ver a seção 6.
## =========================================================================


## =========================================================================
## 1) PACOTES, DIRETORIOS E CHAVES DA RODADA
## =========================================================================

## --- 1.1 pacotes -------------------------------------------------------
if (!requireNamespace("JABBA", quietly = TRUE)) {
  stop("Pacote JABBA ausente. Instale com:\n",
       "  install.packages('devtools')\n",
       "  devtools::install_github('jabbamodel/JABBA')\n",
       "E confirme que o JAGS (programa externo) está instalado — ver seção 0.")
}
library(JABBA)
tem_writexl <- requireNamespace("writexl", quietly = TRUE)

## Checagem explícita do JAGS: falhar AGORA, com mensagem clara, é muito
## melhor do que falhar dentro do primeiro fit_jabba() daqui a 10 minutos.
if (!requireNamespace("rjags", quietly = TRUE))
  stop("Pacote rjags ausente — instale o JAGS (programa externo) e depois install.packages('rjags').")
cat(sprintf("JAGS encontrado: versão %s\n", as.character(rjags::jags.version())))

## --- 1.2 diretórios ----------------------------------------------------
## `DIR_DADOS` é onde estão os CSV da parte 03. Se você colou este bloco no
## fim do Catch_and_index_models.R, o getwd() já é o lugar certo.
DIR_DADOS  <- getwd()
ASSESSMENT <- "macarellus_CV"
## Algumas funções de plotagem "de comparação entre cenários" do JABBA
## (jbplot_summary é a principal, mas jbplot_ensemble tem o mesmo problema)
## vêm do estilo antigo de script de entrada do JABBA — aquele em que você
## definia `assessment` como uma variável solta no ambiente global antes do
## loop de cenários (foi assim no Sp_name_v1.0.R que você me mandou). Essas
## funções ainda procuram por um objeto chamado `assessment` no ambiente
## global em vez de tirar o nome de dentro de cada item da lista de ajustes
## — mesmo quando você já passou `assessment=...` para o build_jabba() de
## cada cenário. Sem essa variável solta, elas quebram com
## "object 'assessment' not found". Criamos ela aqui, redundante com
## ASSESSMENT (maiúsculo, usado no resto do script), só para blindar essas
## chamadas legadas.
assessment <- ASSESSMENT
DIR_SAIDA  <- file.path(DIR_DADOS, "JABBA_macarellus")
dir.create(DIR_SAIDA, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Saídas em: %s\n", DIR_SAIDA))

## --- 1.3 chaves da rodada ----------------------------------------------
## RODADA_RAPIDA = TRUE usa uma cadeia MCMC curta (o `quickmcmc` do JABBA).
## Serve para DESENVOLVER: você descobre em 1 minuto que um cenário está mal
## especificado, em vez de em 20. NUNCA reporte resultado de rodada rápida —
## as caudas das posteriores não estão amostradas direito.
RODADA_RAPIDA <- FALSE

## MCMC da rodada final. O que cada número significa:
##   MCMC_NI = iterações TOTAIS por cadeia
##   MCMC_NB = burn-in: as primeiras iterações, jogadas fora, porque a
##             cadeia ainda está "andando" da posição inicial até a região
##             de alta densidade da posterior
##   MCMC_NT = thinning: guarda 1 a cada NT iterações. Reduz autocorrelação
##             entre amostras guardadas e o tamanho do objeto
##   MCMC_NC = número de cadeias independentes. Precisa ser >= 2 para os
##             diagnósticos de convergência (Gelman-Rubin) fazerem sentido
## Amostras guardadas = (NI - NB)/NT * NC
MCMC_NI <- 30000; MCMC_NB <- 5000; MCMC_NT <- 5; MCMC_NC <- 2
cat(sprintf("MCMC final: %d amostras guardadas por parâmetro\n",
            (MCMC_NI - MCMC_NB) / MCMC_NT * MCMC_NC))

## Rodar os diagnósticos pesados? Retrospectiva e hindcast reajustam o
## modelo uma vez por "peel" — multiplica o tempo por ~6. Deixe FALSE
## enquanto estiver montando os cenários; ligue na rodada final.
RODAR_DIAGNOSTICO <- TRUE
RODAR_PROJECAO    <- TRUE

## --- 1.4 utilitários (só se ainda não existirem) ------------------------
if (!exists("escreve_csv_utf8")) {
  escreve_csv_utf8 <- function(d, caminho) {
    esc <- function(x) {
      x <- ifelse(is.na(x), "", as.character(x))
      pr <- grepl('[",\r\n]', x, useBytes = TRUE)
      x[pr] <- paste0('"', gsub('"', '""', x[pr], fixed = TRUE), '"')
      x
    }
    con <- file(caminho, open = "wb"); on.exit(close(con))
    writeLines(enc2utf8(c(paste(names(d), collapse = ","),
                          do.call(paste, c(lapply(d, esc), sep = ",")))),
               con, useBytes = TRUE)
  }
}
if (!exists("COR_MAC")) { COR_MAC <- "#1F4E79"; COR_AUX <- "#C0501B"; COR_NEU <- "#7F7F7F" }
if (!exists("COR_S2"))  { COR_S2  <- "#2E8B57"; COR_S3  <- "#7030A0"; COR_S2D <- "#00A0B0" }
if (!exists("AA")) AA <- if (.Platform$OS.type == "windows") "cleartype" else "default"


## =========================================================================
## 2) LEITURA DAS ENTRADAS
## =========================================================================
cat("\n===== 2) ENTRADAS =====\n")

le <- function(f) {
  cam <- file.path(DIR_DADOS, f)
  if (!file.exists(cam))
    stop(sprintf("Não achei '%s' em %s.\nRode antes o Catch_and_index_models.R.", f, DIR_DADOS))
  d <- read.csv(cam, stringsAsFactors = FALSE, encoding = "UTF-8")
  ## Arquivos salvos pelo Excel costumam vir com BOM (a marca invisível de
  ## codificação no início do arquivo). O read.csv() a incorpora ao nome da
  ## PRIMEIRA coluna, que vira "X.U.FEFF.Year" em vez de "Year" — e aí
  ## qualquer busca por nome de coluna falha, com uma mensagem que não
  ## sugere BOM nenhum. Limpar aqui resolve para todos os arquivos.
  names(d) <- sub("^X\\.U\\.FEFF\\.", "", names(d))
  names(d) <- sub("^﻿", "", names(d))
  d
}
idx_in  <- le("jabba_indices_macarellus.csv")
cv_in   <- le("jabba_cv_macarellus.csv")
cap_in  <- le("desembarques_macarellus_1989_2025.csv")

cat(sprintf("índices : %d anos x %d cenários (%s)\n", nrow(idx_in), ncol(idx_in) - 1,
            paste(setdiff(names(idx_in), "Yr"), collapse = ", ")))
cat(sprintf("captura : %d anos (%d-%d)\n", nrow(cap_in), min(cap_in$Yr), max(cap_in$Yr)))

## As duas tabelas de índice têm de ser gêmeas — mesmo ano, mesma coluna.
stopifnot(identical(names(idx_in), names(cv_in)), identical(idx_in$Yr, cv_in$Yr))

## Cobertura de cada cenário de índice. Vale olhar: é isso que decide quais
## combinações fazem sentido na seção 6.
cat("\nCobertura de cada cenário de índice:\n")
for (cn in setdiff(names(idx_in), "Yr")) {
  ok <- !is.na(idx_in[[cn]])
  cat(sprintf("  %-28s %2d anos (%d-%d)\n", cn, sum(ok),
              min(idx_in$Yr[ok]), max(idx_in$Yr[ok])))
}


## =========================================================================
## 3) A SERIE DE CAPTURA — a decisão mais consequente do script
## -------------------------------------------------------------------------
## O JABBA trata a captura como REMOÇÃO TOTAL do estoque. Ele não tem como
## saber que existe pescaria que você não contou: se metade das remoções
## fica de fora, o modelo conclui que o estoque aguenta menos do que
## aguenta, e K, MSY e B/Bmsy saem todos deslocados.
##
## >>> PROBLEMA CONHECIDO NESTA SERIE <<<
## A coluna `Catch_macarellus_t` é a captura da REDE DE CERCO INDUSTRIAL —
## foi assim que a parte 03 a montou, e com razão, porque o denominador da
## CPUE tinha de bater com o filtro de arte. Para o índice isso é correto.
## Para a captura do JABBA é uma escolha diferente, e possivelmente errada:
##   - a planilha histórica 1989-2014 tem uma coluna "Total" (cerco + linha
##     de mão) que NAO foi exportada;
##   - a parte artesanal / outras artes de 2015-2025 não está aqui;
##   - a própria FAO (a fonte que você usou no CMSY/DBSRA) tem uma série de
##     desembarques mais abrangente.
##
## TRES CAMINHOS, e o script aceita os três:
##   "cerco"   usa o que está no CSV. Coerente com o índice, mas é UMA arte
##             de UM segmento. Só é defensável se o cerco industrial for
##             de fato quase toda a remoção de cavala — o que precisa ser
##             verificado, não assumido.
##   "externa" você aponta um CSV com a série total (ex.: a da FAO usada no
##             CMSY). É o caminho recomendado se a série existir.
##   "escala"  usa a do cerco multiplicada por um fator fixo, como
##             aproximação declarada. Último recurso, e tem de ir para o
##             texto como premissa.
##
## Stobberup & Erzini (2006), avaliando esta espécie em Cabo Verde, falavam
## em captura sustentável da ordem de milhares de toneladas. Se a sua série
## de cerco somar muito menos que isso, é sinal de que ela NAO representa a
## remoção total — e a opção "cerco" vai enviesar tudo.
## =========================================================================
cat("\n===== 3) SERIE DE CAPTURA =====\n")

## DECISAO REVERTIDA (set/2026): voltamos para "cerco". A série externa
## (Luz & Vieira / FAO) só cobre até 2015, e o diagnóstico direto da CPUE
## (log-CPUE contra captura acumulada) mostrou que o problema de ajuste NAO
## é a série de captura — é que o índice em si não carrega sinal de
## depleção em nenhum dos dois períodos (tendência positiva pré-2015, ligada
## a artefatos de esforço; e pós-2015 o índice está confundido com a
## proporção de cavala na captura da frota, correlação de 0,95). Cortar a
## janela para 1989-2015 perde 10 anos e não resolve nada. "Cerco" cobre a
## série toda (1989-2025) e é o que os cenários de comparação (J1-J12)
## precisam para continuar cobrindo o período inteiro.
FONTE_CAPTURA  <- c("cerco", "externa", "escala")[1]

## Série externa. "Catch_Luz and Vieira.csv" é a compilação de desembarques
## totais (a mesma linhagem de dado que alimentou o CMSY e o DBSRA), com
## colunas "Year" e "Catch". É a opção recomendada, e a razão é específica:
## o JABBA precisa de REMOÇÃO TOTAL no campo de captura, enquanto o índice
## pode vir de UMA frota só (ele ganha um q próprio para isso). Captura de
## todas as artes + índice de uma arte é o arranjo padrão, não um remendo.
##
## O QUE ISSO CUSTA, e tem de ir para o texto:
##   - a avaliação deixa de ser independente do CMSY/DBSRA, porque passa a
##     compartilhar a entrada de captura com eles;
##   - é preciso conferir a resolução taxonômica da série (D. macarellus ou
##     um agregado de Decapterus/Carangidae) e dizer qual é.
ARQ_CAPT_EXT   <- "Catch_Luz and Vieira.csv"   # se FONTE_CAPTURA = "externa"
FATOR_ESCALA   <- 1.0                          # se FONTE_CAPTURA = "escala"

## >>> REGRA INEGOCIAVEL DESTE BLOCO <<<
## A JANELA DO MODELO E A JANELA DA SERIE DE CAPTURA. Nada de estender o
## modelo para além do último ano com captura observada.
##
## Por que isto está escrito em letra grande: a tentação é preencher os anos
## do fim que a compilação externa não cobre (repetindo o último valor,
## escalando a série do cerco, o que for) para que as figuras cheguem até
## 2025. Isso INVENTA exatamente os anos que determinam o status terminal —
## B/Bmsy e F/Fmsy do ano final saem de números que ninguém mediu — e destrói
## a comparabilidade entre os cenários, porque cada um passa a ver um
## pedaço diferente de dado real. Se a captura acaba em 2015, o modelo acaba
## em 2015, e o que se reporta é o estado em 2015.
##
## Buraco NO MEIO da série é outra coisa: ali a interpolação é defensável,
## porque existe observação dos dois lados ancorando. Ponta não tem âncora
## de um dos lados; ponta se corta.
MIN_ANOS_CAPTURA <- 15   # abaixo disso não vale rodar um modelo de produção
## Definido aqui (não só na seção 6, que é onde é mais usado) porque a
## seção 5.8 também precisa dele para decidir se a janela curta 2015-2025
## tem observações suficientes de C5.
MIN_OBS_INDICE   <- 5

capt <- data.frame(Yr = cap_in$Yr, Catch = cap_in$Catch_macarellus_t)

if (FONTE_CAPTURA == "externa") {
  ext <- le(ARQ_CAPT_EXT)
  
  ## Identificação tolerante das colunas: o arquivo pode vir como
  ## Year/Catch, Yr/Catch, ano/captura... Falhar aqui com mensagem clara é
  ## melhor do que seguir com a coluna errada.
  ca <- names(ext)
  col_ano  <- ca[tolower(ca) %in% c("year", "yr", "ano", "anos")][1]
  col_capt <- ca[tolower(ca) %in% c("catch", "captura", "capturas",
                                    "landings", "desembarque", "desembarques")][1]
  if (is.na(col_ano) || is.na(col_capt))
    stop(sprintf("Em '%s' não achei as colunas de ano/captura. Colunas presentes: %s",
                 ARQ_CAPT_EXT, paste(ca, collapse = ", ")))
  ext <- data.frame(Yr = as.integer(ext[[col_ano]]),
                    Catch = as.numeric(ext[[col_capt]]))
  ext <- ext[!is.na(ext$Yr), ]
  ext <- ext[order(ext$Yr), ]
  cat(sprintf("Captura externa: %s | colunas '%s'/'%s' | %d anos (%d-%d)\n",
              ARQ_CAPT_EXT, col_ano, col_capt, nrow(ext), min(ext$Yr), max(ext$Yr)))
  
  ## --- checagem de coerência contra a série do cerco --------------------
  ## Duas coisas saem daqui, e as duas importam:
  ##
  ## (a) A RAZÃO total/cerco mede quanto a série do cerco deixava de fora.
  ##
  ## (b) Um ano em que a série "total" fica ABAIXO da série do cerco é
  ##     LOGICAMENTE IMPOSSIVEL — o cerco é um subconjunto do total. Quando
  ##     isso aparece, quase sempre é ano de cobertura incompleta da
  ##     compilação (o caso clássico: os últimos anos, fechados antes de
  ##     todos os reportes entrarem). Esses anos não podem entrar como se
  ##     fossem remoção observada: com captura subestimada o modelo conclui
  ##     que o estoque se recuperou, e o status terminal sai errado para o
  ##     lado otimista.
  cerco_ref <- data.frame(Yr = cap_in$Yr, cerco = cap_in$Catch_macarellus_t)
  com <- merge(cerco_ref, data.frame(Yr = ext$Yr, total = ext$Catch), by = "Yr")
  com <- com[!is.na(com$cerco) & !is.na(com$total) & com$cerco > 0, ]
  anos_incoerentes <- integer(0)
  if (nrow(com)) {
    com$razao <- com$total / com$cerco
    cat(sprintf("  razão total/cerco: mediana %.2f (%d anos em comum)\n",
                median(com$razao), nrow(com)))
    anos_incoerentes <- com$Yr[com$razao < 1]
    if (length(anos_incoerentes)) {
      cat(sprintf("  [INCOERENTE] total < cerco em %d ano(s): %s\n",
                  length(anos_incoerentes), paste(anos_incoerentes, collapse = ", ")))
      cat("     O cerco é subconjunto do total, então esses anos estão incompletos\n")
      cat("     na compilação externa. Serão CORTADOS (ver CORTAR_INCOERENTES).\n")
    }
  } else {
    cat("  [AVISO] nenhum ano em comum com a série do cerco — confira os anos do arquivo.\n")
  }
  
  ## Cortar os anos incoerentes do FIM da série. Só do fim: um ano
  ## incoerente no meio vira NA e cai na interpolação da seção 3.1, que ali
  ## é defensável. Ponte-se para FALSE se quiser inspecioná-los antes.
  CORTAR_INCOERENTES <- TRUE
  if (CORTAR_INCOERENTES && length(anos_incoerentes)) {
    ult_bom <- max(setdiff(ext$Yr, anos_incoerentes))
    ## só corta a cauda: anos incoerentes que estejam depois do último ano bom
    if (any(anos_incoerentes > ult_bom)) {
      cat(sprintf("  -> cauda cortada: série externa passa a terminar em %d\n", ult_bom))
      ext <- ext[ext$Yr <= ult_bom, ]
    }
    ainda <- intersect(anos_incoerentes, ext$Yr)
    if (length(ainda)) {
      cat(sprintf("  -> %d ano(s) incoerente(s) NO MEIO viram NA (serão interpolados): %s\n",
                  length(ainda), paste(ainda, collapse = ", ")))
      ext$Catch[ext$Yr %in% ainda] <- NA
    }
  }
  
  ## --- a janela do modelo passa a ser a da captura -----------------------
  ## Aqui é onde a regra do topo do bloco é aplicada. Nada de união com a
  ## janela antiga, nada de preencher o que vem depois do último ano
  ## observado: corta-se e pronto.
  obs_ext  <- ext$Yr[!is.na(ext$Catch)]
  anos_mod <- seq(min(obs_ext), max(obs_ext))
  capt <- data.frame(Yr = anos_mod, Catch = ext$Catch[match(anos_mod, ext$Yr)])
  cat(sprintf("  >>> JANELA DO MODELO = %d-%d (era %d-%d com a série do cerco)\n",
              min(anos_mod), max(anos_mod), min(cap_in$Yr), max(cap_in$Yr)))
  cat(sprintf("      tudo daqui para a frente — índices, tabelas, figuras — fica\n"))
  cat(sprintf("      cortado nessa janela. Não há status para depois de %d.\n",
              max(anos_mod)))
  if (length(anos_mod) < MIN_ANOS_CAPTURA)
    stop(sprintf("Só %d anos de captura (mínimo %d). Série curta demais para um modelo de produção.",
                 length(anos_mod), MIN_ANOS_CAPTURA))
  
} else if (FONTE_CAPTURA == "escala") {
  capt$Catch <- capt$Catch * FATOR_ESCALA
  cat(sprintf("Captura do cerco multiplicada por %.2f (PREMISSA — declarar no texto)\n",
              FATOR_ESCALA))
} else {
  cat("Captura = rede de cerco industrial (ver a ressalva acima)\n")
}

## --- 3.1 o JABBA exige captura em TODO ano do modelo -------------------
## Diferente do índice, que pode ter buracos (o JABBA simplesmente pula os
## anos sem observação), a captura entra na equação de biomassa de todo ano.
## Um NA ali quebra o balanço de massa daquele ano e o JAGS não roda.
## PONTA SE CORTA, MEIO SE INTERPOLA. Se sobrou NA no começo ou no fim,
## a janela encolhe até o primeiro/último ano observado — em nenhuma
## hipótese o script extrapola captura para fora do que foi medido.
falta <- is.na(capt$Catch)
if (all(falta)) stop("Nenhum ano com captura observada.")
if (any(falta)) {
  obs <- which(!falta)
  if (min(obs) > 1 || max(obs) < nrow(capt)) {
    cortados <- c(capt$Yr[seq_len(min(obs) - 1)],
                  if (max(obs) < nrow(capt)) capt$Yr[(max(obs) + 1):nrow(capt)])
    cat(sprintf("Sem captura nas pontas -> janela encolhida para %d-%d (fora: %s)\n",
                capt$Yr[min(obs)], capt$Yr[max(obs)], paste(cortados, collapse = ", ")))
    capt  <- capt[min(obs):max(obs), , drop = FALSE]
    falta <- is.na(capt$Catch)
  }
  ## O que restou é buraco interno: tem observação dos dois lados, então a
  ## interpolação linear tem âncora. Continua sendo dado inventado e tem de
  ## constar no texto, mas é uma invenção com suporte dos dois lados.
  if (any(falta)) {
    cat(sprintf("Buraco(s) no MEIO da série: %s\n",
                paste(capt$Yr[falta], collapse = ", ")))
    capt$Catch <- approx(capt$Yr[!falta], capt$Catch[!falta], xout = capt$Yr)$y
    cat("  -> interpolados linearmente entre os vizinhos (declarar no texto)\n")
  }
  if (length(capt$Yr) < MIN_ANOS_CAPTURA)
    stop(sprintf("Restaram %d anos de captura (mínimo %d).",
                 length(capt$Yr), MIN_ANOS_CAPTURA))
}
## Anos com captura sabidamente frágil, herdados da parte 03.
if ("confiavel" %in% names(cap_in)) {
  frag <- cap_in$Yr[!as.logical(cap_in$confiavel)]
  if (length(frag))
    cat(sprintf("Anos de captura marcados frágeis na parte 03 (entram assim mesmo): %s\n",
                paste(frag, collapse = ", ")))
}

ANOS <- capt$Yr
cat(sprintf("\nSérie final: %d-%d | captura média %.1f t | máxima %.1f t (%d) | total %.0f t\n",
            min(ANOS), max(ANOS), mean(capt$Catch), max(capt$Catch),
            capt$Yr[which.max(capt$Catch)], sum(capt$Catch)))


## =========================================================================
## 4) INDICES E ERRO — a matriz que o JABBA espera
## -------------------------------------------------------------------------
## Formato exigido: data.frame(ano, indice1, indice2, ...) com a MESMA
## sequência de anos da captura. Anos sem observação ficam NA.
##
## CV versus SE: o argumento `se` do JABBA quer o ERRO-PADRAO NA ESCALA LOG,
## não o CV. A parte 03 guardou CV, calculado como cv = sqrt(exp(se^2)-1).
## A volta exata é se = sqrt(log(1 + cv^2)). Para CV pequeno os dois quase
## coincidem (CV 0,2 -> se 0,198), mas passar CV como se fosse SE infla o
## erro assumido nos pontos mais incertos, e é justamente ali que o
## deslocamento pesa. Fazemos a conversão explícita.
## =========================================================================
cat("\n===== 4) INDICES E ERRO =====\n")
cv_para_se <- function(cv) sqrt(log(1 + cv^2))

## Alinha as tabelas de índice ao calendário da captura (janela completa,
## usada pela maioria dos cenários).
lin <- match(ANOS, idx_in$Yr)
idx_al <- idx_in[lin, , drop = FALSE]; idx_al$Yr <- ANOS
cv_al  <- cv_in[lin,  , drop = FALSE]; cv_al$Yr  <- ANOS

cat(sprintf("CV -> SE(log): CV 0,20 vira SE %.3f | CV 0,45 vira SE %.3f\n",
            cv_para_se(0.20), cv_para_se(0.45)))

## Monta o par (cpue, se) de um conjunto de colunas de cenário. `rotulos`
## são os nomes que vão aparecer nas figuras do JABBA. `anos` default é a
## janela completa do modelo, mas um cenário pode pedir uma janela mais
## curta (ver `janela_ini` na grade da seção 6, usado pelos cenários
## J13-J17 isolados em 2015-2025) — por isso a função realinha direto de
## idx_in/cv_in em vez de usar idx_al/cv_al, que são fixos na janela cheia.
monta_indices <- function(colunas, rotulos = colunas, anos = ANOS) {
  lin_i <- match(anos, idx_in$Yr)
  cp <- data.frame(Yr = anos)
  sd <- data.frame(Yr = anos)
  for (k in seq_along(colunas)) {
    cp[[rotulos[k]]] <- idx_in[[colunas[k]]][lin_i]
    sd[[rotulos[k]]] <- cv_para_se(cv_in[[colunas[k]]][lin_i])
  }
  list(cpue = cp, se = sd)
}


## =========================================================================
## 5) PRIORS — de onde vem cada número
## -------------------------------------------------------------------------
## Num modelo de produção com dados escassos, o prior não é detalhe técnico:
## é metade do resultado. Cada um abaixo vem com a sua origem.
## =========================================================================
cat("\n===== 5) PRIORS =====\n")

## --- 5.1 r: taxa intrínseca de crescimento populacional ----------------
## `r.dist = "range"` faz o JABBA converter um intervalo (mín, máx) numa
## lognormal. `"lnorm"` esperaria (média, CV).
##
## ORIGEM DO INTERVALO: o FishBase classifica D. macarellus como resiliência
## MEDIA (tempo mínimo de duplicação populacional de 1,4 a 4,4 anos;
## K = 0,28-0,8; tm = 1-2). Na tabela de Froese et al. (2017), usada pelo
## CMSY, resiliência média corresponde a r entre 0,2 e 0,8.
##
## >>> VOCE JA RODOU CMSY NESTE ESTOQUE. Se tem a posterior de r de lá, ela
## é um prior MUITO melhor do que a faixa genérica: troque por
## r.dist = "lnorm" e r.prior = c(mediana, CV). Só não use a posterior do
## CMSY e depois apresente as duas avaliações como independentes — elas
## deixam de ser.
R_DIST  <- "range"
R_PRIOR <- c(0.20, 0.80)

## --- 5.2 K: biomassa não explorada (B0) --------------------------------
## Raciocínio para uma faixa defensável, em vez de um chute:
## no equilíbrio de Schaefer, MSY = r*K/4, logo K = 4*MSY/r. Se a pescaria
## andou em torno do MSY, MSY é da ordem da captura média. Com r na faixa
## acima, K cai entre 4*Cmed/0,8 e 4*Cmed/0,2 — ou seja, de 5 a 20 vezes a
## captura média. Abrimos bem em volta disso para não amarrar o modelo.
C_MED  <- mean(capt$Catch)
K_PLAUS <- c(4 * C_MED / R_PRIOR[2], 4 * C_MED / R_PRIOR[1])
K_DIST  <- "range"
K_PRIOR <- c(max(capt$Catch), 60 * C_MED)
cat(sprintf("Captura média %.1f t -> K plausível por MSY=rK/4: %.0f a %.0f t\n",
            C_MED, K_PLAUS[1], K_PLAUS[2]))
cat(sprintf("Prior de K adotado (bem mais largo): %.0f a %.0f t\n", K_PRIOR[1], K_PRIOR[2]))

## --- 5.2b sensibilidade ao prior de K -----------------------------------
## Por que isto entrou: na primeira rodada, TODOS os oito posteriores de K
## caíram entre 0,74x e 1,28x da mediana do prior. Quando o posterior é o
## prior de volta, e MSY = rK/4, o status que sai é o prior — não o dado.
## O prior de r já tinha os seus dois cenários de sensibilidade (J7/J8); o
## de K não tinha nenhum, e é ele que decide se as capturas históricas
## ficaram acima ou abaixo do MSY.
##
## A faixa base tem mediana em torno de 12x a captura média. Os dois
## cenários abaixo ficam a cerca de 3x para baixo e 3x para cima disso, que
## é o mesmo espaçamento (em escala log) do par J7/J8 para r.
##
## O PISO NÃO É NEGOCIÁVEL: K menor que a maior captura de um único ano é
## fisicamente impossível — o estoque teria sido removido inteiro naquele
## ano. Por isso todo limite inferior passa por max(capt$Catch).
## Extraída como função porque a seção 5.8 repete exatamente este raciocínio
## para a janela curta (2015-2025, cenários J13-J17): o piso físico e o
## espaçamento log-simétrico baixo/alto valem para qualquer sub-período,
## só muda a captura média que ancora a escala.
deriva_k_priors <- function(catch_vec) {
  c_med <- mean(catch_vec)
  piso  <- max(catch_vec)
  list(c_med = c_med, piso = piso,
       base  = c(piso, 60 * c_med),
       baixo = c(piso, max(8 * c_med, piso * 4)),
       alto  = c(max(20 * c_med, piso * 2), 100 * c_med))
}
K_PISO <- max(capt$Catch)
K_PRIOR_BAIXO <- c(K_PISO, max(8 * C_MED, K_PISO * 4))
K_PRIOR_ALTO  <- c(max(20 * C_MED, K_PISO * 2), 100 * C_MED)
.gm <- function(v) sqrt(v[1] * v[2])   # o JABBA converte "range" usando a média geométrica
cat(sprintf("Priors de K para sensibilidade (mediana implícita entre parênteses):\n"))
cat(sprintf("  baixo: %.0f-%.0f t (%.0f) | base: %.0f-%.0f t (%.0f) | alto: %.0f-%.0f t (%.0f)\n",
            K_PRIOR_BAIXO[1], K_PRIOR_BAIXO[2], .gm(K_PRIOR_BAIXO),
            K_PRIOR[1], K_PRIOR[2], .gm(K_PRIOR),
            K_PRIOR_ALTO[1], K_PRIOR_ALTO[2], .gm(K_PRIOR_ALTO)))
cat(sprintf("  MSY implícito em cada um (Schaefer, r = %.2f): %.0f | %.0f | %.0f t/ano\n",
            .gm(R_PRIOR), .gm(R_PRIOR) * .gm(K_PRIOR_BAIXO) / 4,
            .gm(R_PRIOR) * .gm(K_PRIOR) / 4, .gm(R_PRIOR) * .gm(K_PRIOR_ALTO) / 4))
cat(sprintf("  (compare com a captura média de %.0f t/ano: é isso que separa\n", C_MED))
cat("   'a pescaria histórica esteve acima do MSY' de 'nunca chegou perto')\n")

## --- 5.3 psi: depleção inicial B[1]/K ----------------------------------
## Quão explorado já estava o estoque no primeiro ano da série (1989)?
## psi = 1 significa virgem. A pescaria industrial de cerco em Cabo Verde
## já operava antes de 1989, mas a série histórica não sugere um estoque
## colapsado no início. 0,9 com CV 0,25 é um prior frouxo em torno de
## "pouco explorado" — e é o default do JABBA.
## ESTE PRIOR IMPORTA MUITO quando a série de índice não cobre o início,
## que é exatamente o nosso caso nos cenários que só usam C4/C5.
PSI_DIST  <- "lnorm"
PSI_PRIOR <- c(0.90, 0.25)

## --- 5.4 erro de processo ----------------------------------------------
## O erro de processo é a variação da biomassa que a curva de produção não
## explica (recrutamento variável, ambiente). `igamma = c(4, 0.01)` é uma
## inversa-gama fracamente informativa — o default do JABBA. Vale conferir
## o que ela implica antes de aceitar:
IGAMMA <- c(4, 0.01)
set.seed(1); .g <- 1 / rgamma(10000, IGAMMA[1], IGAMMA[2])
cat(sprintf("Prior de erro de processo: sigma médio %.3f (CV %.2f), faixa 10-90%%: %.3f-%.3f\n",
            sqrt(mean(.g)), sd(sqrt(.g)) / mean(sqrt(.g)),
            quantile(sqrt(.g), 0.1), quantile(sqrt(.g), 0.9)))

## --- 5.5 erro de observação --------------------------------------------
## Erro total de observação = sqrt(SE^2 + sigma.est^2 + fixed.obsE^2).
## Como passamos SE de verdade, `fixed.obsE` fica pequeno — ele existe para
## o caso de não haver SE nenhum. `sigma.est = TRUE` deixa o modelo estimar
## variância adicional por índice, o que é o certo quando se desconfia que
## o SE da padronização subestima o erro real (e aqui subestima: ele não
## contém a incerteza da tática inferida).
FIXED_OBSE <- 0.01
SIGMA_EST  <- TRUE

## --- 5.6 erro na captura ------------------------------------------------
## `add.catch.CV = TRUE` deixa o modelo tratar a captura como observada com
## erro em vez de conhecida exatamente. Aqui isso é quase obrigatório: a
## série emenda duas fontes (agregada até 2014, viagem a viagem depois),
## tem um ano traduzido de outro formato (2018) e um ano sem esforço (2013).
## 0,2 é deliberadamente maior que o default de 0,1.
CATCH_CV <- 0.2

## --- 5.8 janela curta 2015-2025: priors dedicados para o C5 isolado -----
## Pedido explícito: testar o índice padronizado (C5) sozinho, isolado no
## período em que ele existe, SEM a âncora da série histórica — e variar r
## e K nele do mesmo jeito que já se faz no caso base (J7/J8, J10/J11).
## Dois ajustes são obrigatórios em relação ao caso base, e os dois mudam o
## que esses cenários conseguem responder:
##
##   (1) K: refeito com deriva_k_priors() usando a captura média DE
##       2015-2025, não a da série toda — a frota mudou de alvo, a captura
##       média do período recente é bem menor, e ancorar K num C_MED da
##       série inteira (1989-2025) puxaria o prior para uma escala que este
##       recorte de dados não tem como sustentar.
##
##   (2) psi (B[1989... ou aqui B[2015]/K): o prior base (0,90, CV 0,25)
##       assume estoque pouco explorado NO PRIMEIRO ANO DO MODELO — razoável
##       para 1989 (início da série), mas FALSO por construção para 2015: a
##       pescaria de cerco já operava havia 26 anos. Não sabemos se em 2015
##       o estoque estava perto do virgem, já bem reduzido, ou em algum
##       ponto no meio — é exatamente essa ignorância que o prior tem de
##       carregar. Por isso PSI_PRIOR_CURTO é deliberadamente vago (média
##       0,5, CV alto), em vez de reciclar o prior de 1989.
##
## CONSEQUENCIA PARA A LEITURA DOS RESULTADOS: J13-J17 não têm como dar um
## status absoluto do estoque — eles respondem só "dado que a depleção em
## 2015 era X (desconhecido, por isso o prior largo), o que os últimos 10
## anos de C5 mudam nisso, e o quanto r/K pesam nessa conta". Não são
## comparáveis por DIC nem pelo painel de status aos cenários J1-J12 (janela
## e prior de psi diferentes) — servem para a pergunta específica que foi
## feita: o problema do ajuste é do período recente sozinho, ou é geral?
ANO_CURTO_INICIO  <- 2015   # primeiro ano com C5 (padronizada com tática)
capt_curto        <- capt[capt$Yr >= ANO_CURTO_INICIO, ]
CURTO_VIAVEL       <- nrow(capt_curto) >= MIN_OBS_INDICE

if (CURTO_VIAVEL) {
  kp_curto <- deriva_k_priors(capt_curto$Catch)
  PSI_PRIOR_CURTO <- c(0.5, 0.6)   # vago de propósito — ver comentário acima
  cat(sprintf("\nJanela curta viável: %d-%d (%d anos, captura média %.1f t)\n",
              min(capt_curto$Yr), max(capt_curto$Yr), nrow(capt_curto), kp_curto$c_med))
  cat(sprintf("  K (base/baixo/alto) para a janela curta: %.0f-%.0f | %.0f-%.0f | %.0f-%.0f t\n",
              kp_curto$base[1], kp_curto$base[2], kp_curto$baixo[1], kp_curto$baixo[2],
              kp_curto$alto[1], kp_curto$alto[2]))
  cat(sprintf("  psi vago adotado: media %.2f, CV %.2f (NAO é o prior de 1989 reciclado)\n",
              PSI_PRIOR_CURTO[1], PSI_PRIOR_CURTO[2]))
} else {
  cat(sprintf("\n[AVISO] Janela curta 2015-2025 inviável nesta rodada (%d anos < mínimo de %d) —\n",
              nrow(capt_curto), MIN_OBS_INDICE))
  cat("        cenários J13-J17 não serão criados.\n")
}


## =========================================================================
## 6) A GRADE DE CENARIOS
## -------------------------------------------------------------------------
## Cada linha é uma rodada do JABBA. O que muda entre elas:
##
##   QUAIS INDICES ENTRAM. Este é o ponto central, e vale entender por quê:
##   C2 (nominal até 2014) e C3/C4/C5 (2015-2025) são períodos com
##   CAPTURABILIDADE DIFERENTE — antes de 2015 a cavala era alvo, depois
##   virou acompanhante. Pôr os dois num índice só (que é o que C1 faz)
##   obriga o modelo a assumir um q constante atravessando justamente a
##   quebra que sabemos existir. Pôr os dois como índices SEPARADOS deixa o
##   JABBA estimar um q para cada um — é assim que se acomoda a mudança de
##   alvo sem jogar fora a série histórica.
##
##   O PREÇO: C2 termina em 2014 e C5 começa em 2015, sem nenhum ano em
##   comum. Índices que nunca se sobrepõem não têm como se calibrar
##   mutuamente — quem liga os dois períodos é a série de CAPTURA e a
##   dinâmica de produção. O modelo funciona, mas os priors de r e K pesam
##   mais. Some isso ao achado de Stobberup & Erzini (2006) de que a
##   posterior de r ficava igual à prior e fica claro por que a
##   sensibilidade de prior (J6/J7) não é enfeite.
##
##   A FORMA DA CURVA DE PRODUCAO (`model.type`):
##     Schaefer  parábola simétrica, Bmsy = K/2. O clássico.
##     Fox       Bmsy = 0,37*K — produção máxima em biomassa mais baixa,
##               costuma ser mais adequado a pelágicos pequenos produtivos.
##     Pella     Pella-Tomlinson com o ponto de inflexão FIXO em `BmsyK`.
##     Pella_m   idem, mas estimando a forma com prior CV = `shape.CV`.
## =========================================================================
cat("\n===== 6) CENARIOS =====\n")

## Edite esta grade à vontade: o resto do script é dirigido por ela.
## `indices` é uma string com as colunas separadas por "+".
cenarios <- data.frame(
  id      = c("J1_C1_ingenuo", "J2_C2C3_nominal", "J3_BASE_C2C5",
              "J4_C2C4_semTat", "J5_soC5_curto", "J6_BASE_Fox",
              "J7_BASE_rBaixo", "J8_BASE_rAlto",
              "J9_soC3_curto", "J10_BASE_kBaixo", "J11_BASE_kAlto",
              "J12_soC2_historico", "J18_BASE_C2C6"),
  indices = c("C1_nominal_total",
              "C2_nominal_pre_alvo+C3_nominal_pos_alvo",
              "C2_nominal_pre_alvo+C5_padronizada_com_tatica",
              "C2_nominal_pre_alvo+C4_padronizada_sem_tatica",
              "C5_padronizada_com_tatica",
              "C2_nominal_pre_alvo+C5_padronizada_com_tatica",
              "C2_nominal_pre_alvo+C5_padronizada_com_tatica",
              "C2_nominal_pre_alvo+C5_padronizada_com_tatica",
              "C3_nominal_pos_alvo",
              "C2_nominal_pre_alvo+C5_padronizada_com_tatica",
              "C2_nominal_pre_alvo+C5_padronizada_com_tatica",
              "C2_nominal_pre_alvo",
              "C2_nominal_pre_alvo+C6_nominal_dirigida"),
  modelo  = c("Schaefer", "Schaefer", "Schaefer", "Schaefer",
              "Schaefer", "Fox", "Schaefer", "Schaefer",
              "Schaefer", "Schaefer", "Schaefer", "Schaefer", "Schaefer"),
  r_min   = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.40,
              0.20, 0.20, 0.20, 0.20, 0.20),
  r_max   = c(0.80, 0.80, 0.80, 0.80, 0.80, 0.80, 0.50, 1.20,
              0.80, 0.80, 0.80, 0.80, 0.80),
  ## K por cenário: NA significa "usa o prior base K_PRIOR". Só os dois
  ## cenários de sensibilidade de K preenchem isto.
  k_min   = c(NA, NA, NA, NA, NA, NA, NA, NA,
              NA, K_PRIOR_BAIXO[1], K_PRIOR_ALTO[1], NA, NA),
  k_max   = c(NA, NA, NA, NA, NA, NA, NA, NA,
              NA, K_PRIOR_BAIXO[2], K_PRIOR_ALTO[2], NA, NA),
  ## janela_ini / psi_mu / psi_cv: NA em todos os cenários J1-J12/J18 (usam
  ## a janela completa e o PSI_PRIOR de 1989). Só os cenários J13-J17/J19
  ## (seção 5.8, anexados abaixo) preenchem isto.
  janela_ini = NA_real_,
  psi_mu     = NA_real_,
  psi_cv     = NA_real_,
  papel   = c("MODELO NULO do teste de capturabilidade: um q só atravessando a troca de alvo",
              "nominal honesto: q separado por período",
              "CASO BASE: nominal antes + padronizada com tática depois",
              "REDUNDANTE com J3 (C4 x C5 correlacionam 0,994) — ver 6.2",
              "índice recente sozinho, mas ainda ancorado na captura histórica",
              "sensibilidade: forma da curva de produção (Fox)",
              "sensibilidade: prior de r mais baixo (menos produtivo)",
              "sensibilidade: prior de r mais alto (mais produtivo)",
              "diagnóstico: só o índice NOMINAL recente (par do J5)",
              "sensibilidade: prior de K mais baixo (estoque menor)",
              "sensibilidade: prior de K mais alto (estoque maior)",
              "CASO BASE da janela histórica: só C2, frota e alvo estáveis",
              "sensibilidade: troca C5 (padronizada) por C6 (nominal, so esforco dirigido) no caso base"),
  ## ATIVO: cenário desligado continua DOCUMENTADO na planilha do racional
  ## (seção 6.2) — não some do registro, só não gasta uma rodada de MCMC.
  ## Para religar, basta trocar para TRUE aqui.
  ativo   = c(TRUE, TRUE, TRUE,
              FALSE,                       # J4 — ver motivo abaixo
              TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE)
cenarios$k_min[is.na(cenarios$k_min)] <- K_PRIOR[1]
cenarios$k_max[is.na(cenarios$k_max)] <- K_PRIOR[2]

## --- 6.1 anexa os cenários da janela curta (2015-2025), se viável -------
## Mesma lógica de J7/J8 (r) e J10/J11 (K), só que aplicada a um índice
## isolado na janela curta em vez de ao caso base na janela completa. Ver
## seção 5.8 para a justificativa de psi vago e de K reancorado na
## captura recente. J19 é o par de J13 com C6 (nominal esforço dirigido)
## no lugar de C5 (padronizada) — mesmo tratamento de janela/psi/K, só
## muda o índice, para isolar o efeito de TROCAR o índice do efeito de
## ANCORAR (ou não) no histórico, do mesmo jeito que o quadrado J2/J9/J3/J5
## já faz para C3/C5.
if (CURTO_VIAVEL) {
  cenarios_curto <- data.frame(
    id      = c("J13_C5curto_base", "J14_C5curto_rBaixo", "J15_C5curto_rAlto",
                "J16_C5curto_kBaixo", "J17_C5curto_kAlto", "J19_C6curto_base"),
    indices = c(rep("C5_padronizada_com_tatica", 5), "C6_nominal_dirigida"),
    modelo  = "Schaefer",
    r_min   = c(R_PRIOR[1], 0.10, 0.40, R_PRIOR[1], R_PRIOR[1], R_PRIOR[1]),
    r_max   = c(R_PRIOR[2], 0.50, 1.20, R_PRIOR[2], R_PRIOR[2], R_PRIOR[2]),
    k_min   = c(kp_curto$base[1], kp_curto$base[1], kp_curto$base[1],
                kp_curto$baixo[1], kp_curto$alto[1], kp_curto$base[1]),
    k_max   = c(kp_curto$base[2], kp_curto$base[2], kp_curto$base[2],
                kp_curto$baixo[2], kp_curto$alto[2], kp_curto$base[2]),
    janela_ini = ANO_CURTO_INICIO,
    psi_mu     = PSI_PRIOR_CURTO[1],
    psi_cv     = PSI_PRIOR_CURTO[2],
    papel   = c("JANELA CURTA 2015-2025: só C5, sem âncora na série histórica",
                "janela curta: prior de r mais baixo",
                "janela curta: prior de r mais alto",
                "janela curta: prior de K mais baixo (ancorado na captura 2015-2025)",
                "janela curta: prior de K mais alto (ancorado na captura 2015-2025)",
                "JANELA CURTA 2015-2025: par do J13 trocando C5 por C6 (esforco dirigido)"),
    ativo   = TRUE,
    stringsAsFactors = FALSE)
  cenarios <- rbind(cenarios, cenarios_curto)
}

## --- 6.2 a planilha do racional dos cenários ----------------------------
## Esta tabela é um PRODUTO, não um log: é ela que um revisor lê para
## entender por que a grade tem o tamanho que tem e o que cada rodada está
## respondendo. Sai em CSV e em XLSX, antes de qualquer ajuste (para
## existir mesmo se a rodada morrer no meio) e de novo na seção 8, aí já
## com os resultados de cada cenário ao lado do racional.
##
## A COLUNA `contraste` É A MAIS IMPORTANTE: cenário de avaliação não tem
## significado sozinho, só contra outro. É ela que diz com quem cada linha
## deve ser lida em par — e é o que evita a leitura errada mais comum, que
## é tratar um cenário isolado como "o resultado".
doc_pergunta <- c(
  J1_C1_ingenuo      = "A capturabilidade mudou na troca de alvo de 2015? J1 e o modelo NULO: um unico q atravessando a quebra.",
  J2_C2C3_nominal    = "O que a CPUE NOMINAL diz quando se admite a quebra de capturabilidade (um q por periodo)?",
  J3_BASE_C2C5       = "CASO BASE: nominal no periodo em que a cavala era alvo + padronizada no periodo recente, um q para cada.",
  J4_C2C4_semTat     = "A covariavel de tatica na padronizacao muda a avaliacao?",
  J5_soC5_curto      = "O indice recente padronizado sozinho, mas com toda a serie de CAPTURA 1989-2025 ainda ancorando a dinamica.",
  J6_BASE_Fox        = "A conclusao depende da forma da curva de producao (Fox: Bmsy/K = 0,37 em vez de 0,5)?",
  J7_BASE_rBaixo     = "O resultado vem do dado ou da prior de r? Lado menos produtivo.",
  J8_BASE_rAlto      = "O resultado vem do dado ou da prior de r? Lado mais produtivo.",
  J9_soC3_curto      = "O indice recente NOMINAL sozinho, com a captura historica ancorando a dinamica.",
  J10_BASE_kBaixo    = "O resultado vem do dado ou da prior de K? Estoque menor.",
  J11_BASE_kAlto     = "O resultado vem do dado ou da prior de K? Estoque maior.",
  J12_soC2_historico = "So o periodo em que a cavala era alvo: o que a serie historica diz sozinha?",
  J18_BASE_C2C6      = "Corrigir o direcionamento por SUBSETTING (C6, so esforco da tatica-cavala) em vez de por PADRONIZACAO (C5) muda a conclusao?",
  J13_C5curto_base   = "Sem nenhuma ancora historica (nem indice nem captura): o que os 10 anos recentes dizem sozinhos?",
  J14_C5curto_rBaixo = "Na janela curta, o resultado vem do dado ou da prior de r? Lado menos produtivo.",
  J15_C5curto_rAlto  = "Na janela curta, o resultado vem do dado ou da prior de r? Lado mais produtivo.",
  J16_C5curto_kBaixo = "Na janela curta, o resultado vem do dado ou da prior de K? Estoque menor.",
  J17_C5curto_kAlto  = "Na janela curta, o resultado vem do dado ou da prior de K? Estoque maior.",
  J19_C6curto_base   = "Na janela curta, trocar o indice padronizado (C5) pelo nominal dirigido (C6) muda a conclusao?")

doc_contraste <- c(
  J1_C1_ingenuo      = "J2. J1 usa OS MESMOS 35 pontos de indice que J2 (C1 e C2 e C3 renormalizados: razao constante, conferida) e esta ANINHADO nele — J2 = J1 + 1 q + 1 variancia. E a UNICA comparacao de DIC formalmente valida da grade.",
  J2_C2C3_nominal    = "J1 (um q so) e J3 (mesma estrutura, indice recente padronizado).",
  J3_BASE_C2C5       = "J2 (recente nominal), J18 (recente dirigido), J6 (forma da curva), J7/J8 e J10/J11 (priors).",
  J4_C2C4_semTat     = "J3. DESLIGADO por redundancia: C4 e C5 correlacionam 0,994 (r em log 0,987) e a propria parte 03 ja imprime esse numero. Rodar os dois no JABBA repete um resultado que a correlacao entre indices ja da, e da de forma mais limpa.",
  J5_soC5_curto      = "J3 (com indice historico) e J13 (mesmo indice, sem ancora nenhuma). O par J5 x J13 separa o efeito de perder o INDICE historico do de perder a CAPTURA historica.",
  J6_BASE_Fox        = "J3 (Schaefer, mesmos dados) — comparacao de DIC valida.",
  J7_BASE_rBaixo     = "J3 e J8. Se a posterior de r seguir a prior nos tres, r nao foi estimado (ver secao 12).",
  J8_BASE_rAlto      = "J3 e J7.",
  J9_soC3_curto      = "J2 (mesmo indice, com o historico) e J5 (mesma posicao, indice padronizado). Fecha o quadrado J2/J9/J3/J5.",
  J10_BASE_kBaixo    = "J3 e J11.",
  J11_BASE_kAlto     = "J3 e J10.",
  J12_soC2_historico = "J9/J5/J13 — e o polo OPOSTO do contraste de ancora: so o periodo antigo contra so o periodo recente.",
  J18_BASE_C2C6      = "J3 (padronizacao) e J19 (mesmo indice C6, sem ancora).",
  J13_C5curto_base   = "J5 (mesmo indice, com a captura historica) e J19 (mesma janela, indice C6).",
  J14_C5curto_rBaixo = "J13 e J15.",
  J15_C5curto_rAlto  = "J13 e J14.",
  J16_C5curto_kBaixo = "J13 e J17.",
  J17_C5curto_kAlto  = "J13 e J16.",
  J19_C6curto_base   = "J13 (mesma janela, indice padronizado) e J18 (mesmo indice, com ancora historica).")

doc_familia <- c(rep("janela completa 1989-2025", 13), rep("janela curta 2015-2025", 6))
names(doc_familia) <- c("J1_C1_ingenuo","J2_C2C3_nominal","J3_BASE_C2C5","J4_C2C4_semTat",
                        "J5_soC5_curto","J6_BASE_Fox","J7_BASE_rBaixo","J8_BASE_rAlto",
                        "J9_soC3_curto","J10_BASE_kBaixo","J11_BASE_kAlto","J12_soC2_historico",
                        "J18_BASE_C2C6","J13_C5curto_base","J14_C5curto_rBaixo",
                        "J15_C5curto_rAlto","J16_C5curto_kBaixo","J17_C5curto_kAlto",
                        "J19_C6curto_base")

## `cenarios_todos` guarda a grade INTEIRA (inclusive o que está desligado),
## porque a planilha do racional tem de registrar também o que ficou de fora
## e por quê — essa é metade da justificativa do desenho.
cenarios_todos <- cenarios
monta_doc <- function(g) {
  data.frame(
    id        = g$id,
    familia   = ifelse(g$id %in% names(doc_familia), doc_familia[g$id], NA),
    ativo     = g$ativo,
    pergunta  = ifelse(g$id %in% names(doc_pergunta),  doc_pergunta[g$id],  NA),
    contraste = ifelse(g$id %in% names(doc_contraste), doc_contraste[g$id], NA),
    papel     = g$papel,
    indices   = g$indices,
    janela    = ifelse(is.na(g$janela_ini), sprintf("%d-%d", min(ANOS), max(ANOS)),
                       sprintf("%d-%d", g$janela_ini, max(ANOS))),
    curva     = g$modelo,
    prior_r   = sprintf("range [%.2f, %.2f] (mediana %.2f)", g$r_min, g$r_max,
                        sqrt(g$r_min * g$r_max)),
    prior_K   = sprintf("range [%.0f, %.0f] t (mediana %.0f)", g$k_min, g$k_max,
                        sqrt(g$k_min * g$k_max)),
    prior_psi = ifelse(is.na(g$psi_mu),
                       sprintf("lnorm(%.2f, CV %.2f) - estoque pouco explorado em %d",
                               PSI_PRIOR[1], PSI_PRIOR[2], min(ANOS)),
                       sprintf("lnorm(%.2f, CV %.2f) - VAGO: nao se sabe a depleção em %d",
                               g$psi_mu, g$psi_cv, ifelse(is.na(g$janela_ini), min(ANOS), g$janela_ini))),
    row.names = NULL, stringsAsFactors = FALSE)
}
doc_cenarios <- monta_doc(cenarios_todos)
escreve_csv_utf8(doc_cenarios, file.path(DIR_SAIDA, "jabba_cenarios_racional.csv"))
if (tem_writexl)
  writexl::write_xlsx(list(racional_cenarios = doc_cenarios),
                      file.path(DIR_SAIDA, "jabba_cenarios_racional.xlsx"))
cat(sprintf("\nRacional dos cenários gravado (CSV + XLSX): %d linhas\n", nrow(doc_cenarios)))

## --- filtro de cenários desligados --------------------------------------
if (any(!cenarios$ativo)) {
  cat("\nCenários DESLIGADOS por redundância (continuam documentados na planilha):\n")
  for (i in which(!cenarios$ativo))
    cat(sprintf("  %-18s %s\n", cenarios$id[i], cenarios$papel[i]))
  cenarios <- cenarios[cenarios$ativo, ]
}

## POR QUE J9 EXISTE. Ele fecha um quadrado que estava com três lados:
##
##                  com histórico (C2 +)      sozinho (2015-2025)
##   nominal        J2                        J9   <- novo
##   padronizado    J3                        J5
##
## Na primeira rodada, J2 (nominal + histórico) deu B/B0 = 0,95 e J3
## (padronizado + histórico) deu 0,055 — mesma captura, mesmos priors, só
## mudando o índice recente. J5 (padronizado sozinho) deu 0,51. Faltava
## saber o que o índice nominal sozinho diz. Com as quatro células dá para
## separar o efeito de PADRONIZAR do efeito de ANCORAR no histórico, que
## hoje estão confundidos.
##
## POR QUE J12 EXISTE, e por que ele é o caso base quando a captura vem da
## compilação externa. Se a série de captura acaba em 2015, os índices de
## 2015-2025 ficam fora da janela e todos os cenários que dependem deles
## caem no filtro logo abaixo. Sobra a janela histórica — e nela C2 é o
## índice certo: mesmo período que a captura, mesma frota, e um período em
## que a cavala era alvo do começo ao fim, sem quebra de capturabilidade
## dentro. É a combinação mais internamente coerente que estes dados
## permitem. O que ele NAO dá é status atual: dá r, K, MSY e o estado do
## estoque no último ano da série de captura.
##
## POR QUE J10 E J11 EXISTEM. São para K o que J7/J8 são para r: mesma
## configuração do caso base, mexendo só no prior de K. Se B/Bmsy andar
## junto com o prior, K não está sendo estimado — está sendo assumido, e o
## status é consequência dessa premissa. Os dois lados importam: com só o
## lado baixo você não distingue "K manda no resultado" de "K baixo manda
## no resultado".

## Um cenário só pode rodar se (a) as colunas que ele pede existirem e
## (b) cada uma delas tiver observação suficiente DENTRO DA JANELA DO
## MODELO. O (b) passou a importar quando a janela virou a da captura: com
## a série externa terminando em 2015, os índices de 2015-2025 sobram com
## um ponto ou nenhum, e um índice de um ponto não informa trajetória
## nenhuma — só serve para o modelo estimar um q que passe por ele. Rodar
## assim produz um resultado que PARECE um cenário e não é.
disp <- setdiff(names(idx_in), "Yr")
n_obs_janela <- vapply(disp, function(cn) {
  v <- idx_in[[cn]][match(ANOS, idx_in$Yr)]
  sum(!is.na(v))
}, integer(1))

cat("Observações de cada índice dentro da janela completa do modelo:\n")
for (cn in disp)
  cat(sprintf("  %-28s %2d %s\n", cn, n_obs_janela[[cn]],
              if (n_obs_janela[[cn]] < MIN_OBS_INDICE) "<<< insuficiente" else ""))

## A checagem é por cenário, não só por índice: um cenário com `janela_ini`
## preenchido (J13-J17) usa uma janela mais curta que a completa, e é contra
## ESSA janela que a cobertura do índice tem de ser conferida.
verifica_cenario <- function(indices_str, janela_ini) {
  anos_i <- if (is.na(janela_ini)) ANOS else ANOS[ANOS >= janela_ini]
  v <- strsplit(indices_str, "+", fixed = TRUE)[[1]]
  if (!all(v %in% disp)) return("coluna inexistente")
  n_obs <- vapply(v, function(cn) sum(!is.na(idx_in[[cn]][match(anos_i, idx_in$Yr)])), integer(1))
  ruins <- v[n_obs < MIN_OBS_INDICE]
  if (length(ruins))
    return(sprintf("índice com < %d obs na janela %d-%d: %s", MIN_OBS_INDICE,
                   min(anos_i), max(anos_i), paste(ruins, collapse = ", ")))
  ""
}
motivo <- mapply(verifica_cenario, cenarios$indices, cenarios$janela_ini)

if (any(motivo != "")) {
  cat("\nCenários DESCARTADOS nesta janela:\n")
  for (i in which(motivo != ""))
    cat(sprintf("  %-16s %s\n", cenarios$id[i], motivo[i]))
  cat("  (não é erro: são cenários que dependem de um período que a série de\n")
  cat("   captura desta rodada não cobre. Para rodá-los é preciso captura\n")
  cat("   naquele período — ver a seção 3.)\n\n")
  cenarios <- cenarios[motivo == "", ]
}
if (!nrow(cenarios))
  stop("Nenhum cenário viável nesta janela. Reveja FONTE_CAPTURA na seção 3.")
cat(sprintf("%d cenários a rodar:\n", nrow(cenarios)))
for (i in seq_len(nrow(cenarios)))
  cat(sprintf("  %-16s %-9s r~[%.2f,%.2f] K~[%.0f,%.0f]  %s\n",
              cenarios$id[i], cenarios$modelo[i], cenarios$r_min[i], cenarios$r_max[i],
              cenarios$k_min[i], cenarios$k_max[i], cenarios$papel[i]))

## Rótulo curto e legível para cada índice nas figuras.
rot_curto <- c(C1_nominal_total = "C1 nom total",
               C2_nominal_pre_alvo = "C2 nom pre",
               C3_nominal_pos_alvo = "C3 nom pos",
               C4_padronizada_sem_tatica = "C4 padr s/tat",
               C5_padronizada_com_tatica = "C5 padr c/tat",
               C6_nominal_dirigida = "C6 nom dirigida")


## =========================================================================
## 7) AJUSTE
## -------------------------------------------------------------------------
## `build_jabba()` monta o objeto de entrada e escreve o modelo em código
## JAGS; `fit_jabba()` roda o MCMC. Separar os dois é útil: o mesmo
## `jbinput` serve para o ajuste, para a retrospectiva e para a projeção.
##
## `sets.q = 1:n` dá um q a cada índice (o default, e o que queremos).
## Para forçar dois índices a COMPARTILHAREM um q — por exemplo, se você
## quisesse testar que a mudança de alvo não alterou a capturabilidade —
## seria sets.q = c(1,1).
## =========================================================================
cat("\n===== 7) AJUSTE =====\n")

fits    <- list()   # os ajustes
jbinputs <- list()  # as entradas correspondentes (a retrospectiva da seção 10 precisa delas)
for (i in seq_len(nrow(cenarios))) {
  cn   <- cenarios[i, ]
  cols <- strsplit(cn$indices, "+", fixed = TRUE)[[1]]
  rots <- ifelse(cols %in% names(rot_curto), rot_curto[cols], cols)
  
  ## Janela deste cenário: a completa, salvo os J13-J17 (janela_ini
  ## preenchido na seção 6.1), que rodam só em 2015-2025.
  anos_i <- if (is.na(cn$janela_ini)) ANOS else ANOS[ANOS >= cn$janela_ini]
  capt_i <- capt[capt$Yr %in% anos_i, , drop = FALSE]
  mi     <- monta_indices(cols, rots, anos = anos_i)
  
  ## psi deste cenário: o prior de 1989 (PSI_PRIOR), salvo se o cenário
  ## pediu um vago próprio (psi_mu/psi_cv — ver seção 5.8).
  psi_prior_i <- c(if (is.na(cn$psi_mu)) PSI_PRIOR[1] else cn$psi_mu,
                   if (is.na(cn$psi_cv)) PSI_PRIOR[2] else cn$psi_cv)
  
  cat(sprintf("\n--- [%d/%d] %s | %s | %d índice(s) | janela %d-%d | r ~ [%.2f, %.2f] | K ~ [%.0f, %.0f] | psi ~ (%.2f, cv %.2f) ---\n",
              i, nrow(cenarios), cn$id, cn$modelo, length(cols),
              min(anos_i), max(anos_i), cn$r_min, cn$r_max,
              cn$k_min, cn$k_max, psi_prior_i[1], psi_prior_i[2]))
  
  jbin <- try(build_jabba(
    catch        = capt_i,
    cpue         = mi$cpue,
    se           = mi$se,
    assessment   = ASSESSMENT,
    scenario     = cn$id,
    model.type   = cn$modelo,
    add.catch.CV = TRUE,
    catch.cv     = CATCH_CV,
    r.dist       = R_DIST,
    r.prior      = c(cn$r_min, cn$r_max),
    K.dist       = K_DIST,
    K.prior      = c(cn$k_min, cn$k_max),
    psi.dist     = PSI_DIST,
    psi.prior    = psi_prior_i,
    sigma.est    = SIGMA_EST,
    fixed.obsE   = FIXED_OBSE,
    sigma.proc   = TRUE,
    igamma       = IGAMMA,
    sets.q       = seq_along(cols),
    sets.var     = seq_along(cols),
    catch.metric = "(t)",
    verbose      = TRUE), silent = TRUE)
  
  if (inherits(jbin, "try-error")) {
    cat("  [ERRO] build_jabba falhou:\n  ", conditionMessage(attr(jbin, "condition")), "\n")
    next
  }
  
  ft <- try(fit_jabba(jbin,
                      ni = MCMC_NI, nb = MCMC_NB, nt = MCMC_NT, nc = MCMC_NC,
                      quickmcmc  = RODADA_RAPIDA,
                      save.jabba = TRUE,
                      save.csvs  = TRUE,
                      output.dir = DIR_SAIDA,
                      verbose    = FALSE), silent = TRUE)
  
  if (inherits(ft, "try-error")) {
    cat("  [ERRO] fit_jabba falhou:\n  ", conditionMessage(attr(ft, "condition")), "\n")
    next
  }
  ## Guardados em listas SEPARADAS de propósito: pendurar o jbinput como
  ## atributo do ajuste funcionaria, mas o objeto `fits` é passado inteiro
  ## para funções do JABBA e é melhor não modificá-lo.
  fits[[cn$id]]     <- ft
  jbinputs[[cn$id]] <- jbin
  cat(sprintf("  ok | MSY = %.1f t | B/Bmsy final = %.2f | F/Fmsy final = %.2f\n",
              ft$refpts$msy[1],
              ft$timeseries[nrow(ft$timeseries), "mu", "BBmsy"],
              ft$timeseries[nrow(ft$timeseries), "mu", "FFmsy"]))
}

if (length(fits) == 0) stop("Nenhum cenário ajustou. Veja as mensagens de erro acima.")
cat(sprintf("\n%d de %d cenários ajustaram.\n", length(fits), nrow(cenarios)))


## =========================================================================
## 8) CONSOLIDACAO — os resultados de todos os cenários numa planilha só
## -------------------------------------------------------------------------
## O objeto que o fit_jabba devolve tem tudo, mas espalhado. Aqui ele vira
## quatro tabelas retangulares, que é o formato em que se escreve artigo:
##   tab_refpts    pontos de referência por cenário (K, r, MSY, Bmsy, Fmsy)
##   tab_status    situação no ano final + probabilidades dos quadrantes Kobe
##   tab_series    as trajetórias completas, formato longo
##   tab_ajuste    qualidade do ajuste e convergência
## =========================================================================
cat("\n===== 8) CONSOLIDACAO =====\n")

## --- 8.1 pontos de referência ------------------------------------------
## `$estimates` traz mu/lci/uci com os parâmetros nas linhas. Puxamos pelo
## nome da linha para não depender da ordem.
pega <- function(est, alvo) {
  i <- match(tolower(alvo), tolower(rownames(est)))
  if (is.na(i)) return(c(NA, NA, NA))
  as.numeric(est[i, c("mu", "lci", "uci")])
}
tab_refpts <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]; e <- f$estimates
  linha <- function(par) { v <- pega(e, par)
  data.frame(cenario = nm, parametro = par, mu = v[1], lci = v[2], uci = v[3],
             row.names = NULL, stringsAsFactors = FALSE) }
  do.call(rbind, lapply(c("K", "r", "psi", "sigma.proc", "m",
                          "Hmsy", "SBmsy", "MSY", "BmsyK"), linha))
}))
tab_refpts <- tab_refpts[!is.na(tab_refpts$mu), ]
cat("\n--- pontos de referência ---\n")
print(transform(tab_refpts, mu = signif(mu, 4), lci = signif(lci, 4), uci = signif(uci, 4)),
      row.names = FALSE)

## --- 8.2 situação final e probabilidades de Kobe -----------------------
## `$kobe` guarda as AMOSTRAS POSTERIORES do ano final (stock = B/Bmsy,
## harvest = F/Fmsy). Contar a fração das amostras em cada quadrante dá
## direto a probabilidade — que é o jeito bayesiano de responder "qual a
## chance de o estoque estar sobrepescado?", e é bem mais informativo do
## que reportar só o ponto estimado.
tab_status <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]; k <- f$kobe; n <- nrow(f$timeseries)
  data.frame(
    cenario     = nm,
    ano_final   = f$yr[length(f$yr)],
    BBmsy       = f$timeseries[n, "mu",  "BBmsy"],
    BBmsy_lci   = f$timeseries[n, "lci", "BBmsy"],
    BBmsy_uci   = f$timeseries[n, "uci", "BBmsy"],
    FFmsy       = f$timeseries[n, "mu",  "FFmsy"],
    FFmsy_lci   = f$timeseries[n, "lci", "FFmsy"],
    FFmsy_uci   = f$timeseries[n, "uci", "FFmsy"],
    BB0         = f$timeseries[n, "mu",  "BB0"],
    ## quadrantes: verde = saudável e sem sobrepesca; vermelho = os dois problemas
    P_verde     = mean(k$stock >= 1 & k$harvest <  1),
    P_amarelo   = mean(k$stock >= 1 & k$harvest >= 1),
    P_laranja   = mean(k$stock <  1 & k$harvest <  1),
    P_vermelho  = mean(k$stock <  1 & k$harvest >= 1),
    P_sobrepescado    = mean(k$stock  <  1),   # B < Bmsy
    P_sofrendo_pesca  = mean(k$harvest >= 1),  # F > Fmsy
    row.names = NULL, stringsAsFactors = FALSE)
}))
cat("\n--- situação no ano final (probabilidades vêm da posterior) ---\n")
print(transform(tab_status,
                BBmsy = round(BBmsy, 2), BBmsy_lci = round(BBmsy_lci, 2),
                BBmsy_uci = round(BBmsy_uci, 2), FFmsy = round(FFmsy, 2),
                FFmsy_lci = round(FFmsy_lci, 2), FFmsy_uci = round(FFmsy_uci, 2),
                BB0 = round(BB0, 2), P_verde = round(P_verde, 3),
                P_amarelo = round(P_amarelo, 3), P_laranja = round(P_laranja, 3),
                P_vermelho = round(P_vermelho, 3),
                P_sobrepescado = round(P_sobrepescado, 3),
                P_sofrendo_pesca = round(P_sofrendo_pesca, 3)), row.names = FALSE)

## --- 8.3 trajetórias completas -----------------------------------------
tab_series <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]; ts <- f$timeseries
  do.call(rbind, lapply(dimnames(ts)[[3]], function(v)
    data.frame(cenario = nm, ano = as.numeric(dimnames(ts)[[1]]), variavel = v,
               mu = ts[, "mu", v], lci = ts[, "lci", v], uci = ts[, "uci", v],
               row.names = NULL, stringsAsFactors = FALSE)))
}))
cat(sprintf("\n--- trajetórias: %d linhas (%s) ---\n", nrow(tab_series),
            paste(unique(tab_series$variavel), collapse = ", ")))

## --- 8.4 ajuste e convergência -----------------------------------------
## `$stats` traz SDNR, RMSE e DIC. (Cuidado: no JABBA a coluna se chama
## literalmente "Stastistic", com o erro de digitação — usamos o nome real.)
##
## COMO LER:
##   RMSE   erro médio do ajuste ao índice, em %. Abaixo de ~30% é bom.
##   SDNR   desvio-padrão dos resíduos normalizados. Perto de 1 significa
##          que o erro assumido bate com o erro observado; muito acima de 1
##          quer dizer que você disse ao modelo que os dados eram mais
##          precisos do que são.
##   DIC    critério de comparação entre modelos com os MESMOS dados. NÃO
##          serve para comparar cenários com índices diferentes — J1 e J3
##          não são comparáveis por DIC, só J3 contra J6 (mesma entrada,
##          curva diferente).
##   Geweke / Heidelberger: diagnósticos de CONVERGENCIA do MCMC, por
##          parâmetro. p pequeno é sinal de cadeia que não convergiu — aí o
##          resultado não vale, independente de quão bonito esteja o ajuste.
tab_ajuste <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  s <- f$stats; v <- setNames(s[[2]], s[[1]])
  pv <- f$pars
  data.frame(cenario = nm,
             N = v[["N"]], RMSE = v[["RMSE"]], SDNR = v[["SDNR"]], DIC = v[["DIC"]],
             min_Geweke_p = if ("Geweke.p" %in% names(pv)) min(pv$Geweke.p, na.rm = TRUE) else NA,
             min_Heidel_p = if ("Heidel.p" %in% names(pv)) min(pv$Heidel.p, na.rm = TRUE) else NA,
             row.names = NULL, stringsAsFactors = FALSE)
}))
cat("\n--- ajuste e convergência ---\n")
print(transform(tab_ajuste, RMSE = round(RMSE, 1), SDNR = round(SDNR, 3),
                DIC = round(DIC, 1), min_Geweke_p = round(min_Geweke_p, 3),
                min_Heidel_p = round(min_Heidel_p, 3)), row.names = FALSE)
if (any(!is.na(tab_ajuste$min_Geweke_p) & tab_ajuste$min_Geweke_p < 0.05))
  cat("[AVISO] há cenário com Geweke p < 0,05: cadeia possivelmente não convergida.\n",
      "        Aumente MCMC_NI/MCMC_NB (e desligue RODADA_RAPIDA) antes de interpretar.\n")

## --- 8.5 resíduos por índice -------------------------------------------
tab_residuos <- do.call(rbind, lapply(names(fits), function(nm) {
  d <- fits[[nm]]$diags
  if (is.null(d)) return(NULL)
  data.frame(cenario = nm, indice = d$name, ano = d$year,
             observado = d$obs, ajustado = d$hat, residuo_log = d$residual,
             row.names = NULL, stringsAsFactors = FALSE)
}))

## --- 8.6 escrita -------------------------------------------------------
## Racional + resultado na MESMA planilha: é assim que ela vira tabela de
## artigo (o leitor vê, na mesma linha, a pergunta que o cenário responde e
## o número que ele produziu) em vez de duas tabelas que ninguém cruza.
doc_cenarios_result <- merge(
  doc_cenarios,
  merge(tab_ajuste[, c("cenario", "N", "RMSE", "SDNR", "DIC")],
        tab_status[, c("cenario", "ano_final", "BBmsy", "BBmsy_lci", "BBmsy_uci",
                       "FFmsy", "BB0", "P_sobrepescado")],
        by = "cenario", all = TRUE),
  by.x = "id", by.y = "cenario", all.x = TRUE)
doc_cenarios_result <- doc_cenarios_result[match(doc_cenarios$id, doc_cenarios_result$id), ]
doc_cenarios_result$rodou <- !is.na(doc_cenarios_result$RMSE)
escreve_csv_utf8(doc_cenarios_result, file.path(DIR_SAIDA, "jabba_cenarios_racional.csv"))
if (tem_writexl)
  writexl::write_xlsx(list(racional_cenarios = doc_cenarios_result),
                      file.path(DIR_SAIDA, "jabba_cenarios_racional.xlsx"))
cat("Racional dos cenários REGRAVADO com os resultados (CSV + XLSX)\n")

escreve_csv_utf8(tab_refpts,   file.path(DIR_SAIDA, "jabba_pontos_referencia.csv"))
escreve_csv_utf8(tab_status,   file.path(DIR_SAIDA, "jabba_status_final.csv"))
escreve_csv_utf8(tab_series,   file.path(DIR_SAIDA, "jabba_trajetorias.csv"))
escreve_csv_utf8(tab_ajuste,   file.path(DIR_SAIDA, "jabba_ajuste_convergencia.csv"))
if (!is.null(tab_residuos))
  escreve_csv_utf8(tab_residuos, file.path(DIR_SAIDA, "jabba_residuos.csv"))
if (tem_writexl) {
  abas <- list(racional_cenarios = doc_cenarios_result,
               pontos_referencia = tab_refpts, status_final = tab_status,
               trajetorias = tab_series, ajuste = tab_ajuste, cenarios = cenarios_todos)
  if (!is.null(tab_residuos)) abas$residuos <- tab_residuos
  writexl::write_xlsx(abas, file.path(DIR_SAIDA, "JABBA_macarellus_resultados.xlsx"))
}
cat(sprintf("\nPlanilhas gravadas em %s\n", DIR_SAIDA))


## =========================================================================
## 9) FIGURAS
## -------------------------------------------------------------------------
## Duas famílias: as do próprio JABBA (uma pasta por cenário, completas) e
## um painel comparativo construído aqui a partir das tabelas consolidadas —
## esse último é o que costuma ir para o artigo.
## =========================================================================
cat("\n===== 9) FIGURAS =====\n")

## --- 9.1 o pacote completo, por cenário ---------------------------------
## jabba_plots() escreve de uma vez: ajuste aos índices, resíduos, teste de
## runs, distribuições prior-posterior, trajetórias, fase de produção, Kobe.
for (nm in names(fits)) {
  d <- file.path(DIR_SAIDA, nm); dir.create(d, showWarnings = FALSE)
  try(jabba_plots(jabba = fits[[nm]], output.dir = d), silent = TRUE)
}
cat("Figuras padrão do JABBA: uma subpasta por cenário.\n")

## --- 9.2 comparação entre cenários (função do próprio JABBA) ------------
try({
  jbplot_summary(fits, prefix = "comparacao_cenarios", as.png = TRUE,
                 save.summary = TRUE, output.dir = DIR_SAIDA)
  cat("PNG salvo: comparacao_cenarios (jbplot_summary)\n")
}, silent = TRUE)

## --- 9.3 painel comparativo próprio -------------------------------------
## Quatro quadros: captura, B/Bmsy, F/Fmsy e o plano de Kobe com o ponto
## final de cada cenário. É a figura que responde "a conclusão muda conforme
## o cenário?" — que é a pergunta real quando se roda um conjunto deles.
## Cores geradas dinamicamente para o número de cenários que de fato
## ajustaram (a grade cresceu de 8 para até 17 com J13-J17): uma paleta fixa
## e curta recicla cor quando a grade passa do seu tamanho, e duas linhas da
## mesma cor no painel comparativo tornariam a figura ilegível justamente
## onde ela precisa ser lida.
cor_cen <- setNames(grDevices::hcl.colors(length(fits), palette = "Dark 3"),
                    names(fits))

png(file.path(DIR_SAIDA, "jabba_painel_cenarios.png"),
    width = 26, height = 22, res = 300, antialias = AA, units = "cm")
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1), bty = "l",
          cex.main = 0.95, cex = 0.8)

## (A) captura. Mostra a série completa; os cenários J13-J17 (janela curta,
## 2015-2025) usam só o trecho final dela — ver seção 5.8.
plot(capt$Yr, capt$Catch, type = "h", lwd = 4, lend = 1, col = COR_NEU,
     ylim = c(0, max(capt$Catch) * 1.05),
     xlab = "Ano", ylab = "Captura (t)",
     main = "A. Remocoes usadas no modelo (J13-J17 usam so 2015-2025)")

## (B) e (C) trajetórias
quadro <- function(v, ylab, titulo, ref = 1) {
  s <- tab_series[tab_series$variavel == v, ]
  plot(NA, xlim = range(s$ano), ylim = c(0, max(s$uci, na.rm = TRUE) * 1.05),
       xlab = "Ano", ylab = ylab, main = titulo)
  abline(h = ref, lty = 2, col = "grey40")
  for (nm in names(fits)) {
    d <- s[s$cenario == nm, ]
    lines(d$ano, d$mu, lwd = 2.2, col = cor_cen[nm])
  }
}
quadro("BBmsy", "B / Bmsy", "B. Biomassa relativa")
quadro("FFmsy", "F / Fmsy", "C. Mortalidade por pesca relativa")

## (D) Kobe — só o ponto final de cada cenário, sobre os quadrantes
s_max <- max(c(tab_status$BBmsy_uci, 2), na.rm = TRUE)
f_max <- max(c(tab_status$FFmsy_uci, 2), na.rm = TRUE)
plot(NA, xlim = c(0, s_max), ylim = c(0, f_max), xlab = "B / Bmsy", ylab = "F / Fmsy",
     main = "D. Kobe - ano final por cenario")
rect(0, 1, 1, f_max, col = adjustcolor("#D7301F", 0.18), border = NA)  # vermelho
rect(1, 0, s_max, 1, col = adjustcolor("#4DAF4A", 0.18), border = NA)  # verde
rect(0, 0, 1, 1, col = adjustcolor("#FF7F00", 0.15), border = NA)      # laranja
rect(1, 1, s_max, f_max, col = adjustcolor("#FFD92F", 0.20), border = NA) # amarelo
abline(h = 1, v = 1, lty = 2, col = "grey30")
for (nm in names(fits)) {
  r <- tab_status[tab_status$cenario == nm, ]
  segments(r$BBmsy_lci, r$FFmsy, r$BBmsy_uci, r$FFmsy, col = cor_cen[nm], lwd = 1.4)
  segments(r$BBmsy, r$FFmsy_lci, r$BBmsy, r$FFmsy_uci, col = cor_cen[nm], lwd = 1.4)
  points(r$BBmsy, r$FFmsy, pch = 21, bg = cor_cen[nm], col = "white", cex = 1.8, lwd = 1.5)
}
## Legenda com fundo: sem ele os rótulos somem dentro do quadrante colorido.
legend("topright", names(fits), col = cor_cen, pch = 19, pt.cex = 1.2, cex = 0.58,
       bg = adjustcolor("white", 0.82), box.col = NA)
par(op); dev.off()
cat("PNG salvo: jabba_painel_cenarios.png\n")


## =========================================================================
## 10) DIAGNOSTICO — retrospectiva e validacao por hindcast
## -------------------------------------------------------------------------
## Ajuste bom não é o mesmo que modelo confiável. Dois testes separam as
## duas coisas, e ambos são exigidos em avaliação séria hoje (Carvalho et
## al. 2021, "A cookbook for using model diagnostics in integrated stock
## assessments", Fish. Res. 240:105959):
##
##   RETROSPECTIVA: reajusta o modelo removendo 1, 2, ... anos do fim e vê
##   se a estimativa do mesmo ano MUDA quando chegam dados novos. Um padrão
##   sistemático (sempre revisando a biomassa para baixo, por exemplo) é
##   viés retrospectivo. Resume-se no rho de Mohn; a referência usual é
##   |rho| dentro de cerca de 0,15-0,20.
##
##   HINDCAST CROSS-VALIDATION: mesma poda, mas em vez de olhar o parâmetro
##   olha a PREVISÃO do índice nos anos removidos, e compara com o
##   observado. O resumo é o MASE: MASE < 1 significa que o modelo prevê
##   melhor do que a regra ingênua de "ano que vem é igual a este ano".
##   MASE > 1 significa que ele não tem poder preditivo — e aí o índice
##   provavelmente não está informando o modelo.
##
## Rode isto no CASO BASE e nos cenários que você pretende apresentar, não
## em todos: cada um custa (n_peels + 1) ajustes.
## =========================================================================
if (RODAR_DIAGNOSTICO) {
  cat("\n===== 10) DIAGNOSTICO =====\n")
  CEN_DIAG <- intersect(c("J3_BASE_C2C5", "J2_C2C3_nominal"), names(fits))
  diag_dir <- file.path(DIR_SAIDA, "diagnostico")
  dir.create(diag_dir, showWarnings = FALSE)
  resumo_diag <- list()
  
  for (nm in CEN_DIAG) {
    cat(sprintf("\n--- %s ---\n", nm))
    jbin <- jbinputs[[nm]]
    hc <- try(hindcast_jabba(jbinput = jbin, fit = fits[[nm]], peels = 1:5), silent = TRUE)
    if (inherits(hc, "try-error")) {
      cat("  hindcast falhou: ", conditionMessage(attr(hc, "condition")), "\n"); next
    }
    rho  <- try(jbplot_retro(hc,  as.png = TRUE, single.plots = FALSE, output.dir = diag_dir),
                silent = TRUE)
    mase <- try(jbplot_hcxval(hc, as.png = TRUE, single.plots = FALSE, output.dir = diag_dir),
                silent = TRUE)
    if (!inherits(rho, "try-error"))  { cat("  rho de Mohn:\n"); print(rho) }
    if (!inherits(mase, "try-error")) { cat("  MASE:\n");        print(mase) }
    resumo_diag[[nm]] <- list(rho = rho, mase = mase)
  }
  
  ## Teste de runs nos resíduos: procura SEQUENCIAS de resíduos do mesmo
  ## sinal. Resíduo que fica anos seguidos acima e depois anos seguidos
  ## abaixo indica que o modelo não está acompanhando a tendência do índice
  ## — falha estrutural, que o RMSE sozinho não pega.
  for (nm in names(fits)) {
    rt <- try(jbplot_runstest(fits[[nm]], as.png = TRUE, output.dir = diag_dir), silent = TRUE)
    if (!inherits(rt, "try-error") && !is.null(rt)) {
      cat(sprintf("\nteste de runs — %s\n", nm)); print(rt)
    }
  }
  saveRDS(resumo_diag, file.path(diag_dir, "resumo_diagnostico.rds"))
} else {
  cat("\n[10] Diagnóstico desligado (RODAR_DIAGNOSTICO = FALSE).\n",
      "     Ligue antes de fechar os resultados: sem retrospectiva e hindcast\n",
      "     a avaliação não tem como ser defendida.\n")
}


## =========================================================================
## 11) PROJECOES — o que acontece com o estoque sob capturas alternativas
## -------------------------------------------------------------------------
## `TACs` é o vetor de capturas anuais constantes a testar; o modelo projeta
## a biomassa para cada uma. É o insumo direto de recomendação de manejo:
## a saída responde "com captura de X toneladas, qual a probabilidade de
## B > Bmsy em N anos?".
##   TACint  captura do ano-ponte (entre o fim dos dados e a implementação)
##   imp.yr  primeiro ano em que a TAC vale
##   pyrs    quantos anos projetar
## =========================================================================
if (RODAR_PROJECAO) {
  cat("\n===== 11) PROJECOES =====\n")
  ## Preferência: o caso base da janela que de fato rodou. Com a captura
  ## externa (janela histórica) o J3 nem existe — quem manda é o J12. Os
  ## cenários J13-J17 (janela curta, psi vago) ficam de fora desta escolha
  ## de propósito: projetar a partir de um psi deliberadamente vago não dá
  ## uma base defensável para recomendação de manejo.
  CEN_PROJ_CANDIDATOS <- setdiff(names(fits), grep("^J1[3-7]_", names(fits), value = TRUE))
  CEN_PROJ <- c("J3_BASE_C2C5", "J12_soC2_historico", CEN_PROJ_CANDIDATOS[1])
  CEN_PROJ <- CEN_PROJ[CEN_PROJ %in% names(fits)][1]
  cn <- cenarios[cenarios$id == CEN_PROJ, ]
  cols <- strsplit(cn$indices, "+", fixed = TRUE)[[1]]
  rots <- ifelse(cols %in% names(rot_curto), rot_curto[cols], cols)
  anos_proj <- if (is.na(cn$janela_ini)) ANOS else ANOS[ANOS >= cn$janela_ini]
  capt_proj <- capt[capt$Yr %in% anos_proj, , drop = FALSE]
  mi <- monta_indices(cols, rots, anos = anos_proj)
  psi_prior_proj <- c(if (is.na(cn$psi_mu)) PSI_PRIOR[1] else cn$psi_mu,
                      if (is.na(cn$psi_cv)) PSI_PRIOR[2] else cn$psi_cv)
  
  ## Grade de TAC em torno da captura recente: nada de números redondos
  ## arbitrários — ancorar no que a pescaria de fato retira torna a leitura
  ## da figura imediata.
  C_REC <- mean(tail(capt_proj$Catch, 3))
  TACS  <- round(seq(0.4 * C_REC, 1.6 * C_REC, length.out = 7))
  cat(sprintf("Cenário projetado: %s | Captura média dos últimos 3 anos: %.1f t | TACs testadas: %s\n",
              CEN_PROJ, C_REC, paste(TACS, collapse = ", ")))
  
  jbp <- build_jabba(catch = capt_proj, cpue = mi$cpue, se = mi$se,
                     assessment = ASSESSMENT, scenario = paste0(CEN_PROJ, "_proj"),
                     model.type = cn$modelo, add.catch.CV = TRUE, catch.cv = CATCH_CV,
                     r.dist = R_DIST, r.prior = c(cn$r_min, cn$r_max),
                     K.dist = K_DIST, K.prior = c(cn$k_min, cn$k_max),
                     psi.dist = PSI_DIST, psi.prior = psi_prior_proj,
                     sigma.est = SIGMA_EST, fixed.obsE = FIXED_OBSE,
                     sigma.proc = TRUE, igamma = IGAMMA,
                     sets.q = seq_along(cols), sets.var = seq_along(cols),
                     projection = TRUE, TACs = TACS,
                     TACint = C_REC, imp.yr = max(anos_proj) + 1, pyrs = 10,
                     catch.metric = "(t)", verbose = FALSE)
  
  fit_proj <- fit_jabba(jbp, ni = MCMC_NI, nb = MCMC_NB, nt = MCMC_NT, nc = MCMC_NC,
                        quickmcmc = RODADA_RAPIDA, save.csvs = TRUE,
                        output.dir = DIR_SAIDA, verbose = FALSE)
  
  ## jbplot_prj() é outra função "antiga" do JABBA (mesma família da
  ## jbplot_summary() que deu problema na seção 9): ela tenta ler
  ## fit_proj$projections, um campo de array que a versão atual do pacote —
  ## a que usa build_jabba()/fit_jabba() — nunca chega a criar. Conferi o
  ## código-fonte de fit_jabba() no repositório: o bloco de projeção
  ## ("compiling Future Projections under fixed quota") escreve as
  ## trajetórias projetadas dentro de fit_proj$kbtrj (a mesma tabela usada
  ## pelas figuras do painel da seção 9), marcando type="prj" e usando a
  ## coluna `run` para guardar o rótulo da TAC (formato "C<valor>"). Como
  ## fit_proj$projections nunca é preenchido, jbplot_prj() indexa um array
  ## vazio e o max()/min() internos dele reclamam "no non-missing arguments"
  ## — o aviso não é um problema do seu ajuste, é a função de plot que está
  ## desatualizada em relação ao objeto que fit_jabba() retorna hoje.
  ##
  ## Solução: montar o gráfico direto de fit_proj$kbtrj, do mesmo jeito que
  ## o JABBA monta internamente, só que sem depender da função quebrada.
  kbtrj_proj <- fit_proj$kbtrj[fit_proj$kbtrj$run %in% paste0("C", TACS), ]
  kbtrj_proj$run <- factor(as.character(kbtrj_proj$run),
                           levels = paste0("C", TACS))
  
  ## mediana e IC95% por ano e por TAC, pra cada variável de interesse
  resume_kb <- function(kb, var) {
    do.call(rbind, lapply(split(kb, list(kb$year, kb$run), drop = TRUE), function(d) {
      data.frame(ano = d$year[1], TAC = as.character(d$run[1]),
                 mediana = median(d[[var]], na.rm = TRUE),
                 li = as.numeric(quantile(d[[var]], 0.025, na.rm = TRUE)),
                 ls = as.numeric(quantile(d[[var]], 0.975, na.rm = TRUE)))
    }))
  }
  res_stock   <- resume_kb(kbtrj_proj, "stock")   # B/Bmsy
  res_harvest <- resume_kb(kbtrj_proj, "harvest") # F/Fmsy
  
  ## rampa de cor: verde (TAC mais conservadora) -> cinza -> laranja (TAC
  ## mais agressiva), na mesma paleta usada no resto do script
  cores_tac <- colorRampPalette(c(COR_S2, COR_NEU, COR_AUX))(length(TACS))
  names(cores_tac) <- paste0("C", TACS)
  
  plot_proj <- function(res, ylab) {
    plot(res$ano, res$mediana, type = "n",
         xlim = range(res$ano), ylim = range(0, res$li, res$ls, na.rm = TRUE),
         xlab = "Ano", ylab = ylab, bty = "l")
    abline(h = 1, lty = 2, col = COR_NEU)
    for (tac in levels(kbtrj_proj$run)) {
      d <- res[res$TAC == tac, ]
      d <- d[order(d$ano), ]
      if (nrow(d) < 2) next
      polygon(c(d$ano, rev(d$ano)), c(d$li, rev(d$ls)),
              col = adjustcolor(cores_tac[tac], 0.15), border = NA)
      lines(d$ano, d$mediana, col = cores_tac[tac], lwd = 2)
    }
    legend("topleft", legend = TACS, col = cores_tac, lwd = 2, bty = "n",
           title = "TAC (t)", cex = 0.7, ncol = 2)
  }
  
  png(file.path(DIR_SAIDA, "jabba_projecoes.png"), width = 26, height = 11,
      res = 300, antialias = AA, units = "cm")
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1), bty = "l", cex = 0.8)
  plot_proj(res_stock,   "B/Bmsy projetado")
  plot_proj(res_harvest, "F/Fmsy projetado")
  par(op); dev.off()
  cat("PNG salvo: jabba_projecoes.png\n")
} else {
  cat("\n[11] Projeções desligadas (RODAR_PROJECAO = FALSE).\n")
}


## =========================================================================
## 12) PRIORI x POSTERIORI, POSTERIORES COMPARADAS, RESIDUOS E TORNADO
## -------------------------------------------------------------------------
## Esta seção existe para responder, com figura e não com adjetivo, a
## pergunta central desta avaliação: O MODELO APRENDEU ALGUMA COISA COM OS
## DADOS, ou devolveu a prior de volta? Stobberup & Erzini (2006) já haviam
## reportado, neste mesmo estoque, que "a posterior marginal de r era quase
## idêntica à prior". Aqui isso deixa de ser citação e vira medida.
##
## DE ONDE VEM A PRIOR DESENHADA NAS FIGURAS: não é recalculada a partir da
## grade de cenários — é lida de dentro do próprio objeto ajustado, em
## `fit$settings`, onde o JABBA guarda os parâmetros JÁ CONVERTIDOS para a
## parametrização que o JAGS usou:
##   r.pr   = c(mediana, desvio-padrão em log)  <- de r.dist="range", a
##            mediana é sqrt(min*max) e o sd em log é log(max/min)/4
##   K.pr   = idem para K
##   psi.pr = c(mediana, sd em log)  <- o CV informado já vem convertido
##            por sqrt(log(1+CV^2))
##   igamma = c(forma, taxa) da inversa-gama do erro de processo: o JAGS
##            amostra 1/sigma2 ~ dgamma(igamma[1], igamma[2])
## Ler do objeto em vez de reconstruir elimina a classe de erro mais chata
## aqui: desenhar uma prior que não é a que o modelo usou.
##
## O q FICA DE FORA DO TESTE PRIORI x POSTERIORI, e o motivo está no próprio
## objeto: `settings$q1.pr` é c(1e-30, 1e+20) — cinquenta ordens de
## grandeza, ou seja, a prior de q é deliberadamente não-informativa. Dizer
## que "a posterior de q atualizou em relação à prior" seria verdadeiro e
## vazio. O que INTERESSA no q é outra coisa, e essa sim vai para a figura:
## o valor de cada q e, nos cenários com dois índices, a RAZÃO entre eles —
## que é a medida direta do salto de capturabilidade de 2015.
## =========================================================================
cat("\n===== 12) PRIORI x POSTERIORI E FIGURAS COMPARATIVAS =====\n")

DIR_FIG <- file.path(DIR_SAIDA, "figuras_posteriores")
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

## --- 12.1 utilitários ----------------------------------------------------

## Nome de arquivo seguro (o id do cenário pode ganhar caracteres estranhos
## se a grade for editada à mão).
arquivo_seguro <- function(x) {
  x <- gsub("[/\\\\]", "-", x); x <- gsub("[()]", "", x)
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

## Amostra da PRIOR de cada parâmetro, na mesma parametrização do JAGS.
## `n` alto porque a densidade da prior de K/MSY é muito espalhada e com
## poucas amostras o desenho fica serrilhado.
amostra_priori <- function(fit, n = 20000) {
  s <- fit$settings; out <- list()
  if (!is.null(s$r.pr))   out$r   <- rlnorm(n, log(s$r.pr[1]),   s$r.pr[2])
  if (!is.null(s$K.pr))   out$K   <- rlnorm(n, log(s$K.pr[1]),   s$K.pr[2])
  ## psi só entra se a prior for lognormal (é o nosso caso sempre; se algum
  ## dia virar "beta" esta amostragem não valeria e o painel some sozinho).
  if (!is.null(s$psi.pr) && identical(s$psi.dist, "lnorm"))
    out$psi <- rlnorm(n, log(s$psi.pr[1]), s$psi.pr[2])
  ## erro de processo: só faz sentido desenhar a prior se ele foi ESTIMADO.
  if (isTRUE(s$sigma.proc) && !is.null(s$igamma))
    out$sigma2 <- 1 / rgamma(n, s$igamma[1], s$igamma[2])
  ## MSY implícito pela prior: a mesma fórmula de Pella-Tomlinson que o
  ## JABBA usa, com o m do cenário. Schaefer (m=2) cai em rK/4; Fox (m->1)
  ## cai em rK/e. É isto que torna o painel de MSY comparável ao do CMSY.
  m  <- if (!is.null(s$mu.m)) s$mu.m else 2
  bk <- if (abs(m - 1) < 1e-6) exp(-1) else (1 / m)^(1 / (m - 1))
  if (!is.null(out$r) && !is.null(out$K)) out$MSY <- out$r * out$K * bk / m
  out
}

## Amostras da POSTERIOR, já reunidas por cenário.
draws_cen <- function(nm) {
  f  <- fits[[nm]]
  pp <- f$pars_posterior
  rp <- f$refpts_posterior
  qs <- pp[, grep("^q($|\\.)", names(pp)), drop = FALSE]
  ## B/K do último ano: vem das trajetórias por iteração (kbtrj), que é
  ## onde o JABBA guarda a posterior completa ano a ano. Usar a coluna BB0
  ## do último ano ajustado (type == "fit") dá o análogo exato do "Bt/K do
  ## último ano" que o CMSY++ reporta.
  bb0 <- NULL
  if (!is.null(f$kbtrj)) {
    kbf <- f$kbtrj[f$kbtrj$type == "fit", ]
    if (nrow(kbf)) bb0 <- kbf$BB0[kbf$year == max(kbf$year)]
  }
  list(r = pp$r, K = pp$K, MSY = rp$MSY, psi = pp$psi,
       sigma2 = if ("sigma2" %in% names(pp)) pp$sigma2 else NULL,
       BB0 = bb0, q = qs)
}

## Densidade normalizada ao pico, com a opção de trabalhar em escala log.
## POR QUE A ESCALA LOG PARA K, MSY E sigma2: a prior desses parâmetros é
## larga de propósito (K varia mais de uma ordem de grandeza dentro da
## própria prior). Desenhada em escala linear, a posterior vira um risco
## vertical e a figura não mostra nada. Em log, as duas curvas ficam
## legíveis e a COMPARAÇÃO — que é o objetivo — aparece. A densidade é
## calculada em log e convertida de volta pela mudança de variável
## f_X(x) = f_U(u)/x, com u = log(x), como no script do CMSY++.
dens_pico <- function(x, log_x = FALSE, from = NULL, to = NULL, adjust = 1.2) {
  x <- x[is.finite(x)]
  if (log_x) x <- x[x > 0]
  if (length(x) < 10) return(NULL)
  if (log_x) {
    d  <- density(log(x), adjust = adjust)
    xs <- exp(d$x); ys <- d$y / xs
  } else {
    d  <- if (!is.null(from) && !is.null(to))
      density(x, adjust = adjust, from = from, to = to) else density(x, adjust = adjust)
    xs <- d$x; ys <- d$y
  }
  list(x = xs, y = ys / max(ys))
}

## Um painel priori x posteriori.
## A JANELA DO EIXO X sai dos quantis 0,5%-99,5% das DUAS amostras juntas,
## não do intervalo total: uma única amostra extrema da prior de K (que é
## lognormal e tem cauda pesada) esticaria o eixo até deixar tudo ilegível.
## A prior continua existindo fora da janela; a legenda avisa.
painel_pri_post <- function(prior, post, xlab, col_post, main, log_x = FALSE,
                            lim01 = FALSE) {
  if (is.null(post) || !length(post)) { plot.new(); title(main = main, cex.main = 0.85)
    text(0.5, 0.5, "sem amostras", cex = 0.8, col = COR_NEU); return(invisible(NULL)) }
  vals <- c(if (!is.null(prior)) prior, post)
  vals <- vals[is.finite(vals)]
  xl <- as.numeric(quantile(vals, c(0.005, 0.995), na.rm = TRUE))
  if (lim01) xl <- c(0, 1)
  dp  <- if (!is.null(prior)) dens_pico(prior, log_x) else NULL
  dpo <- dens_pico(post, log_x, from = if (lim01) 0 else NULL, to = if (lim01) 1 else NULL)
  plot(NA, xlim = xl, ylim = c(0, 1.08), xlab = xlab,
       ylab = "densidade (normalizada ao pico)", main = main,
       bty = "l", log = if (log_x) "x" else "", cex.main = 0.85)
  if (!is.null(dp))  { polygon(c(dp$x, rev(dp$x)), c(dp$y, rep(0, length(dp$y))),
                               col = adjustcolor("grey45", 0.12), border = NA)
    lines(dp$x, dp$y, col = "grey45", lwd = 2.6, lty = 2) }
  lines(dpo$x, dpo$y, col = col_post, lwd = 2.6)
  ## Medianas marcadas: é delas que sai a razão posterior/prior do teste.
  if (!is.null(prior)) abline(v = median(prior), col = "grey45", lty = 3)
  abline(v = median(post), col = col_post, lty = 3)
  legend("topright", c(if (!is.null(prior)) "Priori", "Posteriori"),
         col = c(if (!is.null(prior)) "grey45", col_post),
         lty = c(if (!is.null(prior)) 2, 1), lwd = 2.2, bty = "n", cex = 0.72)
}

## Cores por cenário, geradas para o número de cenários que de fato rodou.
cor_cen_fig <- setNames(grDevices::hcl.colors(length(fits), palette = "Dark 3"),
                        names(fits))

## --- 12.2 painel priori x posteriori, um PNG por cenário -----------------
## Layout dinâmico: 6 painéis fixos (r, K, MSY, psi, sigma2, B/K final) +
## um painel por q do cenário. Cenário de 1 índice fica com 7; de 2 índices,
## com 8.
##
## COMO LER ESTA FIGURA, que é o ponto de tê-la: se a curva colorida
## (posterior) tiver praticamente a mesma forma e a mesma mediana da curva
## cinza (prior), aquele parâmetro NÃO foi estimado — foi assumido. Aí o
## MSY que sai dele, e o status que sai do MSY, são consequência da premissa
## e não do dado. A tabela `tab_atualizacao` logo abaixo põe número nisso.
cat("\n--- 12.2 priori x posteriori por cenário ---\n")

tab_atualizacao <- NULL
for (nm in names(fits)) {
  f  <- fits[[nm]]
  pr <- amostra_priori(f)
  po <- draws_cen(nm)
  nq <- ncol(po$q)
  
  paineis <- list(
    list(pri = pr$r,      pos = po$r,      lab = "r",              cor = COR_MAC, log = FALSE),
    list(pri = pr$K,      pos = po$K,      lab = "K (t)",          cor = COR_AUX, log = TRUE),
    list(pri = pr$MSY,    pos = po$MSY,    lab = "MSY (t/ano)",    cor = COR_S2,  log = TRUE),
    list(pri = pr$psi,    pos = po$psi,    lab = "psi (B inicial/K)", cor = COR_S3, log = FALSE),
    list(pri = pr$sigma2, pos = po$sigma2, lab = "sigma2 (erro de processo)", cor = COR_S2D, log = TRUE),
    list(pri = NULL,      pos = po$BB0,    lab = "B/K no ultimo ano", cor = "#B03060", log = FALSE, lim01 = TRUE))
  for (j in seq_len(nq))
    paineis[[length(paineis) + 1]] <- list(pri = NULL, pos = po$q[[j]],
                                           lab = sprintf("%s (capturabilidade)", names(po$q)[j]),
                                           cor = "#8C6D31", log = TRUE)
  
  nc <- 4; nr <- ceiling(length(paineis) / nc)
  png(file.path(DIR_FIG, sprintf("priori_posteriori_%s.png", arquivo_seguro(nm))),
      width = 9 * nc, height = 8.2 * nr, res = 300, antialias = AA, units = "cm")
  op <- par(mfrow = c(nr, nc), mar = c(4.2, 4.4, 2.8, 1), oma = c(0, 0, 2.2, 0),
            bty = "l", cex = 0.72)
  for (p in paineis)
    painel_pri_post(p$pri, p$pos, p$lab, p$cor,
                    main = sprintf("%s - %s", sub(" .*", "", p$lab), nm),
                    log_x = isTRUE(p$log), lim01 = isTRUE(p$lim01))
  mtext(sprintf("Priori x posteriori - %s", nm), outer = TRUE, cex = 0.85, font = 2)
  par(op); dev.off()
  
  ## --- o teste em número, não em figura ---------------------------------
  ## razao = mediana(posterior) / mediana(prior). Perto de 1 = o dado não
  ## moveu o parâmetro. `sobrepos` é a fração da posterior que cai dentro
  ## do intervalo de 95% da prior: 1 significa que a posterior é um
  ## subconjunto da prior (nenhuma informação nova sobre a região), e
  ## valores baixos significam que o dado empurrou o parâmetro para fora.
  linha <- function(par, pri, pos) {
    if (is.null(pri) || is.null(pos)) return(NULL)
    ic <- quantile(pri, c(0.025, 0.975), na.rm = TRUE)
    data.frame(cenario = nm, parametro = par,
               mediana_priori = median(pri, na.rm = TRUE),
               mediana_posterior = median(pos, na.rm = TRUE),
               razao_post_priori = median(pos, na.rm = TRUE) / median(pri, na.rm = TRUE),
               reducao_IC95_pct = 100 * (1 - (diff(quantile(pos, c(0.025, 0.975), na.rm = TRUE)) /
                                                diff(ic))),
               fracao_dentro_IC95_priori = mean(pos >= ic[1] & pos <= ic[2], na.rm = TRUE),
               row.names = NULL, stringsAsFactors = FALSE)
  }
  tab_atualizacao <- rbind(tab_atualizacao,
                           linha("r", pr$r, po$r), linha("K", pr$K, po$K),
                           linha("MSY", pr$MSY, po$MSY), linha("psi", pr$psi, po$psi),
                           linha("sigma2", pr$sigma2, po$sigma2))
}
cat(sprintf("PNGs salvos: priori_posteriori_<cenario>.png (%d cenários) em %s\n",
            length(fits), basename(DIR_FIG)))

cat("\n--- quanto o dado moveu cada parâmetro (razão posterior/priori) ---\n")
cat("razão ~1 + redução de IC ~0 = parâmetro NÃO estimado, apenas assumido\n")
print(transform(tab_atualizacao,
                mediana_priori = signif(mediana_priori, 4),
                mediana_posterior = signif(mediana_posterior, 4),
                razao_post_priori = round(razao_post_priori, 3),
                reducao_IC95_pct = round(reducao_IC95_pct, 1),
                fracao_dentro_IC95_priori = round(fracao_dentro_IC95_priori, 3)),
      row.names = FALSE)
escreve_csv_utf8(tab_atualizacao, file.path(DIR_SAIDA, "jabba_priori_posteriori.csv"))
if (tem_writexl)
  writexl::write_xlsx(list(priori_posteriori = tab_atualizacao),
                      file.path(DIR_SAIDA, "jabba_priori_posteriori.xlsx"))

## --- 12.2a a posterior de B/K é unimodal? --------------------------------
## POR QUE ESTA CHECAGEM ENTROU: no painel do caso base a posterior de B/K
## do último ano saiu com DOIS PICOS — um perto de 0,05 (estoque colapsado)
## e outro perto de 0,95 (estoque quase virgem) — e quase nenhuma massa no
## meio. Isso tem uma consequência direta que nenhuma tabela de mediana
## mostra: A MEDIANA CAI NO VALE. Reportar "B/Bmsy = 0,92" para um cenário
## assim é reportar justamente o valor em que o modelo acredita MENOS; o
## que ele está dizendo é "ou o estoque colapsou, ou está intacto, e os
## dados não escolhem". Um intervalo de credibilidade de 0,05 a 2,15 não
## avisa disso — parece apenas incerteza larga.
##
## O teste é deliberadamente simples e sem pacote extra: mede a massa da
## posterior em três faixas. Se as duas pontas têm massa relevante e o meio
## não, a distribuição é bimodal e a mediana não resume nada.
tab_bimodal <- do.call(rbind, lapply(names(fits), function(nm) {
  b <- draws_cen(nm)$BB0
  if (is.null(b) || !length(b)) return(NULL)
  p_baixo <- mean(b < 0.2); p_meio <- mean(b >= 0.2 & b <= 0.8); p_alto <- mean(b > 0.8)
  data.frame(cenario = nm,
             P_BK_menor_0.2 = p_baixo, P_BK_entre = p_meio, P_BK_maior_0.8 = p_alto,
             mediana_BK = median(b),
             bimodal = (p_baixo > 0.15 & p_alto > 0.15 & p_meio < 0.5),
             row.names = NULL, stringsAsFactors = FALSE)
}))
cat("\n--- a posterior de B/K final é unimodal? ---\n")
cat("bimodal = TRUE: a MEDIANA cai no vale entre dois picos e NAO resume o resultado\n")
print(transform(tab_bimodal, P_BK_menor_0.2 = round(P_BK_menor_0.2, 3),
                P_BK_entre = round(P_BK_entre, 3), P_BK_maior_0.8 = round(P_BK_maior_0.8, 3),
                mediana_BK = round(mediana_BK, 3)), row.names = FALSE)
if (any(tab_bimodal$bimodal)) {
  cat(sprintf("[ATENCAO] %d cenário(s) com posterior BIMODAL de B/K: %s\n",
              sum(tab_bimodal$bimodal),
              paste(tab_bimodal$cenario[tab_bimodal$bimodal], collapse = ", ")))
  cat("          Nesses, reportar as DUAS soluções e a probabilidade de cada uma,\n")
  cat("          não a mediana nem o IC de 95% (que atravessa o vale).\n")
}
escreve_csv_utf8(tab_bimodal, file.path(DIR_SAIDA, "jabba_bimodalidade_BK.csv"))

## --- 12.2b capturabilidade: o q e o salto de 2015 ------------------------
## Nos cenários de dois índices, q.1 é do índice histórico e q.2 do recente.
## A razão q.2/q.1 é a medida direta de quanto a capturabilidade mudou na
## troca de alvo — e é um número que o texto precisa, porque é ele que
## justifica ter separado os q em vez de emendar a série.
tab_q <- do.call(rbind, lapply(names(fits), function(nm) {
  q <- draws_cen(nm)$q
  if (!ncol(q)) return(NULL)
  ## a coluna da razão existe sempre (NA quando o cenário tem um índice só),
  ## senão o rbind entre cenários de 1 e de 2 índices quebra
  d <- data.frame(cenario = nm, q = names(q),
                  mediana = sapply(q, median),
                  lci = sapply(q, function(z) quantile(z, 0.025)),
                  uci = sapply(q, function(z) quantile(z, 0.975)),
                  razao_q2_q1 = NA_real_,
                  row.names = NULL, stringsAsFactors = FALSE)
  if (ncol(q) == 2) d$razao_q2_q1[2] <- median(q[[2]] / q[[1]])
  d
}))
cat("\n--- capturabilidade (q) por cenário ---\n")
print(transform(tab_q, mediana = signif(mediana, 3), lci = signif(lci, 3),
                uci = signif(uci, 3)), row.names = FALSE)
escreve_csv_utf8(tab_q, file.path(DIR_SAIDA, "jabba_capturabilidade_q.csv"))

## --- 12.3 boxplots das posteriores, todos os cenários lado a lado --------
## Mesma figura que existe no CMSY++ e no DB-SRA, para o conjunto ficar
## comparável entre os três métodos.
cat("\n--- 12.3 boxplots das posteriores por cenário ---\n")
vars_box <- list(r = "r", K = "K (t)", MSY = "MSY (t/ano)", BB0 = "B/K no ultimo ano")
draws_todos <- setNames(lapply(names(fits), draws_cen), names(fits))

png(file.path(DIR_FIG, "jabba_boxplots_posteriores.png"),
    width = 30, height = 24, res = 300, antialias = AA, units = "cm")
op <- par(mfrow = c(2, 2), mar = c(7.5, 5, 3, 1), oma = c(5, 0, 0, 0),
          bty = "l", cex = 0.78, cex.main = 0.9)
for (v in names(vars_box)) {
  lst <- lapply(names(fits), function(nm) {
    z <- draws_todos[[nm]][[v]]; if (is.null(z)) NA_real_ else z[is.finite(z)] })
  ## escala log em K e MSY pelo mesmo motivo do painel individual
  usa_log <- v %in% c("K", "MSY")
  boxplot(lst, names = names(fits), main = vars_box[[v]], col = "#8FAADC",
          outline = FALSE, ylab = vars_box[[v]], xaxt = "n", bty = "l",
          log = if (usa_log) "y" else "")
  if (v == "BB0") abline(h = 0.5, col = "firebrick", lty = 2)
  axis(1, at = seq_along(fits), labels = FALSE)
  text(x = seq_along(fits), y = par("usr")[3], labels = names(fits),
       srt = 45, adj = 1, xpd = NA, cex = 0.6)
}
par(op); dev.off()
cat("PNG salvo: jabba_boxplots_posteriores.png\n")

## --- 12.4 densidade conjunta de MSY e de B/K final -----------------------
## ATENÇÃO METODOLOGICA, e ela precisa ir para a legenda da figura: juntar
## num pool único as posteriores de cenários que discordam ESTRUTURALMENTE
## (aqui o B/Bmsy final varia de 0,07 a 2,6 conforme o cenário) não produz
## "a" distribuição do MSY — produz uma mistura, cuja largura é o desacordo
## entre cenários, não a incerteza de um modelo. Ela é útil exatamente como
## isso: uma medida visual do tamanho do desacordo. Por isso o pool exclui
## por padrão a família de janela curta, que responde a outra pergunta e
## nem sequer cobre os mesmos anos.
CEN_POOL <- names(fits)[is.na(cenarios$janela_ini[match(names(fits), cenarios$id)])]
if (!length(CEN_POOL)) CEN_POOL <- names(fits)
cat(sprintf("\n--- 12.4 densidades conjuntas (pool de %d cenários de janela completa) ---\n",
            length(CEN_POOL)))

dens_conjunta <- function(x, log_x, xlab, titulo, arquivo, lim01 = FALSE) {
  x <- x[is.finite(x)]; if (log_x) x <- x[x > 0]
  q <- quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))
  d <- dens_pico(x, log_x, from = if (lim01) 0 else NULL, to = if (lim01) 1 else NULL, adjust = 2)
  png(file.path(DIR_FIG, arquivo), width = 25, height = 16, res = 300,
      antialias = AA, units = "cm")
  op <- par(mar = c(4.8, 5, 3.4, 1), bty = "l", cex.main = 0.9)
  plot(d$x, d$y, type = "n", log = if (log_x) "x" else "",
       xlim = if (lim01) c(0, 1) else range(d$x),
       xlab = xlab, ylab = "densidade (normalizada ao pico)", main = titulo)
  faixa <- d$x >= q[1] & d$x <= q[5]
  polygon(c(d$x[faixa], rev(d$x[faixa])), c(d$y[faixa], rep(0, sum(faixa))),
          col = adjustcolor(COR_MAC, 0.18), border = NA)
  lines(d$x, d$y, col = COR_MAC, lwd = 2.4)
  abline(v = q[c(1, 5)], col = COR_AUX, lty = 2, lwd = 1.6)
  abline(v = q[3], col = COR_MAC, lty = 1, lwd = 1.8)
  mtext("pool de cenários: a largura mede o DESACORDO entre cenários, não a incerteza de um modelo",
        side = 3, line = 0.2, cex = 0.62, col = "grey30")
  legend("topright", c("densidade conjunta", "faixa de 95%", "mediana"),
         col = c(COR_MAC, adjustcolor(COR_MAC, 0.4), COR_MAC),
         lwd = c(2.4, 8, 1.8), bty = "n", cex = 0.72)
  par(op); dev.off()
  data.frame(quantil = c("2.5%", "25%", "mediana (50%)", "75%", "97.5%"),
             valor = as.numeric(q), row.names = NULL)
}

msy_pool <- unlist(lapply(CEN_POOL, function(nm) draws_todos[[nm]]$MSY))
btk_pool <- unlist(lapply(CEN_POOL, function(nm) draws_todos[[nm]]$BB0))
q_msy <- dens_conjunta(msy_pool, TRUE, "MSY (t/ano)",
                       "Distribuicao conjunta do MSY - cenarios de janela completa (JABBA)",
                       "jabba_msy_densidade_conjunta.png")
q_btk <- dens_conjunta(btk_pool, FALSE, "B/K no ultimo ano",
                       "Distribuicao conjunta de B/K final - cenarios de janela completa (JABBA)",
                       "jabba_btk_densidade_conjunta.png", lim01 = TRUE)
names(q_msy)[2] <- "MSY"; names(q_btk)[2] <- "BK_final"
cat("\nQuantis do MSY (pool):\n");     print(transform(q_msy, MSY = round(MSY, 1)), row.names = FALSE)
cat("\nQuantis do B/K final (pool):\n"); print(transform(q_btk, BK_final = round(BK_final, 3)), row.names = FALSE)
escreve_csv_utf8(q_msy, file.path(DIR_SAIDA, "jabba_quantis_msy_conjunto.csv"))
escreve_csv_utf8(q_btk, file.path(DIR_SAIDA, "jabba_quantis_btk_conjunto.csv"))
if (tem_writexl)
  writexl::write_xlsx(list(msy_conjunto = q_msy, btk_conjunto = q_btk,
                           capturabilidade_q = tab_q),
                      file.path(DIR_SAIDA, "jabba_posteriores_conjuntas.xlsx"))

## --- 12.5 trajetórias de biomassa de todos os cenários -------------------
## Dois painéis, como no CMSY++: relativa (B/Bmsy, comparável entre
## cenários) e absoluta (em toneladas, que é o que o manejo lê). As faixas
## são o IC de 95% — e aqui elas são o assunto da figura, não o enfeite:
## é a sobreposição delas que mostra que os cenários não se distinguem.
cat("\n--- 12.5 trajetórias de biomassa ---\n")
png(file.path(DIR_FIG, "jabba_trajetorias_biomassa.png"),
    width = 32, height = 15, res = 300, antialias = AA, units = "cm")
op <- par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3, 1), bty = "l", cex = 0.8, cex.main = 0.9)
for (v in c("BBmsy", "B")) {
  s <- tab_series[tab_series$variavel == v, ]
  plot(NA, xlim = range(s$ano), ylim = c(0, quantile(s$uci, 0.995, na.rm = TRUE) * 1.05),
       xlab = "Ano", ylab = if (v == "BBmsy") "B / Bmsy" else "Biomassa (t)",
       main = if (v == "BBmsy") "Biomassa relativa (B/Bmsy)" else "Biomassa absoluta (t)")
  for (nm in names(fits)) {
    d <- s[s$cenario == nm, ]; d <- d[order(d$ano), ]
    if (!nrow(d)) next
    polygon(c(d$ano, rev(d$ano)), c(d$lci, rev(d$uci)),
            col = adjustcolor(cor_cen_fig[nm], 0.10), border = NA)
  }
  for (nm in names(fits)) {
    d <- s[s$cenario == nm, ]; d <- d[order(d$ano), ]
    if (nrow(d)) lines(d$ano, d$mu, col = cor_cen_fig[nm], lwd = 2.6)
  }
  if (v == "BBmsy") abline(h = 1, col = "firebrick", lty = 2)
  legend("topright", names(fits), col = cor_cen_fig, lwd = 2.4, bty = "n", cex = 0.5)
}
par(op); dev.off()
cat("PNG salvo: jabba_trajetorias_biomassa.png\n")

## --- 12.6 resíduos de todos os cenários numa figura só -------------------
## O QUE ESTA FIGURA RESPONDE: o ajuste ruim é culpa da CONFIGURAÇÃO do
## modelo (índice escolhido, prior, forma da curva) ou é do DADO? Se cada
## cenário errasse em anos diferentes, seria configuração. Se todos erram
## nos MESMOS anos, na mesma direção, o problema está no índice — nenhuma
## reconfiguração vai resolver. As caixas cinza por ano mostram a dispersão
## ENTRE cenários; as linhas coloridas, cada cenário.
cat("\n--- 12.6 resíduos por índice + RMSE ---\n")
if (!is.null(tab_residuos) && nrow(tab_residuos)) {
  idxs <- unique(tab_residuos$indice)
  nc <- min(2, length(idxs)); nr <- ceiling(length(idxs) / nc)
  lay <- matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE)
  lay[lay > length(idxs)] <- 0
  lay <- rbind(lay, matrix(max(lay) + 1, nrow = 1, ncol = nc))
  
  png(file.path(DIR_FIG, "jabba_residuos_cenarios.png"),
      width = 13 * nc, height = 8 * nr + 9, res = 300, antialias = AA, units = "cm")
  op <- par(mar = c(4, 4.6, 2.8, 1), bty = "l", cex = 0.72)
  layout(lay, heights = c(rep(1, nr), 1.15))
  
  for (ix in idxs) {
    d <- tab_residuos[tab_residuos$indice == ix, ]
    anos_ix <- sort(unique(d$ano))
    plot(NA, xlim = range(anos_ix), ylim = range(d$residuo_log, na.rm = TRUE) * 1.1,
         xlab = "Ano", ylab = "Residuo (log)",
         main = sprintf("%s  (n = %d anos, %d cenarios)", ix, length(anos_ix),
                        length(unique(d$cenario))))
    ## dispersão entre cenários, ano a ano
    if (length(unique(d$cenario)) > 2)
      boxplot(residuo_log ~ ano, data = d, add = TRUE, at = anos_ix, boxwex = 0.55,
              col = adjustcolor("grey60", 0.25), border = "grey55", outline = FALSE,
              axes = FALSE)
    abline(h = 0, col = "firebrick", lty = 2, lwd = 1.4)
    for (nm in unique(d$cenario)) {
      dd <- d[d$cenario == nm, ]; dd <- dd[order(dd$ano), ]
      lines(dd$ano, dd$residuo_log, col = adjustcolor(cor_cen_fig[nm], 0.9), lwd = 1.5)
      points(dd$ano, dd$residuo_log, col = cor_cen_fig[nm], pch = 16, cex = 0.55)
    }
    ## sequência de resíduos do mesmo sinal = o modelo não acompanha a
    ## tendência do índice (é o que o teste de runs formaliza na seção 10)
    med_ano <- tapply(d$residuo_log, d$ano, median)
    lines(as.numeric(names(med_ano)), med_ano, col = "grey20", lwd = 2.4, lty = 1)
  }
  
  ## painel inferior: RMSE por cenário, com o SDNR anotado
  par(mar = c(4.2, 11, 3, 2))
  ta <- tab_ajuste[order(tab_ajuste$RMSE), ]
  bp <- barplot(ta$RMSE, horiz = TRUE, names.arg = ta$cenario, las = 1,
                col = cor_cen_fig[ta$cenario], border = NA, xlim = c(0, max(ta$RMSE) * 1.25),
                xlab = "RMSE (%)", main = "Erro de ajuste ao indice por cenario (RMSE) e SDNR",
                cex.names = 0.62)
  abline(v = 30, col = "forestgreen", lty = 2, lwd = 1.6)
  text(max(ta$RMSE) * 1.02, bp, sprintf("SDNR %.2f", ta$SDNR), cex = 0.6, adj = 0, xpd = NA)
  mtext("linha verde = 30%, referencia usual de bom ajuste | SDNR ~1 = o erro assumido bate com o observado",
        side = 3, line = 0.15, cex = 0.58, col = "grey30", adj = 0)
  par(op); layout(1); dev.off()
  cat("PNG salvo: jabba_residuos_cenarios.png\n")
} else {
  cat("[aviso] sem tabela de resíduos — figura de resíduos não gerada.\n")
}

## --- 12.7 tornado: o que mais move o resultado ---------------------------
## A LOGICA, e ela é diferente da do CMSY++: lá os dois eixos de cenário
## (método de r, método de bk) são colunas explícitas. Aqui os cenários
## diferem do caso base em campos diferentes da grade, então o fator é
## DEDUZIDO: para cada cenário, compara-se a linha dele com a do caso base
## e vê-se quais campos mudaram. Só entram no tornado os cenários que
## mudaram EXATAMENTE UM campo — que é a condição para o efeito ser
## atribuível àquele fator. Cenário que muda dois campos ao mesmo tempo
## (a família de janela curta muda janela, psi e K juntos) não é
## interpretável num tornado e sai, com aviso.
cat("\n--- 12.7 tornado de sensibilidade ---\n")

campos_fator <- c(indices = "Indice usado", modelo = "Forma da curva de producao",
                  r_min = "Prior de r", r_max = "Prior de r",
                  k_min = "Prior de K", k_max = "Prior de K",
                  janela_ini = "Janela do modelo", psi_mu = "Prior de psi",
                  psi_cv = "Prior de psi")

tornado <- function(base_id, metrica, metrica_nome, arquivo) {
  if (!base_id %in% names(fits)) { cat(sprintf("  [pulado] base %s não ajustou\n", base_id)); return(NULL) }
  cb <- cenarios[cenarios$id == base_id, ]
  val <- function(nm) {
    d <- draws_todos[[nm]]
    switch(metrica, MSY = median(d$MSY), BB0 = median(d$BB0),
           r = median(d$r), K = median(d$K),
           BBmsy = tab_status$BBmsy[tab_status$cenario == nm])
  }
  v_base <- val(base_id)
  linhas <- NULL
  for (nm in setdiff(names(fits), base_id)) {
    cc <- cenarios[cenarios$id == nm, ]
    difs <- names(campos_fator)[sapply(names(campos_fator), function(k) {
      a <- cb[[k]]; b <- cc[[k]]
      !isTRUE(all.equal(a, b)) && !(is.na(a) && is.na(b)) })]
    fatores <- unique(campos_fator[difs])
    if (length(fatores) != 1) next          # muda mais de um fator: não entra
    linhas <- rbind(linhas, data.frame(
      fator = fatores, nivel = nm, cenario_id = nm, valor = val(nm),
      delta_pct = 100 * (val(nm) - v_base) / v_base,
      row.names = NULL, stringsAsFactors = FALSE))
  }
  if (is.null(linhas)) { cat("  [pulado] nenhum cenário de fator único\n"); return(NULL) }
  
  ## empilha os níveis de cada fator dos dois lados do zero (convenção de
  ## leiaute, igual à do CMSY++: não representa soma de efeitos)
  empilhar <- function(df) {
    df <- df[order(abs(df$delta_pct)), ]
    neg <- df[df$delta_pct < 0, , drop = FALSE]; pos <- df[df$delta_pct >= 0, , drop = FALSE]
    if (nrow(neg)) { cum <- 0; for (i in seq_len(nrow(neg))) {
      neg$xmax[i] <- cum; cum <- cum + neg$delta_pct[i]; neg$xmin[i] <- cum } }
    if (nrow(pos)) { cum <- 0; for (i in seq_len(nrow(pos))) {
      pos$xmin[i] <- cum; cum <- cum + pos$delta_pct[i]; pos$xmax[i] <- cum } }
    rbind(neg, pos)
  }
  te <- do.call(rbind, lapply(split(linhas, linhas$fator), empilhar))
  amp <- sapply(split(te, te$fator), function(d) max(d$xmax) - min(d$xmin))
  ordem_f <- names(sort(amp, decreasing = TRUE))
  te$y <- match(te$fator, rev(ordem_f))
  
  cores_n <- setNames(grDevices::hcl.colors(length(unique(te$nivel)), palette = "Dynamic"),
                      unique(te$nivel))
  xl <- range(c(te$xmin, te$xmax, 0)) * 1.15
  
  png(file.path(DIR_FIG, arquivo), width = 30, height = 4 + 3.2 * length(ordem_f),
      res = 300, antialias = AA, units = "cm")
  op <- par(mar = c(4.6, 13, 5, 13), bty = "l")
  plot(NA, xlim = xl, ylim = c(0.5, length(ordem_f) + 0.5), yaxt = "n", ylab = "",
       xlab = sprintf("Variacao da mediana de %s em relacao ao cenario base (%%)", metrica_nome))
  mtext(sprintf("Tornado - sensibilidade de %s (JABBA)", metrica_nome),
        side = 3, line = 2.6, cex = 1.05, font = 2, adj = 0)
  mtext(sprintf("Base: %s | mediana de %s = %.3g", base_id, metrica_nome, v_base),
        side = 3, line = 1.1, cex = 0.8, adj = 0)
  abline(v = pretty(xl), col = "grey92"); abline(v = 0, col = "black", lwd = 1.4)
  for (i in seq_len(nrow(te))) rect(te$xmin[i], te$y[i] - 0.32, te$xmax[i], te$y[i] + 0.32,
                                    col = cores_n[te$nivel[i]], border = "white")
  axis(2, at = seq_along(ordem_f), labels = rev(ordem_f), las = 1, tick = FALSE, cex.axis = 0.82)
  legend(x = xl[2] * 1.08, y = length(ordem_f) + 0.5, xpd = NA, legend = names(cores_n),
         fill = cores_n, bty = "n", cex = 0.72, title = "Cenario alternativo", xjust = 0)
  par(op); dev.off()
  cat(sprintf("PNG salvo: %s\n", arquivo))
  cbind(base = base_id, metrica = metrica_nome, te[, c("fator", "nivel", "valor", "delta_pct")])
}

BASE_TORNADO       <- if ("J3_BASE_C2C5" %in% names(fits)) "J3_BASE_C2C5" else names(fits)[1]
BASE_TORNADO_CURTO <- if ("J13_C5curto_base" %in% names(fits)) "J13_C5curto_base" else NA

tab_tornado <- rbind(
  tornado(BASE_TORNADO, "MSY",   "MSY",     "jabba_tornado_msy.png"),
  tornado(BASE_TORNADO, "BBmsy", "B/Bmsy final", "jabba_tornado_bbmsy.png"),
  if (!is.na(BASE_TORNADO_CURTO))
    tornado(BASE_TORNADO_CURTO, "MSY", "MSY (janela curta)", "jabba_tornado_msy_curto.png"))

if (!is.null(tab_tornado)) {
  cat("\n--- variação em relação ao cenário base ---\n")
  print(transform(tab_tornado, valor = signif(valor, 4), delta_pct = round(delta_pct, 1)),
        row.names = FALSE)
  escreve_csv_utf8(tab_tornado, file.path(DIR_SAIDA, "jabba_tornado_sensibilidade.csv"))
  if (tem_writexl)
    writexl::write_xlsx(list(tornado = tab_tornado),
                        file.path(DIR_SAIDA, "jabba_tornado_sensibilidade.xlsx"))
}

cat(sprintf("\nFiguras da seção 12 em: %s\n", DIR_FIG))


## =========================================================================
## 13) COMO REPORTAR
## =========================================================================
cat("\n===== COMO REPORTAR =====\n")
cat("1. Apresentar os cenários como ALTERNATIVAS, nunca escolher um e\n")
cat("   esconder os outros. A pergunta que o leitor faz é 'a conclusão\n")
cat("   depende de qual cenário?' — a figura do painel responde isso.\n")
cat("2. O caso base (J3) usa índices com q SEPARADO por período. Explicar\n")
cat("   que isso é o que acomoda a mudança de alvo de 2015, e que os dois\n")
cat("   índices NAO se sobrepõem em nenhum ano — quem liga os períodos é a\n")
cat("   série de captura e a dinâmica, não os dados de abundância.\n")
cat("3. DECLARAR o que a série de captura cobre. Se ficou só a rede de\n")
cat("   cerco industrial, dizer isso e dizer o que sobrou de fora; K e MSY\n")
cat("   valem para as remoções que entraram, não para o estoque inteiro.\n")
cat("   Se usou a série externa (Luz & Vieira), dizer a resolução taxonômica\n")
cat("   dela e que a avaliação deixa de ser independente do CMSY/DBSRA.\n")
cat("4. Reportar a SENSIBILIDADE AOS PRIORS DE r (J7/J8) E DE K (J10/J11)\n")
cat("   junto com o caso base. Stobberup & Erzini (2006) mostraram, neste\n")
cat("   mesmo estoque, que a posterior de r saía praticamente igual à prior.\n")
cat("   O TESTE: comparar a mediana de cada posterior com a mediana do prior\n")
cat("   que aquele cenário usou (raiz de r_min*r_max, idem para K). Se a\n")
cat("   razão for ~1 em todos, o modelo não atualizou o parâmetro — MSY =\n")
cat("   rK/4 vira premissa, e o status com ela. Isso NAO invalida o\n")
cat("   trabalho: é um resultado, e é o argumento técnico para pedir ao IMar\n")
cat("   a série de desembarques completa.\n")
cat("4b. O quadrado J2/J9/J3/J5 separa duas coisas que estavam confundidas:\n")
cat("   o efeito de PADRONIZAR o índice recente e o de ANCORAR no histórico.\n")
cat("   Reportar as quatro células juntas, não só a que deu o caso base.\n")
cat("4c. J13-J17 (janela curta 2015-2025, só C5, psi vago) respondem uma\n")
cat("   pergunta diferente dos demais: 'o problema de ajuste é do período\n")
cat("   recente sozinho, ou aparece em qualquer recorte?'. NÃO reportar\n")
cat("   B/Bmsy ou status do estoque a partir deles — o psi vago (seção 5.8)\n")
cat("   torna o nível absoluto de biomassa não identificável; só a REACAO\n")
cat("   do ajuste às mudanças de prior de r/K é informativa aqui, e só é\n")
cat("   comparável entre J13-J17/J19, nunca contra J1-J12/J18 (janela e\n")
cat("   prior de psi diferentes tornam DIC e status incomparáveis).\n")
cat("4d. J18 (BASE trocando C5 por C6) e J19 (par de J13 trocando C5 por C6\n")
cat("   na janela curta) testam se o índice muda de conclusão quando o\n")
cat("   viés de tática é corrigido por SUBSETTING (C6: só esforço da tática\n")
cat("   mais associada à cavala) em vez de por PADRONIZAÇÃO (C5: GLM com a\n")
cat("   tática como covariável, toda viagem dentro). Ver a parte 03, seção\n")
cat("   'C6', para a ressalva: nenhuma tática desta frota é de fato\n")
cat("   DOMINADA pela cavala, então C6 mede 'esforço da tática que mais\n")
cat("   encontra cavala', não 'esforço dirigido' no sentido estrito.\n")
cat("5. Reportar convergência (Geweke/Heidelberger), ajuste (RMSE, SDNR) e\n")
cat("   diagnóstico (rho de Mohn, MASE, teste de runs). Um modelo com MASE\n")
cat("   acima de 1 não prevê o índice melhor que uma regra ingênua.\n")
cat("6. Reportar as probabilidades de Kobe, não só o ponto estimado: é a\n")
cat("   vantagem de ser bayesiano e é o que o manejo consegue usar.\n")
cat("7. Declarar as premissas de preenchimento: anos de captura interpolados\n")
cat("   e anos marcados frágeis na parte 03 (2013 sem esforço, 2018 de\n")
cat("   levantamento diferente) entraram na série mesmo assim.\n")
cat("8. DIC só compara cenários com os MESMOS dados de entrada. Serve para\n")
cat("   Schaefer contra Fox (J3 x J6); não serve para J1 contra J3.\n")
cat("8b. A COMPARACAO MAIS LIMPA DA GRADE E J1 x J2, e ela merece estar no\n")
cat("   texto: C1 e literalmente C2 e C3 renormalizados (a razao C1/C2 e\n")
cat("   C1/C3 e constante ano a ano), entao J1 e J2 ajustam OS MESMOS 35\n")
cat("   pontos e J1 (um q) esta ANINHADO em J2 (dois q). A diferenca de DIC\n")
cat("   entre eles e o teste formal de que a capturabilidade mudou em 2015 —\n")
cat("   e o que justifica separar os q em vez de emendar a serie.\n")
cat("9. Confirmar que os resultados finais NAO saíram de RODADA_RAPIDA.\n")
cat("10. REPORTAR O TESTE PRIORI x POSTERIORI (secao 12) como resultado, nao\n")
cat("   como diagnostico interno. jabba_priori_posteriori.csv traz, por\n")
cat("   cenario e parametro, a razao entre as medianas e a reducao do IC de\n")
cat("   95%. Razao perto de 1 COM reducao de IC perto de zero = parametro\n")
cat("   assumido, nao estimado — e ai o MSY que sai dele e premissa. E a\n")
cat("   versao quantitativa do achado de Stobberup & Erzini (2006) neste\n")
cat("   mesmo estoque.\n")
cat("11. CHECAR A BIMODALIDADE (jabba_bimodalidade_BK.csv) ANTES de citar\n")
cat("   qualquer mediana de B/K ou B/Bmsy. Quando a posterior tem um pico\n")
cat("   perto de zero e outro perto de um, a mediana cai no VALE: e o valor\n")
cat("   em que o modelo acredita MENOS. Nesses cenarios reporte as duas\n")
cat("   solucoes e a probabilidade de cada uma, nunca o ponto nem o IC.\n")
cat("12. Reportar a RAZAO ENTRE OS q (jabba_capturabilidade_q.csv) nos\n")
cat("   cenarios de dois indices: e a medida direta do salto de\n")
cat("   capturabilidade na troca de alvo e da suporte numerico ao item 8b.\n")
cat("13. A planilha jabba_cenarios_racional.csv/.xlsx traz, na mesma linha,\n")
cat("   a pergunta que cada cenario responde, com quem ele deve ser lido em\n")
cat("   par (coluna `contraste`) e o resultado que ele produziu. Cenario de\n")
cat("   avaliacao nao tem significado sozinho, so contra outro.\n")

cat("\n===== FIM =====\n")


# =====================================================================
# Gráfico da série de capturas de Decapterus macarellus (Cabo Verde)
# Versão revisada: destaca o ano de 2018 (corrigido por média móvel)
# =====================================================================
# Depende de: ct já lido e corrigido, ou seja, já ter rodado:
#
#   ct <- read.csv2("Desembarques cavala_1989-2023.csv", sep = ";", dec = ".")
#   names(ct) <- c("Year", "Catch")
#   ct <- ct[ct$Year <= 2023, ]
#   ct$Catch[ct$Year == 2018] <- mean(ct$Catch[ct$Year %in% c(2015, 2016, 2017, 2019, 2020, 2021)])
#
# library(ggplot2)  # (já deve estar carregado no seu script principal)

# ---------------------------------------------------------------------
# 1) marca, dentro dos dados, qual ponto foi corrigido -> isso permite
#    colorir/rotular o ano de 2018 de forma diferente dos demais no
#    gráfico, sem precisar de um data.frame separado
# ---------------------------------------------------------------------
ct$Status <- ifelse(ct$Year == 2018, "Corrigido (média móvel)", "Observado")
ct$Status <- factor(ct$Status, levels = c("Observado", "Corrigido (média móvel)"))

# valor original de 2018, só para exibir na anotação do gráfico
valor_original_2018 <- 2362
valor_corrigido_2018 <- ct$Catch[ct$Year == 2018]

# ---------------------------------------------------------------------
# 2) monta o gráfico
# ---------------------------------------------------------------------
p_ct <- ggplot(ct, aes(x = Year, y = Catch)) +

  # tendência suavizada (loess) ao fundo, em cinza -- serve de contexto,
  # sem competir visualmente com a série observada/corrigida
  geom_smooth(
    method = "loess", formula = y ~ x, se = TRUE,
    color = "#898781", fill = "#898781", alpha = 0.15, linewidth = 0.9) +

  # linha conectando os pontos da série -- fica sempre azul, contínua;
  # só os PONTOS mudam de cor (não a linha), pra não fragmentar o traço
  geom_line(color = "#2a78d6", linewidth = 1) +

  # pontos: azul = observado; laranja e maior = ano corrigido (2018)
  geom_point(aes(color = Status, size = Status), shape = 16) +
  scale_color_manual(
    values = c("Observado" = "#2a78d6", "Corrigido (média móvel)" = "#eb6834")) +
  scale_size_manual(
    values = c("Observado" = 1.8, "Corrigido (média móvel)" = 3.4), guide = "none") +

  # marca a possível quebra metodológica visível na série (queda abrupta
  # a partir de 2014) -- ajuste o ano ou remova se não fizer mais sentido
  geom_vline(xintercept = 2014, linetype = "dashed", linewidth = 0.6, color = "#c3c2b7") +
  annotate("text", x = 2014.3, y = max(ct$Catch, na.rm = TRUE) * 0.97,
           label = "possível mudança\nde metodologia", hjust = 0,
           size = 3.2, color = "#898781", lineheight = 0.9) +

  # anota a correção do outlier de 2018, mostrando o valor original
  # descartado ao lado do valor que entrou de fato na análise
  annotate("text", x = 2001, y = valor_corrigido_2018 + 250,
           label = paste0("2018: ", valor_original_2018, " t (observado) →\n",
                           round(valor_corrigido_2018), " t (média móvel 2015-17 e 2019-21)"),
           size = 3.1, color = "#eb6834", hjust = 0, lineheight = 0.9) +

  scale_x_continuous(breaks = seq(1990, 2023, 5)) +
  scale_y_continuous(labels = scales::comma) +

  labs(
    title = "Capturas de Decapterus macarellus em Cabo Verde",
    subtitle = "Série 1989–2023 (todas as artes) · ano de 2018 corrigido por média móvel",
    x = "Ano", y = "Captura (t)", color = NULL) +

  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "#52514e", size = 11, margin = margin(b = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#e1e0d9", linewidth = 0.4),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_blank(),
    plot.margin = margin(10, 15, 5, 5))

p_ct

# ---------------------------------------------------------------------
# 3) salva em alta resolução
# ---------------------------------------------------------------------
ggsave(
  filename = "capturas_Decapterus_macarellus.png",
  plot = p_ct,
  device = "png",
  units = "cm",
  width = 25,
  height = 15,
  dpi = 300,
  bg = "white")

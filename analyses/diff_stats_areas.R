library(PopGenBounds)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

plot(sapply(2:1000, function(K) mean(Fup(K,seq(0.01,1-0.01,0.01)))))

K_range <- 2:1000
F_area <- sapply(K_range, function(K)
  mean(Fup(K, seq(0.01, 1 - 0.01, 0.01)))
)

G_area <- sapply(K_range, function(K)
  mean(Gpup(K, seq(0.01, 1 - 0.01, 0.01)))
)

D_area <- sapply(K_range, function(K)
  mean(Dup(K, seq(0.01, 1 - 0.01, 0.01)))
)

df <- data.frame(
  K = K_range,
  F_area = F_area,
  G_area = G_area,
  D_area = D_area
)

df_long <- pivot_longer(
  df,
  cols = c(F_area, G_area, D_area),
  names_to = "curve",
  values_to = "value"
)
df_long$curve <- factor(df_long$curve,
                        levels = c("F_area", "G_area", "D_area"))

# K values where points should appear
K_points <- c(2, 3, 6, 40)

# Subset for points
df_points <- df_long %>%
  filter(K %in% K_points)

p <- ggplot(df_long, aes(x = K, y = value, color = curve)) +
  geom_line(size = 1) +
  geom_point(
    data = df_points,
    size = 2
  ) +
  scale_x_log10(breaks = c(2, 3, 6, 40, 200, 1000)) +
  ylim(0, 1) +
  labs(
    x = expression(italic(K)),
    y = "Area",
    color = expression(bold("Statistic"))
  ) +
  scale_color_manual(
    values = c(
      "F_area" = "#C93842",
      "G_area" = "#2678B3",
      "D_area" = "#65A856"
    ),
  labels = c(
    expression(italic(F[ST])),
    expression(italic(G*"'"[ST])),
    expression(italic(D))
    )
  )+
  theme_light()+
  theme(
    axis.title.x = element_text(face = "italic"),
    legend.title = element_text(face = "bold")
  )

p

filepath <- paste0("../plots/", "stats_areas.svg")

ggsave(filepath, p, width = 3, height = 2, dpi = 300)



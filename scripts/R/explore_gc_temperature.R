# ============================================================
# Exploratory analysis: GC content vs temperature in Archaea
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Read data
# ------------------------------------------------------------

traits <- read.table(
  "data/archaea_traits.nex",
  skip = 7,
  nrows = 33,
  header = FALSE
)

colnames(traits) <- c("species", "GC_content", "OGT")

# Convert GC from proportion to percentage
traits$GC_percent <- traits$GC_content * 100


# ------------------------------------------------------------
# 2. Pearson correlation
# ------------------------------------------------------------

correlation <- cor.test(
  traits$OGT,
  traits$GC_percent,
  method = "pearson"
)

r <- unname(correlation$estimate)
p <- correlation$p.value


# ------------------------------------------------------------
# 3. Linear regression
# ------------------------------------------------------------

model <- lm(GC_percent ~ OGT, data = traits)

r_squared <- summary(model)$r.squared


# ------------------------------------------------------------
# 4. Statistical label
# ------------------------------------------------------------

stats_label <- sprintf(
  "Pearson r = %.3f\nR² = %.3f\np < 2.2 × 10⁻¹⁶\nn = %d",
  r,
  r_squared,
  nrow(traits)
)


# ------------------------------------------------------------
# 5. Publication-style figure
# ------------------------------------------------------------

p1 <- ggplot(
  traits,
  aes(x = OGT, y = GC_percent)
) +

  # Individual species
  geom_point(
    size = 3.5,
    colour = "#2878B5",
    alpha = 0.85
  ) +

  # Linear regression + 95% confidence interval
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "#D9534F",
    fill = "#D9534F",
    alpha = 0.15,
    linewidth = 1.2
  ) +

  # Statistical results
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = stats_label,
    hjust = 1.15,
    vjust = -0.7,
    size = 4.5
  ) +

  labs(
    title = "GC content increases with optimal growth temperature",
    subtitle = "Archaea (n = 33)",
    x = "Optimal growth temperature (°C)",
    y = "GC content (%)"
  ) +

  theme_classic(base_size = 14) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    axis.title = element_text(
      face = "bold"
    )
  )


# ------------------------------------------------------------
# 6. Save figure
# ------------------------------------------------------------

ggsave(
  "figures/gc_vs_temperature.png",
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/gc_vs_temperature.pdf",
  plot = p1,
  width = 8,
  height = 6
)


# ------------------------------------------------------------
# 7. Print statistical results
# ------------------------------------------------------------

cat("\n============================\n")
cat("PEARSON CORRELATION\n")
cat("============================\n")

print(correlation)

cat("\n============================\n")
cat("LINEAR REGRESSION\n")
cat("============================\n")

print(summary(model))
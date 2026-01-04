# 确保安装了 ggplot2 包
library(ggplot2)
library(latex2exp)

draw_straight_line <- function(models, filename) {
  # 1. 定义数据生成函数
  x_vals <- seq(1, 5, length.out = 100)

  # 2. 构建绘图数据框
  df_list <- list()
  for (name in names(models)) {
    a <- models[[name]]$a
    b <- models[[name]]$b
    y_vals <- a * exp(b * x_vals)
    df_list[[name]] <- data.frame(
      x = x_vals,
      y = y_vals,
      Group = name
    )
  }
  plot_data <- do.call(rbind, df_list)

  # 3. 准备标签文本
  label_data <- data.frame(
    Group = names(models),
    # 调整标签位置，稍微错开
    x_pos = c(3.5, 3.5, 3.5),
    y_pos = c(
      models[[1]]$a * exp(models[[1]]$b * 3.5),
      models[[2]]$a * exp(models[[2]]$b * 3.5),
      models[[3]]$a * exp(models[[3]]$b * 3.5)
    ),
    label_text = sapply(names(models), function(n) {
      a <- models[[n]]$a
      b <- models[[n]]$b
      intercept <- round(log(a), 2)
      slope <- b
      return(paste0(n, ": ln(y) = ", intercept, " ", slope, "x"))
    })
  )

  color_values <- c(models[[1]]$color, models[[2]]$color, models[[3]]$color)
  names(color_values) <- names(models)
  # 4. 绘图
  ggfig <- ggplot(plot_data, aes(x = x, y = y, color = Group)) +
    geom_line(linewidth = 1.2) +

    # Y轴设置：对数坐标，范围从 0.01 到 0.5
    scale_y_continuous(
      trans = "log10",
      limits = c(10^(-2), 10^(-0.5)),
      breaks = c(10^(-2), 10^(-1.5), 10^(-1.0), 10^(-0.5)),
      labels = function(x) TeX(sprintf("$10^{%.1f}$", log10(x))),
      expand = c(0, 0)
    ) +

    # X轴设置
    scale_x_continuous(expand = c(0, 0), limits = c(1, 5.1)) +

    geom_text(
      data = label_data,
      aes(x = x_pos, y = y_pos, label = label_text),
      vjust = -0.5,
      hjust = 0,
      size = 4,
      show.legend = FALSE
    ) +

    # 设置颜色 (对应新的组名)
    scale_color_manual(
      values = color_values
    ) +

    labs(
      title = "Linearized Deletion Frequency Curves (New Variants)",
      subtitle = "Log-transformed y-axis",
      x = "Microhomology size (bp)",
      y = "Deletion frequency (%)"
    ) +

    # 样式：黑色坐标轴，无网格
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 12),
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top"
    )

  ggsave(filename, ggfig)
}

draw_straight_line(
  models = list(
    "PolQ-ΔHel" = list(a = 0.36, b = -0.61, color = "#56B4E9"), # 蓝色
    "PolQ-ΔPol" = list(a = 0.33, b = -0.60, color = "#E15759"), # 红色
    "ΔPolQ" = list(a = 0.36, b = -0.60, color = "#F28E2B") # 橙色
  ),
  filename = "result/PolQ.pdf"
)

draw_straight_line(
  models = list(
    "SpyCas9" = list(a = 0.29, b = -0.54, color = "#56B4E9"), # 蓝色
    "SpyMac" = list(a = 0.4, b = -0.81, color = "#E15759"), # 红色
    "iSpyMac" = list(a = 0.41, b = -0.74, color = "#F28E2B") # 橙色
  ),
  filename = "result/mix.pdf"
)

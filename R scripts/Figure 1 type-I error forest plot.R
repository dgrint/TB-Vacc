
## =============================================================================
## TB infection precision simulation study
## Forest plot of type-I error
## =============================================================================

library(pacman)

p_load(tidyverse,
       gridExtra,
       ggsci,
       esci)

# type-I error

theme_set(theme_bw())
theme_update(
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank()
)


data <- tibble(order = 18:1,
               prec = rep(c("A","B","C","D","E","F"), each = 3),
               scen = rep(c(0.02, 0.05, 0.08), times = 6),
               est = c(2.2, 2.5, 2.0,
                       97.8, 51.4, 26.0,
                       2.2, 2.6, 2.1,
                       96.8, 53.2, 27.8,
                       51.6, 15.7, 8.8,
                       100, 100, 94.8),
               mc_err = c(0.3, 0.3, 0.3,
                          0.1, 1.1, 0.9,
                          0.3, 0.4, 0.3,
                          0.4, 1.1, 1.0,
                          1.1, 0.8, 0.6,
                          0.07, 0.07, 0.5),
               group = c(1, 1, 1,
                         0, 0, 0,
                         1, 1, 1,
                         0, 0, 0,
                         0, 0, 0,
                         0, 0, 0)
               )

data$low <- data$est-(1.96*data$mc_err)
data$high <- data$est+(1.96*data$mc_err)
data$high <- if_else(data$high > 100, 100, data$high)


plot <- ggplot(data, aes(est, order)) + 
  geom_point(size=4, shape = 18, aes(color = factor(group)), show.legend = FALSE) +
  geom_errorbarh(aes(xmax = high, xmin = low, color = factor(group)), height = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = 2.5, linetype = "longdash") +
  geom_vline(xintercept = c(0, 100), linetype = "solid") +
  scale_x_continuous(breaks = seq(0, 100, 10), labels = seq(0, 100, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(1, 18), clip = "off", xlim = c(0, 100)) +
  annotate("text", x = 86, y = 19.5, label = "Non-inferior by default") +
  annotate("text", x = 10, y = 19.5, label = "Intended type-I error") +
  labs(x="Type-I error (%, Monte Carlo error bounds)", y="")

plot
plot1 <- plot + theme(plot.margin = unit(c(0.75,0.3,0.1,0), "cm")) + scale_color_nejm()
plot1


theme_set(theme_bw())
theme_update(
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
)
  

labs <- tibble(order = rep(18:1, times = 3),
               xord = rep(c(1, 1.3, 1.6), each = 18),
               lab = c("", "100%", "",
                        "", "100%", "",
                        "", "95%", "",
                        "", "95%", "",
                        "", "95%", "",
                        "", "64%", "",
                        "", "100%", "",
                        "", "95%", "",
                        "", "100%", "",
                        "", "95%", "",
                        "", "98%", "",
                        "", "85%", "",
                        "2%", "5%", "8%",
                        "2%", "5%", "8%",
                        "2%", "5%", "8%",
                        "2%", "5%", "8%",
                        "2%", "5%", "8%",
                        "2%", "5%", "8%"))
  
data_table <- ggplot(labs, aes(x = xord, y = order, label = format(lab, nsmall = 1))) +
  geom_text(size = 4) + theme_bw() +
  geom_hline(aes(yintercept=c(18.5))) + 
  annotate("text", x = 1.04, y = 19.5, label = "Sensitivity") +
  annotate("text", x = 1.3, y = 19.5, label = "Specificity") +
  annotate("text", x = 1.55, y = 19.5, label = "TB Risk") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(colour="white"),#element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_line(colour="white"),#element_blank(),
        plot.margin = unit(c(0.75,0.3,0.1,0), "cm")) +
  labs(x="",y="") +
  coord_cartesian(ylim = c(1, 18), clip = "off", xlim = c(1, 1.6))
  
data_table


fig1 <- grid.arrange(data_table, plot1, ncol=2, widths = c(1, 2))

ggsave("./Output/Figure 1.pdf",
       fig1,
       height = 5,
       width = 9,
       units = "in")


# type-II error

theme_set(theme_bw())
theme_update(
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank()
)


data <- data.frame(order = 18:1,
               prec = rep(c("A","B","C","D","E","F"), each = 3),
               scen = rep(c(0.02, 0.05, 0.08), times = 6),
               est = c(79.5, 79.2, 80.0,
                       100, 98.2, 95.0,
                       78.1, 76.8, 78.3,
                       100, 97.4, 93.9,
                       97.8, 90.6, 87.7,
                       100, 99.8, 99.3),
               mc_err = c(0.9, 0.9, 0.9,
                          0.07, 0.3, 0.5,
                          0.9, 0.9, 0.9,
                          0.07, 0.4, 0.5,
                          0.3, 0.6, 0.7,
                          0.07, 0.1, 0.2),
               group = c(1, 1, 1,
                         0, 0, 0,
                         0, 0, 1,
                         0, 0, 0,
                         0, 0, 0,
                         0, 0, 0)
)

data$low <- data$est - (qnorm(0.975) * data$mc_err)
data$high <- data$est + (qnorm(0.975) * data$mc_err)
data$high <- if_else(data$high > 100, 100, data$high)


plot <- ggplot(data, aes(est, order)) + 
  geom_point(size=4, shape = 18, aes(color = factor(group)), show.legend = FALSE) +
  geom_errorbarh(aes(xmax = high, xmin = low, color = factor(group)), height = 0.5, show.legend = FALSE) +
  geom_vline(xintercept = 80, linetype = "longdash") +
  geom_vline(xintercept = c(50, 100), linetype = "solid") +
  scale_x_continuous(breaks = seq(50, 100, 10), labels = seq(50, 100, 10), expand = c(0, 0)) +
  coord_cartesian(ylim = c(1, 18), clip = "off") +
  annotate("text", x = 100, y = 19.5, label = "Overwhelming power", hjust = 1) +
  annotate("text", x = 80, y = 19.5, label = "Intended type-II error", hjust = 1) +
  labs(x="Type-II error (%, Monte Carlo error bounds)", y="")

plot
plot1 <- plot + theme(plot.margin = unit(c(0.75,0.3,0.1,0), "cm")) + scale_color_nejm()
plot1


theme_set(theme_bw())
theme_update(
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
)


labs <- tibble(order = rep(18:1, times = 3),
               xord = rep(c(1, 1.3, 1.6), each = 18),
               lab = c("", "100%", "",
                       "", "100%", "",
                       "", "95%", "",
                       "", "95%", "",
                       "", "95%", "",
                       "", "64%", "",
                       "", "100%", "",
                       "", "95%", "",
                       "", "100%", "",
                       "", "95%", "",
                       "", "98%", "",
                       "", "85%", "",
                       "2%", "5%", "8%",
                       "2%", "5%", "8%",
                       "2%", "5%", "8%",
                       "2%", "5%", "8%",
                       "2%", "5%", "8%",
                       "2%", "5%", "8%"))

data_table <- ggplot(labs, aes(x = xord, y = order, label = format(lab, nsmall = 1))) +
  geom_text(size = 4) + theme_bw() +
  geom_hline(aes(yintercept=c(18.5))) + 
  annotate("text", x = 1.04, y = 19.5, label = "Sensitivity") +
  annotate("text", x = 1.3, y = 19.5, label = "Specificity") +
  annotate("text", x = 1.55, y = 19.5, label = "TB Risk") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none",
        panel.border = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(colour="white"),#element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_line(colour="white"),#element_blank(),
        plot.margin = unit(c(0.75,0.3,0.1,0), "cm")) +
  labs(x="",y="") +
  coord_cartesian(ylim = c(1, 18), clip = "off", xlim = c(1, 1.6))

data_table


fig1 <- grid.arrange(data_table, plot1, ncol=2, widths = c(1, 2))

ggsave("./Output/Figure 2.pdf",
       fig1,
       height = 5,
       width = 9,
       units = "in")

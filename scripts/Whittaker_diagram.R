

library(ggplot2)
library(plotbiomes)
library(gridExtra)
library(grid)

whittaker_plot <- whittaker_base_plot() +
  guides(color = guide_legend(title = NULL)) 
whittaker_plot

whittaker_plot <- whittaker_base_plot() +
  theme(legend.position = c(0.15, 0.84),
        legend.text = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 30, colour = 'black'),
        axis.title.x = element_text(size = 30, colour = 'black'),
        axis.text.x = element_text(size = 30, colour = 'black'),
        axis.text.y = element_text(size = 30, colour = 'black'))
whittaker_plot 

stat_clean$pre_mean_cm = stat_clean$pre_mean/10
stat_clean$tmp
names(stat_clean)
beta_diagram <- subset(stat_clean,pre_mean_cm>0 ) ## these points where pre_mean =0 are in fact NA
nrow(beta_diagram) # 4176
names(beta_diagram)
beta_diagram$Economy

myco_color_palette <- c("Scavenging" = "magenta", "Mining" = "blue", "NM" = "black")
myco_shape_palette <- c("Scavenging" = 16,  # Filled circle
                        "Mining" = 17,  # Filled triangle
                        "NM" = 21)  # Filled diamond

whittaker_beta <- whittaker_plot +
  geom_point(data = beta_diagram, aes(x = tmp, y = pre_mean_cm, color=Economy, shape=Economy),size = 2.5)+
  scale_color_manual(values = myco_color_palette)+
  scale_shape_manual(values = myco_shape_palette)+
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 1, override.aes = list(size = 5))) +
  theme(legend.position = c(0.15, 0.70),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.title.y = element_text(size = 30, colour = 'black'),
        axis.title.x = element_text(size = 30, colour = 'black'),
        axis.text.x = element_text(size = 30, colour = 'black'),
        axis.text.y = element_text(size = 30, colour = 'black'))+
  guides(fill = guide_legend(
    title= NULL,
    title.position = "top",
    override.aes = list(size = 8),
    keywidth = unit(1, "cm"),
    keyheight = unit(1, "cm"),
    #nrow = 20,
    byrow = TRUE))
whittaker_beta


library(dplyr)
library(ggplot2)
library(scales)
library(patchwork) 

rm(list = ls())
setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Objs_3y4")

sup <- read.csv("Superficies_MSPA.csv")

#Limpieza de datos
sup_clean <- sup %>%
  mutate(across(everything(), ~ as.numeric(gsub(",", "", .))))

#Modelos lineales del periodo
lm_core <- lm(Core_ha ~ year, data = sup_clean)
lm_core
r2_core <- summary(lm_core)$r.squared
r2_core
p_core <- summary(lm_core)$coefficients[2, 4]
p_core

p_core_str <- ifelse(p_core < 0.001, "p < 0.001", sprintf("p = %.3f", p_core))
label_core <- paste0("R² = ", round(r2_core, 2), "\n", p_core_str)

lm_con <- lm(Conectores_ha ~ year, data = sup_clean)
r2_con <- summary(lm_con)$r.squared
r2_con
p_con <- summary(lm_con)$coefficients[2, 4]
p_con

p_con_str <- ifelse(p_con < 0.001, "p < 0.001", sprintf("p = %.3f", p_con))
label_con <- paste0("R² = ", round(r2_con, 2), "\n", p_con_str)

tema_cuadrado <- theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1, 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10, color = "black"), # Activamos los títulos centrados y en negrita
    axis.text = element_text(color = "black") # Fuerza los números de los ejes a negro puro
  )

#Panel A
plot_A <- ggplot(sup_clean, aes(x = year, y = Core_ha)) +
  geom_line(color = "#44693e", linewidth = 1.2) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 0.8) +
  # Forzar el límite inferior a 0 y dejar un 10% de margen superior para las etiquetas
  scale_y_continuous(labels = label_number(big.mark = ""), limits = c(0, max(sup_clean$Core_ha) * 1.1)) + 
  labs(title = "Núcleos de hábitat natural", x = "Año", y = "Superficie (ha)") +
  annotate("text", x = Inf, y = Inf, label = label_core, 
           hjust = 1.1, vjust = 5, size = 4, color = "black") +
  tema_cuadrado

#Panel B
plot_B <- ggplot(sup_clean, aes(x = year, y = Conectores_ha)) +
  geom_line(color = "#800000", linewidth = 1.2) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, linewidth = 0.8) +
  # Forzar el límite inferior a 0 y dejar un 10% de margen superior para las etiquetas
  scale_y_continuous(labels = label_number(big.mark = ""), limits = c(0, max(sup_clean$Conectores_ha) * 1.1)) +
  labs(title = "Conectores estructurales o corredores", x = "Año", y = "Superficie (ha)") +
  annotate("text", x = Inf, y = Inf, label = label_con, 
           hjust = 1.1, vjust = 5, size = 4, color = "black") +
  tema_cuadrado

#Ensamblaje de figura
figura_final <- plot_A + plot_B +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

print(figura_final)


#ismej_theme.R
#author: "Mike Connelly"
#date: "03/22/2022"

# Create ggplot2 theme for journal submission
# ISME figure guidelines
# 1-column width: 85 mm, 2-column width: 175 mm
# All text should be sans-serif typeface, preferably Helvetica or Arial.
# Maximum text size is 7pt. Minimum text size is 5pt.
theme_ismej <- function(base_size = 10) {
  (theme_foundation(base_size=base_size)
   + theme(
     # 
     plot.background = element_rect(colour = NA),
     panel.background = element_blank(),
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(), 
     panel.border = element_rect(color = "black", fill = NA),
     # 
     plot.title = element_text(face = "plain", size = rel(1.2), hjust = 0),
     plot.subtitle = element_text(face = "plain", size = rel(0.8)),
     # 
     axis.title = element_text(face = "plain",size = rel(1)),
     axis.title.x = element_text(vjust = 0),
     axis.title.y = element_text(angle = 90, vjust = 2),
     text = element_text(),
     axis.text = element_text(size = rel(0.8)), 
     axis.text.x = element_text(angle = 0),
     axis.text.y = element_text(), 
     #axis.line = element_line(colour = "black"),
     axis.ticks = element_line(),
     # 
     legend.title = element_text(size = rel(1)),
     legend.text = element_text(size = rel(0.8), margin = margin(t = 2, b = 2, unit = "mm")),
     legend.key = element_rect(color = NA),
     legend.background = element_rect(fill = NA, colour = NA),
     legend.position = "right",
     legend.direction = "vertical",
     legend.spacing.x = unit(2, "mm"),
     legend.spacing.y = unit(0, "mm"),
     legend.key.size = unit(2, "mm"),
     # 
     plot.margin = unit(c(2,2,2,2), "mm"),
     # 
     strip.background = element_rect(color = "black", fill = "grey"),
     strip.text = element_text(size = rel(1), face = "plain")
   ))
}
###
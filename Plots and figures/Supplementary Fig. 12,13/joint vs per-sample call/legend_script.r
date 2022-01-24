library(ComplexHeatmap)


lgd = Legend(labels = c("Full coverage joint call", "Full coverage per sample call"), type = "lines", 
             legend_gp = gpar(col=c("red", "blue")), grid_width = unit(1, "cm"))

lgd2 = Legend(labels = c("Downsampled joint call", "Downsampled per sample call"), type = "lines", 
              legend_gp = gpar(col=c("red", "blue"), lty=2), grid_width = unit(1, "cm"))

pd = packLegend(lgd, lgd2, direction="vertical")


draw(pd, x = unit(4.2, "cm"), y = unit(3.6, "cm"), just = c("left"))

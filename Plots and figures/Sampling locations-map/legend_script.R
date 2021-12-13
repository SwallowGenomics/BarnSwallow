library(ComplexHeatmap)



lgd = Legend(labels = c("ds1 (N=5)", "ds2 (N=135)", "ds4 (N=533)", "ds5 (N=216)", "ds6 (N=161)", "ds2, ds3 (N=44)", "ds2, ds3, ds6 (N=179)", "ds2, ds6 (N=53)"), type = "points", pch = c(18,3,18,18,18,18,18,18), 
             legend_gp = gpar(col=c("blue", "darkorange", "black", "red", "green4", "purple", "#ed21ff", "turquoise1")), background="white", grid_width = unit(1, "cm"))

             
draw(lgd, x = unit(4.2, "cm"), y = unit(3.6, "cm"), just = c("left"))


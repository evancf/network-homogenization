# Need some other packages for this
ipak(c("gridExtra", "base2grob", "gplots"))


# First metanetworks with and without introduced interactions ------------------

node.size <- 0.25

# Native metanetwork
p1.labs <- as.data.frame(c("Americas", "Europe", "Africa", "Asia",
                           "New Zealand", "Hawaii", "Madagascar", "Australia"))
colnames(p1.labs) <- "reg"
p1.labs$x <- c(-0.75,-0.35,-0.375,-0.15,
               -0.4,-0.52,-0.16, -0.06)
p1.labs$y <- c(0.55,0.58,0.27,0.67,
               0,0.99,0.21,0.45)
p1.labs$size <- c(rep(3,4),rep(2.8,4))

gnet.nat$col <- ifelse(gnet.nat$vertex.names %in% colnames(net.all.nat),
                       rgb(90,180,172, maxColorValue = 255), rgb(216,179,101, maxColorValue = 255))

p1 <- ggplot() +
  geom_edges(data=gnet.nat,
             aes(x=-x, y=y, xend=-xend, yend=yend),
             color="grey50", curvature=0.1, size=0.1, alpha=1/2) +
  geom_nodes(data=gnet.nat,
             size = node.size,
             aes(x=-x, y=y),
             color=rep(c(plant.rgb(), animal.rgb()), times = dim(net.all.nat)),
             alpha=0.9) +
  theme_void() +
  theme(legend.position="none") + 
  geom_text(data = p1.labs, aes(x,y, label = reg), size = p1.labs$size) +
  annotate("text", x = -1, y = 0.09, label = "Native\ninteractions\nonly", hjust = 0, size = 3.5, fontface = "italic", lineheight = 0.9) +
  annotate("text", x = -1, y = 1, label = "a", hjust = 0, size = 4, fontface = "bold")
  




# All metanetwork
gnet.all$col <- ifelse(gnet.all$vertex.names %in% colnames(net.all),
                   rgb(90,180,172, maxColorValue = 255), rgb(216,179,101, maxColorValue = 255))

p2 <- ggplot() +
  geom_edges(data=gnet.all,
             aes(x=-x, y=y, xend=-xend, yend=yend),
             color="grey50", curvature=0.1, size=0.1, alpha=1/2) +
  geom_nodes(data=gnet.all,
             size = node.size,
             aes(x=-x, y=y),
             color=rep(c(plant.rgb(), animal.rgb()), times = dim(net.all)),
             alpha=0.9) +
  theme_void() +
  theme(legend.position="none") +
  annotate("text", x = -1, y = 0.09, label = "Including\nintroduced\ninteractions", hjust = 0, size = 3.5, fontface = "italic", lineheight = 0.9) +
  annotate("text", x = -1, y = 1, label = "b", hjust = 0, size = 4, fontface = "bold") +
  
  annotate("point", x = -0.16, y = 0.995, size = 2, color = plant.rgb()) + 
  annotate("text", x = -0.08, y = 0.996, label = "plants", size = 3) +
  annotate("point", x = -0.16, y = 0.945, size = 2, color = animal.rgb()) + 
  annotate("text", x = -0.066, y = 0.946, label = "animals", size = 3)







# Make a figure showing spatial distribution of networks ----------------------

metanet$reg.for.col <- ifelse(metanet$oceanic.island == "yes", 
                              "Oceanic Islands",
                              metanet$reg) %>% factor()
metanet$reg.for.col <- relevel(metanet$reg.for.col, ref = "Oceanic Islands")

set.seed(123) 
reg.cols <- brewer.pal(12, "Paired")[sample(1:length(levels(metanet$reg.for.col)),
                                            length(levels(metanet$reg.for.col)), 
                                            replace = F)]

cex.pt.map <- 0.5
cex.lab <- 0.8
cex.inner.text <- 0.65
cex.axis <- 0.65

p3.fun <- function(){
  
  
  blank.map(col = alpha(reg.cols[as.numeric(metanet$reg.for.col)], 0.5),
            cex = cex.pt.map,
            add.box = F)
  
  legend(20, -55, xpd = T,
         bty = "n",
         pch = 16,
         pt.cex = cex.pt.map * 1.2,
         cex = cex.axis,
         col = reg.cols,
         legend = levels(metanet$reg.for.col),
         ncol = 3,
         #x.intersp = 4,
         y.intersp = 1.2,
         xjust = 0.5,
         text.width = 110)
  
  text(x = -180, y = 87.5, "c", font = 2)
}
p3.null <- base2grob(plot.new)
p3 <- base2grob(p3.fun)



# A figure showing the distribution of links separating species ----------------
my.tck <- -0.01
my.tick.label.line <- 0.5


p4.fun <- function(){
  
  # op <- par()
  # par(pty = "s")
  set.seed(4)
  plot(close.dat$closeness.nat.ig,
       close.dat$closeness.all.ig,
       xlim = round(range(c(close.dat$closeness.nat.ig, close.dat$closeness.all.ig), na.rm=T), 2),
       ylim = round(range(c(close.dat$closeness.nat.ig, close.dat$closeness.all.ig), na.rm=T), 2),
       xlab = "",
       ylab = "",
       pch = 16,
       cex = 0.3,
       asp = 1,
       las = 1,
       frame = F,
       cex.lab = cex.lab,
       col = rgb(0.6,0.6,0.6),
       # col = ifelse(close.dat$node.type == "animal",
       #              animal.rgb(190),
       #              plant.rgb(190)),
       xaxt = "n",
       yaxt = "n"
  )
  
  
  mtext(side = 2, line = 1.5, "Closeness (observed)", cex = cex.lab)
  mtext(side = 1, line = 1, "Closeness (native only)", adj = 0.7, cex = cex.lab)
  lab1 <- c(0.1, 0.2, 0.3)
  axis(1, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, tck = my.tck)
  axis(1, at = lab1, lwd = 0, lwd.ticks = 0, line = -0.9, cex.axis = cex.axis)
  
  axis(2, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, las = 1, tck = my.tck)
  axis(2, at = lab1, lwd = 0, lwd.ticks = 0, line = -0.5, cex.axis = cex.axis, las = 1)
  
  
  
  curve(x*1, add = T, lty = 2, xpd = F)
  
  curve(coef(close.sma.mod)[1] + x * coef(close.sma.mod)[2], 
        add = T, lwd = 1, col = 1,
        from = min(close.dat$closeness.nat.ig, na.rm = T), xpd = F)
  
  x <- seq(min(close.dat$closeness.nat.ig, na.rm = T),
           max(close.dat$closeness.all.ig, na.rm = T), length.out = 100)
  y1 <- close.sma.mod$groupsummary$Int_lowCI[1] + x * close.sma.mod$groupsummary$Slope_lowCI[1]
  y2 <- close.sma.mod$groupsummary$Int_highCI[1] + x * close.sma.mod$groupsummary$Slope_highCI[1]
  
  xx <- c(x, rev(x))
  yy <- c(y1, rev(y2))
  yy <- ifelse(yy < 0, 0, yy)
  yy <- ifelse(yy > 1, 1, yy)
  
  polygon(xx, yy, col = rgb(0,0,0,0.3), border = F, xpd = F)
  
  #par(op)
  
  
  text(x = 0.1 - 0.2*.27, y = 0.3 + (0.2)*.2, "e", font = 2)
  
  text(x = 0.29, y = 0.25, "1:1\nline", font = 3, cex = cex.inner.text)
  
}

p4 <- base2grob(p4.fun)



# A figure showing how nodes are distributed within clusters -------------------

p5.fun <- function(){
  
  plot(NA,
       xlim = c(0,17),
       ylim = c(0,0.3),
       cex.lab = cex.lab,
       las = 1,
       xaxt = "n",
       yaxt = "n",
       #ann = F,
       xlab = "", #Degrees of separation
       ylab = "", # Proportion #"Portion of species pairs"
       frame = F)
  #axis(1, at = seq(0, 15, by = 5), cex.axis = cex.axis, tck = my.tck)
  #axis(2, at = seq(0, 0.3, by = 0.1), cex.axis = cex.axis, las = 1, tck = my.tck)
  
  
  mtext(side = 2, line = 1.5, "Proportion", cex = cex.lab)
  mtext(side = 1, line = 1, "Degrees of separation", adj = 0.7, cex = cex.lab)
  lab1 <- seq(0, 15, by = 5)
  lab2 <- seq(0, 0.3, by = 0.1)
  axis(1, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, tck = my.tck)
  axis(1, at = lab1, lwd = 0, lwd.ticks = 0, line = -0.9, cex.axis = cex.axis)
  
  axis(2, at = lab2, labels = rep("", length(lab1)), cex.axis = cex.axis, las = 1, tck = my.tck)
  axis(2, at = lab2, lwd = 0, lwd.ticks = 0, line = -0.5, cex.axis = cex.axis, las = 1)
  
  dens.bw <- 0.6
  nat.dens <- density(dist.net.nat.ig, bw = dens.bw, from = 0, to = net.nat.diam)
  polygon(nat.dens, col = native.rgb(75), border = native.rgb(190), lwd = 2)
  
  
  all.dens <- density(dist.net.ig, bw = dens.bw, from = 0, to = net.diam)
  polygon(all.dens, col = all.rgb(75), border = all.rgb(190), xpd = F, lwd = 2)
  
  text(9,.27, "Including \nintroduced \ninteractions", pos = 4, cex = cex.inner.text, font = 3)
  segments(x0 = 7.5, x1 = 9, y0 = 0.28, lwd = 2, col = all.rgb(190))
  
  text(9,.155, "Native \ninteractions \nonly", pos = 4, cex = cex.inner.text, font = 3)
  segments(x0 = 7.5, x1 = 9, y0 = 0.166, lwd = 2, col = native.rgb(190))
  
  
  text(x = 0-17*.27, y = 0.3*1.2, "d", font = 2)
  
  
}

p5 <- base2grob(p5.fun)




width.fig1 <- 7.25
height.fig1 <- 6

if(make.pdf){
  setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
  pdf(file = "Figure 1.pdf", width = width.fig1, height = height.fig1)
}


grid.arrange(
  p1,p2,p3.null,p5,p4,
  widths = c(1,1,1,1),
  layout_matrix = matrix(c(1,1,2,2,
                           1,1,2,2,
                           1,1,2,2,
                           3,3,4,5,
                           3,3,4,5), 
                         ncol=4, byrow = T))


vp <- grid::viewport(x=0.22,y=0.215, width = 0.5, height = 0.75)
grid::pushViewport(vp)
grid::grid.draw(p3)





if(make.pdf){
  dev.off()
}







if(make.pdf){
  setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
  pdf(file = "Extended Data Figure 1.pdf", width = 6, height = 4.5)
}



q.nat.boot.dens <- density(q.nat.boot)
q.all.boot.dens <- density(q.all.boot)
q.reduced.null.dens <- density(q.reduced.null)
q.biome.null.dens <- density(q.biome.null)

par(mfrow=c(1,1))
plot(NA, xlim = c(0,0.7), ylim = c(0,100),
     xlab = "Modularity",
     ylab = "Probability density",
     las = 1,
     frame = F)
polygon(q.nat.boot.dens, col = native.rgb(100), lwd = 2, border = F)
polygon(q.all.boot.dens, col = all.rgb(100), lwd = 2, border = F)
lines(q.reduced.null.dens, col = all.rgb(175), lty = 2, lwd = 2)
lines(q.biome.null.dens, col = rgb(0,0,0,0.7), lty = 2, lwd = 2)

legend(x = 0.2, y = 105, 
       legend = c("Native interactions only", "Including introduced interactions",
                  "Null: reduced", "Null: randomized by biome"),
       lty = c(1,1,2,2),
       col = c(native.rgb(100),
               all.rgb(100),
               all.rgb(175),
               rgb(0,0,0,0.7)),
       lwd = c(6,6,2,2),
       cex = 0.8,
       bty = "n",
       text.font = 3)

if(make.pdf){
  dev.off()
}


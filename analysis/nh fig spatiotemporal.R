


setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
pdf("Figure 4.pdf", width = 4.75, height = 3)


layout(mat = matrix(c(1,3,3,
                      2,3,3), nrow=2, byrow=T))

cex.names <- 0.8
cex.labels <- 0.8
cex.axis <- 0.9
cex.panel.lab <- 1.6

my.tck <- -0.05
my.tick.label.line <- 0.5

metanet.reg <- metanet.reg[order(metanet.reg$prop.intro.int),]

ocean.reg <- metanet.reg$reg[which(metanet.reg$oceanic.island=="yes")]
reg.inds <- which(!rownames(metanet.reg) %in% ocean.reg)
ocean.inds <- which(rownames(metanet.reg) %in% ocean.reg)

intervals <- c(1,1,10)
reg.x.vals <- cumsum(c(0,rep(intervals, length(reg.inds))))
ocean.x.vals <- cumsum(c(max(reg.x.vals)+12, rep(intervals, length(ocean.inds))))

reg.2 <- rep(1:c(length(reg.inds)+length(ocean.inds)), each = 3)




par(mar = c(4.05,5,0.5,0))
plot(NA,
     xlim = c(-60, 20),
     ylim = c(0, 0.4),
     las = 1,
     xlab = "", # "Year"
     xaxt = "n",
     yaxt = "n",
     ylab = "", # Proportion introduced
     frame = F)
mtext("Year", side = 1, line = 2, xpd = T, cex = cex.labels)

mtext("Proportion introduced", side = 2, line = 3, xpd = T, adj = 2.1, cex = cex.labels)

lab1 <- seq(-60,20, by = 40)
lab2 <- seq(0, 0.4, by = 0.1)
axis(1, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, tck = my.tck)
axis(1, at = lab1, labels = seq(1940, 2020, by = 40), lwd = 0, lwd.ticks = 0, line = -0.4, cex.axis = cex.axis)

axis(2, at = lab2, labels = rep("", length(lab2)), cex.axis = cex.axis, las = 1, tck = my.tck)
axis(2, at = lab2, lwd = 0, lwd.ticks = 0, line = -0.3, cex.axis = cex.axis, las = 1)


lines(year.pred$year2000$x, year.pred$year2000$predicted,
      lwd = 2, lend = "butt")
polygon(c(year.pred$year2000$x,rev(year.pred$year2000$x)),
        c(year.pred$year2000$conf.low, rev(year.pred$year2000$conf.high)),
        col = int.rgb(50),
        border = F)
text(label = "a", x = -52, y = 0.395, xpd = T, cex = cex.panel.lab, font = 2)



plot(NA,
     xlim = c(0, 1),
     ylim = c(0, 0.8),
     las = 1,
     xlab = "", #gHM index
     xaxt = "n",
     yaxt = "n",
     ylab = "", #Proportion introduced
     frame = F)
mtext("Human modification", side = 1, line = 2, xpd = T, cex = cex.labels)

lab1 <- seq(0,1, by = 0.5)
lab2 <- seq(0, 0.8, by = 0.2)
axis(1, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, tck = my.tck)
axis(1, at = lab1, lwd = 0, lwd.ticks = 0, line = -0.4, cex.axis = cex.axis)

axis(2, at = lab2, labels = rep("", length(lab2)), cex.axis = cex.axis, las = 1, tck = my.tck)
axis(2, at = lab2, lwd = 0, lwd.ticks = 0, line = -0.3, cex.axis = cex.axis, las = 1)

lines(ghm.pred$ghm$x, ghm.pred$ghm$predicted,
      lwd = 2, lend = "butt")
polygon(c(ghm.pred$ghm$x,rev(ghm.pred$ghm$x)),
        c(ghm.pred$ghm$conf.low, rev(ghm.pred$ghm$conf.high)),
        col = int.rgb(50),
        border = F)
text(label = "b", x = 0.1, y = 0.79, xpd = T, cex = cex.panel.lab, font = 2)



par(mar = c(3.4,4.5,1,0))

plot(NA, 
     xlim = c(range(c(reg.x.vals, ocean.x.vals))),
     ylim = c(0.001,1),
     las = 1,
     frame = F,
     xaxt = "n",
     #log = "y",
     xlab = "",
     ylab = "") #Proportion introduced


text(label = "c", x = 1, y = 1.05, xpd = T, cex = cex.panel.lab, font = 2)

x.counter <- 1

for(i in reg.inds){
  
  yy <- metanet.reg[i, c("prop.intro.animal", 
                         "prop.intro.plant", 
                         "prop.intro.int")]
  
  xx <- reg.x.vals[which(reg.2==x.counter)]
  points(xx, 
         yy,
         pch = ifelse(yy==0,1,16),
         col = c(animal.rgb(150),
                 plant.rgb(150),
                 rgb(0,0,0,0.6)),
         cex = 1.5)
  
  
  x.counter <- x.counter + 1
  if(x.counter %% 2 == 0){
    polygon(rep(c(min(xx) - 5,
                  max(xx + 5)),each = 2),
            #c(1.02,-0.02,-0.02,1.02),  
            c(1,0,0,1),
            col = rgb(0,0,0,0.05),
            border = NA)
  }
  
  text(max(xx)+5, -0.05, metanet.reg$reg[i], 
       srt = 40, xpd = T, pos = 2, cex = cex.names)
  
}

x.counter <- 1

abline(v = mean(c(max(reg.x.vals[-length(reg.x.vals)]), 
                  min(ocean.x.vals))), lty = 2)

text(mean(reg.x.vals[-length(reg.x.vals)]), 1.05,
     "Continental",
     xpd = T,
     font = 3)

text(mean(ocean.x.vals[-length(ocean.x.vals)]), 1.05,
     "Oceanic",
     xpd = T,
     font = 3)

for(i in ocean.inds){
  
  yy <- metanet.reg[i, c("prop.intro.animal", 
                         "prop.intro.plant", 
                         "prop.intro.int")]
  
  xx <- ocean.x.vals[which(reg.2==x.counter)]
  points(xx, 
         yy,
         pch = ifelse(yy==0,1,16),
         col = c(animal.rgb(150),
                 plant.rgb(150),
                 rgb(0,0,0,0.6)),
         cex = 1.5)
  
  
  x.counter <- x.counter + 1
  if(x.counter %% 2 == 0){
    polygon(rep(c(min(xx) - 5,
                  max(xx + 5)),each = 2),
            #c(1.02,-0.02,-0.02,1.02),  
            c(1,0,0,1),
            col = rgb(0,0,0,0.05),
            border = NA)
  }
  
  text(max(xx)+5, -0.05, metanet.reg$reg[i], 
       srt = 40, xpd = T, pos = 2, cex = cex.names)
  
  
}

legend(4,1, pch = 16, pt.cex = 1.5,
       col = c(animal.rgb(150),
               plant.rgb(150),
               rgb(0,0,0,0.6)),
       legend = c("Animals", "Plants", "Interactions"),
       bty = "n")




dev.off()




# # Figures showing results using different subsets of networks
#
# pdf("Figure 4a alternative.pdf", width = 3, height = 4)
# 
# plot(NA,
#      xlim = c(-40, 20),
#      ylim = c(0, 0.4),
#      las = 1,
#      xlab = "Year", # 
#      xaxt = "n",
#      yaxt = "n",
#      ylab = "Proportion introduced", 
#      frame = F)
# 
# lab1 <- seq(-40,20, by = 20)
# lab2 <- seq(0, 0.4, by = 0.1)
# axis(1, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, tck = my.tck)
# axis(1, at = lab1, labels = seq(1960, 2020, by = 20), lwd = 0, lwd.ticks = 0, line = -0.4, cex.axis = cex.axis)
# 
# axis(2, at = lab2, labels = rep("", length(lab2)), cex.axis = cex.axis, las = 1, tck = my.tck)
# axis(2, at = lab2, lwd = 0, lwd.ticks = 0, line = -0.3, cex.axis = cex.axis, las = 1)
# 
# 
# lines(year.predX1$year2000$x, year.predX1$year2000$predicted,
#       lwd = 2, lend = "butt")
# polygon(c(year.predX1$year2000$x,rev(year.predX1$year2000$x)),
#         c(year.predX1$year2000$conf.low, rev(year.predX1$year2000$conf.high)),
#         col = int.rgb(50),
#         border = F)
# dev.off()
# 
# 
# 
# 
# 
# 
# pdf("Figure 4a alternatives.pdf", width = 6, height = 4)
# 
# par(mfrow=c(1,2))
# 
# 
# 
# plot(NA,
#      xlim = c(0, 1),
#      ylim = c(0, 0.8),
#      las = 1,
#      xlab = "Human modification", #gHM index
#      xaxt = "n",
#      yaxt = "n",
#      ylab = "Proportion introduced", #
#      frame = F)
# 
# lab1 <- seq(0,1, by = 0.5)
# lab2 <- seq(0, 0.8, by = 0.2)
# axis(1, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, tck = my.tck)
# axis(1, at = lab1, lwd = 0, lwd.ticks = 0, line = -0.4, cex.axis = cex.axis)
# 
# axis(2, at = lab2, labels = rep("", length(lab2)), cex.axis = cex.axis, las = 1, tck = my.tck)
# axis(2, at = lab2, lwd = 0, lwd.ticks = 0, line = -0.3, cex.axis = cex.axis, las = 1)
# 
# lines(ghm.pred.2000$ghm$x, ghm.pred.2000$ghm$predicted,
#       lwd = 2, lend = "butt")
# polygon(c(ghm.pred.2000$ghm$x,rev(ghm.pred.2000$ghm$x)),
#         c(ghm.pred.2000$ghm$conf.low, rev(ghm.pred.2000$ghm$conf.high)),
#         col = int.rgb(50),
#         border = F)
# 
# text(label = "a", x = 0.1, y = 0.79, xpd = T, cex = 1, font = 2)
# 
# 
# 
# plot(NA,
#      xlim = c(0, 1),
#      ylim = c(0, 0.8),
#      las = 1,
#      xlab = "Human modification", #gHM index
#      xaxt = "n",
#      yaxt = "n",
#      ylab = "Proportion introduced", #
#      frame = F)
# 
# lab1 <- seq(0,1, by = 0.5)
# lab2 <- seq(0, 0.8, by = 0.2)
# axis(1, at = lab1, labels = rep("", length(lab1)), cex.axis = cex.axis, tck = my.tck)
# axis(1, at = lab1, lwd = 0, lwd.ticks = 0, line = -0.4, cex.axis = cex.axis)
# 
# axis(2, at = lab2, labels = rep("", length(lab2)), cex.axis = cex.axis, las = 1, tck = my.tck)
# axis(2, at = lab2, lwd = 0, lwd.ticks = 0, line = -0.3, cex.axis = cex.axis, las = 1)
# 
# 
# lines(ghm.pred.2010$ghm$x, ghm.pred.2010$ghm$predicted,
#       lwd = 2, lend = "butt")
# polygon(c(ghm.pred.2010$ghm$x,rev(ghm.pred.2010$ghm$x)),
#         c(ghm.pred.2010$ghm$conf.low, rev(ghm.pred.2010$ghm$conf.high)),
#         col = int.rgb(50),
#         border = F)
# 
# text(label = "b", x = 0.1, y = 0.79, xpd = T, cex = 1, font = 2)
# 
# dev.off()


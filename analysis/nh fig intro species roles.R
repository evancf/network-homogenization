# Figure 3 ---------------------------------------------------------------------

# Version with supressed intercept for plotting ease for plants
plant.intro.prop.mod1 <- glmer(cbind(n.int, n.nat) ~ plant.native.status - 1 + (1|net.id) + (1|reg),
                               family = "binomial",
                               data = subset(plant.intro.prop.df, net.id %in% net.ids.with.intro.animals))


# Version with supressed intercept for plotting ease for animals
animal.intro.prop.mod1 <- glmer(cbind(n.int, n.nat) ~ animal.native.status - 1 + (1|net.id) + (1|reg),
                                family = "binomial",
                                data = subset(animal.intro.prop.df, net.id %in% net.ids.with.intro.plants))


animal.ests <- summary(animal.intro.prop.mod1)$coefficients
animal.mean <- plogis(animal.ests[,1])
animal.lo <- plogis(animal.ests[,1] + animal.ests[,2])
animal.hi <- plogis(animal.ests[,1] - animal.ests[,2])

plant.ests <- summary(plant.intro.prop.mod1)$coefficients
plant.mean <- plogis(plant.ests[,1])
plant.lo <- plogis(plant.ests[,1] + plant.ests[,2])
plant.hi <- plogis(plant.ests[,1] - plant.ests[,2])

if(make.pdf){
  setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
  pdf(file = "Figure 3.pdf", width = 3.5, height = 2)
}

ylimz <- c(0,0.5)
xvals <- c(1,2,4,5)
cexx <- 0.8
tckxy <- -0.075
par(mar = c(2,3.5,1,0.25),
    mfrow = c(1,2))

plot(c(animal.mean,plant.mean)  ~ xvals,
     frame.plot = F,
     pch = 16,
     xaxt = "n",
     yaxt = "n",
     type = "h",
     xlab = "",
     ylab = "",
     lend = "butt",
     col = rep(c(native.rgb(100), all.rgb(100)),2),
     lwd = 15,
     las = 1,
     xlim = c(0.25,5.75),
     ylim = c(-0.01,0.5))#ylimz)
#mtext(side = 2, line = 2, "Proportion\nintroduced partners", cex = cexx)
mtext(side = 2, line = 1.8, "Prob. of interaction with\nintroduced partner", cex = cexx)
axis(2, at = c(0,0.1,0.2,0.3,0.4, 0.5), las = 1, cex.axis = cexx*0.9, hadj = 0.5, tck = tckxy)
arrows(x0 = xvals,
       y0 = c(animal.lo, plant.lo),
       y1 = c(animal.hi, plant.hi),
       angle = 90,
       length = 0.05,
       code = 3,
       xpd = T,
       col = rep(c(native.rgb(), all.rgb()),2))

yvals <- c(0.42,0.47)
xvals <- c(0.5,1.12)
segments(x0=xvals[1], x1 = xvals[2], y0 = yvals,
         col = c(all.rgb(100),native.rgb(100)), lwd = 3, lend = "butt")
text(x = max(xvals)-0.15, y = yvals, c("introduced", "native"), cex = cexx*0.79, font = 3, pos = 4)

text(-0.2, 0.5*1.1, "a", font=2, xpd = T)

text(1.5,-0.087, "Animals", xpd = T, cex = cexx * 0.9)
text(4.5,-0.087, "Plants", xpd = T, cex = cexx * 0.9)

mtext(side = 1, line = 1, "Trophic level ", cex = cexx)


par(mar = c(2,2.75,1,1))
set.seed(4)
plot(jitter(modules, 0.75)  ~ jitter(prop.intro.int, 100),
     data = metanet,
     frame.plot = F,
     pch = 16,
     cex = 1,
     xaxt = "n",
     yaxt = "n",
     xlab = "",
     ylab = "",
     col = int.rgb(75),
     las = 1,
     xpd = T,
     xlim = c(0,1),
     ylim = c(0.9,6))
mtext(side = 2, line = 1.5, "Number of modules", cex = cexx)
mtext(side = 1, line = 1, "Proportion introduced", cex = cexx)
axis(1, at = c(0, 0.5, 1), las = 1, cex.axis = cexx*0.9, padj = -1.9, tck = tckxy)
axis(2, at = c(1:6), las = 1, cex.axis = cexx*0.9, hadj = 0.3, tck = tckxy)

text(-0.075, 6*1.092, "b", font=2, xpd = T)


if(make.pdf){
  dev.off()
}





# Extended Data Figure 3 --------------------------------------------------------------------


if(make.pdf){
  setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
  pdf(file = "Extended Data Figure 3.pdf", width = 3, height = 4)
}

par(mfcol=c(2,2), mar = c(4,4,0.5,0))

plant.sp.pred <- get_model_data(plant.sp.mod, type = "pred")
animal.sp.pred <- get_model_data(animal.sp.mod, type = "pred")


ylimz <- c(0,0.6)

plot(predicted ~ x, data = plant.sp.pred[[1]],
     frame.plot = F,
     pch = 16,
     xaxt = "n",
     yaxt = "n",
     type = "h",
     xlab = "",
     ylab = "",
     lend = "butt",
     col = c(native.rgb(100), all.rgb(100)),
     lwd = 20,
     las = 1,
     xlim = c(0.5,2.5),
     ylim = ylimz)
axis(2, at = c(0,0.2,0.4,0.6), las = 1)
arrows(x0 = plant.sp.pred[[1]]$x,
       y0 = plant.sp.pred[[1]]$conf.low,
       y1 = plant.sp.pred[[1]]$conf.high,
       angle = 90,
       length = 0.08,
       code = 3,
       col = c(native.rgb(), all.rgb()))

text(2.5,0.565, "plants", pos = 2, 
     font = 3, xpd = T)

text(c(1.35,2.35), -0.05, srt = 45,
     labels = c("Native", "Introduced"), 
     pos = 2, xpd = T)

text(label = "a", x = 0.75, y = 0.57, xpd = T, font = 2, cex = 1.3)



plot(predicted ~ x, data = animal.sp.pred[[1]],
     frame.plot = F,
     pch = 16,
     xaxt = "n",
     yaxt = "n",
     type = "h",
     xlab = "",
     ylab = "",
     lend = "butt",
     col = c(native.rgb(100), all.rgb(100)),
     lwd = 20,
     las = 1,
     xlim = c(0.5,2.5),
     ylim = ylimz)
axis(2, at = c(0,0.2,0.4,0.6), las = 1)
arrows(x0 = animal.sp.pred[[1]]$x,
       y0 = animal.sp.pred[[1]]$conf.low,
       y1 = animal.sp.pred[[1]]$conf.high,
       angle = 90,
       length = 0.08,
       code = 3,
       col = c(native.rgb(), all.rgb()))

text(2.5,0.565, "animals", pos = 2, 
     font = 3, xpd = T)

mtext("Normalized degree in native range",
      side = 2,
      line = 3, xpd = T, adj = -0.28, cex = 1.00)

text(c(1.35,2.35), -0.05, srt = 45,
     labels = c("Native", "Introduced"), 
     pos = 2, xpd = T)

text(label = "b", x = 0.75, y = 0.57, xpd = T, font = 2, cex = 1.3)






plant.sp.db.pred <- get_model_data(plant.sp.db.mod, type = "pred")
animal.sp.db.pred <- get_model_data(animal.sp.db.mod, type = "pred")


ylimz <- c(0,0.6)

plot(predicted ~ x, data = plant.sp.db.pred[[1]],
     frame.plot = F,
     pch = 16,
     xaxt = "n",
     yaxt = "n",
     type = "h",
     xlab = "",
     ylab = "",
     lend = "butt",
     col = c(native.rgb(100), all.rgb(100)),
     lwd = 20,
     las = 1,
     xlim = c(0.5,2.5),
     ylim = ylimz)
axis(2, at = c(0,0.2,0.4,0.6), las = 1)
arrows(x0 = plant.sp.db.pred[[1]]$x,
       y0 = plant.sp.db.pred[[1]]$conf.low,
       y1 = plant.sp.db.pred[[1]]$conf.high,
       angle = 90,
       length = 0.08,
       code = 3,
       col = c(native.rgb(), all.rgb()))

text(2.5,0.565, "plants", pos = 2, 
     font = 3, xpd = T)

text(c(1.35,2.35), -0.05, srt = 45,
     labels = c("Native", "Introduced"), 
     pos = 2, xpd = T)

text(label = "c", x = 0.75, y = 0.57, xpd = T, font = 2, cex = 1.3)


plot(predicted ~ x, data = animal.sp.db.pred[[1]],
     frame.plot = F,
     pch = 16,
     xaxt = "n",
     yaxt = "n",
     type = "h",
     xlab = "",
     ylab = "",
     lend = "butt",
     col = c(native.rgb(100), all.rgb(100)),
     lwd = 20,
     las = 1,
     xlim = c(0.5,2.5),
     ylim = ylimz)
axis(2, at = c(0,0.2,0.4,0.6), las = 1)
arrows(x0 = animal.sp.db.pred[[1]]$x,
       y0 = animal.sp.db.pred[[1]]$conf.low,
       y1 = animal.sp.db.pred[[1]]$conf.high,
       angle = 90,
       length = 0.08,
       code = 3,
       col = c(native.rgb(), all.rgb()))

text(2.5,0.565, "animals", pos = 2, 
     font = 3, xpd = T)


text(c(1.35,2.35), -0.05, srt = 45,
     labels = c("Native", "Introduced"), 
     pos = 2, xpd = T)

text(label = "d", x = 0.75, y = 0.57, xpd = T, font = 2, cex = 1.3)


if(make.pdf){
  dev.off()
}










if(make.pdf){
  setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
  pdf(file = "Extended Data Figure 4.pdf", width = 4, height = 6)
}

set.seed(4)
plot(c(Modularity - Modularity.null) ~ jitter(prop.intro.int,50),
     data = metanet,
     pch = 16,
     col = int.rgb(75),
     xpd = T,
     xlim = c(0, 1),
     #ylim = c(-0.08, 0.15),
     ylim = c(-0.1, 0.25),
     las = 1,
     xlab = "Proportion introduced interactions",
     ylab = "Null model-corrected modularity",
     frame = F)


xx <- c(0,1)
polygon(c(0,1,1,0),
        c(coef(fit.lqmm.mod)[1,1],
          coef(fit.lqmm.mod)[1,1] + coef(fit.lqmm.mod)[2,1],
          coef(fit.lqmm.mod)[1,3] + coef(fit.lqmm.mod)[2,3],
          coef(fit.lqmm.mod)[1,3]),
        col = int.rgb(50),
        border = F)
curve(coef(fit.lqmm.mod)[1,2]  + coef(fit.lqmm.mod)[2,2]*(x), add = TRUE,
      lwd = 2, lend = "butt")

text(0.8, 0.07, "95th", cex = 0.9, font = 3, pos = 4)
text(0.8, -0.0, "50th", cex = 0.9, font = 3)
text(0.8, -0.07, "5th", cex = 0.9, font = 3, pos = 2)


if(make.pdf){
  dev.off()
}


all.plant.mat.reg <- all.plant.mat.reg.num / all.plant.mat.reg.denom
nat.plant.mat.reg <- nat.plant.mat.reg.num / nat.plant.mat.reg.denom

all.animal.mat.reg <- all.animal.mat.reg.num / all.animal.mat.reg.denom
nat.animal.mat.reg <- nat.animal.mat.reg.num / nat.animal.mat.reg.denom

all.int.mat.reg <- all.int.mat.reg.num / all.int.mat.reg.denom
nat.int.mat.reg <- nat.int.mat.reg.num / nat.int.mat.reg.denom



if(make.pdf){
  setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
  pdf("Figure 2.pdf", 
      width = 7.25, height = 5,
      useDingbats = F)
}

par(mfrow=c(3,2), mar = c(0,0,0,0), oma = c(0,2,2,0), xpd = T)

gc.text.cex <- 0.9



shared.lines <- 0.5

blank.map(meta = metanet.reg, cex = log10(metanet.reg$rich.plant * (1 - metanet.reg$prop.intro.plant)))
add.leg.vals(col = plant.rgb(200))
mtext("Native species only", side = 3, line = 1, cex = gc.text.cex, font = 3)
mtext("Shared plant species", side = 2, line = shared.lines, cex = gc.text.cex, font = 3)
reg.gci(nat.plant.mat.reg, col = plant.rgb(200))


blank.map(meta = metanet.reg, cex = log10(metanet.reg$rich.plant))
mtext("Including introduced species", side = 3, line = 1, cex = gc.text.cex, font = 3)
reg.gci(all.plant.mat.reg, col = plant.rgb(200))

blank.map(meta = metanet.reg, cex = log10(metanet.reg$rich.animal * (1 - metanet.reg$prop.intro.animal)))
add.leg.vals(col = animal.rgb(200))
mtext("Shared animal species", side = 2, line = shared.lines, cex = gc.text.cex, font = 3)
reg.gci(nat.animal.mat.reg, col = animal.rgb(200))
blank.map(meta = metanet.reg, cex = log10(metanet.reg$rich.animal))
reg.gci(all.animal.mat.reg, col = animal.rgb(200))

blank.map(meta = metanet.reg, cex = log10(metanet.reg$rich.int * (1 - metanet.reg$prop.intro.int)))
add.leg.vals(lwd.factor = 10^3, col = rgb(0,0,0,0.8), leg.vals = c(0.0001, 0.001, 0.01),
             dot.vals = c(10,200,4000), custom.decimals = T)
mtext("Shared interactions", side = 2, line = shared.lines, cex = gc.text.cex, font = 3)
reg.gci(nat.int.mat.reg, lwd.factor = 10^3, col = rgb(0,0,0,0.8))
blank.map(meta = metanet.reg, cex = log10(metanet.reg$rich.int))
reg.gci(all.int.mat.reg, lwd.factor = 10^3, col = rgb(0,0,0,0.8))

if(make.pdf){
  dev.off()
}


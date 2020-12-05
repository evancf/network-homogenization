# Need some other packages for this
ipak(c("gridExtra", "base2grob", "gplots"))


# Make the plots

a.groups <- c("birds", "bats", "primates", "mamm.carns", "mamm.herbs")

p1 <- ggplot(data = subset(beta.nat, animal.group %in% a.groups), 
             aes(x = dist, y = (1-S))) +
  ylim(0,1) +
  theme_classic() +
  geom_point(color = native.rgb(10)) +
  geom_smooth(method = "loess", 
              color = int.rgb(200), size = 1.5) +
  labs(y = expression(paste("1 - ", beta["S"])), 
       x = "Distance (km)") +
  annotate("text", x = 19000, y = 1, label = "a", size = 4, fontface = "bold")

p2 <- ggplot(data = subset(beta.nat, animal.group %in% a.groups), 
             aes(x = dist, y = (1-WN))) +
  ylim(0,1) +
  theme_classic() +
  geom_point(color = native.rgb(10)) +
  geom_smooth(method = "loess", 
              color = int.rgb(200), size = 1.5) +
  labs(y = expression(paste("1 - ", beta["WN"])), 
       x = "Distance (km)") +
  annotate("text", x = 19000, y = 1, label = "b", size = 4, fontface = "bold")

p3 <- ggplot(data = subset(beta.all, animal.group %in% a.groups), 
             aes(x = dist, y = (1-S))) +
  ylim(0,1) +
  theme_classic() +
  geom_point(color = all.rgb(10)) +
  geom_smooth(method = "loess", 
              color = int.rgb(200), size = 1.5) +
  labs(y = expression(paste("1 - ", beta["S"])), 
       x = "Distance (km)") +
  annotate("text", x = 19000, y = 1, label = "c", size = 4, fontface = "bold")

p4 <- ggplot(data = subset(beta.all, animal.group %in% a.groups), 
             aes(x = dist, y = (1-WN))) +
  ylim(0,1) +
  theme_classic() +
  geom_point(color = all.rgb(10)) +
  geom_smooth(method = "loess", 
              color = int.rgb(200), size = 1.5) +
  labs(y = expression(paste("1 - ", beta["WN"])), 
       x = "Distance (km)") +
  annotate("text", x = 19000, y = 1, label = "d", size = 4, fontface = "bold")


setwd(paste(top.wd, "analysis", "homogenization figures", sep = "/"))
pdf("Extended Data Figure 2.pdf", 
    width = 6, height = 6,
    useDingbats = F)

grid.arrange(p1,p2,p3,p4,
             layout_matrix = matrix(1:4, nrow=2, byrow=T))

dev.off()
# Load data and packages -------------------------------------------------------

top.wd <- getwd()

load(file = "./analysis/homogenization.RData")

packages <- c("tidyverse", "plyr",
              "ggmap", "maps", "mapdata", "leaflet", 
              "RColorBrewer",
              "sp", "vegan",
              "bipartite", "igraph", "network", "ggnetwork", "ggrepel", 
              "geosphere","sf","rnaturalearth",
              "smatr", "lme4", "lqmm", "sjPlot", 
              "betalink")
ipak(packages) # This will download any not yet downloaded packages


# A 'long' version of the dataset (1 row per local network by species 
# combination) is called net.long0

# net.long is the subset of this with just interactions among species (not unidentified taxa)

# Get a version of net.long with only native species, call it net.long.nat
net.long.nat <- subset(net.long, plant.native.status == "native" & animal.native.status == "native")


# Run Extended Data Figure scripts and print to PDF? ---------------------------

run.figures <- T
make.pdf <- T



# By the numbers ---------------------------------------------------------------

# Total unique pairwise interactions
total.int <- sum(net.long$value)
total.int

# Total unique pairwise interactions only involving native species
total.int.nat <- sum(net.long.nat$value)
total.int.nat

# Percent of interactions with an introduced species
1 - (total.int.nat / total.int)

# Number of spatially/temporally distinct networks
l.unique(net.long$net.id)

# Number of studies
l.unique(net.long$study.id)

# Number of taxa
tax.counts <- with(net.long, 
                   data.frame("taxa" = c(l.unique(animal.id), l.unique(plant.id)),
                              "accepted.species" = c(l.unique(animal.accepted.species), l.unique(plant.accepted.species)),
                              "genera" = c(l.unique(animal.genus), l.unique(plant.genus)),
                              "families" = c(l.unique(animal.family), l.unique(plant.family))))
rownames(tax.counts) <- c("animals", "plants")
tax.counts

# 
table(net.long[which(net.long$value == 1), c("plant.native.status", "animal.native.status")])


# Global metanetwork -----------------------------------------------------------

# Make global adjacency matrices
net.all <- make.global.bin.net(net.long)
net.all.nat <- make.global.bin.net(net.long.nat)





# See how the network fits modular structure imposed by bioregion --------------

# Set up for the first "reduced" null model to account for the lower number
# of interactions in the "native" meta-network
n.int.all <- sum(unlist(net.all)) # Total interactions incl. introduced
n.int.nat <- sum(unlist(net.all.nat)) # Total interaction w/o introduced

n.row.all <- nrow(net.all)
n.col.all <- ncol(net.all)
n.all <- n.row.all * n.col.all

n.row.nat <- nrow(net.all.nat)
n.col.nat <- ncol(net.all.nat)
n.nat <- n.row.nat * n.col.nat

p.int.diff <- 1 - (n.int.nat / n.int.all)
p.int.nat <- n.int.nat / n.nat





# Set up to output modularity values

q.nat.boot <- NA
q.all.boot <- NA
q.reduced.null <- NA
q.biome.null <- NA

# Need to have biome on net.long for this
net.long <- left_join(net.long, metanet[,c("net.id", "biome")])

n.boot <- 200 #
set.seed(4)
for(i in 1:n.boot){
  
  nl <- net.long[sample(1:nrow(net.long), replace = T), ]
  nl.nat <- subset(nl, plant.native.status == "native" & animal.native.status == "native")
  
  nl.net.boot <- make.global.bin.net(nl)
  q.all.boot[i] <- bioreg.mod(ig = graph_from_incidence_matrix(nl.net.boot),
                              nat.df.a = a.native.reg, 
                              nat.df.p = p.native.reg)
  
  nl.net.nat.boot <- make.global.bin.net(nl.nat)
  q.nat.boot[i] <- bioreg.mod(ig = graph_from_incidence_matrix(nl.net.nat.boot),
                              nat.df.a = a.native.reg, 
                              nat.df.p = p.native.reg)
  
  
  # The "reduced" null
  
  # Want to randomly remove interactions so that there are on average the
  # number of interactions in n.int.nat, so this first diff.mat will have 
  # mostly 1s and some zeros. By multiplying, the zeros turn existing 1s to 0s
  nl.net.dim <- dim(nl.net.boot)
  diff.mat <- (1 - rbinom(prod(nl.net.dim),
                          size = 1,
                          prob = p.int.diff)) %>%
    matrix(nrow = nl.net.dim[1], 
           ncol = nl.net.dim[2])
  reduced.mat <- diff.mat * nl.net.boot
  reduced.mat <- reduced.mat[which(rowSums(reduced.mat) > 0),
                             which(colSums(reduced.mat) > 0)]
  
  q.reduced.null[i] <- bioreg.mod(ig = graph_from_incidence_matrix(reduced.mat),
                                  nat.df.a = a.native.reg, 
                                  nat.df.p = p.native.reg)
  
  
  # The "randomized biome" null
  
  nl.intx <- subset(nl, value > 0)
  
  # Use this loop to randomize species identities with respect to interactions
  for(j in 1:length(unique(nl.intx$biome))){
    bio.inds <- which(nl.intx$biome == unique(nl.intx$biome)[j])
    
    # randomize species
    
    nl.intx[bio.inds,
            "plant.accepted.species"] <- nl.intx[sample(bio.inds),
                                                 "plant.accepted.species"]
    nl.intx[bio.inds,
            "animal.accepted.species"] <- nl.intx[sample(bio.inds),
                                                  "animal.accepted.species"]
  }
  
  biome.mat <- make.global.bin.net(nl.intx)
  q.biome.null[i] <- bioreg.mod(ig = graph_from_incidence_matrix(biome.mat),
                                nat.df.a = a.native.reg, 
                                nat.df.p = p.native.reg)
  
  print(i)
  
}

mean(q.nat.boot - q.all.boot) / mean(q.nat.boot - q.biome.null)

z.to.p <- function(z, sided = 2){
  sided * pnorm(-abs(z))
}

z.to.p(mean(q.nat.boot-q.all.boot) / sd(q.nat.boot-q.all.boot))

z.to.p(mean(q.all.boot-q.reduced.null) / sd(q.all.boot-q.reduced.null))

z.to.p(mean(q.nat.boot-q.reduced.null) / sd(q.nat.boot-q.reduced.null))



# Diameter, degrees of separation, proportion unconnected

net.ig <- make.global.bin.net(net.long) %>% graph_from_incidence_matrix()
net.nat.ig <- make.global.bin.net(net.long.nat) %>% graph_from_incidence_matrix()

net.diam <- diameter(net.ig, unconnected = T)
net.nat.diam <- diameter(net.nat.ig, unconnected = T)

dist.net.ig <- distances(net.ig)
mean_distance(net.ig)

dist.net.nat.ig <- distances(net.nat.ig)
mean_distance(net.nat.ig)

net.csize <- components(net.ig)$csize
net.nat.csize <- components(net.nat.ig)$csize

net.unconnected <- sum(net.csize[2:length(net.csize)]) / net.csize[1]
net.nat.unconnected <- sum(net.nat.csize[2:length(net.nat.csize)]) / net.nat.csize[1]



# Looking at centrality

net.ig.main <- induced_subgraph(
  net.ig, V(net.ig)[components(net.ig)$membership == which.max(components(net.ig)$csize)]
)

net.nat.ig.main <- induced_subgraph(
  net.nat.ig, V(net.nat.ig)[components(net.nat.ig)$membership == which.max(components(net.nat.ig)$csize)]
)


closeness.all.ig <- closeness(net.ig.main, normalized = T)
closeness.nat.ig <- closeness(net.nat.ig.main, normalized = T)


# Prep a dataframe for this analysis
intro.spp <- c(subset(net.long, animal.native.status == "introduced")$animal.accepted.species,
               subset(net.long, plant.native.status == "introduced")$plant.accepted.species) %>% unique()

close.dat <- as.data.frame(closeness.all.ig)
close.dat.nat <- as.data.frame(closeness.nat.ig)
close.dat.nat$spp <- rownames(close.dat.nat)
close.dat$spp <- rownames(close.dat)

close.dat <- left_join(close.dat, close.dat.nat)

close.dat$node.type <- ifelse(close.dat$spp %in% net.long$animal.accepted.species,
                              "animal", "plant")
close.dat$intro <- ifelse(close.dat$spp %in% intro.spp,
                          "yes", "no")


# Do the closeness analysis
close.sma.mod <- sma(closeness.all.ig ~ closeness.nat.ig, data = close.dat, slope = 1)
summary(close.sma.mod)





# Get network geometries for plots
set.seed(4)
gnet.nat <- ggnetwork(net.all.nat, layout="fruchtermanreingold",
                      arrow.gap=0, cell.jitter=0)

set.seed(4)
gnet.all <- ggnetwork(net.all, layout="fruchtermanreingold",
                      arrow.gap=0, cell.jitter=0)


# Figure 1
if(run.figures == T){
  source(paste(top.wd, "analysis", "nh fig global metanetwork.R", sep = "/"))
}





# Shared species map -----------------------------------------------------------

# Note that the same could be done with networks by replacing reg with net.id in
# the following chunk (but it would take a long time)

# Make a matrix with regs as rows and columns
regs <- unique(net.long$reg)
reg.mat <- matrix(NA, 
                  nrow=length(regs),
                  ncol=length(regs)) %>% as.data.frame()
rownames(reg.mat) <- regs
colnames(reg.mat) <- regs

all.plant.mat.reg.num <- reg.mat
all.animal.mat.reg.num <- reg.mat
nat.plant.mat.reg.num <- reg.mat
nat.animal.mat.reg.num <- reg.mat
all.int.mat.reg.num <- reg.mat
nat.int.mat.reg.num <- reg.mat

all.plant.mat.reg.denom <- reg.mat
all.animal.mat.reg.denom <- reg.mat
nat.plant.mat.reg.denom <- reg.mat
nat.animal.mat.reg.denom <- reg.mat
all.int.mat.reg.denom <- reg.mat
nat.int.mat.reg.denom <- reg.mat


skip.these <- NA

for(i in regs){
  focal.all <- subset(net.long, reg == i)
  focal.nat <- subset(net.long.nat, reg == i)
  
  focal.all.plant <- unique(focal.all$plant.accepted.species)
  focal.all.animal <- unique(focal.all$animal.accepted.species)
  focal.nat.plant <- unique(focal.nat$plant.accepted.species)
  focal.nat.animal <- unique(focal.nat$animal.accepted.species)
  
  focal.all.int <- unique(paste(focal.all$plant.accepted.species, focal.all$animal.accepted.species))
  focal.nat.int <- unique(paste(focal.nat$plant.accepted.species, focal.nat$animal.accepted.species))
  
  
  # Will only fill out the upper triangle of the matrix and not the diagonal
  # which is always 1
  skip.these <- c(skip.these, i)
  
  for(j in regs[!regs %in% skip.these]){
    compare.all <- subset(net.long, reg == j)
    compare.nat <- subset(net.long.nat, reg == j)
    
    compare.all.plant <- unique(compare.all$plant.accepted.species)
    compare.all.animal <- unique(compare.all$animal.accepted.species)
    compare.nat.plant <- unique(compare.nat$plant.accepted.species)
    compare.nat.animal <- unique(compare.nat$animal.accepted.species)
    
    compare.all.int <- unique(paste(compare.all$plant.accepted.species, compare.all$animal.accepted.species))
    compare.nat.int <- unique(paste(compare.nat$plant.accepted.species, compare.nat$animal.accepted.species))
    
    # For this, the goal is to look at the proportions, but we will make 
    # matrices reflecting the numerator (the intersection) and the denominator
    # so that we can perform binomial GLMs
    all.plant.mat.reg.num[i,j] <- length(intersect(compare.all.plant, focal.all.plant))
    all.plant.mat.reg.denom[i,j] <- length(union(compare.all.plant, focal.all.plant))
    
    nat.plant.mat.reg.num[i,j] <- length(intersect(compare.nat.plant, focal.nat.plant))
    nat.plant.mat.reg.denom[i,j] <- length(union(compare.nat.plant, focal.nat.plant))
    
    
    all.animal.mat.reg.num[i,j] <- length(intersect(compare.all.animal, focal.all.animal))
    all.animal.mat.reg.denom[i,j] <- length(union(compare.all.animal, focal.all.animal))
    
    nat.animal.mat.reg.num[i,j] <- length(intersect(compare.nat.animal, focal.nat.animal))
    nat.animal.mat.reg.denom[i,j] <- length(union(compare.nat.animal, focal.nat.animal))
    
    
    all.int.mat.reg.num[i,j] <- length(intersect(compare.all.int, focal.all.int))
    all.int.mat.reg.denom[i,j] <- length(union(compare.all.int, focal.all.int))
    
    nat.int.mat.reg.num[i,j] <- length(intersect(compare.nat.int, focal.nat.int))
    nat.int.mat.reg.denom[i,j] <- length(union(compare.nat.int, focal.nat.int))
    

  }
  
}



share.any.fun <- function(all.mat, nat.mat){
  
  share.any.df <- rbind(ifelse(unlist(all.mat) > 0, 1, 0) %>% cbind(., "all"),
                        ifelse(unlist(nat.mat) > 0, 1, 0) %>% cbind(., "nat")) %>% as.data.frame()
  share.any.df[,1] <- as.numeric(as.character(share.any.df[,1]))
  colnames(share.any.df) <- c("share.any", "all.or.nat")
  share.any.df$reg.combo <- gsub(".1", "", rownames(share.any.df))
  share.any.df <- share.any.df[complete.cases(share.any.df),]
  share.any.mod <- glmer(as.numeric(share.any) ~ all.or.nat + (1|reg.combo), data = share.any.df, family = "binomial")
  share.any.mod0 <- glmer(as.numeric(share.any) ~ 1 + (1|reg.combo), data = share.any.df, family = "binomial")
  
  out.without <- mean(unlist(nat.mat) > 0, na.rm = T)
  out.with <- mean(unlist(all.mat) > 0, na.rm = T)
  out.p <- anova(share.any.mod, share.any.mod0)$"Pr(>Chisq)"[2]
  
  cbind(out.without, out.with, out.p)
  
}

prop.share.fun <- function(all.num.mat, all.denom.mat,
                           nat.num.mat, nat.denom.mat){
  
  prop.share.df <- rbind(cbind(unlist(all.num.mat), unlist(all.denom.mat), "all"),
                         cbind(unlist(nat.num.mat), unlist(nat.denom.mat), "nat")) %>% as.data.frame()
  prop.share.df[,1] <- as.numeric(as.character(prop.share.df[,1]))
  prop.share.df[,2] <- as.numeric(as.character(prop.share.df[,2]))
  colnames(prop.share.df) <- c("num", "denom", "all.or.nat")
  prop.share.df$reg.combo <- gsub(".1", "", rownames(prop.share.df))
  prop.share.df <- prop.share.df[complete.cases(prop.share.df),]
  prop.share.mod <- glmer(cbind(num, denom-num) ~ all.or.nat + (1|reg.combo), data = prop.share.df, family = "binomial")
  prop.share.mod0 <- glmer(cbind(num, denom-num) ~ 1 + (1|reg.combo), data = prop.share.df, family = "binomial")
  
  out.without <- mean(unlist(nat.num.mat / nat.denom.mat), na.rm = T)
  out.with <- mean(unlist(all.num.mat / all.denom.mat), na.rm = T)
  out.p <- anova(prop.share.mod, prop.share.mod0)$"Pr(>Chisq)"[2]
  
  cbind(out.without, out.with, out.p)
  
}

shared.table <- matrix(NA, ncol = 6, nrow = 4) %>% as.data.frame()
rownames(shared.table) <- c("", "Plants", "Animals", "Interactions")
colnames(shared.table)[1] <- "Proportion of regions that that share..."
colnames(shared.table)[4] <- "Average proportion shared..."
shared.table[1,] <- c("native only", "with introd.", "P")
shared.table[2,1:3] <- share.any.fun(all.plant.mat.reg.num, nat.plant.mat.reg.num) %>% signif(2)
shared.table[3,1:3] <- share.any.fun(all.animal.mat.reg.num, nat.animal.mat.reg.num) %>% signif(2)
shared.table[4,1:3] <- share.any.fun(all.int.mat.reg.num, nat.int.mat.reg.num) %>% signif(2)

shared.table[2,4:6] <- prop.share.fun(all.plant.mat.reg.num, all.plant.mat.reg.denom, 
                                      nat.plant.mat.reg.num, nat.plant.mat.reg.denom) %>% signif(2)
shared.table[3,4:6] <- prop.share.fun(all.animal.mat.reg.num, all.animal.mat.reg.denom, 
                                      nat.animal.mat.reg.num, nat.animal.mat.reg.denom) %>% signif(2)
shared.table[4,4:6] <- prop.share.fun(all.int.mat.reg.num, all.int.mat.reg.denom, 
                                      nat.int.mat.reg.num, nat.int.mat.reg.denom) %>% signif(2)

shared.table





# Make a figure
if(run.figures == T){
  source(paste(top.wd, "analysis", "nh fig great circle.R", sep = "/"))
}



# Network Beta diversity -------------------------------------------------------

if(run.figures == T){
  source(paste(top.wd, "analysis", "nh betalink.R", sep = "/"))
}



# Local network scale analyses -------------------------------------------------

# Figure 3 analyses

# This chunk builds dataframes saying how many native and introduced
# partners that individual focal species have.

# First, plants
# Determine if a plant species appears as an introduced species in the database
net.id.plants.sp <- paste(net.long0$net.id, net.long0$plant.accepted.species, sep = "&")

unique.net.id.plants.sp <- unique(net.id.plants.sp)

plant.intro.prop.df <- matrix(nrow  = length(unique.net.id.plants.sp),
                              ncol = 5) %>%  as.data.frame()
colnames(plant.intro.prop.df) <- c("net.id", "plant.accepted.species", "plant.native.status", "n.int", "n.nat")
plant.intro.prop.df$net.id <- word(unique.net.id.plants.sp, 1, sep = "&")
plant.intro.prop.df$plant.accepted.species <- word(unique.net.id.plants.sp, 2, sep = "&")

for(i in 1:length(unique.net.id.plants.sp)){
  inds <- which(net.id.plants.sp == unique.net.id.plants.sp[i])
  
  statuses <- subset(net.long0[inds,], value>0)$animal.native.status
  
  plant.intro.prop.df$plant.native.status[i] <- as.character(net.long0$plant.native.status[inds[1]])
  
  plant.intro.prop.df[i, c("n.int", "n.nat")] <- table(statuses)
  
}

plant.intro.prop.df <- left_join(plant.intro.prop.df, metanet[,c("net.id", "reg")])


# remove unidentified plants for analysis
plant.intro.prop.df <- plant.intro.prop.df[complete.cases(plant.intro.prop.df),]

# Also need to keep ones that interact with at least 1 accepted species
plant.intro.prop.df <- plant.intro.prop.df[(plant.intro.prop.df$n.int + plant.intro.prop.df$n.nat) != 0,]

# Find the net.ids that have potential introduced partners for these plants
net.ids.with.intro.animals <- ddply(plant.intro.prop.df, "net.id", summarize, sum(n.int))
net.ids.with.intro.animals <- net.ids.with.intro.animals[net.ids.with.intro.animals[,2]>0,1]

plant.intro.prop.df$plant.native.status <- factor(plant.intro.prop.df$plant.native.status, levels = c("native", "introduced"))
plant.intro.prop.mod <- glmer(cbind(n.int, n.nat) ~ plant.native.status + (1|net.id) + (1|reg),
                              family = "binomial",
                              data = subset(plant.intro.prop.df, net.id %in% net.ids.with.intro.animals))
plant.intro.prop.mod0 <- glmer(cbind(n.int, n.nat) ~ 1 + (1|net.id) + (1|reg),
                               family = "binomial",
                               data = subset(plant.intro.prop.df, net.id %in% net.ids.with.intro.animals))
summary(plant.intro.prop.mod)
anova(plant.intro.prop.mod, plant.intro.prop.mod0)




# Second, animals
# Determine if a animal species appears as an introduced species in the database
net.id.animals.sp <- paste(net.long0$net.id, net.long0$animal.accepted.species, sep = "&")

unique.net.id.animals.sp <- unique(net.id.animals.sp)

animal.intro.prop.df <- matrix(nrow  = length(unique.net.id.animals.sp),
                               ncol = 5) %>%  as.data.frame()
colnames(animal.intro.prop.df) <- c("net.id", "animal.accepted.species", "animal.native.status", "n.int", "n.nat")
animal.intro.prop.df$net.id <- word(unique.net.id.animals.sp, 1, sep = "&")
animal.intro.prop.df$animal.accepted.species <- word(unique.net.id.animals.sp, 2, sep = "&")

for(i in 1:length(unique.net.id.animals.sp)){
  inds <- which(net.id.animals.sp == unique.net.id.animals.sp[i])
  
  statuses <- subset(net.long0[inds,], value>0)$plant.native.status
  
  animal.intro.prop.df$animal.native.status[i] <- as.character(net.long0$animal.native.status[inds[1]])
  
  animal.intro.prop.df[i, c("n.int", "n.nat")] <- table(statuses)
  
}

animal.intro.prop.df <- left_join(animal.intro.prop.df, metanet[,c("net.id", "reg")])

# remove unidentified animals for analysis
animal.intro.prop.df <- animal.intro.prop.df[complete.cases(animal.intro.prop.df),]

# Also need to keep ones that interact with at least 1 accepted species
animal.intro.prop.df <- animal.intro.prop.df[(animal.intro.prop.df$n.int + animal.intro.prop.df$n.nat) != 0,]


# Find the net.ids that have potential introduced partners for these animals
net.ids.with.intro.plants <- ddply(animal.intro.prop.df, "net.id", summarize, sum(n.int))
net.ids.with.intro.plants <- net.ids.with.intro.plants[net.ids.with.intro.plants[,2]>0,1]

animal.intro.prop.df$animal.native.status <- factor(animal.intro.prop.df$animal.native.status, levels = c("native", "introduced"))
animal.intro.prop.mod <- glmer(cbind(n.int, n.nat) ~ animal.native.status + (1|net.id) + (1|reg),
                               family = "binomial",
                               data = subset(animal.intro.prop.df, net.id %in% net.ids.with.intro.plants))
animal.intro.prop.mod0 <- glmer(cbind(n.int, n.nat) ~ 1 + (1|net.id) + (1|reg),
                                family = "binomial",
                                data = subset(animal.intro.prop.df, net.id %in% net.ids.with.intro.plants))
summary(animal.intro.prop.mod)
anova(animal.intro.prop.mod,animal.intro.prop.mod0, test = "LRT")






### Extended Data Figure 3 analyses

# How interactive are introduced plant species in their native range?
# VERSION considering plants that appear ANYWHERE (as determined by global databases) as introduced
plant.sp.mod <- glmer(cbind(degree, 
                            partner.richness - degree) ~ plant.introduced.anywhere + (1|net.id),
                      data = subset(nd.plant.analysis, plant.native.status == "native" & !is.na(plant.introduced.anywhere)), 
                      family = "binomial")
plant.sp.mod0 <- glmer(cbind(degree, 
                             partner.richness - degree) ~ 1 + (1|net.id),
                       data = subset(nd.plant.analysis, plant.native.status == "native" & !is.na(plant.introduced.anywhere)), 
                       family = "binomial")
summary(plant.sp.mod)
anova(plant.sp.mod, plant.sp.mod0, test = "LRT")


# How interactive are introduced plant species in their native range?
# VERSION considering only plants that appear in the DATABASE as introduced
plant.sp.db.mod <- glmer(cbind(degree,
                               partner.richness - degree) ~ appears.in.db.as.intro + (1|net.id),
                         data = subset(nd.plant.analysis, plant.native.status == "native"),
                         family = "binomial")
plant.sp.db.mod0 <- glmer(cbind(degree,
                                partner.richness - degree) ~ 1 + (1|net.id),
                          data = subset(nd.plant.analysis, plant.native.status == "native"),
                          family = "binomial")
summary(plant.sp.db.mod)
anova(plant.sp.db.mod, plant.sp.db.mod0, test = "LRT")





# How interactive are introduced animal species in their native range?
# VERSION considering animals that appear ANYWHERE (as determined by global databases) as introduced
animal.sp.mod <- glmer(cbind(degree, 
                             partner.richness - degree) ~ animal.introduced.anywhere + (1|net.id),
                       data = subset(nd.animal.analysis, animal.native.status == "native"), 
                       family = "binomial")
animal.sp.mod0 <- glmer(cbind(degree, 
                              partner.richness - degree) ~ 1 + (1|net.id),
                        data = subset(nd.animal.analysis, animal.native.status == "native"), 
                        family = "binomial")
summary(animal.sp.mod)
anova(animal.sp.mod, animal.sp.mod0, test = "LRT")


# How interactive are introduced animal species in their native range?
# VERSION considering only animals that appear in the DATABASE as introduced
animal.sp.db.mod <- glmer(cbind(degree,
                                partner.richness - degree) ~ appears.in.db.as.intro + (1|net.id),
                          data = subset(nd.animal.analysis, animal.native.status == "native"),
                          family = "binomial")
animal.sp.db.mod0 <- glmer(cbind(degree,
                                 partner.richness - degree) ~ 1 + (1|net.id),
                           data = subset(nd.animal.analysis, animal.native.status == "native"),
                           family = "binomial")
summary(animal.sp.db.mod)
anova(animal.sp.db.mod, animal.sp.db.mod0, test = "LRT")




# Local network modularity analyses
metanet$latlon.id <- paste(metanet$latitude, metanet$longitude)
metanet.mod.set <- metanet[complete.cases(metanet[,c("Modularity", "Modularity.null", "prop.intro.int", "latlon.id")]),]
modularity.mod <- lmer(c(Modularity - Modularity.null) ~ prop.intro.int + (1|latlon.id),
     data = metanet.mod.set)
summary(modularity.mod)
modularity.mod0 <- lmer(c(Modularity - Modularity.null) ~ 1 + (1|latlon.id),
     data = metanet.mod.set)
anova(modularity.mod, modularity.mod0, test = "LRT")

modularity.mod.full <- lmer(c(Modularity - Modularity.null) ~ prop.intro.int + n.target + connectance + oceanic.island + (1|latlon.id),
                        data = metanet.mod.set, na.action = "na.fail", REML = F)
modularity.dredge <- MuMIn::dredge(modularity.mod.full, rank = "AIC")
modularity.dredge
modularity.table <- modularity.dredge[,c("prop.intro.int", "n.target", "connectance","oceanic.island",
                                         "AIC", "delta")]
modularity.table[,c("AIC", "delta")] <- modularity.table[,c("AIC", "delta")] %>% round(digits = 1)
modularity.table

fit.lqmm.mod <- lqmm(fixed = c(Modularity - Modularity.null) ~ prop.intro.int, random = ~ 1, group = latlon.id, 
                     data = metanet.mod.set,
                     na.action = na.omit, type = "normal",
                     nK = 3,
                     tau = c(0.05, 0.5, 0.95),
                     control = list(LP_max_iter=100000))
summary(fit.lqmm.mod)

# Make the figures

if(run.figures == T){
  source(paste(top.wd, "analysis", "nh fig intro species roles.R", sep = "/"))
}






# Spatiotemporal patterns ------------------------------------------------------

# Analyses
metanet$int.intro <- metanet$prop.intro.int * metanet$rich.int
metanet$int.native <- metanet$rich.int - metanet$int.intro
ghm.mod <- glmer(cbind(int.intro, int.native) ~ ghm + (1|reg), data = metanet, family = "binomial")
mod0 <- glmer(cbind(int.intro, int.native) ~ 1 + (1|reg), data = metanet, family = "binomial")
summary(ghm.mod)
anova(ghm.mod, mod0, test = "LRT")
ghm.pred <- get_model_data(ghm.mod, type = "pred")

metanet$year2000 <- metanet$year - 2000
year.mod <- glmer(cbind(int.intro, int.native) ~ year2000 + (1|reg), data = metanet, family = "binomial")
summary(year.mod)
anova(year.mod, mod0, test = "LRT")
year.pred <- get_model_data(year.mod, type = "pred")


# Supplementary analysis

# Using a more conservative set of local networks where all potential interactions were observed shows a qualitatively identical relationship
year.modX1 <- glmer(cbind(int.intro, int.native) ~ year2000 + (1|reg), 
                    data = subset(metanet, plant.coverage == "comp" & target.birds == "comp" | target.bats == "comp" | target.primates == "comp"), 
                    family = "binomial")
year.modX0 <- glmer(cbind(int.intro, int.native) ~ 1 + (1|reg), 
                    data = subset(metanet, plant.coverage == "comp" & target.birds == "comp" | target.bats == "comp" | target.primates == "comp"), 
                    family = "binomial")
anova(year.modX1, year.modX0, test = "LRT")
year.predX1 <- get_model_data(year.modX1, type = "pred")


# We found no relationship between study year and the study location’s human modification values 
# Compare AIC for models involving these two predictors and their interaction
year.ghm.int.mod <- glmer(cbind(int.intro, int.native) ~ year2000 * ghm + (1|reg), data = metanet, family = "binomial")
year.ghm.mod <- glmer(cbind(int.intro, int.native) ~ year2000 + ghm + (1|reg), data = metanet, family = "binomial")

anova(year.ghm.int.mod, year.ghm.mod, ghm.mod, year.mod, mod0)

# Reanalyzing using only local networks collected after 2000 or 2010 shows qualitatively identical relationships

ghm.mod.2000 <- glmer(cbind(int.intro, int.native) ~ ghm + (1|reg), 
                      data = subset(metanet[!is.na(metanet$ghm),], year >= 2000), 
                      family = "binomial")
ghm.mod0.2000 <- glmer(cbind(int.intro, int.native) ~ 1 + (1|reg), 
                       data = subset(metanet[!is.na(metanet$ghm),], year >= 2000), 
                       family = "binomial")
summary(ghm.mod.2000)
anova(ghm.mod.2000, ghm.mod0.2000, test = "LRT")
ghm.pred.2000 <- get_model_data(ghm.mod.2000, type = "pred")

ghm.mod.2010 <- glmer(cbind(int.intro, int.native) ~ ghm + (1|reg), 
                      data = subset(metanet[!is.na(metanet$ghm),], year >= 2010), 
                      family = "binomial")
ghm.mod0.2010 <- glmer(cbind(int.intro, int.native) ~ 1 + (1|reg), 
                       data = subset(metanet[!is.na(metanet$ghm),], year >= 2010), 
                       family = "binomial")
summary(ghm.mod.2010)
anova(ghm.mod.2010, ghm.mod0.2010, test = "LRT")
ghm.pred.2010 <- get_model_data(ghm.mod.2010, type = "pred")


# Make the figures (and run some supplementary analyses)
if(run.figures == T){
  source(paste(top.wd, "analysis", "nh fig spatiotemporal.R", sep = "/"))
}


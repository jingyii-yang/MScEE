rm(list=ls())
setwd("~/R/# ProJ")
        
        library(readxl)
        library(tidyr)
        library(dplyr)
        library(stringr)

    library(ape)
    library(caper)
    library(phylolm)
    library(rr2)

        library(phylobase)
        library(ggplot2)
        library(ggtree)
        library(cowplot)

## (The data is currently unavailable for public access. Use upon request.)
a <- read.csv('Analysis.csv')

# Remove the least certain data & NA.
str_which(colnames(a), 'certain')
a <- subset(a, a[, 7] != 4 & a[,9] != 4 & a[,11] != 4)  
a <- na.omit(a)

str(a)                          # 1136 species.   

        

################################# 1. DATA PREPARATION ###################################

## Check Xs.

# Rare levels in each categorical variable? 

table(a$Incubation.roles)         # 0:2 = 549:587 (species)
table(a$SS_Scores)                # 0:1:2:3 = 986:19:21:110. Should unite 1-2(Monogamous), 3-4(Polygynous).
table(a$Nest.invisibility)        # 0:1 = 517:619
table(a$Territory)                # 0:1 = 122:1014

# Distributions of the continuous vairable?

hist(abs(a$Latitude_median))
  hist(sqrt(abs(a$Latitude_median)))   # Need square-root transformation.

  
##======================== Data Transformation ===========================##
  
# Assign intepretable names to predictors. Set reference levels.
  
a$inc <- ifelse(a$Incubation.roles == 2, '1.Biparental', '2.Female only')
a$ss  <- ifelse(a$SS_Scores %in% c(0, 1), '1.Monogamous', '2.Polygynous')
a$nest_inv <- ifelse(a$Nest.invisibility == 1, '1.invisible', '2.visible')
a$terr <- ifelse(a$Territory %in% c('none'), '1.no_terr', '2.yes_terr')

# Transform continuous variables for modelling.

a$lat <- sqrt(abs(a$Latitude_median))
    

# Get the partial dichromatism scores. 

a$dichro_vis <- a$Head_dichro + a$Upper_dichro + a$Tail_dichro 
a$dichro_inv <- a$Under_dichro + a$Wing_dichro


# Rare combinations? (for using interactions in models)

a %>% count(Incubation.roles, Nest.invisibility) # 4 combinations even.(305:244:212:375)


## Check Ys.

# Distributions of 8x dichromatism data.

par(mfrow = c(3,3), mar = c(3,3,2,1))
for (i in c(15:20, 26, 27)){
  hist(a[,i], main = colnames(a)[i])     
}            


## Relationships between each X and Y. (see below: Visualisation -> Figure S1.)


##========= Use high certainty data (certainty = 1-2) for models. ==========##

a2 <- subset(a, a$Incubation.certainty %in% c(1, 2) &         # 879 species.
                 a$SS_cert %in% c(1, 2) &
                 a$Nest.certainty.scores %in% c(1, 2))

## Collinearity?

a2$ss_n <- as.numeric( as.factor(a2$ss) )
a2$terr_n <- as.numeric( as.factor(a2$terr))

# VIF:

VIF <- function(i){
    v <- c('Incubation.roles', 'Nest.invisibility', 'ss_n', 'terr_n', 'lat')
    vi <- which(colnames(a2) %in% v)
    fml <- formula(paste(v[i], '~ . -', v[i]))
    r2 <- summary(lm(fml, data = a2[, vi]))$adj.r.squared
    1 / (1 - r2)
}

VIF(1)    #  1.29    
VIF(2)    #  1.10          
VIF(3)    #  2.67          
VIF(4)    #  2.42          
VIF(5)    #  1.02         

# mean VIFs: 1.70


##======================== Merge GBD & the tree ==========================##

# Load the tree and the file for tip names matching.

t <- ape::read.tree('T400F_AOS_sppnames.tre')
n <- read_xlsx('tip_names.xlsx')

# Function for merging GBD and the tree.
tGBD <- function(d){
    an <- merge(n, d, by = 'Unique_Scientific_Name', all.y = TRUE)
    # remove the duplicated tip names
    an <- distinct(an, AOS.Howard.Moore.Name, .keep_all = TRUE)   
    # for name matching
    an$AOS.Howard.Moore.Name <- gsub(' ', '_', an$AOS.Howard.Moore.Name) 
    an <-na.omit(an)
    # Merge the tree with GBD
    caper::comparative.data(data = an, phy = t,           
                 names.col = 'AOS.Howard.Moore.Name',
                 vcv = TRUE, warn.dropped = TRUE)
}

# Use this data for the main analysis.

ant_2 <- tGBD(a2)

length(ant_2$phy$tip.label)                   # 879 species. 



##################################### 2. ANALYSIS ######################################

# PGLS models and variance partition.

Mod.Var <- function(bp){
  
  ## prepare formulas.
  
      fml.full <- formula(paste(bp, '~', 'inc', '*', 'nest_inv', '+', 'ss', 
                                '+', 'terr', '+', 'lat'))
      fml.reduce <- formula(paste(bp, '~', 1))

  ## PGLS (lambda) models & get the three types of R2 values.
      
  r2 <- data.frame(Body_part = bp, Full.model.R2 = NA, Phylogeny.R2 = NA, Predictors.R2 = NA)

    # 1) Full model -- to get the full model R2.
      bm <- phylolm(fml.full,
                    data = ant_2$data, phy = ant_2$phy, model = 'lambda')
      r2$Full.model.R2 <- R2.lik(bm)          
      
    # 2) Remove phylo structure (lambda = 0) -- to get phylogeny R2.
      nop <- data.frame(lambda = 0)
      bmP <- phylolm(fml.full,
                     data = ant_2$data, phy = ant_2$phy, model = 'lambda',
                     starting.value = nop, lower.bound = nop$lambda, upper.bound = nop$lambda)
      r2$Phylogeny.R2 <- R2.lik(bm, bmP)    
      
     # 3) Intercept-only model (no predictors) -- to get predictors R2.
      nov <- data.frame(lambda = bm$optpar)
      bmV <- phylolm(fml.reduce, 
                     data = ant_2$data, phy = ant_2$phy, model = 'lambda', 
                     starting.value = nov, lower.bound = nov$lambda, upper.bound = nov$lambda) 
      r2$Predictors.R2 <- R2.lik(bm, bmV)   
      
  ## save the model results and R2 values.
  res <- list(bm, r2) # bm: Table S2-4.  r2: Table S5.
  return(res)

}

# Model results.

bp1 <- Mod.Var('Dichro_total')      
  bp2 <- Mod.Var('dichro_vis')          
    bp2.1 <- Mod.Var('Head_dichro')        
    bp2.2 <- Mod.Var('Upper_dichro')       
    bp2.3 <- Mod.Var('Tail_dichro')        
  bp3 <- Mod.Var('dichro_inv')          
    bp3.1 <- Mod.Var('Under_dichro')        
    bp3.2 <- Mod.Var('Wing_dichro')        
    
    
### Table S2.
summary(bp1[[1]])

  ### Table S3.
  summary(bp2[[1]])
  summary(bp3[[1]])
  
    ### Table S4.
    summary(bp2.1[[1]])
    summary(bp2.2[[1]])
    summary(bp2.3[[1]])
    summary(bp3.1[[1]])
    summary(bp3.2[[1]])



###################################### 3. VISUALISATION #######################################
    
    
##### 1) Figure 2. The tree (dichroamtism) and heatmaps (natural+sexual selection) #####
    
# Use this data for non-analysis plots.

ant <- tGBD(a)     # 1136 sp.
      
## Tree plot:
    
# Select the total dichromatism data (need tip names later)
    
d4p <- data.frame(dich = ant$data$Dichro_total, sp = row.names(ant$data))

# Replace the dicrhomatism scores (0:11) with 11 gradient colours (navy:gold).

cc <- scales::seq_gradient_pal(low = 'navy', high = 'gold')(seq(0, 1, length.out = 11))
cc <- data.frame(colors = cc, dich = 0:10)
d4p <- merge(cc, d4p, by = 'dich')
row.names(d4p) <- d4p$sp

# Assemble the tree and dichromatism data (shown as branch colours).
d4pp <- as.matrix(d4p)[,2]
tp <- phylo4d(as(ant$phy,'phylo4'), d4pp) 

# Add node colours to the tree.
nodeData(tp) <- data.frame(dt = rep('grey42', times = nNodes(tp)), row.names = nodeId(tp, 'internal'))


## The central tree:
tp1 <- ggtree(tp, aes(col = I(dt)), layout = 'circular')


## Add heatmap-1 (inner circle): incubation roles.
i4p <- as.data.frame(ant$data$Incubation.roles)    
i4p[i4p == '0'] <- 'Female only'
i4p[i4p == '2'] <- 'Biparental'
rownames(i4p) <- rownames(ant$data)

tp2 <- gheatmap(tp1, i4p, 
                offset = 0.1, width = 0.05, colnames = FALSE) 


## Add heatmap-2 (middle circle): visibility.
v4p <- as.data.frame(ant$data$Nest.invisibility)    
v4p[v4p == '0'] <- 'Visible'
v4p[v4p == '1'] <- 'Invisible'
rownames(v4p) <- rownames(ant$data)

tp3 <- gheatmap(tp2, v4p,
                offset = 3.6, width = 0.05, colnames = FALSE) +
        scale_fill_manual(values = c('lightblue1', 'skyblue3', 'coral2', 'peachpuff1')) +
        labs(fill = 'Natural Selection')  # need its legend.


## Add heatmap-3 (outer circle): mating system.
s4p <- as.data.frame(ant$data$ss)
s4p[s4p == '1.Monogamous'] <- 'Monogamous'
s4p[s4p == '2.Polygynous'] <- 'Polygynous'
rownames(s4p) <- rownames(ant$data)

tp4 <- gheatmap(tp1, s4p,
                offset = 3.6, width = 0.05, colnames = FALSE) +
        scale_fill_manual(values = c('olivedrab', 'olivedrab1')) + 
        labs(fill = 'Sexual Selection')  # need its legend.


### Final plot with the tree and 3 heatmaps (Figure 2):

gheatmap (tp3, 
          s4p, 
          legend_title = '',
          offset = 7.2, width = 0.05, colnames = FALSE) +
  scale_fill_manual(values = c('lightblue1', 'skyblue3', 'coral2',
                               'olivedrab', 'olivedrab1', 'peachpuff1')) +
  theme(legend.position = "none")
 

### retrieve legends of previous layers.

plot(get_legend(tp3))      # << heatmap 1-2 (Natural Selection).

plot(get_legend(tp4))      # << heatmap 3   (Sexual Selection).

d4p <- as.data.frame(ant$data$Dichro_total)   
rownames(d4p) <- rownames(ant$data)
gheatmap (tp1,               # << the tree (total Dichromatism).
          d4p, 
          low = 'navy', high = 'gold',
          legend_title = 'Total Dichromatism',
          offset = 7.2, width = 0.05, colnames = FALSE) 
 
  

############## 2) Figure 3: variance partition (bar plot) ################
  
### Table S5.
  
order <- c('TOTAL', 'Visible body parts', 'Head', 'Upper Part', 'Tail', 'Invisible body parts', 'Underparts', 'Wings')
bp.R2 <- rbind(bp1[[2]], bp2[[2]], bp2.1[[2]], bp2.2[[2]], 
               bp2.3[[2]], bp3[[2]], bp3.1[[2]], bp3.2[[2]])
  bp.R2$Body_part <- order
  
  
# To use full model r2 as the background (layer 1).
bp.R2.full <- data.frame(Body_part = order, R2_value = bp.R2$Full.model.R2, R2.types = 'Full.model.R2')

# To stack phylogeny and predictors r2 in the forground (layer 2).
bp.R2.redu <- tidyr::gather(bp.R2, key = 'R2.types', value = 'R2_value', c(3, 4))

### Figure 3.

ggplot(data = bp.R2.full, aes(factor(Body_part, level = rev(order)), R2_value, fill = R2.types)) + 
  geom_col(alpha = 0.5) +
  geom_col(data = bp.R2.redu) + 
  scale_fill_manual(values=c('grey60', 'darkolivegreen3', 'slateblue1')) + 
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs(x = 'Body Parts', y = 'Partial R2 Values') + ggtitle('Variance Partition') +
  coord_flip() 



############## 3) Figure 4 & 5: parameter estimates (forest plots) ##################

figs <- function(x, A){
    s <- as.data.frame(summary(x[[1]])$coefficients)
    s$Predictors <- c('(intercept)','Female only incubation (invisible)', 'Visible (biparental incubation)', 
                     'Polygynous', 'Territorial', 'Latitude', 'Female only incubation : Visible')  
    s$lower <- s$Estimate - 1.96 * s$StdErr  
    s$upper <- s$Estimate + 1.96 * s$StdErr
    s$Significance <- ifelse(s$lower * s$upper > 0, 'Significant', 'Non-significant')
    s$Effect <- ifelse(s$Estimate > 0, 'Positive', 'Negative')
    
    ggplot(s[-1, ], aes(Predictors, Estimate)) + 
              geom_pointrange(aes(ymin = lower, ymax = upper, shape = Effect, col = Significance)) +
              geom_hline(yintercept = 0, lty = 2) + 
              ylim(-A, A) + labs(y = 'Estimates') +
              coord_flip() 

}

p1 <- figs(bp1, 3) + ggtitle('Total Dichromatism')   
p2 <- figs(bp2, 2) + ggtitle('Visible Body Parts')
p3 <- figs(bp3, 2) + ggtitle('Invisible Body Parts')


### Figure 4.
    
p1 
    
### Figure 5.
    
pvis <- p2 +  theme(legend.position = "none") 

pinv <- p3 + theme(axis.title.y=element_blank(), axis.text.y=element_blank())

cowplot::plot_grid(pvis, pinv, ncol = 2, rel_widths = c(1.15, 1))



########## 4) Figure S1: dichromatism ~ predictors in raw data (boxplots) ################

a$NS[a$Incubation.roles == 0 & a$Nest.invisibility == 0] <- 'Female:Visible'
a$NS[a$Incubation.roles == 0 & a$Nest.invisibility == 1] <- 'Female:Invisible'
a$NS[a$Incubation.roles == 2 & a$Nest.invisibility == 0] <- 'Biparental:Visible'
a$NS[a$Incubation.roles == 2 & a$Nest.invisibility == 1] <- 'Biparental:Invisible'

a$ss_names[a$SS_Scores == 0] <- 'Strictly Monogamous'
a$ss_names[a$SS_Scores == 1] <- 'EPP <5%'
a$ss_names[a$SS_Scores == 2] <- 'EPP 5-20%'
a$ss_names[a$SS_Scores == 3] <- 'Polygynous'

a$terr_names[a$terr == '1.no_terr'] <- 'Non-territorial'
a$terr_names[a$terr == '2.yes_terr'] <- 'Territorial'


ps <- ggplot(a, aes(factor(ss_names, levels = c('Strictly Monogamous', 'EPP <5%', 'EPP 5-20%', 
                                                'Polygynous')), Dichro_total), group = ss_names) + 
            geom_boxplot(fill = 'grey90') +
            labs(x = 'Sexual Selection', y = 'Total Dichromatism') +
            theme_minimal()

pn <- ggplot(a, aes(NS, Dichro_total), group = NS) + 
            geom_boxplot(fill = 'grey90') +
            labs(x = 'Natural Selection', y = '') + 
            theme_minimal()

pt <- ggplot(a, aes(terr_names, Dichro_total), group = terr_names) + 
            geom_boxplot(fill = 'grey90') +
            labs(x = 'Social Selection', y = 'Total Dichromatism') +
            theme_minimal()

pl <- ggplot(a, aes(Latitude_median, jitter(Dichro_total)) ) + 
            geom_point(size = 0.7) +
            labs(x = 'Median Latitude', y = '')


### Figure S1.

cowplot::plot_grid(ps, pn, pt, pl, ncol = 2, scale = 0.95)




 

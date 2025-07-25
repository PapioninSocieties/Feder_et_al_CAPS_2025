## DATA AND CODE FOR FEDER ET AL. 2025
## COMPARATIVE ANALYSIS OF PAPION SOCIETIES
## Disparate social structures are underpinned by distinct social rules
## across a primate radiation

## LOAD NECESSARY PACKAGES
library(lubridate)
library(dplyr)
library(ggplot2)
library(igraph)
library(brms)
library(metafor)
library(orchaRd)
library(phytools)
library(ape)
library(treeio)
library(tidybayes)
library(performance)
library(modelr)
library(ggbeeswarm)
library(patchwork)
theme_set(theme_classic(base_size = 16))

## SETWD
setwd("~/Desktop/SPRF/PNAS/GitHub")

## GET GROUP-LEVEL DATA
caps_groups<-read.csv("CAPS_Groups.csv")

## LOAD TREE FILE FROM KUDERNA ET AL. 2023
treefile<-ape::read.tree("Kuderna_Tree2.nex.tree")
phylo1<-treefile
phylo1$edge.length

## TRIM TO INCLUDE SAMPLED TAXA
## AND CONGENERIC COGNATES
keep<-unique(caps_groups$phylo)
phylo1<-keep.tip(phylo1, keep)
phylo.check<-phylo1
plot.phylo(phylo.check,  edge.width = 2,cex=1, label.offset = 0)
A <- ape::vcv.phylo(phylo1, corr=T)

colnames(caps_groups)

## PART 1
## Papionin grooming networks vary as a function of network size and social system

## DENSITY MODEL 1 -- CLASSIC SOCIAL CATEGORIES
brm_density1<-brm(data=caps_groups, data2=list(A=A),
                       N.Groomed|trials(N.Dyads) ~ sqrt(GroupSize) + Cat.Classic + scale(Obs.Effort.Per.Capita) +
                         (1|gr(phylo, cov=A)) +
                         (1|species/population_id/group_id),
                       prior = c(
                         prior(normal(0, 2), "b")),
                       family="beta_binomial", iter=4000, chains=4, cores=2,
                       control=list(adapt_delta=0.95, max_treedepth=15))
summary(brm_density1, prob=0.89)
conditional_effects(brm_density1, "Cat.Classic")
conditional_effects(brm_density1, "GroupSize:Cat.Classic")
conditional_effects(brm_density1, "Obs.Effort.Per.Capita")

## DENSITY MODEL 2 -- NEW SOCIAL CATEGORIES
brm_density2<-brm(data=caps_groups, data2=list(A=A),
                       N.Groomed|trials(N.Dyads) ~ sqrt(GroupSize) + Cat.New + scale(Obs.Effort.Per.Capita) +
                         (1|gr(phylo, cov=A)) +
                         (1|species/population_id/group_id),
                       prior = c(
                         prior(normal(0, 2), "b")),
                       family="beta_binomial", iter=4000, chains=4, cores=2,
                       control=list(adapt_delta=0.95, max_treedepth=15))
summary(brm_density2, prob=0.89)
conditional_effects(brm_density2, "Cat.New")
conditional_effects(brm_density2, "GroupSize:Cat.New")

## DENSITY MODEL 3 -- GROUP SIZE ALONE
brm_density3<-brm(data=caps_groups, data2=list(A=A),
                       N.Groomed|trials(N.Dyads)  ~ sqrt(GroupSize) + scale(Obs.Effort.Per.Capita) +
                         (1|gr(phylo, cov=A)) +
                         (1|species/population_id/group_id),
                       prior = c(
                         prior(normal(0, 2), "b")),
                       family="beta_binomial", iter=4000, chains=4, cores=2,
                       control=list(adapt_delta=0.95, max_treedepth=15))
summary(brm_density3, prob=0.89)
conditional_effects(brm_density3, "GroupSize")

## DENSITY -- COMPARE THREE
loo1<-brms::loo(brm_density1, pointwise = T)
loo2<-brms::loo(brm_density2, pointwise = T)
loo3<-brms::loo(brm_density3, pointwise = T)
loo_compare(loo1,loo2, loo3)
model_weights(brm_density1, brm_density2, brm_density3, weights="stacking")

## MODULARITY MODEL 1 -- CLASSIC SOCIAL CATEGORIES
caps_groups$Modularity[caps_groups$Modularity==0.01]<-0
caps_groups$Cat.New<-factor(caps_groups$Cat.New, levels=c("Cohesive","Cliquish", "Multi-level"))
brm_modularity1<-brm(data=caps_groups, data2=list(A=A),
                          Modularity ~ sqrt(GroupSize) + Cat.Classic + scale(Obs.Effort.Per.Capita) +
                            (1|gr(phylo, cov=A)) +  
                            (1|species/population_id/group_id),
                          prior = c(
                            prior(normal(0, 2), "b")),
                          family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                          control=list(adapt_delta=0.99, max_treedepth=10))
summary(brm_modularity1, prob=0.89)

## MODULARITY MODEL 2 -- NEW SOCIAL CATEGORIES
brm_modularity2<-brm(data=caps_groups, data2=list(A=A),
                          Modularity ~ sqrt(GroupSize) + Cat.New + scale(Obs.Effort.Per.Capita) +
                            (1|gr(phylo, cov=A)) +  
                            (1|species/population_id/group_id),
                          prior = c(
                            prior(normal(0, 2), "b")),
                          family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                          control=list(adapt_delta=0.99, max_treedepth=10))
summary(brm_modularity2,prob=0.89)
conditional_effects(brm_modularity2, "Cat.New", prob=0.89)

## MODULARITY MODEL 3 -- GROUP SIZE ALONE
brm_modularity3<-brm(data=caps_groups, data2=list(A=A),
                          Modularity ~ sqrt(GroupSize) + scale(Obs.Effort.Per.Capita) +
                            (1|gr(phylo, cov=A)) +  
                            (1|species/population_id/group_id),
                          prior = c(
                            prior(normal(0, 2), "b")),
                          family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                          control=list(adapt_delta=0.99, max_treedepth=10))
summary(brm_modularity3, prob=0.89)

## MODULARITY -- COMPARE THREE
loo1<-brms::loo(brm_modularity1)
loo2<-brms::loo(brm_modularity2)
loo3<-brms::loo(brm_modularity3)
loo_compare(loo1,loo2, loo3)
model_weights(brm_modularity1,brm_modularity2,brm_modularity3, weights="stacking")

## MAKE TIDY PLOTS
grid = caps_groups %>%
  data_grid(GroupSize=seq(min(caps_groups$GroupSize), max(caps_groups$GroupSize)), N.Dyads=1, Cat.Classic, Cat.New, Obs.Effort.Per.Capita=15)

min.SLS<-min(caps_groups$GroupSize[caps_groups$Cat.New=="Cohesive"])
max.SLS<-max(caps_groups$GroupSize[caps_groups$Cat.New=="Cohesive"])
min.MLS<-min(caps_groups$GroupSize[caps_groups$Cat.New=="Multi-level"])
max.MLS<-max(caps_groups$GroupSize[caps_groups$Cat.New=="Multi-level"])
min.SemiLS<-min(caps_groups$GroupSize[caps_groups$Cat.New=="Cliquish"])
max.SemiLS<-max(caps_groups$GroupSize[caps_groups$Cat.New=="Cliquish"])
grid<-grid[!(grid$Cat.New=="Cohesive" & grid$GroupSize>max.SLS),]
grid<-grid[!(grid$Cat.New=="Multi-level" & grid$GroupSize<min.MLS),]
grid<-grid[!(grid$Cat.New=="Cliquish" & grid$GroupSize<min.SemiLS),]
grid<-grid[!(grid$Cat.New=="Cliquish" & grid$GroupSize>max.SemiLS),]
grid<-grid[!(grid$Cat.New=="Multi-level" & grid$Cat.Classic=="Single-level"),]
grid<-grid[!(grid$Cat.New=="Cliquish" & grid$Cat.Classic=="Multi-level"),]
grid<-grid[!(grid$Cat.New=="Cohesive" & grid$Cat.Classic=="Multi-level"),]

## MAKE DENSITY PLOT
means = grid %>%
  add_epred_draws(brm_density1, ndraws=100, allow_new_levels=T, re_formula = NA)

check<-means %>%
  group_by(GroupSize, Cat.Classic, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

range1<-means %>%
  group_by(Cat.Classic, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-means %>%
  group_by(GroupSize, Cat.Classic) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

check$type.draw<-paste(check$Cat.Classic, check$.draw, sep="_")

check$Cat.Classic<-factor(check$Cat.Classic, levels=c("Single-level","Multi-level"))
caps_groups$Cat.Classic<-factor(caps_groups$Cat.Classic, levels=c("Single-level","Multi-level"))

branded_colors <- c("#4D759C","#9c68b3", "#66a182")
branded_colors2<-c("#7373A4", "#66a182")
plot2a<-ggplot() +  
  scale_color_manual(values=branded_colors2) + labs(size="Observation hours\nper individidual") +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=type.draw, color=Cat.Classic), data = check, alpha = 0.1, lwd=1.2) +
  geom_point(data=caps_groups, aes(y=Density,x=GroupSize,size=Sqrt.Effort, color=Cat.Classic), alpha=0.35) +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=Cat.Classic), data = overall, alpha = 1, lwd=1.2, color = "black") +
  ylab("Density") + xlab("") + theme(legend.position ="none") + ylim(0,1) + xlim(0,80)
plot2a

## PLOT WITH SLOPES -- MODULARITY
means = grid %>%
  add_epred_draws(brm_modularity2, ndraws=100, allow_new_levels=T, re_formula = NA)

check<-means %>%
  group_by(GroupSize, Cat.New, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

range2<-means %>%
  group_by(Cat.New, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall<-means %>%
  group_by(GroupSize, Cat.New) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

check$type.draw<-paste(check$Cat.New, check$.draw, sep="_")

RColorBrewer::brewer.pal(8, "PRGn")

plot2c<-ggplot() + 
  scale_color_manual(values=branded_colors) + labs(size="Observation hours\nper individidual") +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=type.draw, color=Cat.New), data = check, alpha = 0.1, lwd=1.2) +
  geom_point(data=caps_groups, aes(y=Modularity,x=GroupSize,size=Sqrt.Effort, color=Cat.New),alpha=0.35) +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=Cat.New), data = overall, alpha = 1, lwd=1.2, color = "black") +
  ylab("Modularity") + xlab("Network size") + theme(legend.position ="none") + ylim(0,1) + xlim(0,80)
plot2c

## QUASIRANDOM PLOTS
## DENSITY
summary.range1<-range1 %>%
  group_by(Cat.Classic) %>%
  dplyr::summarise(Median=median(Mean.Pred), LL=quantile(Mean.Pred,0.055), UL=quantile(Mean.Pred,0.945), LL.75=quantile(Mean.Pred,0.125), UL.75=quantile(Mean.Pred,0.875))

summary.range1$Cat.Classic<-as.factor(summary.range1$Cat.Classic)
summary.range1$Cat.Classic<-factor(summary.range1$Cat.Classic, levels=c("Single-level","Multi-level"))
caps_groups$Cat.Classic<-factor(caps_groups$Cat.Classic, levels=c("Single-level","Multi-level"))

plot2b<-ggplot() + 
  scale_color_manual(values=branded_colors2) + labs(size="Observation hours\nper individidual") +
  geom_quasirandom(data=caps_groups, aes(y=Density,x=Cat.Classic,size=Sqrt.Effort, color=Cat.Classic),alpha=0.35) +
  geom_point(data=summary.range1,aes(x=Cat.Classic,y=Median)) +
  geom_errorbar(data=summary.range1,aes(x=Cat.Classic,ymin=LL,ymax=UL), width=0.08, lwd=1) +
  ylab("") + xlab("") + theme(legend.position ="none") + ylim(0,1)
plot2b

## MODULARITY
summary.range2<-range2 %>%
  group_by(Cat.New) %>%
  dplyr::summarise(Median=median(Mean.Pred), LL=quantile(Mean.Pred,0.055), UL=quantile(Mean.Pred,0.945), LL.75=quantile(Mean.Pred,0.125), UL.75=quantile(Mean.Pred,0.875))

plot2d<-ggplot() + 
  scale_color_manual(values=branded_colors) + labs(size="Observation hours\nper individidual") +
  geom_quasirandom(data=caps_groups, aes(y=Modularity,x=Cat.New,size=Sqrt.Effort, color=Cat.New),alpha=0.35) +
  geom_point(data=summary.range2,aes(x=Cat.New,y=Median)) +
  geom_errorbar(data=summary.range2,aes(x=Cat.New,ymin=LL,ymax=UL), width=0.08, lwd=1) +
  ylab("") + xlab("Social structure") + theme(legend.position ="none") + ylim(0,1)
plot2d

## PLOT COMPLETE VERSION
design<-"ABCC
         DEFF"

theme_set(theme_classic(base_size = 16))
plot2a + plot2b + plot2c + plot2d + plot_annotation(tag_levels = "a") 
ggsave(file="Figure_2.jpg", units="cm", width=20, height=16, dpi=300)

## PART 2
## Female-female grooming ties were variably shaped by kinship, rank, and shared male ties

## FEMALE-FEMALE META-ANALYSES
## CLASSIC METHODS
## USING METAFOR
caps_groups$Scale.Group<-as.numeric(scale(sqrt(caps_groups$GroupSize)))

## KINSHIP
meta1=rma.mv(yi=ff.kinship.r, V=ff.kinship.se^2, mods= ~ Cat.New + Philopatry + Scale.Group -1,
             random = list (~1|species/population_id/group_id/group.year,
                            ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps_groups)
summary(meta1)

## RANK
meta2=rma.mv(yi=ff.rank.r, V=ff.rank.se^2, mods= ~ Cat.New + Philopatry + Scale.Group -1,
             random = list (~1|species/population_id/group_id/group.year,
                            ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps_groups)
summary(meta2)

meta3=rma.mv(ff.male.r, V=ff.male.se^2, mods= ~ Cat.New + Philopatry + Scale.Group -1,
             random = list (~1|species/population_id/group_id/group.year,
                            ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps_groups)
summary(meta3)

meta4=rma.mv(ff.groomee.r, V=ff.groomee.se^2, mods= ~ Cat.New + Philopatry + Scale.Group-1,
             random = list (~1|species/group_id/group.year,
                            ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps_groups)
summary(meta4)

## CLASSIC VISUALIZATIONS
plot3a<-orchaRd::orchard_plot(meta1, mod = "Cat.New", group = "species", 
                              xlab = "Kinship effect", k=F,g=F,legend.pos="none", twig.size = 0)
plot3a<-plot3a + ylim(-2,8) +scale_fill_manual(values=branded_colors) + scale_color_manual(values=branded_colors)
plot3a

plot3b<-orchaRd::orchard_plot(meta2, mod = "Cat.New", group = "species", 
                              xlab = "Rank similarity effect",k=F,g=F, legend.pos="none",twig.size = 0)
plot3b<-plot3b + scale_fill_manual(values=branded_colors) + scale_color_manual(values=branded_colors) + ylim(-1,4) 
plot3b 

plot3c<-orchaRd::orchard_plot(meta3, mod = "Cat.New", group = "species", 
                              xlab = "Shared male effect",k=F,g=F, twig.size=0, legend.pos="none")
plot3c<-plot3c + scale_fill_manual(values=branded_colors) + scale_color_manual(values=branded_colors) + ylim(-2,8) 
plot3c

## BAYESIAN
## META-ANALYSES
## BRMS

## RELATEDNESS META-ANALYTIC MODEL
caps_groups$Philopatry<-factor(caps_groups$Philopatry, levels=c("Female philopatric", "Female dispersal"))
relatedness.meta<-brm(data = caps_groups, family = "student", data2=list(A=A),
                      ff.kinship.r|se(ff.kinship.se) ~ 0 + Cat.New + scale(sqrt(GroupSize)) + Philopatry + (1|gr(phylo, cov=A)) +
                        (1|species/population_id/group_id/group.year),
                      prior = c(prior(normal(0, 2), class = b),
                                prior(cauchy(0, 0.3), class = sd)),
                      iter = 4000, warmup = 1000, cores = 2, chains = 4,
                      control=list(adapt_delta=0.99, max_treedepth=15))
summary(relatedness.meta, prob=0.89)
plot(conditional_effects(relatedness.meta, "Cat.New", prob=0.89), points=T)
plot(conditional_effects(relatedness.meta, "GroupSize:Cat.New", prob=0.89))
conditional_effects(relatedness.meta, "Cat.New:Philopatry", prob=0.89)

## RANK META-ANALYTIC MODEL
rank.meta<-brm(data = caps_groups, family = "student", data2=list(A=A),
               ff.rank.r| se(ff.rank.se) ~  0 + Cat.New + scale(sqrt(GroupSize)) + Philopatry + (1|gr(phylo, cov=A)) +
                 (1|species/population_id/group_id/group.year),
               prior = c(prior(normal(0, 2), class = b),
                         prior(cauchy(0, 0.3), class = sd)),
               iter = 4000, warmup = 1000, cores = 2, chains = 4,
               control=list(adapt_delta=0.99, max_treedepth=15))
summary(rank.meta, prob=0.89)
conditional_effects(rank.meta, "Cat.New", prob=0.89)
conditional_effects(rank.meta, "GroupSize", prob=0.89)

## GROOMEEE RANK META-ANALYTIC MODEL
groomee.meta<-brm(data = caps_groups, family = "student", data2=list(A=A),
                  ff.groomee.r| se(ff.groomee.se) ~  0 + Cat.New + scale(sqrt(GroupSize)) + Philopatry + (1|gr(phylo, cov=A)) +
                    (1|species/population_id/group_id/group.year),
                  prior = c(prior(normal(0, 2), class = b),
                            prior(cauchy(0, 0.3), class = sd)),
                  iter = 4000, warmup = 1000, cores = 2, chains = 4,
                  control=list(adapt_delta=0.99, max_treedepth=15))
summary(groomee.meta, prob=0.89)
conditional_effects(groomee.meta, "Cat.New", prob=0.89)
conditional_effects(groomee.meta, "GroupSize")
conditional_effects(groomee.meta, "GroupSize:Cat.New")
conditional_effects(groomee.meta, "Cat.New:Philopatry")

## SHARED MALE META-ANALYTIC MODEL
male.meta<-brm(data = caps_groups, family = "student", data2=list(A=A),
               ff.male.r|se(ff.male.se) ~  0 + Cat.New + scale(sqrt(GroupSize)) + Philopatry + (1|gr(phylo, cov=A)) +
                 (1|species/population_id/group_id/group.year),
               prior = c(prior(normal(0, 2), class = b),
                         prior(cauchy(0, 0.3), class = sd)),
               iter = 4000, warmup = 1000, cores = 2, chains = 4,
               control=list(adapt_delta=0.99, max_treedepth=15))
summary(male.meta, prob=0.89)
conditional_effects(male.meta, "Cat.New")
conditional_effects(male.meta, "GroupSize")
conditional_effects(male.meta, "GroupSize:Cat.New")

## PART 3
## Female-male grooming ties were variably shaped by social dominance

## FEMALE-MALE META-ANALYSES
## CLASSIC METHODS
## USING METAFOR
meta5=rma.mv(yi=fm.rank.male.r, V=fm.rank.male.se^2, mods= ~ Cat.New + Scale.Group - 1,
             random = list(~1|species/population_id/group_id,
                           ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps_groups)
summary(meta5)

meta6=rma.mv(fm.rank.fem.r, V=fm.rank.fem.se^2, mods= ~ Cat.New + Scale.Group - 1,
             random = list (~1|species/population_id/group_id/group.year,
                            ~1|phylo),
             R=list(phylo=A),
             data = caps_groups)
summary(meta6)

meta7=rma.mv(fm.rank.intx.r, V=fm.rank.intx.se^2, mods= ~ Cat.New + Scale.Group - 1 ,
             random = list (~1|species/population_id/group_id,
                            ~1|phylo),
             R=list(phylo=A),
             data = caps_groups)
summary(meta7)

meta8=rma.mv(fm.sex.r, V=fm.sex.se^2, mods= ~ Cat.New + Scale.Group - 1,
             random = list (~1|species/population_id/group_id/group.year,
                            ~1|phylo),
             R=list(phylo=A),
             data = caps_groups)
summary(meta8)

## PLOT IT
model_results <- orchaRd::mod_results(meta5, mod = "Cat.New", at = NULL, group = "species")
plot3d<-orchaRd::orchard_plot(meta5, mod = "Cat.New", group = "species", twig.size=0, 
                              xlab = "Male rank effect", k=F,g=F, legend.pos="none")

plot3d<-plot3d+scale_fill_manual(values=branded_colors) + scale_color_manual(values=branded_colors) + ylim(-1,4)
plot3d
plot3e<-orchaRd::orchard_plot(meta6, mod = "Cat.New", group = "species", twig.size=0,
                              xlab = "Female rank effect",k=F,g=F, legend.pos="none")

plot3e<-plot3e+scale_fill_manual(values=branded_colors) + scale_color_manual(values=branded_colors) + ylim(-1,4)
plot3e

plot3f<-orchaRd::orchard_plot(meta7, mod = "Cat.New", group = "species", twig.size = 0,
                              xlab = "Rank interaction effect",k=F,g=F, legend.pos="none")

plot3f<-plot3f+scale_fill_manual(values=branded_colors) + scale_color_manual(values=branded_colors) + ylim(-1,4)
plot3f

layout1<-"AD
          BE
          CF"

theme_set(theme_test(base_size = 16))

plot3a+plot3b+plot3c+plot3d+plot3e+plot3f + plot_annotation(tag_levels = "a") + plot_layout(design=layout1)
ggsave(file="Figure_S6.jpg", units="cm", width=22, height=22, dpi=300)

## MALE RANK FM META-ANALYTIC MODEL
male.effect.fm<-brm(data = caps_groups, data2=list(A=A),
                    fm.rank.male.r|se(fm.rank.male.se) ~ 0 + Cat.New + scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                      (1|species/population_id/group_id/group.year),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(cauchy(0, 0.3), class = sd)),    
                    family=student,
                    iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                    control=list(adapt_delta=0.99, max_treedepth=15))
summary(male.effect.fm, prob=0.89)
plot(conditional_effects(male.effect.fm, "Cat.New", prob=0.89), points=T)

## FEMALE RANK FM META-ANALYTIC MODEL
caps_groups$ASR<-caps_groups$N.Fem/caps_groups$N.Male
fem.effect.fm<-brm(data = caps_groups, family = student, data2=list(A=A),
                   fm.rank.fem.r| se(fm.rank.fem.se) ~  0 + Cat.New + scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                     (1|species/population_id/group_id/group.year),
                   prior = c(prior(normal(0, 2), class = b),
                             prior(cauchy(0, 0.3), class = sd)),    
                   iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                   control=list(adapt_delta=0.99, max_treedepth=15))
summary(fem.effect.fm, prob=0.89)
conditional_effects(fem.effect.fm, "Cat.New", prob=0.89)
conditional_effects(fem.effect.fm, "GroupSize", prob=0.89)

## INTERACTION RANK FM META-ANALYTIC MODEL
intx.effect.fm<-brm(data = caps_groups, family = student, data2=list(A=A),
                    fm.rank.intx.r|se(fm.rank.intx.se) ~ 0 + Cat.New + scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                      (1|species/population_id/group_id/group.year),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(cauchy(0, 0.5), class = sd)),
                    iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                    control=list(adapt_delta=0.99, max_treedepth=15))
summary(intx.effect.fm, prob=0.89)
conditional_effects(intx.effect.fm, "Cat.New", prob=0.89)
conditional_effects(intx.effect.fm, "GroupSize")

## SEX EFFECT
sex.effect.fm<-brm(data = caps_groups, family = student, data2=list(A=A),
                   fm.sex.r| se(fm.sex.se) ~  0 + Cat.New +  scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                     (1|species/population_id/group_id/group.year),
                   prior = c(prior(normal(0, 2), class = b),
                             prior(cauchy(0, 0.3), class = sd)),    
                   iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                   control=list(adapt_delta=0.99, max_treedepth=15))
summary(sex.effect.fm)
plot(conditional_effects(sex.effect.fm, "Cat.New"), points=T)
conditional_effects(sex.effect.fm, "GroupSize")

## EXTRACT PLOTS - FF DYADS
effects1<-as.data.frame(conditional_effects(relatedness.meta, "Cat.New", prob=0.89)$Cat.New)
effects2<-as.data.frame(conditional_effects(rank.meta, "Cat.New", prob=0.89)$Cat.New)
effects3<-as.data.frame(conditional_effects(male.meta, "Cat.New", prob=0.89)$Cat.New)
effects4<-as.data.frame(conditional_effects(groomee.meta, "Cat.New", prob=0.89)$Cat.New)

## ADD PHILOPATRY INFO
effects.intx<-as.data.frame(conditional_effects(relatedness.meta, "Cat.New:Philopatry", prob=0.89)$Cat.New)
effects.intx<-effects.intx[!(effects.intx$Cat.New=="Cliquish" & effects.intx$Philopatry=="Female dispersal"),]
effects.intx<-effects.intx[!(effects.intx$Cat.New=="Cohesive" & effects.intx$Philopatry=="Female dispersal"),]
effects.intx2<-as.data.frame(conditional_effects(rank.meta, "Cat.New:Philopatry", prob=0.89)$Cat.New)
effects.intx2<-effects.intx2[!(effects.intx2$Cat.New=="Cliquish" & effects.intx2$Philopatry=="Female dispersal"),]
effects.intx2<-effects.intx2[!(effects.intx2$Cat.New=="Cohesive" & effects.intx2$Philopatry=="Female dispersal"),]

## CATEGORY PLOTS
library(ggbeeswarm)
plot3a<-ggplot() + 
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=ff.kinship.r, size=1/(ff.kinship.se), color=Cat.New, shape=Philopatry), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  geom_point(data=effects.intx, aes(x=Cat.New, y=estimate__,shape=Philopatry),
             position=position_dodge(width=0.5), size=3) +
  geom_errorbar(data=effects.intx, aes(x=Cat.New, group=Philopatry,ymin=lower__, ymax=upper__), 
                position=position_dodge(width=0.5), width=0, lwd=1) +
  ggtitle("Female-female dyads") +
  xlab("") + ylab("Kinship effect") + theme(legend.position="none",plot.title=element_text(hjust=0.5, face="bold")) + 
  scale_y_continuous(limits=c(-2,6), breaks=seq(-2,6, by=2)) + coord_flip()

plot3b<-ggplot() + 
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=ff.rank.r, size=1/(ff.rank.se), color=Cat.New,shape=Philopatry), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  geom_point(data=effects2, aes(x=Cat.New, y=estimate__), size = 3) +
  geom_errorbar(data=effects2, aes(x=Cat.New, ymax=lower__, ymin=upper__), width=0, lwd=1) +
  xlab("") + ylab("Rank similarity effect") + theme(legend.position="none") + ylim(-1,3) + coord_flip()

plot.groomee<-ggplot() + 
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=ff.groomee.r, size=1/(ff.groomee.se), color=Cat.New,shape=Philopatry), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  geom_point(data=effects4, aes(x=Cat.New, y=estimate__),size = 3) +
  geom_errorbar(data=effects4, aes(x=Cat.New, ymax=lower__, ymin=upper__), width=0, lwd=1) +
  xlab("") + ylab("Partner rank effect") + theme(legend.position="none") + ylim(-1,3) + coord_flip()

plot3c<-ggplot() + 
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=ff.male.r, size=1/(ff.male.se), color=Cat.New,shape=Philopatry), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  geom_point(data=effects3, aes(x=Cat.New, y=estimate__), size=3) +
  geom_errorbar(data=effects3, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  xlab("") + ylab("Shared male effect") + theme(legend.position="none") + 
  scale_y_continuous(limits=c(-2.1,6.3), breaks=seq(-2,6, by=2)) +coord_flip()

## EXTRACT PLOTS - FF DYADS
effects1<-as.data.frame(conditional_effects(male.effect.fm, "Cat.New", prob=0.89)$Cat.New)
effects2<-as.data.frame(conditional_effects(fem.effect.fm, "Cat.New", prob=0.89)$Cat.New)
effects3<-as.data.frame(conditional_effects(intx.effect.fm, "Cat.New", prob=0.89)$Cat.New)
effects4<-as.data.frame(conditional_effects(sex.effect.fm, "Cat.New", prob=0.89)$Cat.New)

plot3d<-ggplot() + 
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=fm.rank.male.r, size=1/(fm.rank.male.se), color=Cat.New,shape=Philopatry), stroke=0, alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  geom_point(data=effects1, aes(x=Cat.New, y=estimate__), size=3) + 
  geom_errorbar(data=effects1, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  ggtitle("Female-male dyads") + ylim(-1,3) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("") + ylab("Male rank effect") + theme(legend.position="none",plot.title=element_text(hjust=0.5, face="bold")) + coord_flip()
plot3d

plot3e<-ggplot() +
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=fm.rank.fem.r, size=1/(fm.rank.fem.se),shape=Philopatry, color=Cat.New), stroke=0, alpha=0.5) + 
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  scale_size(range=c(1,4)) +
  geom_point(data=effects2, aes(x=Cat.New, y=estimate__), size=3) +
  geom_errorbar(data=effects2, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("") + ylab("Female rank effect") + theme(legend.position="none") + ylim(-1,3) + coord_flip()

plot3f<-ggplot() + 
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=fm.rank.intx.r, shape=Philopatry,size=1/(fm.rank.intx.se), color=Cat.New), stroke=0, alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  geom_point(data=effects3, aes(x=Cat.New, y=estimate__), size=3) +
  geom_errorbar(data=effects3, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("") + ylab("Rank interaction effect") + theme(legend.position="none") + ylim(-1,3) + coord_flip()

plot.sex<-ggplot() + 
  geom_quasirandom(data=caps_groups, aes(x=Cat.New,y=-fm.sex.r, size=(1/fm.sex.se),shape=Philopatry, color=Cat.New), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=branded_colors) +
  geom_point(data=effects4, aes(x=Cat.New, y=-estimate__)) +
  geom_errorbar(data=effects4, aes(x=Cat.New, ymin=-upper__, ymax=-lower__), width=0, lwd=1) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  xlab("") + ylab("Sex effect") + theme(legend.position="none") + ylim(-1,3) + coord_flip()

layout1<-"AD
          BE
          CF"
plot3a + plot3b + plot3c +
  plot3d + plot3e + plot3f + plot_annotation(tag_levels = "a") + plot_layout(design = layout1)
ggsave(file="Figure_3.jpg", units="cm", width=25, height=25, dpi=300)

## PART 4
## Nepotistic biases and shared male effects are independently linked with network structure

## DYADIC BIASES AND NETWORK DENSITY
## FF EFFECTS
brm_density_phylo_eff<-brm(data=caps_groups, data2=list(A=A),
                           N.Groomed|trials(N.Dyads) ~ Cat.Classic + sqrt(GroupSize) + me(ff.kinship.r, ff.kinship.se) + me(ff.rank.r,ff.rank.se) + me(ff.groomee.r,ff.groomee.se) + me(ff.male.r,ff.male.se)  + 
                             scale(Obs.Effort.Per.Capita) +
                             (1|gr(phylo, cov=A)) +
                             (1|species/population_id/group_id),
                           prior = c(
                             prior(normal(0, 2), "b")),
                           family="beta_binomial", iter=4000, chains=4, cores=2,
                           control=list(adapt_delta=0.99, max_treedepth=15))
summary(brm_density_phylo_eff, prob=0.89)

plots2<-conditional_effects(brm_density_phylo_eff, "ff.rank.r", prob=0.89)
plots3<-conditional_effects(brm_density_phylo_eff, "ff.male.r", prob=0.89)
plots1<-conditional_effects(brm_density_phylo_eff, "ff.kinship.r", prob = 0.89)

plot.related<-plot(plots1,plot = FALSE, prob=0.89, line_args = list(color="black"),)[[1]] + ylim(0,1) + xlab("Kinship effect") + ylab("Density") + scale_fill_manual(values="black")
plot.rank<-plot(plots2, plot = FALSE, prob=0.89, line_args = list(color="black"))[[1]] + ylim(0,1) + xlab("Rank similarity effect") + ylab("")
plot.male<-plot(plots3, plot = FALSE, prob=0.89, line_args = list(color="black"))[[1]] + ylim(0,1) + xlab("Shared male effect") + ylab("")

library(patchwork)
plot.related + plot.rank + plot.male + plot_annotation(tag_levels = "a")
ggsave(file="Figure_S2.jpg", units="cm", width=28, height=10, dpi=300)

## DYADIC BIASES AND NETWORK MODULARITY
## FF EFFECTS
brm_modularity_phylo_eff<-brm(data=caps_groups, data2=list(A=A),
                              Modularity ~  scale(sqrt(GroupSize)) + Cat.Classic + me(ff.kinship.r, ff.kinship.se) + me(ff.groomee.r, ff.groomee.se) + me(ff.rank.r,ff.rank.se) + me(ff.male.r,ff.male.se) + scale(Obs.Effort.Per.Capita) +
                                (1|gr(phylo, cov=A)) +  
                                (1|species/population_id/group_id),
                              prior = c(
                                prior(normal(0, 2), "b")),
                              family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                              control=list(adapt_delta=0.99, max_treedepth=12))
summary(brm_modularity_phylo_eff, prob=0.89)
plot2<-conditional_effects(brm_modularity_phylo_eff, "ff.rank.r:Cat.Classic", prob=0.89)
plot3<-conditional_effects(brm_modularity_phylo_eff, "ff.groomee.r:Cat.Classic", prob=0.89)
plot4<-conditional_effects(brm_modularity_phylo_eff, "ff.male.r:Cat.Classic", prob=0.89)
plot1<-conditional_effects(brm_modularity_phylo_eff, "ff.kinship.r:Cat.Classic", prob = 0.89)

plot.related<-plot(plot1,plot = FALSE, prob=0.89, line_args = list(color="black"))[[1]] + ylim(0,1) + xlab("Kinship effect")
plot.related<-plot.related + geom_point(data=caps_groups, aes(x=ff.kinship.r, y=Modularity, color=Cat.Classic, size=(1/related.se)), inherit.aes = F, alpha=0.7) +
  scale_color_manual(values=branded_colors2) + scale_fill_manual(values=branded_colors2) + theme(legend.position = "none") + xlim(-1,6) + scale_size(range=c(1,4)) 
plot.rank<-plot(plot2, plot = FALSE, prob=0.89, line_args = list(color="black"))[[1]] + ylim(0,1) + xlab("Rank similarity effect") + ylab("")
plot.rank<-plot.rank + geom_point(data=caps_groups, aes(x=ff.rank.r, y=Modularity, color=Cat.Classic, size=(1/rank.se)), inherit.aes = F, alpha=0.7) +
  scale_color_manual(values=branded_colors2) + scale_fill_manual(values=branded_colors2) + theme(legend.position = "none") + scale_size(range=c(1,4)) 
plot.groomee<-plot(plot3, plot = FALSE, prob=0.89, line_args = list(color="black"))[[1]] + ylim(0,1) + xlab("Partner rank effect") + ylab("")
plot.groomee<-plot.groomee + geom_point(data=caps_groups, aes(x=ff.groomee.r, y=Modularity, color=Cat.Classic, size=(1/groomee.se)), inherit.aes = F, alpha=0.7) +
  scale_color_manual(values=branded_colors2) + scale_fill_manual(values=branded_colors2) + theme(legend.position = "none") +  scale_size(range=c(1,4)) 
plot.male<-plot(plot4, plot = FALSE, prob=0.89, line_args = list(color="black"))[[1]] + ylim(0,1) + xlab("Shared male effect") + ylab("")
plot.male<-plot.male + geom_point(data=caps_groups, aes(x=ff.male.r, y=Modularity, color=Cat.Classic, size=(1/male.se)), inherit.aes = F, alpha=0.7) +
  scale_color_manual(values=branded_colors2) + scale_fill_manual(values=branded_colors2) + theme(legend.position = "none") +  scale_size(range=c(1,4)) 

## RELATED PLOT TRUNCATED
plot.related<-plot.related$data
related.min.ML<-min(caps_groups$ff.kinship.r[caps_groups$Cat.Classic=="Multi-level"], na.rm = T)
related.min.SL<-min(caps_groups$ff.kinship.r[caps_groups$Cat.Classic=="Single-level"], na.rm=T)
related.max.ML<-max(caps_groups$ff.kinship.r[caps_groups$Cat.Classic=="Multi-level"], na.rm = T)
related.max.SL<-max(caps_groups$ff.kinship.r[caps_groups$Cat.Classic=="Single-level"],na.rm=T)
plot.related<-plot.related[!(plot.related$ff.kinship.r>related.max.SL & plot.related$Cat.Classic=="Single-level"),]
plot.related<-plot.related[!(plot.related$ff.kinship.r<related.min.SL & plot.related$Cat.Classic=="Single-level"),]
plot.related<-plot.related[!(plot.related$ff.kinship.r>related.max.ML & plot.related$Cat.Classic=="Multi-level"),]
plot.related<-plot.related[!(plot.related$ff.kinship.r<related.min.ML & plot.related$Cat.Classic=="Multi-level"),]

plot.related<-ggplot() + geom_point(data=caps_groups, aes(x=ff.kinship.r, y=Modularity,size=(1/ff.kinship.se), color=Cat.Classic, alpha=0.7)) +
  scale_color_manual(values=branded_colors2) + scale_fill_manual(values=branded_colors2) + theme(legend.position = "none") +  scale_size(range=c(1,4)) +
  geom_ribbon(data=plot.related, aes(x=ff.kinship.r, ymin=lower__, ymax=upper__, fill=Cat.Classic, alpha=0.7)) + ylim(0,1) + xlab("Kinship effect") + ylab("Modularity") +
  geom_line(data=plot.related, aes(x=ff.kinship.r, y=estimate__, group=Cat.Classic, alpha=0.7)) 

## RANK PLOT TRUNCATED
plot.rank<-plot.rank$data
rank.min.ML<-min(caps_groups$ff.rank.r[caps_groups$Cat.Classic=="Multi-level"],na.rm=T)
rank.min.SL<-min(caps_groups$ff.rank.r[caps_groups$Cat.Classic=="Single-level"],na.rm=T)
rank.max.ML<-max(caps_groups$ff.rank.r[caps_groups$Cat.Classic=="Multi-level"],na.rm=T)
rank.max.SL<-max(caps_groups$ff.rank.r[caps_groups$Cat.Classic=="Single-level"],na.rm=T)
plot.rank<-plot.rank[!(plot.rank$ff.rank.r>rank.max.SL & plot.rank$Cat.Classic=="Single-level"),]
plot.rank<-plot.rank[!(plot.rank$ff.rank.r<rank.min.SL & plot.rank$Cat.Classic=="Single-level"),]
plot.rank<-plot.rank[!(plot.rank$ff.rank.r>rank.max.ML & plot.rank$Cat.Classic=="Multi-level"),]
plot.rank<-plot.rank[!(plot.rank$ff.rank.r<rank.min.ML & plot.rank$Cat.Classic=="Multi-level"),]

plot.rank<-ggplot() + geom_point(data=caps_groups, aes(x=ff.rank.r, y=Modularity,size=(1/ff.rank.se), color=Cat.Classic, alpha=0.7)) +
  scale_color_manual(values=branded_colors2) + scale_fill_manual(values=branded_colors2) + theme(legend.position = "none") +  scale_size(range=c(1,4)) +
  geom_ribbon(data=plot.rank, aes(x=ff.rank.r, ymin=lower__, ymax=upper__, fill=Cat.Classic, alpha=0.7)) + ylim(0,1) + xlab("Rank similarity effect") + ylab("") +
  geom_line(data=plot.rank, aes(x=ff.rank.r, y=estimate__, group=Cat.Classic, alpha=0.7)) 

## MALE PLOT TRUNCATED
plot.male<-plot.male$data
male.min.ML<-min(caps_groups$ff.male.r[caps_groups$Cat.Classic=="Multi-level"], na.rm=T)
male.min.SL<-min(caps_groups$ff.male.r[caps_groups$Cat.Classic=="Single-level"], na.rm=T)
male.max.ML<-max(caps_groups$ff.male.r[caps_groups$Cat.Classic=="Multi-level"], na.rm=T)
male.max.SL<-max(caps_groups$ff.male.r[caps_groups$Cat.Classic=="Single-level"], na.rm=T)
plot.male<-plot.male[!(plot.male$ff.male.r>male.max.SL & plot.male$Cat.Classic=="Single-level"),]
plot.male<-plot.male[!(plot.male$ff.male.r<male.min.SL & plot.male$Cat.Classic=="Single-level"),]
plot.male<-plot.male[!(plot.male$ff.male.r>male.max.ML & plot.male$Cat.Classic=="Multi-level"),]
plot.male<-plot.male[!(plot.male$ff.male.r<male.min.ML & plot.male$Cat.Classic=="Multi-level"),]

plot.male<-ggplot() + geom_point(data=caps_groups, aes(x=ff.male.r, y=Modularity,size=(1/ff.male.se), color=Cat.Classic, alpha=0.7)) +
  scale_color_manual(values=branded_colors2) + scale_fill_manual(values=branded_colors2) + theme(legend.position = "none") +  scale_size(range=c(1,4)) +
  geom_ribbon(data=plot.male, aes(x=ff.male.r, ymin=lower__, ymax=upper__, fill=Cat.Classic, alpha=0.7)) + ylim(0,1) + xlab("Shared male effect") + ylab("") +
  geom_line(data=plot.male, aes(x=ff.male.r, y=estimate__, group=Cat.Classic, alpha=0.7))

plot.related + plot.rank + plot.male + plot_annotation(tag_levels = "a")
ggsave(file="Figure_4.jpg", units="cm", width=27, height=9, dpi=300)

## DYADIC BIASES AND NETWORK DENSITY
## FM EFFECTS
brm_density_phylo_eff2<-brm(data=caps_groups, data2=list(A=A),
                            N.Groomed|trials(N.Dyads) ~ Cat.Classic + sqrt(GroupSize) + me(fm.rank.fem.r, fm.rank.fem.se) + me(fm.rank.male.r, fm.rank.male.se) + me(fm.rank.intx.r, fm.rank.intx.se) +
                              me(fm.sex.r, fm.sex.se) + scale(Obs.Effort.Per.Capita) +
                              (1|gr(phylo, cov=A)) +  
                              (1|species/population_id/group_id),
                            prior = c(
                              prior(normal(0, 2), "b")),
                            family="beta_binomial", iter=4000, chains=4, cores=2,
                            control=list(adapt_delta=0.99, max_treedepth=12))
summary(brm_density_phylo_eff2)

## DYADIC BIASES AND NETWORK MODULARITY
## FM EFFECTS
brm_modularity_phylo_eff2<-brm(data=caps_groups, data2=list(A=A),
                               Modularity ~ Cat.Classic + sqrt(GroupSize) + me(fm.rank.fem.r, fm.rank.fem.se) + me(fm.rank.male.r, fm.rank.male.se) + me(fm.rank.intx.r, fm.rank.intx.se) +
                                 me(fm.sex.r, fm.sex.se) + scale(Obs.Effort.Per.Capita) +
                                 (1|gr(phylo, cov=A)) +  
                                 (1|species/population_id/group_id),
                               prior = c(
                                 prior(normal(0, 2), "b")),
                               family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                               control=list(adapt_delta=0.99, max_treedepth=12))
summary(brm_modularity_phylo_eff2)

## GET INTRACLASS CORRELATIONS
## FF DYADS
phylo.relat<-as.data.frame(icc(relatedness.meta, by_group = T))
phylo.relat$Effect<-"Kinship"
phylo.rank<-as.data.frame(icc(rank.meta, by_group = T))
phylo.rank$Effect<-"Rank similarity"
phylo.groomee<-as.data.frame(icc(groomee.meta, by_group = T))
phylo.groomee$Effect<-"Partner's rank"
phylo.male<-as.data.frame(performance::icc(male.meta, by_group = T))
phylo.male$Effect<-"Shared male"

phylo.ests<-rbind(phylo.relat, phylo.male, phylo.groomee, phylo.rank)
phylo.ests$Effect<-as.factor(phylo.ests$Effect)
phylo.ests$Effect<-factor(phylo.ests$Effect, levels=c("Kinship","Rank similarity","Partner's rank", "Shared male"))
phylo.ests$Group<-as.factor(phylo.ests$Group)
levels(phylo.ests$Group)<-c("Phylogeny","Species","Population","Group","Residual")

## FM DYADS
phylo.fm.male<-as.data.frame(icc(male.effect.fm, by_group = T))
phylo.fm.male$Effect<-"Male rank"
phylo.fm.fem<-as.data.frame(icc(fem.effect.fm, by_group = T))
phylo.fm.fem$Effect<-"Female rank"
phylo.intx.fem<-as.data.frame(icc(intx.effect.fm, by_group = T))
phylo.intx.fem$Effect<-"Rank interaction"
phylo.fm.sex<-as.data.frame(icc(sex.effect.fm, by_group = T, robust=T))
phylo.fm.sex$Effect<-"Sex"

phylo.ests2<-rbind(phylo.fm.male, phylo.fm.fem, phylo.intx.fem, phylo.fm.sex)
phylo.ests2$Effect<-as.factor(phylo.ests2$Effect)
phylo.ests2$Effect<-factor(phylo.ests2$Effect, levels=c("Male rank","Female rank","Rank interaction", "Sex"))
phylo.ests2$Group<-as.factor(phylo.ests2$Group)
unique(phylo.ests2$Group)
levels(phylo.ests2$Group)<-c("Phylogeny","Species","Population","Group","Residual")

phylo.plot1<-ggplot(data=phylo.ests, aes(x=Effect, y=ICC,fill=Group)) + geom_col(alpha=0.8, col="black") + theme(legend.position = "none") +
  scale_fill_brewer(palette="Accent") + ylab("Intraclass correlation coefficient") + xlab("Female-female effects")

phylo.plot2<-ggplot(data=phylo.ests2, aes(x=Effect, y=ICC,fill=Group)) + geom_col(alpha=0.8, col="black") + 
  scale_fill_brewer(palette="Accent") + ylab("Intraclass correlation coefficient") + xlab("Female-male effects")

layout1<-"AAAA
          BBBB"
phylo.plot1+phylo.plot2 + plot_annotation(tag_levels = "a") + plot_layout(design = layout1)
ggsave(file="Figure_S3.jpg", units="cm", width=20, height=25, dpi=300)


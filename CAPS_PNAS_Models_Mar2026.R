## DATA AND CODE FOR FEDER ET AL. 2025
## COMPARATIVE ANALYSIS OF PAPION SOCIETIES
## Disparate social structures are underpinned by distinct social rules
## across a primate radiation

## LOAD NECESSARY PACKAGES
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
library(ggtree)
theme_set(theme_classic(base_size = 16))

## SETWD
setwd("~/Desktop/SPRF/PNAS/GitHub/NEW DATA")

## SET PALETTES
pal1 <- c("#4D759C","#9c68b3", "#66a182")
pal2<-c("#7373A4", "#66a182")

## GET GROUP-LEVEL DATA
caps.groups<-readRDS("CAPS_Groups.rds")

## ADD FACTOR FOR SAMPLING STYLE
i<-1
for (i in 1:length(caps.groups)) {
  caps.groups[[i]]$SamplingType<-"Focal"
  caps.groups[[i]]$SamplingType[caps.groups[[i]]$population_id=="AMBOSELI"|caps.groups[[i]]$population_id=="GASHAKA"]<-"Representative"
}

## FIX LEVELS
i<-1
for (i in 1:length(caps.groups)) {
  caps.groups[[i]]$Cat.Classic<-factor(caps.groups[[i]]$Cat.Classic, levels=c("Single-level","Multi-level"))
  caps.groups[[i]]$Cat.New<-factor(caps.groups[[i]]$Cat.New, levels=c("Cohesive","Cliquish", "Multi-level"))
}

## LOAD TREE FILE FROM KUDERNA ET AL. 2023
treefile<-ape::read.tree("Kuderna_Tree2.nex.tree")
phylo1<-treefile
phylo1$edge.length

## TRIM TO INCLUDE SAMPLED TAXA
## AND CONGENERIC COGNATES
keep<-unique(caps.groups[[1]]$phylo)
phylo1<-keep.tip(phylo1, keep)
phylo.check<-phylo1
plot.phylo(phylo.check,  edge.width = 2,cex=1, label.offset = 0)
A <- ape::vcv.phylo(phylo1, corr=T)
colnames(caps.groups)

## PART 1
## Papionin grooming networks vary as a function of network size and social system

## MAKE PLOT WITH TREE
tips<-as.character(unique(caps.groups[[1]]$phylo))
tips.kept<-keep.tip(phylo1, tips)
tips.kept$tip.label
tips.kept$tip.label<-c("Gelada","Gray-cheeked mangabey", "Hamadryas baboon", "Olive baboon","Guinea baboon",
                       "Yellow baboon", "Chacma baboon","Kinda baboon", "Mandrill", "Sanje mangabey", "Sooty mangabey")
tips.kept$Cat.New<-c("Multi-level","Cohesive", "Multi-level", "Cohesive","Multi-level",
                  "Cohesive","Cliquish","Cliquish","Cliquish","Cohesive","Cohesive")

pal1[2]
gg_tr <- ggtree(tips.kept)  + geom_tiplab(color=c(pal1[3],pal1[1],pal1[3],
                                                  pal1[1],pal1[3],pal1[1],
                                                  pal1[2],pal1[2],pal1[2],
                                                  pal1[1],pal1[1]),
                                          align=T, size=5)# make more room for the labels
gg_tr<- gg_tr + xlim(0,28)
gg_tr

## DENSITY MODEL 1 -- NEW SOCIAL CATEGORIES
repeated_list <- rep(list(list(A=A)),20)
repeated_list
brm_density1<-brm_multiple(data=caps.groups, data2=repeated_list,
                        Density ~ scale(sqrt(GroupSize)) + Cat.New +
                         (1|gr(phylo, cov=A)) +
                         (1|Species/population_id/group_id/group.year),
                        prior = c(
                          prior(normal(0,0.5), "Intercept"),
                          prior(student_t(3,0,10), "sd"),
                          prior(gamma(0.1, 0.1), "phi"),
                          prior(normal(0, 2), "b")),    
                       family="zero_one_inflated_beta", iter=4000, chains=4, cores=2,
                       control=list(adapt_delta=0.99, max_treedepth=15))
summary(brm_density1, prob=0.89)
conditional_effects(brm_density1, "Cat.New")
conditional_effects(brm_density1, "GroupSize:Cat.New")

## DENSITY MODEL 2 -- OLD SOCIAL CATEGORIES
brm_density2<-brm_multiple(data=caps.groups, data2=repeated_list,
                    Density ~ scale(sqrt(GroupSize)) + Cat.Classic +
                        (1|gr(phylo, cov=A)) +
                        (1|Species/population_id/group_id/group.year),
                  prior = c(
                    prior(normal(0,0.5), "Intercept"),
                    prior(student_t(3,0,10), "sd"),
                    prior(gamma(0.1, 0.1), "phi"),
                    prior(normal(0, 2), "b")),    
                       family="zero_one_inflated_beta", iter=4000, chains=4, cores=2,
                       control=list(adapt_delta=0.99, max_treedepth=15))
summary(brm_density2, prob=0.89)
conditional_effects(brm_density2, "Cat.Classic")
conditional_effects(brm_density2, "GroupSize:Cat.Classic")

## DENSITY MODEL 3 -- GROUP SIZE ALONE
brm_density3<-brm_multiple(data=caps.groups, data2=repeated_list,
                       Density ~ scale(sqrt(GroupSize))+
                       (1|gr(phylo, cov=A)) +
                       (1|Species/population_id/group_id/group.year),
                       prior = c(
                       prior(normal(0,0.5), "Intercept"),
                       prior(student_t(3,0,10), "sd"),
                       prior(gamma(0.1, 0.1), "phi"),
                       prior(normal(0, 2), "b")
                       ),                      
                       family="zero_one_inflated_beta", iter=4000, chains=4, cores=2,
                       control=list(adapt_delta=0.99, max_treedepth=15))
summary(brm_density3, prob=0.89)
conditional_effects(brm_density3, "GroupSize")

## DENSITY -- COMPARE THREE
loo1<-brms::loo(brm_density1, pointwise = T)
loo2<-brms::loo(brm_density2, pointwise = T)
loo_compare(loo1,loo2)
waic(brm_density1)
waic(brm_density2)
model_weights(brm_density1, brm_density2, weights="loo")

as.data.frame(icc(brm_density3, by_group = T))

## MODULARITY MODEL 1 -- CLASSIC SOCIAL CATEGORIES
brm_modularity1<-brm_multiple(data=caps.groups, data2=repeated_list,
                          Modularity ~ scale(sqrt(GroupSize)) + Cat.New +
                            (1|gr(phylo, cov=A)) +  
                            (1|Species/population_id/group_id/group.year),
                       prior = c(
                         prior(normal(0, 2), "b"),
                         prior(gamma(0.1, 0.1), "phi")),
                          family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                          control=list(adapt_delta=0.99, max_treedepth=15))
summary(brm_modularity1, prob=0.89)
conditional_effects(brm_modularity1, "Cat.New")

## MODULARITY MODEL 2 -- NEW SOCIAL CATEGORIES
brm_modularity2<-brm_multiple(data=caps.groups, data2=repeated_list,
                          Modularity ~ scale(sqrt(GroupSize)) + Cat.Classic +
                            (1|gr(phylo, cov=A)) +  
                            (1|Species/population_id/group_id/group.year),
                        prior = c(
                        prior(normal(0, 2), "b"),
                        prior(gamma(0.1, 0.1), "phi")),
                          family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                          control=list(adapt_delta=0.99, max_treedepth=15))
summary(brm_modularity2,prob=0.89)
conditional_effects(brm_modularity2, "Cat.Classic", prob=0.89)

## MODULARITY MODEL 3 -- GROUP SIZE ALONE
brm_modularity3<-brm_multiple(data=caps.groups, data2=repeated_list,
                          Modularity ~ scale(sqrt(GroupSize)) +
                            (1|gr(phylo, cov=A)) +  
                            (1|Species/population_id/group_id/group.year),
                          prior = c(
                          prior(normal(0, 2), "b"),
                          prior(gamma(0.1, 0.1), "phi")),
                          family="zero_inflated_beta", iter=4000, chains=4, cores=2,
                          control=list(adapt_delta=0.99, max_treedepth=15))
summary(brm_modularity3, prob=0.89)

## MODULARITY -- COMPARE THREE
loo1<-brms::loo(brm_modularity1, pointwise = T)
loo2<-brms::loo(brm_modularity2, pointwise = T)
loo_compare(loo1,loo2)
waic(brm_modularity1)
waic(brm_modularity2)
model_weights(brm_modularity1, brm_modularity2, weights="loo")
as.data.frame(icc(brm_modularity3, by_group = T))

## ANALYZE PREDICTIVE ACCURACY
correct.perc<-c()
library(e1071)
for (i in 1:length(caps.groups)) {
  bayes.data<-caps.groups[[i]]
  colnames(caps.groups[[i]])
  bayes.data<-bayes.data[,c(1,12,14:16)]
  bayes.data<-bayes.data[bayes.data$Cat.New!="Multi-level",]
  boxplot(data=bayes.data, Modularity~Cat.New)
  
  # You would typically use a package like 'caret' for robust splitting, but a simple method works for demonstration:
  set.seed(123)
  train_indices <- sample(1:nrow(bayes.data), 0.5 * nrow(bayes.data))
  train_data <- bayes.data[train_indices, ]
  test_data <- bayes.data[-train_indices, ]
  
  # Train the Naive Bayes model
  # The formula 'Species ~ .' means predict 'Species' using all other variables
  model <- naiveBayes(Cat.New ~ ., data = train_data[,-1])
  
  # Make predictions on the test data
  predictions <- predict(model, test_data[,-1])
  
  # Evaluate the model (e.g., using a confusion matrix)
  confusion_matrix <- table(predictions, test_data$Cat.New)
  
  print(confusion_matrix)
  test_data$Prediction<-predictions
  test_data$correct<-0
  test_data$correct[test_data$Prediction==test_data$Cat.New]<-1
  correct.perc[i]<-sum(test_data$correct[test_data$Cat.New!="Multi-level"])/nrow(test_data[test_data$Cat.New!="Multi-level",])
  binom.test(sum(test_data$correct), nrow(test_data))
}

range(correct.perc)
mean(correct.perc)

## AGGREGATE DATA FOR VISUALIZATION
library(data.table)
caps.groups[[1]]$group_id
group_summary_long <- rbindlist(caps.groups, fill = TRUE, idcol = "source_id")
group_summary_avg<-group_summary_long %>%
  group_by(group.year, group_id, Species, phylo, Cat.New, GroupSize, population_id) %>%
  dplyr::summarise(Density.Avg=mean(Density), Density.Size=sqrt(1/sd(Density)), Modularity.Avg=mean(Modularity), Modularity.Size=sqrt(1/sd(Modularity)),
                   CV.Avg=mean(CV), Density.Size=sqrt(1/sd(CV)))
group_summary_avg$Modularity.Size[is.infinite(group_summary_avg$Modularity.Size)]<-5

group_summary_avg$Species<-factor(group_summary_avg$Species, levels=c("Mandrillus_sphinx","Cercocebus_sanjei", "Cercocebus_atys_atys","Theropithecus_gelada", "Lophocebus_albigena", "Papio_hamadryas", "Papio_anubis","Papio_papio", "Papio_cynocephalus", "Papio_ursinus", "Papio_kindae"))
group_summary_avg$Cat.New<-as.factor(group_summary_avg$Cat.New)

## SUPPLEMENTAL FIGURE
## PLOT BY SPECIES
library(tidybayes)
library(modelr)
data.for.grid<-data.frame(population_id=group_summary_avg$population_id, group_id=group_summary_avg$group_id,
                          year=group_summary_long$year, GroupSize=group_summary_long$GroupSize, Species=group_summary_long$Species, N.Dyads=group_summary_long$N.Dyads,
                          Cat.New=group_summary_avg$Cat.New, Cat.Classic=group_summary_long$Cat.Classic)

grid = data.for.grid %>%
  data_grid(GroupSize=seq(min(group_summary_avg$GroupSize), max(group_summary_avg$GroupSize)), N.Dyads=1, Cat.New, Obs.Effort.Per.Capita=15)

group_summary_avg$Combo<-paste(group_summary_avg$group_id, group_summary_avg$population_id, sep = "_")
min.SLS<-min(group_summary_avg$GroupSize[group_summary_avg$Cat.New=="Cohesive"])
max.SLS<-max(group_summary_avg$GroupSize[group_summary_avg$Cat.New=="Cohesive"])
min.MLS<-min(group_summary_avg$GroupSize[group_summary_avg$Cat.New=="Multi-level"])
max.MLS<-max(group_summary_avg$GroupSize[group_summary_avg$Cat.New=="Multi-level"])
min.SemiLS<-min(group_summary_avg$GroupSize[group_summary_avg$Cat.New=="Cliquish"])
max.SemiLS<-max(group_summary_avg$GroupSize[group_summary_avg$Cat.New=="Cliquish"])
grid<-grid[!(grid$Cat.New=="Cohesive" & grid$GroupSize>max.SLS),]
grid<-grid[!(grid$Cat.New=="Multi-level" & grid$GroupSize<min.MLS),]
grid<-grid[!(grid$Cat.New=="Cliquish" & grid$GroupSize<min.SemiLS),]
grid<-grid[!(grid$Cat.New=="Cliquish" & grid$GroupSize>max.SemiLS),]

## PLOT WITH SLOPES -- DENSITY
## CREATE AVERAGE ACROSS IMPUTED DATASETS
library(tidybayes)
means = grid %>%
  add_epred_draws(brm_density1, ndraws=250, allow_new_levels=T, re_formula = NA)

check<-means %>%
  group_by(GroupSize, Cat.New, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

range1<-means %>%
  group_by(Cat.New, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall.dens<-means %>%
  group_by(GroupSize, Cat.New) %>%
  dplyr::summarise(Mean.Pred=mean(.epred), LL.Pred=quantile(.epred,0.055), UL.Pred=quantile(.epred,0.945))

check$type.draw<-paste(check$Cat.New, check$.draw, sep="_")

#check$Type2<-factor(check$Type2, levels=c("Single-level","Multi-level"))
#group.summary$Type2<-factor(group.summary$Type2, levels=c("Single-level","Multi-level"))

pal1 <- c("#4D759C","#9c68b3", "#66a182")
pal12<-c("#7373A4", "#66a182")
plot1<-ggplot() +  
  scale_color_manual(values=pal1) + labs(size="Observation hours\nper individidual") +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=type.draw, color=Cat.New), data = check, alpha = 0.1, lwd=1.2) +
  geom_point(data=group_summary_avg, aes(y=Density.Avg,x=GroupSize,size=Density.Size, color=Cat.New), alpha=0.35) +
  scale_size_continuous(range = c(1, 4)) +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=Cat.New), data = overall.dens, alpha = 1, lwd=1.2, color = "black") +
  ylab("Density") + xlab("") + theme(legend.position ="none") + ylim(0,1) + xlim(0,80)
plot1

## PLOT WITH SLOPES -- MODULARITY
means = grid %>%
  add_epred_draws(brm_modularity1, ndraws=250, allow_new_levels=T, re_formula = NA)

check<-means %>%
  group_by(GroupSize, Cat.New, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

range2<-means %>%
  group_by(Cat.New, .draw) %>%
  dplyr::summarise(Mean.Pred=mean(.epred))

overall.mod<-means %>%
  group_by(GroupSize, Cat.New) %>%
  dplyr::summarise(Mean.Pred=mean(.epred), LL.Pred=quantile(.epred,0.055), UL.Pred=quantile(.epred,0.945))

check$type.draw<-paste(check$Cat.New, check$.draw, sep="_")

RColorBrewer::brewer.pal(8, "PRGn")

plot2<-ggplot() + 
  scale_color_manual(values=pal1) + labs(size="Observation hours\nper individidual") +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=type.draw, color=Cat.New), data = check, alpha = 0.1, lwd=1.2) +
  geom_point(data=group_summary_avg, aes(y=Modularity.Avg,x=GroupSize,size=Modularity.Size, color=Cat.New),alpha=0.35) +
  scale_size_continuous(range = c(1, 4)) +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=Cat.New), data = overall.mod, alpha = 1, lwd=1.2, color = "black") +
  ylab("Modularity") + xlab("Network size") + theme(legend.position ="none") + ylim(0,1) + xlim(0,80)
plot2

## QUASIRANDOM PLOTS
library(ggbeeswarm)
summary.range1<-range1 %>%
  group_by(Cat.New) %>%
  dplyr::summarise(Median=median(Mean.Pred), LL=quantile(Mean.Pred,0.055), UL=quantile(Mean.Pred,0.945), LL.75=quantile(Mean.Pred,0.125), UL.75=quantile(Mean.Pred,0.875))

summary.range1$Cat.New<-as.factor(summary.range1$Cat.New)
summary.range1$Cat.New<-factor(summary.range1$Cat.New, levels=c("Cohesive","Cliquish","Multi-level"))

plot1b<-ggplot() + 
  scale_color_manual(values=pal1) +
  geom_quasirandom(data=group_summary_avg, aes(y=Density.Avg,x=Cat.New,size=Density.Size, color=Cat.New),alpha=0.35) +
  scale_size_continuous(range = c(1, 4)) +
  geom_point(data=summary.range1,aes(x=Cat.New,y=Median)) +
  geom_errorbar(data=summary.range1,aes(x=Cat.New,ymin=LL,ymax=UL), width=0.08, lwd=1) +
  ylab("") + xlab("") + theme(legend.position ="none") + ylim(0,1)
plot1b

summary.range2<-range2 %>%
  group_by(Cat.New) %>%
  dplyr::summarise(Median=median(Mean.Pred), LL=quantile(Mean.Pred,0.055), UL=quantile(Mean.Pred,0.945), LL.75=quantile(Mean.Pred,0.125), UL.75=quantile(Mean.Pred,0.875))

plot2b<-ggplot() + 
  scale_color_manual(values=pal1) +
  geom_quasirandom(data=group_summary_avg, aes(y=Modularity.Avg,x=Cat.New,size=Modularity.Size, color=Cat.New),alpha=0.35) +
  scale_size_continuous(range = c(1, 4)) +
  geom_point(data=summary.range2,aes(x=Cat.New,y=Median)) +
  geom_errorbar(data=summary.range2,aes(x=Cat.New,ymin=LL,ymax=UL), width=0.08, lwd=1) +
  ylab("") + xlab("Social structure") + theme(legend.position ="none") + ylim(0,1)
plot2b

## PLOT COMPLETE VERSION
design<-"ABCC
         DEFF"

library(patchwork)
theme_set(theme_classic(base_size = 16))
plot1 + plot1b + plot2 + plot2b + plot_annotation(tag_levels = "a") 
ggsave(file="Figure_S2.jpg", units="cm", width=20, height=16, dpi=300)

## FIGURE 1b
## PLOT DENSITY
plot1b.dens<-group_summary_avg %>%
  group_by(Species) %>%
  add_epred_draws(brm_density1, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + xlim(0,1) + scale_fill_manual(values=pal1) + xlab("density") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=group_summary_avg,
             aes(x=Density.Avg,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) +
  ylab("") + xlab("Density") + theme(axis.title.y=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.ticks.y=element_blank())
## PLOT MODULARITY
plot1b.mod<-group_summary_avg %>%
  group_by(Species) %>%
  add_epred_draws(brm_modularity1, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + xlim(0,1) + scale_fill_manual(values=pal1) + xlab("modularity") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=group_summary_avg,
             aes(x=Modularity.Avg,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) +
  ylab("") + xlab("Modularity") + theme(axis.title.y=element_blank(),
                                        axis.text.y=element_blank(),
                                        axis.ticks.y=element_blank())

layout.tree<-"AAABBCC"
gg_tr+ plot1b.dens + plot1b.mod + plot_layout(design = layout.tree)& theme(plot.tag = element_text(size = 18), axis.title = element_text(size=25)) 
ggsave(file="Figure_1b.jpg", units="cm", width=28.3, height=12, dpi=300)

## SYNTHESIZE FOR FIGURE 2
data.by.group <- group_summary_long %>%
  group_by(Cat.New, Species)  %>%
  summarize(N.Years=length(unique(group.year)), Point.Size=sqrt(sum(N.Years)),
            GroupMean=mean(GroupSize), LL.Group=mean(GroupSize)-2*sd(GroupSize)/sqrt(N.Years), UL.Group=mean(GroupSize)+2*sd(GroupSize)/sqrt(N.Years),
            Mean.Density=mean(Density), LL.Density=mean(Density)-2*sd(Density)/sqrt(N.Years), UL.Density=mean(Density)+2*sd(Density)/sqrt(N.Years),
            Mean.Mod=mean(Modularity), LL.Mod=mean(Modularity)-2*sd(Modularity)/sqrt(N.Years), UL.Mod=mean(Modularity)+2*sd(Modularity)/sqrt(N.Years),
            Mean.CV=mean(CV), LL.CV=mean(CV)-2*sd(CV)/sqrt(N.Years), UL.CV=mean(CV)+2*sd(CV)/sqrt(N.Years))


data.by.group$Species<-factor(data.by.group$Species,levels =c("Lophocebus_albigena","Cercocebus_atys_atys","Cercocebus_sanjei",
                                                              "Papio_anubis", "Papio_cynocephalus",  
                                                              "Papio_kindae","Papio_ursinus", "Mandrillus_sphinx",
                                                              "Papio_papio", "Papio_hamadryas", "Theropithecus_gelada"))
levels(data.by.group$Species)<-c("Gray mangabey","Sooty mangabey","Sanje mangabey",
                                 "Olive baboon","Yellow baboon","Kinda baboon", "Chacma baboon","Mandrill",
                                 "Guinea baboon","Hamadryas","Gelada")
colors.new<-c("#6BAED6","#4292C6","#2171B5","#08519C","#08306B",
              "#807DBA","#6A51A3","#3F007D",
              "#31A354","#006D2C", "#00441B")

library(ggnewscale)
plot2a<-ggplot() +  
  geom_ribbon(aes(ymin=LL.Pred,ymax = UL.Pred, x= GroupSize, fill=Cat.New), data = overall.dens, alpha = 0.5) +
  geom_line(aes(y = Mean.Pred, x= GroupSize, group=Cat.New, color=Cat.New), data = overall.dens, alpha = 1, lwd=1.2) +
  scale_color_manual(values=pal1) +
  scale_fill_manual(values=pal1) +
  scale_shape_manual(values=c(21:25,21:24,21:23)) +
  new_scale_color() +
  geom_errorbar(data=data.by.group, aes(ymin=LL.Density, ymax=UL.Density,x=GroupMean, color=Species), alpha=1,width=0) +
  geom_errorbarh(data=data.by.group, aes(y=Mean.Density, xmin=LL.Group, xmax=UL.Group, color=Species), alpha=1,width=0) +
  geom_point(data=data.by.group, aes(y=Mean.Density,x=GroupMean, color=Species, shape = Species), alpha=1,size=3, stroke=1) +
  scale_color_manual(values=c(rep(pal1[1],5),
                              rep(pal1[2],3),
                              rep(pal1[3],3))) +
  ylab("Density") + xlab("Network size") + theme(legend.position ="none") + ylim(0,1) + xlim(0,80)
plot2a

plot2b<-ggplot() +  
  geom_ribbon(aes(ymin=LL.Pred,ymax = UL.Pred, x= GroupSize, fill=Cat.New), data = overall.mod, alpha = 0.5) +
  geom_line(aes(y = Mean.Pred, x= GroupSize, color=Cat.New), data = overall.mod, alpha = 1, lwd=1.2) +
  scale_color_manual(values=pal1) +
  guides(fill="none", color="none", shape = guide_legend(override.aes = list(lineCat.New = 0))) +
  new_scale_color() +
  scale_color_manual(values=c(rep(pal1[1],5),
                              rep(pal1[2],3),
                              rep(pal1[3],3))) +
  geom_errorbar(data=data.by.group, aes(ymin=LL.Mod, ymax=UL.Mod,x=GroupMean, color=Species), alpha=1, width=0) +
  geom_errorbarh(data=data.by.group, aes(y=Mean.Mod, xmin=LL.Group, xmax=UL.Group, color=Species), alpha=1,width=0) +
  geom_point(data=data.by.group, aes(y=Mean.Mod,x=GroupMean, color=Species, shape=Species), alpha=1, size=3, stroke=1) +
  scale_shape_manual(values=c(21:25,21:24,21:23)) +
  scale_fill_manual(values=pal1) +
  ylab("Modularity") + xlab("Network size") + ylim(0,1) + xlim(0,80)
plot2b

getwd()
library(patchwork)
plot2a+plot2b + plot_annotation(tag_levels = "a")
ggsave(file="Figure_2.jpg", units="cm", width=25, height=10, dpi=300)

## PART 2
## Female-female grooming ties were variably shaped by kinship, rank, and shared male ties

## FEMALE-FEMALE META-ANALYSES
## CLASSIC METHODS
## USING METAFOR
## GET GROUP-LEVEL DATA
caps.meta<-read.csv("CAPS_Meta_Groups.csv")
caps.meta$Scale.Group<-as.numeric(scale(sqrt(caps.meta$GroupSize)))
caps.meta$Cat.New<-factor(caps.meta$Cat.New, levels=c("Cohesive","Cliquish","Multi-level"))
unique(caps.meta$Philopatry)
caps.meta$Philopatry<-factor(caps.meta$Philopatry, levels=c("Female philopatry","Female dispersal"))

## KINSHIP
meta1=rma.mv(yi=ff.related.r, V=ff.related.se^2, mods= ~ Cat.New + Philopatry + Scale.Group -1,
             random = list (~1|Species/population_id/group_id/group.year,
                            ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps.meta)
summary(meta1)

## RANK
meta2=rma.mv(yi=ff.rank.r, V=ff.rank.se^2, mods= ~ Cat.New + Philopatry + Scale.Group -1,
             random = list (~1|Species/population_id/group_id/group.year,
                            ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps.meta)
summary(meta2)

meta3=rma.mv(ff.male.r, V=ff.male.se^2, mods= ~ Cat.New + Philopatry + Scale.Group -1,
             random = list (~1|Species/population_id/group_id,
                            ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps.meta)
summary(meta3)

## CLASSIC VISUALIZATIONS
plot3a<-orchaRd::orchard_plot(meta1, mod = "Cat.New", group = "Species", 
                              xlab = "Kinship effect", k=F,g=F,legend.pos="none", twig.size = 0)
plot3a<-plot3a + ylim(-1.5,4.5) +scale_fill_manual(values=pal1) + scale_color_manual(values=pal1)
plot3a

plot3b<-orchaRd::orchard_plot(meta2, mod = "Cat.New", group = "Species", 
                              xlab = "Rank similarity effect",k=F,g=F, legend.pos="none",twig.size = 0)
plot3b<-plot3b + scale_fill_manual(values=pal1) + scale_color_manual(values=pal1) + ylim(-0.5,1.5) 
plot3b 

plot3c<-orchaRd::orchard_plot(meta3, mod = "Cat.New", group = "Species", 
                              xlab = "Shared male effect",k=F,g=F, twig.size=0, legend.pos="none")
plot3c<-plot3c + scale_fill_manual(values=pal1) + scale_color_manual(values=pal1) + ylim(-1.5,4.5) 
plot3c

## BAYESIAN
## META-ANALYSES
## BRMS

## RELATEDNESS META-ANALYTIC MODEL
relatedness.meta<-brm(data = caps.meta, family = "student", data2=list(A=A),
                      ff.related.r|se(ff.related.se) ~ 0 + Cat.New + scale(sqrt(GroupSize)) + Philopatry + (1|gr(phylo, cov=A)) +
                        (1|Species/population_id/group_id/group.year),
                      prior = c(prior(normal(0, 2), class = b),
                                prior(cauchy(0, 0.3), class = sd)),
                      iter = 4000, warmup = 1000, cores = 2, chains = 4,
                      control=list(adapt_delta=0.99, max_treedepth=15))
summary(relatedness.meta, prob=0.89)
plot(conditional_effects(relatedness.meta, "Cat.New", prob=0.89), points=T)
plot(conditional_effects(relatedness.meta, "GroupSize:Cat.New", prob=0.89))
conditional_effects(relatedness.meta, "Cat.New:Philopatry", prob=0.89)

## RANK META-ANALYTIC MODEL
rank.meta<-brm(data = caps.meta, family = "student", data2=list(A=A),
               ff.rank.r| se(ff.rank.se) ~  0 + Cat.New + scale(sqrt(GroupSize)) + Philopatry + (1|gr(phylo, cov=A)) +
                 (1|Species/population_id/group_id/group.year),
               prior = c(prior(normal(0, 2), class = b),
                         prior(cauchy(0, 0.3), class = sd)),
               iter = 4000, warmup = 1000, cores = 2, chains = 4,
               control=list(adapt_delta=0.99, max_treedepth=15))
summary(rank.meta, prob=0.89)
conditional_effects(rank.meta, "Cat.New", prob=0.89)
conditional_effects(rank.meta, "GroupSize", prob=0.89)

## SHARED MALE META-ANALYTIC MODEL
male.meta<-brm(data = caps.meta, family = "student", data2=list(A=A),
               ff.male.r|se(ff.male.se) ~  0 + Cat.New + scale(sqrt(GroupSize)) + Philopatry + (1|gr(phylo, cov=A)) +
                 (1|Species/population_id/group_id/group.year),
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
meta4=rma.mv(yi=fm.rank.male.r, V=fm.rank.male.se^2, mods= ~ Cat.New + Scale.Group - 1,
             random = list(~1|Species/population_id/group_id,
                           ~1|phylo), method="REML",
             R=list(phylo=A),
             data = caps.meta)
summary(meta4)

meta5=rma.mv(fm.rank.fem.r, V=fm.rank.fem.se^2, mods= ~ Cat.New + Scale.Group - 1,
             random = list (~1|Species/population_id/group_id/group.year,
                            ~1|phylo),
             R=list(phylo=A),
             data = caps.meta)
summary(meta5)

meta6=rma.mv(fm.rank.intx.r, V=fm.rank.intx.se^2, mods= ~ Cat.New + Scale.Group - 1 ,
             random = list (~1|Species/population_id/group_id,
                            ~1|phylo),
             R=list(phylo=A),
             data = caps.meta)
summary(meta6)

## PLOT IT
model_results <- orchaRd::mod_results(meta4, mod = "Cat.New", at = NULL, group = "Species")
plot3d<-orchaRd::orchard_plot(meta4, mod = "Cat.New", group = "Species", twig.size=0, 
                              xlab = "Male rank effect", k=F,g=F, legend.pos="none")

plot3d<-plot3d+scale_fill_manual(values=pal1) + scale_color_manual(values=pal1) + ylim(-2/3,2)
plot3d

plot3e<-orchaRd::orchard_plot(meta5, mod = "Cat.New", group = "Species", twig.size=0,
                              xlab = "Female rank effect",k=F,g=F, legend.pos="none")

plot3e<-plot3e+scale_fill_manual(values=pal1) + scale_color_manual(values=pal1) + ylim(-2/3,2)
plot3e

plot3f<-orchaRd::orchard_plot(meta6, mod = "Cat.New", group = "Species", twig.size = 0,
                              xlab = "Rank interaction effect",k=F,g=F, legend.pos="none")

plot3f<-plot3f+scale_fill_manual(values=pal1) + scale_color_manual(values=pal1) + ylim(-2/3,2)
plot3f

layout1<-"AD
          BE
          CF"

theme_set(theme_test(base_size = 16))

plot3a+plot3b+plot3c+plot3d+plot3e+plot3f + plot_annotation(tag_levels = "a") + plot_layout(design=layout1)
ggsave(file="Figure_S3.jpg", units="cm", width=22, height=22, dpi=300)

## MALE RANK FM META-ANALYTIC MODEL
male.effect.fm<-brm(data = caps.meta, data2=list(A=A),
                    fm.rank.male.r|se(fm.rank.male.se) ~ 0 + Cat.New + scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                      (1|Species/population_id/group_id/group.year),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(cauchy(0, 0.3), class = sd)),    
                    family=student,
                    iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                    control=list(adapt_delta=0.99, max_treedepth=15))
summary(male.effect.fm, prob=0.89)
plot(conditional_effects(male.effect.fm, "Cat.New", prob=0.89), points=T)

## FEMALE RANK FM META-ANALYTIC MODEL
caps.meta$ASR<-caps.meta$N.Fem/caps.meta$N.Male
fem.effect.fm<-brm(data = caps.meta, family = student, data2=list(A=A),
                   fm.rank.fem.r| se(fm.rank.fem.se) ~  0 + Cat.New + scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                     (1|Species/population_id/group_id/group.year),
                   prior = c(prior(normal(0, 2), class = b),
                             prior(cauchy(0, 0.3), class = sd)),    
                   iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                   control=list(adapt_delta=0.99, max_treedepth=15))
summary(fem.effect.fm, prob=0.89)
conditional_effects(fem.effect.fm, "Cat.New", prob=0.89)
conditional_effects(fem.effect.fm, "GroupSize", prob=0.89)

## INTERACTION RANK FM META-ANALYTIC MODEL
intx.effect.fm<-brm(data = caps.meta, family = student, data2=list(A=A),
                    fm.rank.intx.r|se(fm.rank.intx.se) ~ 0 + Cat.New + scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                      (1|Species/population_id/group_id/group.year),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(cauchy(0, 0.5), class = sd)),
                    iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                    control=list(adapt_delta=0.99, max_treedepth=15))
summary(intx.effect.fm, prob=0.89)
conditional_effects(intx.effect.fm, "Cat.New", prob=0.89)
conditional_effects(intx.effect.fm, "GroupSize")

## EXTRACT PLOTS - FF DYADS
effects1<-as.data.frame(conditional_effects(relatedness.meta, "Cat.New", prob=0.89)$Cat.New)
effects2<-as.data.frame(conditional_effects(rank.meta, "Cat.New", prob=0.89)$Cat.New)
effects3<-as.data.frame(conditional_effects(male.meta, "Cat.New", prob=0.89)$Cat.New)

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
  geom_quasirandom(data=caps.meta, aes(x=Cat.New,y=ff.related.r, size=1/(ff.related.se), color=Cat.New, shape=Philopatry), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=pal1) +
  geom_point(data=effects.intx, aes(x=Cat.New, y=estimate__,shape=Philopatry),
             position=position_dodge(width=0.5), size=3) +
  geom_errorbar(data=effects.intx, aes(x=Cat.New, group=Philopatry,ymin=lower__, ymax=upper__), 
                position=position_dodge(width=0.5), width=0, lwd=1) +
  ggtitle("Female-female dyads") +
  xlab("") + ylab("Kinship effect") + theme(legend.position="none",plot.title=element_text(hjust=0.5, face="bold")) + 
  scale_y_continuous(limits=c(-1.5,4.5), breaks=seq(-1.5,4.5, by=1.5)) + coord_flip()

plot3b<-ggplot() + 
  geom_quasirandom(data=caps.meta[!caps.meta$population_id=="FILOHA",], aes(x=Cat.New,y=ff.rank.r, size=1/(ff.rank.se), color=Cat.New,shape=Philopatry), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=pal1) +
  geom_point(data=effects2, aes(x=Cat.New, y=estimate__), size = 3) +
  geom_errorbar(data=effects2, aes(x=Cat.New, ymax=lower__, ymin=upper__), width=0, lwd=1) +
  xlab("") + ylab("Rank similarity effect") + theme(legend.position="none") + ylim(-0.5,1.5) + coord_flip()

plot3c<-ggplot() + 
  geom_quasirandom(data=caps.meta, aes(x=Cat.New,y=ff.male.r, size=1/(ff.male.se), color=Cat.New,shape=Philopatry), alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=pal1) +
  geom_point(data=effects3, aes(x=Cat.New, y=estimate__), size=3) +
  geom_errorbar(data=effects3, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  xlab("") + ylab("Shared male effect") + theme(legend.position="none") + 
  scale_y_continuous(limits=c(-1.5,4.5), breaks=seq(-1.5,4.5, by=1.5)) + coord_flip()

## EXTRACT PLOTS - FM DYADS
effects1<-as.data.frame(conditional_effects(male.effect.fm, "Cat.New", prob=0.89)$Cat.New)
effects2<-as.data.frame(conditional_effects(fem.effect.fm, "Cat.New", prob=0.89)$Cat.New)
effects3<-as.data.frame(conditional_effects(intx.effect.fm, "Cat.New", prob=0.89)$Cat.New)

plot3d<-ggplot() + 
  geom_quasirandom(data=caps.meta, aes(x=Cat.New,y=fm.rank.male.r, size=1/(fm.rank.male.se), color=Cat.New,shape=Philopatry), stroke=0, alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=pal1) +
  geom_point(data=effects1, aes(x=Cat.New, y=estimate__), size=3) + 
  geom_errorbar(data=effects1, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  ggtitle("Female-male dyads") + ylim(-1,2) +
  xlab("") + ylab("Male rank effect") + theme(legend.position="none",plot.title=element_text(hjust=0.5, face="bold")) + coord_flip()
plot3d

plot3e<-ggplot() +
  geom_quasirandom(data=caps.meta[!caps.meta$population_id=="FILOHA",], aes(x=Cat.New,y=fm.rank.fem.r, size=1/(fm.rank.fem.se),shape=Philopatry, color=Cat.New), stroke=0, alpha=0.5) + 
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=pal1) +
  scale_size(range=c(1,4)) +
  geom_point(data=effects2, aes(x=Cat.New, y=estimate__), size=3) +
  geom_errorbar(data=effects2, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  xlab("") + ylab("Female rank effect") + theme(legend.position="none") + ylim(-1,2) + coord_flip()

plot3f<-ggplot() + 
  geom_quasirandom(data=caps.meta[!caps.meta$population_id=="FILOHA",], aes(x=Cat.New,y=fm.rank.intx.r, shape=Philopatry,size=1/(fm.rank.intx.se), color=Cat.New), stroke=0, alpha=0.5) + 
  scale_size(range=c(1,4)) +
  geom_hline(yintercept=0, lty=2) + scale_color_manual(values=pal1) +
  geom_point(data=effects3, aes(x=Cat.New, y=estimate__), size=3) +
  geom_errorbar(data=effects3, aes(x=Cat.New, ymin=lower__, ymax=upper__), width=0, lwd=1) +
  xlab("") + ylab("Rank interaction effect") + theme(legend.position="none") + ylim(-1,2) + coord_flip()

layout1<-"AD
          BE
          CF"
plot3a + plot3b + plot3c +
  plot3d + plot3e + plot3f + plot_annotation(tag_levels = "a") + plot_layout(design = layout1)
ggsave(file="Figure_3.jpg", units="cm", width=25, height=25, dpi=300)

## GROUP SIZE ONLY MODELS
## BAYESIAN STYLE
## RELATEDNESS META-ANALYTIC MODEL
relatedness.meta.null<-brm(data = caps.meta, family = "student", data2=list(A=A),
                      related.r| se(related.se) ~ scale(sqrt(GroupSize)) +
                        (1|gr(phylo, cov=A)) +
                        (1|Species/population_id/group_id/group.year),
                      prior = c(prior(normal(0, 2), class = "b"),
                                prior(cauchy(0, 0.3), class = "sd")),
                      iter = 4000, warmup = 1000, cores = 2, chains = 4,
                      control=list(adapt_delta=0.99, max_treedepth=15))

## RANK META-ANALYTIC MODEL
rank.meta.null<-brm(data = caps.meta, family = "student", data2=list(A=A),
               rank.r| se(rank.se) ~ scale(sqrt(GroupSize))  + 
                 (1|gr(phylo, cov=A)) + (1|Species/population_id/group_id/group.year),
               prior = c(prior(normal(0, 2), class = b),
                         prior(cauchy(0, 0.3), class = sd)),
               iter = 4000, warmup = 1000, cores = 2, chains = 4,
               control=list(adapt_delta=0.99, max_treedepth=15))

## MALE META-ANALYTIC MODEL
male.meta.null<-brm(data = caps.meta, family = "student", data2=list(A=A),
               male.r| se(male.se) ~  scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                 (1|Species/population_id/group_id/group.year),
               prior = c(prior(normal(0, 2), class = b),
                         prior(cauchy(0, 0.3), class = sd)),
               iter = 4000, warmup = 1000, cores = 2, chains = 4,
               control=list(adapt_delta=0.99, max_treedepth=15))

## VARIANCE DECOMPOSITION
## F-F EFFECTS
phylo.relat<-as.data.frame(icc(relatedness.meta.null, by_group = T))
phylo.relat$Effect<-"Kinship"

phylo.rank<-as.data.frame(icc(rank.meta.null, by_group = T))
phylo.rank$Effect<-"Rank similarity"

phylo.male<-as.data.frame(performance::icc(male.meta.null, by_group = T))
phylo.male$Effect<-"Shared male"

phylo.ests<-rbind(phylo.relat, phylo.male, phylo.rank)
phylo.ests$Effect<-as.factor(phylo.ests$Effect)
phylo.ests$Effect<-factor(phylo.ests$Effect, levels=c("Kinship","Rank similarity", "Shared male"))
phylo.ests$Group<-as.factor(phylo.ests$Group)
levels(phylo.ests$Group)<-c("Phylogeny","Species","Population","Group","Residual")

ggplot(data=phylo.ests, aes(x=Effect, y=ICC,fill=Group)) + geom_col(alpha=0.8, col="black") + theme(legend.position = "none") +
  scale_fill_brewer(palette="Set2") + ylab("Intraclass correlation coefficient") + xlab("Female-female effects")

## GET GROUP SIZE ONLY MODELS
## BAYESIAN MODELS
## MALE RANK FM META-ANALYTIC MODEL
male.effect.fm.null<-brm(data = caps.meta, data2=list(A=A),
                    fm.rank.male.r|se(fm.rank.male.se) ~ scale(sqrt(GroupSize)) + 
                      (1|gr(phylo, cov=A)) +
                      (1|Species/population_id/group_id/group.year),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(cauchy(0, 0.3), class = sd)),    
                    family=student,
                    iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                    control=list(adapt_delta=0.99, max_treedepth=15))

## FEMALE RANK FM META-ANALYTIC MODEL
fem.effect.fm.null<-brm(data = caps.meta, family = student, data2=list(A=A),
                   fm.rank.fem.r| se(fm.rank.fem.se) ~  scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                     (1|Species/population_id/group_id/group.year),
                   prior = c(prior(normal(0, 2), class = b),
                             prior(cauchy(0, 0.3), class = sd)),    
                   iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                   control=list(adapt_delta=0.99, max_treedepth=15))

## INTERACTION RANK FM META-ANALYTIC MODEL
intx.effect.fm.null<-brm(data = caps.meta, family = student, data2=list(A=A),
                    fm.rank.intx.r|se(fm.rank.intx.se) ~ scale(sqrt(GroupSize)) + (1|gr(phylo, cov=A)) +
                      (1|Species/population_id/group_id/group.year),
                    prior = c(prior(normal(0, 2), class = b),
                              prior(cauchy(0, 0.5), class = sd)),
                    iter = 4000, warmup = 1000, cores = 2, chains = 4, sample_prior = T,
                    control=list(adapt_delta=0.99, max_treedepth=15))

## F-M EFFECTS
phylo.fm.male<-as.data.frame(icc(male.effect.fm.null, by_group = T))
phylo.fm.male$Effect<-"Male rank"

phylo.fm.fem<-as.data.frame(icc(fem.effect.fm.null, by_group = T))
phylo.fm.fem$Effect<-"Female rank"

phylo.intx.fem<-as.data.frame(icc(intx.effect.fm.null, by_group = T))
phylo.intx.fem$Effect<-"Rank interaction"

phylo.ests2<-rbind(phylo.fm.male, phylo.fm.fem, phylo.intx.fem)
phylo.ests2$Effect<-as.factor(phylo.ests2$Effect)
phylo.ests2$Effect<-factor(phylo.ests2$Effect, levels=c("Male rank","Female rank","Rank interaction"))
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
ggsave(file="Figure_S5.jpg", units="cm", width=20, height=20, dpi=300)

## NEW SUPPLEMENTARY FIGURES 
caps.meta$Species<-as.character(caps.meta$Species)
caps.meta$Species[caps.meta$phylo=="Papio_cynocephalus"]<-"Yellow\nbaboon"
caps.meta$Species[caps.meta$phylo=="Papio_hamadryas"]<-"Hamadryas\nbaboon"
caps.meta$Species[caps.meta$phylo=="Papio_kindae"]<-"Kinda\nbaboon"
caps.meta$Species[caps.meta$phylo=="Papio_ursinus"]<-"Chacma\nbaboon"
caps.meta$Species[caps.meta$phylo=="Papio_papio"]<-"Guinea\nbaboon"
caps.meta$Species[caps.meta$phylo=="Papio_anubis"]<-"Olive\nbaboon"
caps.meta$Species[caps.meta$phylo=="Lophocebus_aterrimus"]<-"Gray-cheeked\nmangabey"
caps.meta$Species[caps.meta$phylo=="Cercocebus_chrysogaster"]<-"Sanje\nmangabey"
caps.meta$Species[caps.meta$phylo=="Cercocebus_lunulatus"]<-"Sooty\nmangabey"

## FEMALE-FEMALE EFFECTS BY SPECIES
## PLOT RELATED EFFECT BY SPECIES
caps.meta.short<-caps.meta[!caps.meta$phylo=="Papio_hamadryas",]
caps.meta.short$phylo<-as.factor(caps.meta.short$phylo)
caps.meta.short<-droplevels(caps.meta.short)
levels(caps.meta.short$phylo)
caps.meta$Species<-factor(caps.meta$Species,levels=c("Mandrill","Sanje\nmangabey", "Sooty\nmangabey","Gelada", "Gray-cheeked\nmangabey", "Hamadryas\nbaboon", "Olive\nbaboon","Guinea\nbaboon", "Yellow\nbaboon", "Chacma\nbaboon", "Kinda\nbaboon"))

plot.ff.related<-caps.meta.short %>%
  group_by(Species) %>%
  add_epred_draws(relatedness.meta, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + scale_fill_manual(values=pal1) + xlab("modularity") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=caps.meta.short,
             aes(x=ff.related.r,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) +
  ggtitle("Female-female effects") + theme(plot.title = element_text(face="bold",hjust=0.5)) +
  ylab("Species") + xlab("Kinship effect") + xlim(-1.5,4.5) + geom_vline(xintercept=0,lty=2) +
  scale_y_discrete(labels=c('Mandrill',"Sanje","Sooty","Gelada","Gray","Olive","Guinea","Yellow", "Chacma","Kinda"))

## PLOT RANK EFFECT BY SPECIES
plot.ff.ranks<-caps.meta %>%
  group_by(Species) %>%
  add_epred_draws(rank.meta, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + scale_fill_manual(values=pal1) + xlab("modularity") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=caps.meta,
             aes(x=ff.rank.r,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) +
  ylab("Species") + xlab("Rank effect") + xlim(-0.5,1.5) + geom_vline(xintercept=0,lty=2) +
  scale_y_discrete(labels=c('Mandrill',"Sanje","Sooty","Gelada","Gray","Hamadryas","Olive","Guinea","Yellow", "Chacma","Kinda"))

## PLOT RANK EFFECT BY SPECIES
plot.ff.male<-caps.meta %>%
  group_by(Species) %>%
  add_epred_draws(male.meta, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + scale_fill_manual(values=pal1) + xlab("modularity") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=caps.meta,
             aes(x=ff.male.r,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) +
  ylab("Species") + xlab("Shared male effect")  + geom_vline(xintercept=0,lty=2) + xlim(-1.5,4.5) +
  scale_y_discrete(labels=c('Mandrill',"Sanje","Sooty","Gelada","Gray","Hamadryas","Olive","Guinea","Yellow", "Chacma","Kinda"))

## FEMALE-MALE EFFECTS BY SPECIES
## PLOT MALE RANK EFFECT BY SPECIES
caps.meta.short<-caps.meta[!is.na(caps.meta$fm.rank.fem.r),]
plot.fm.male<-caps.meta.short %>%
  group_by(Species) %>%
  add_epred_draws(male.effect.fm, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + scale_fill_manual(values=pal1) + xlab("modularity") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=caps.meta.short,
             aes(x=fm.rank.male.r,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) + 
  ggtitle("Female-male effects") + theme(plot.title = element_text(face="bold",hjust=0.5)) +
  ylab("") + xlab("Male rank effect") + xlim(-1,2) + geom_vline(xintercept=0,lty=2) +
  scale_y_discrete(labels=c('Mandrill',"Sanje","Sooty","Gelada","Gray","Hamadryas","Olive","Guinea","Yellow", "Chacma","Kinda"))

## PLOT RANK EFFECT BY SPECIES
plot.fm.fem<-caps.meta.short %>%
  group_by(Species) %>%
  add_epred_draws(fem.effect.fm, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + scale_fill_manual(values=pal1) + xlab("modularity") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=caps.meta.short,
             aes(x=fm.rank.fem.r,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) + 
  ylab("") + xlab("Female rank effect") + xlim(-1,2) + geom_vline(xintercept=0,lty=2) +
  scale_y_discrete(labels=c('Mandrill',"Sanje","Sooty","Gelada","Gray","Hamadryas","Olive","Guinea","Yellow", "Chacma","Kinda"))

## PLOT RANK EFFECT BY SPECIES
plot.fm.intx<-caps.meta.short %>%
  group_by(Species) %>%
  add_epred_draws(intx.effect.fm, ndraws=1000) %>%
  ggplot(aes(x = .epred, y = Species, fill=Cat.New)) + scale_fill_manual(values=pal1) + xlab("modularity") +
  stat_halfeye(aes(x = .epred), scale = 0.6, adjust=5,  alpha=0.8,position = position_nudge(y = 0.175)) + theme(legend.position = "none") +
  geom_point(data=caps.meta.short,
             aes(x=fm.rank.intx.r,y=Species, color=Cat.New), shape="|", size=2, position=position_nudge(y=-0.1),
             alpha = 1) + scale_color_manual(values=pal1) + 
  ylab("") + xlab("Rank interaction effect") + xlim(-1,2) + geom_vline(xintercept=0,lty=2) +
  scale_y_discrete(labels=c('Mandrill',"Sanje","Sooty","Gelada","Gray","Hamadryas","Olive","Guinea","Yellow", "Chacma","Kinda"))

layout1<-"AD
          BE
          CF"
plot.ff.related + plot.ff.ranks + plot.ff.male +
plot.fm.male + plot.fm.fem + plot.fm.intx + plot_layout(design = layout1) + plot_annotation(tag_levels="a")
ggsave(file="Figure_S4.jpg", units="cm", width=30, height=30, dpi=300)

## PRINT MODEL OUTPUTS
sink("CAPS_models.txt")
summary(brm_density1, prob=0.89)
summary(brm_density2, prob=0.89)
summary(brm_density3, prob=0.89)
summary(brm_modularity1, prob=0.89)
summary(brm_modularity2, prob=0.89)
summary(brm_modularity3, prob=0.89)
summary(relatedness.meta, prob=0.89)
summary(rank.meta, prob=0.89)
summary(male.meta, prob=0.89)
summary(male.effect.fm, prob=0.89)
summary(fem.effect.fm, prob=0.89)
summary(intx.effect.fm, prob=0.89)
summary(meta1)
summary(meta2)
summary(meta3)
summary(meta4)
summary(meta5)
summary(meta6)
sink()

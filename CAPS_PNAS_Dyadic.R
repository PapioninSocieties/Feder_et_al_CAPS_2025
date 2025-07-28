## LOOPS TO GENERATE DYADIC REGRESSION MODELS FOR EACH GROUP-YEAR WITHIN THE DATASET
## DRAWING FROM 5 IMPUTED DATASETS FOR EACH GROUP-YEAR

## FEMALE-FEMALE DYADS
## RUN LOOOP
models.to.check<-c()
models.to.check2<-c()
i<-1
for(i in 1:nrow(caps.groups)) {
  group.year<-caps.groups$group.year[i]
  edges <- lapply(FF_Dyads, \(x) {
    x <- x[x$group.year == group.year, ]
  })
  boxplot(data=edges[[5]], groom_z~related, main=group.year)
  if(caps.groups$population_id[i]=="FILOHA") {
    model<-brm_multiple(data=edges, Grooming~ rank.z + Higher.Target + male.combo + offset(log(Obs.effort)) + 
                          (1|mm(ID1,ID2)), family="zero_inflated_negbinomial", chains=2, cores=2,
                        prior = c(
                          prior(normal(0, 2), "b"),
                          prior(student_t(3, 0, 5), "sd")),
                        control=list(adapt_delta=0.9))
    models.to.check[i]<-max(rhat(model))
    models.to.check2[i]<-min(summary(model)$fixed$Bulk_ESS)
    caps.groups$ff.rank.r[i]<-fixef(model)[2,1]
    caps.groups$ff.rank.se[i]<-fixef(model)[2,2]
    caps.groups$ff.groomee.r[i]<-fixef(model)[3,1]
    caps.groups$ff.groomee.se[i]<-fixef(model)[3,2]
    caps.groups$ff.male.r[i]<-fixef(model)[4,1]
    caps.groups$ff.male.se[i]<-fixef(model)[4,2]
  } else{
    model<-brm_multiple(data=edges, Grooming~related + rank.z + Higher.Target + male.combo + offset(log(Obs.effort)) + 
                          (1|mm(ID1,ID2)), family="zero_inflated_negbinomial", chains=2, cores=2,
                        prior = c(
                          prior(normal(0, 2), "b"),
                          prior(student_t(3, 0, 5), "sd")),
                        control=list(adapt_delta=0.95))
    models.to.check[i]<-max(rhat(model))
    models.to.check2[i]<-min(summary(model)$fixed$Bulk_ESS)
    caps.groups$ff.kinship.r[i]<-fixef(model)[2,1]
    caps.groups$ff.kinship.se[i]<-fixef(model)[2,2]
    caps.groups$ff.rank.r[i]<-fixef(model)[3,1]
    caps.groups$ff.rank.se[i]<-fixef(model)[3,2]
    caps.groups$ff.groomee.r[i]<-fixef(model)[4,1]
    caps.groups$ff.groomee.se[i]<-fixef(model)[4,2]
    caps.groups$ff.male.r[i]<-fixef(model)[5,1]
    caps.groups$ff.male.se[i]<-fixef(model)[5,2]
  }
}

## FEMALE-MALE DYADS
## RUN LOOOP
models.to.check<-c()
models.to.check2<-c()
for(i in 1:nrow(caps.groups)) {
  group.year<-caps.groups$group.year[i]
  edges <- lapply(male.edge.list, \(x) {
    x <- x[x$group.year == group.year, ]
  }) 
  check<-edges[[1]]
  boxplot(data=check, groom_z~rank.male.z, main=group.year)
  model<-brm_multiple(data=edges, Grooming~rank.male.z*rank.fem.z + subject_sex1 + offset(log(Obs.effort)) + 
                        (1|mm(ID1,ID2)), family="zero_inflated_negbinomial", chains=2, cores=2,
                      prior = c(
                        prior(normal(0, 2), "b"),
                        prior(student_t(3, 0, 5), "sd")),
                      control=list(adapt_delta=0.95))
  models.to.check[i]<-max(rhat(model))
  models.to.check2[i]<-min(summary(model)$fixed$Bulk_ESS)
  caps.groups$fm.rank.male.r[i]<-fixef(model)[2,1]
  caps.groups$fm.rank.male.se[i]<-fixef(model)[2,2]
  caps.groups$fm.rank.fem.r[i]<-fixef(model)[3,1]
  caps.groups$fm.rank.fem.se[i]<-fixef(model)[3,2]
  caps.groups$fm.rank.intx.r[i]<-fixef(model)[5,1]
  caps.groups$fm.rank.intx.se[i]<-fixef(model)[5,2]
  caps.groups$fm.sex.r[i]<-fixef(model)[4,1]
  caps.groups$fm.sex.se[i]<-fixef(model)[4,2]
}
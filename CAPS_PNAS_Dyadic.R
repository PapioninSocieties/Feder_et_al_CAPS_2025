## LOOPS TO GENERATE DYADIC REGRESSION MODELS FOR EACH GROUP-YEAR WITHIN THE DATASET
## DRAWING FROM 20 IMPUTED DATASETS FOR EACH GROUP-YEAR

## FEMALE-FEMALE DYADS
## RUN LOOOP
for(i in 1:nrow(meta.group.ff)) {
  group.year<-meta.group.ff$group.year[i]
  edges <- lapply(imputed.ff.new, \(x) {
    x <- x[x$group.year == group.year, ]
  })
  if(meta.group.ff$population_id[i]=="FILOHA") {
    model<-brm_multiple(data=edges, Groom.Index  ~ scale(RankDiff) + male.combo +  (1|mm(ID1,ID2)), 
                        family="Gamma", chains=4, cores=4,
                        prior = c(
                          prior(normal(0, 2), "b"),
                          prior(student_t(3, 0, 5), "sd")),
                        control=list(adapt_delta=0.95))
    meta.group.ff$rank.r[i]<-fixef(model)[2,1]
    meta.group.ff$rank.se[i]<-fixef(model)[2,2]
    meta.group.ff$male.r[i]<-fixef(model)[3,1]
    meta.group.ff$male.se[i]<-fixef(model)[3,2]
  } else{
    model<-brm_multiple(data=edges, Groom.Index ~ related + scale(RankDiff) + male.combo + (1|mm(ID1,ID2)),
                        family="gamma", chains=4, cores=4,
                        prior = c(
                          prior(normal(0, 2), "b"),
                          prior(student_t(3, 0, 5), "sd")),
                        control=list(adapt_delta=0.95))
    meta.group.ff$related.r[i]<-fixef(model)[2,1]
    meta.group.ff$related.se[i]<-fixef(model)[2,2]
    meta.group.ff$rank.r[i]<-fixef(model)[3,1]
    meta.group.ff$rank.se[i]<-fixef(model)[3,2]
    meta.group.ff$male.r[i]<-fixef(model)[4,1]
    meta.group.ff$male.se[i]<-fixef(model)[4,2]
  }
    print(group.year)
}

## FEMALE-MALE DYADS
## RUN LOOOP
for(i in 1:nrow(meta.group.fm)) {
  group.year<-meta.group.fm$group.year[i]
  edges <- lapply(imputed.fm.new, \(x) {
    x <- x[x$group.year == group.year, ]
  }) 
  model<-brm_multiple(data=edges, Groom.Index~scale(Rank.Male)*scale(Rank.Fem)+
                        (1|mm(ID1,ID2)), family="gamma", chains=4, cores=4,
                      prior = c(
                        prior(normal(0, 2), "b"),
                        prior(student_t(3, 0, 5), "sd")),
                      control=list(adapt_delta=0.95))
  meta.group.fm$fm.rank.male.r[i]<-fixef(model)[2,1]
  meta.group.fm$fm.rank.male.se[i]<-fixef(model)[2,2]
  meta.group.fm$fm.rank.fem.r[i]<-fixef(model)[3,1]
  meta.group.fm$fm.rank.fem.se[i]<-fixef(model)[3,2]
  meta.group.fm$fm.rank.intx.r[i]<-fixef(model)[4,1]
  meta.group.fm$fm.rank.intx.se[i]<-fixef(model)[4,2]
  print(group.year)
}

load(file = "./data/llik_NL1.Rdata")
llik_NL1<-llik
load(file = "./data/llik_NL2.Rdata")
llik_NL2<-llik
load(file = "./data/llik.Rdata")
llik_NL <-list()
llik_NL[[4]] <-c(llik_NL1[[4]][llik_NL1[[4]]<llik_NL2[[4]]],llik_NL2[[4]][llik_NL2[[4]]<=llik_NL1[[4]]])
llik_NL[[4]] <- llik_NL[[4]][names(llik[[4]])]
print(llik_NL[[4]]<llik[[4]])
print(llik_NL[[4]])
print(llik[[4]])
# plot2
llik3<-llik[[4]][llik[[4]]<200 & llik[[4]]>0]
llik4<-llik_NL[[4]][llik[[4]]<200 & llik[[4]]>0]
print(llik4)
llik3_grand<-c(llik3[llik3>llik4],0*llik3[llik3<=llik4])
llik3_petit<-c(llik3[llik3<=llik4],0*llik3[llik3>llik4])
llik4_grand<-c(llik4[llik3<=llik4],0*llik4[llik3>llik4])
llik4_petit<-c(llik4[llik3>llik4],0*llik4[llik3<=llik4])
llik3_grand<-llik3_grand[names(llik4_petit)]
llik3_petit<-llik3_petit[names(llik4_petit)]
llik4_grand<-llik4_grand[names(llik4_petit)]
print(llik4_grand)
print(llik4_petit)
barplot(llik4_grand,main="Compared likelihood of model 4 and model 4 combined non-linear",ylab = "likelihood",ylim=c(0,200),cex.names=0.5,cex.main=0.75)
barplot(llik3_grand,col="cadetblue1",add=TRUE,ylim=c(0,200),cex.names=0.5)
barplot(llik4_petit,add=TRUE,ylim=c(0,200),cex.names=0.5)
barplot(llik3_petit,col="cadetblue1",add = TRUE,ylim=c(0,200),cex.names=0.5)
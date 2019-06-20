load(file = "./data/AIC_NL1.Rdata")
AIC_NL1<-AIC
load(file = "./data/AIC_NL2.Rdata")
AIC_NL2<-AIC
load(file = "./data/AIC.Rdata")
AIC_NL <-list()
AIC_NL[[4]] <-c(AIC_NL1[[4]][AIC_NL1[[4]]<AIC_NL2[[4]]],AIC_NL2[[4]][AIC_NL2[[4]]<=AIC_NL1[[4]]])
AIC_NL[[4]] <- AIC_NL[[4]][names(AIC[[4]])]
print(AIC_NL[[4]]<AIC[[3]])
print(sum(AIC_NL[[4]]<AIC[[3]]))

# plot2
AIC3<-AIC[[4]][AIC[[4]]<200 & AIC[[4]]>0]
AIC4<-AIC_NL[[4]][AIC[[4]]<200 & AIC[[4]]>0]
AIC3_grand<-c(AIC3[AIC3>AIC4],0*AIC3[AIC3<=AIC4])
AIC3_petit<-c(AIC3[AIC3<=AIC4],0*AIC3[AIC3>AIC4])
AIC4_grand<-c(AIC4[AIC3<=AIC4],0*AIC4[AIC3>AIC4])
AIC4_petit<-c(AIC4[AIC3>AIC4],0*AIC4[AIC3<=AIC4])
AIC3_grand<-AIC3_grand[names(AIC4_petit)]
AIC3_petit<-AIC3_petit[names(AIC4_petit)]
AIC4_grand<-AIC4_grand[names(AIC4_petit)]

barplot(AIC4_grand,main="Compared AIC of model 4 and model 4 combined non-linear",ylab = "likelihood",ylim=c(0,200),cex.names=0.5,cex.main=0.75)
barplot(AIC3_grand,col="cadetblue1",add=TRUE,ylim=c(0,200),cex.names=0.5)
barplot(AIC4_petit,add=TRUE,ylim=c(0,200),cex.names=0.5)
barplot(AIC3_petit,col="cadetblue1",add = TRUE,ylim=c(0,200),cex.names=0.5)


diversity<-read.csv("diversity.csv",header = TRUE)
fix(diversity)
attach(diversity)
Nlevelf<-factor(diversity$Nlevel)
Blockf<-factor(diversity$Block)
table(Nadd,Nlevel)
fit<-aov(S~Nadd*Nlevelf+Blockf)
summary(fit)

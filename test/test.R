#! /software/additional/bin/Rscript

#setwd("../")
require("devtools")
load_all()

data(mouse)


###################################################
### code chunk number 3: internal
###################################################
express <- mouse[,c("M1","M2","M3","NC1","NC2","NC3")]
rownames(express) <- mouse$ID

intern <- clValid(express, 2:16, clMethods=c("hierarchical","kmeans","pam"), ncore, validation="internal")
summary(intern)

stab <- clValid(express, 2:16, clMethods=c("hierarchical","kmeans","pam"), ncore, validation="stability")
optimalScores(stab)

fc <- tapply(rownames(express),mouse$FC, c)
fc <- fc[!names(fc)%in%c("EST","Unknown")]
bio <- clValid(express, 2:16, clMethods=c("hierarchical","kmeans","pam"), ncore, validation="biological", annotation=fc)
optimalScores(bio)

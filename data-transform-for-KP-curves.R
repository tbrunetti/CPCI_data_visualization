library(survival)
#makes a new matrix of patients that have KRAS and TP53 mutations
tempTable<-matrix(ncol=87, nrow=11)

numPos=0
for (i in 1:82){
  if (decisionTree_input_UC$KRAS[i]=="+" && decisionTree_input_UC$TP53[i]=="+"){
    numPos<-numPos+1
    tempTable[numPos,]<-as.matrix(decisionTree_input_UC[i,])
  }
}

colnames(tempTable)<-names(decisionTree_input_UC)
##delete columns
tempTable<-tempTable[,-2]
tempTable<-tempTable[,-6]
decisionTree_input_UC<-decisionTree_input_UC[,-2]
decisionTree_input_UC<-decisionTree_input_UC[,-6]

varlistToIterate<-names(decisionTree_input_UC)[2:83]  
dataFrameTemp<-data.frame(tempTable)
dataFrameTemp<-transform(dataFrameTemp, days=as.numeric(as.character(days)), outcome=as.numeric(as.character(outcome)))
#test<-survfit(Surv(days, outcome) ~BAGE2, data=dataFrameTemp)
survCalc<-survdiff(Surv(days, outcome) ~CIITA, data=dataFrameTemp)



model<-lapply(varlistToIterate, function(x){
  summary<-survfit(substitute(Surv(days, outcome) ~i, list(i=as.name(x))), data=dataFrameTemp)
})


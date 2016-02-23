##################################################################################
###�˩w��լ����t��
##Version 1.0
##Known Issue: Use too much loop and temp variable, should use apply family in the future 
##2015.12.11
##Author:GregTsai   Email:gregorytsai@gmail.com
##Developed on R Ver.3.2.2
##MBESS 3.3.3
##################################################################################

##�Ĥ@������άO�ɯ�R���w�ˮM��
install.packages("Hmisc") 
install.packages("psych") 
install.packages("cocor") # for CI
##�}��R����]�Ĥ@���ݸ��J�M��
library("Hmisc")
library("psych") 
library("cocor") 


####�ϥΪ̻ݫ��w����
##Ū�����
##�۰ʿ�ܰ��榹��
data=read.csv(file.choose(), header=T) 
##��ʫ��w���榹��
data=read.csv("C://Users/Greg/OneDrive/���פ���/Raw data/BioPsychoSocial.csv", header=T) 
data=read.csv("C://Users/Grego/OneDrive/���פ���/Raw data/CorrelationTest.csv", header=T) 
data=read.csv("E://OneDrive/���פ���/Raw data/BioPsychoSocial.csv", header=T) 
data1=read.csv("C://Users/Grego/OneDrive/���פ���/Raw data/103-2-SenzaMedi.csv", header=T)  #%~3!<K%;



dataTemp=data1[,c(42,3,4,15,43,18,51,55,52,56,59,54,57,60,71,72)]

data = dataTemp[complete.cases(dataTemp),]

data=dataDiversity[,c(2,1,3,4,5,14:ncol(dataDiversity))]


##�g�J���G�ANot Yet Finished�A�٤��ब�ʦ�����ܼg�J�a�I�M�ɮצW��
write.csv("Cor-CESD-Nomed.csv", x=result, row.names = T)

##�]�w�n�]�����d��
#�w���ܶ�
xStart=2
xEnd=ncol(data)
#�էO�ܶ��Ҧb���
groupCol=1
#�էO��CODING,��r���[�W�޸�""
groupA=0
groupB=1
#���G�n��ܪ��էO�W��
nameGroupA="Low CESD"
nameGroupB="High CESD"



#######################�H�U�����ϥΪ̤������#######################################
if(require("RevoUtilsMath")){ setMKLthreads(2) } 


nVar=xEnd-xStart+1


#���G�[�W�P�����
addStar <- function(betaValue,pvalue)
{
  if (pvalue<0.001) {betaStar=paste0(betaValue,"***")}
  else if (pvalue<0.01) {betaStar=paste0(betaValue,"**")}
  else if (pvalue<=0.05) {betaStar=paste0(betaValue,"*")}
  else if (pvalue<0.10) {betaStar=paste0(betaValue,".")}
  else {betaStar=betaValue}
  
  return(betaStar)
}




# Correlations with significance levels
# rcorr(as.matrix(mydata[,c(3,31,26,8,12,9,13,16,11,14,17)]), type="pearson") 

corMatrixA=rcorr(as.matrix(data[which(data[,groupCol]==groupA),xStart:xEnd]), type="pearson") # type can be pearson or spearman
corMatrixB=rcorr(as.matrix(data[which(data[,groupCol]==groupB),xStart:xEnd]), type="pearson") # type can be pearson or spearman

corA=corMatrixA$r
corB=corMatrixB$r

numA=corMatrixA$n
numB=corMatrixB$n

proA=corMatrixA$P
proB=corMatrixB$P

proA[is.na(proA)]<-1
proB[is.na(proB)]<-1

resultA=corA
resultB=corB

compareTemp=r.test(n=numA,corA,corB,n2=numB)
resultCompare=matrix(nrow=nVar,ncol=nVar)
proCom=compareTemp$p


for (i in 1:nrow(corA))
  {
  for (j in 1:ncol(corA))
    {
    resultA[i,j]=paste0(addStar(round(corA[i,j],2),proA[i,j])," (",numA[i,j],")")
    resultB[i,j]=paste0(addStar(round(corB[i,j],2),proB[i,j])," (",numB[i,j],")")
    if (is.na(compareTemp$z[i,j])==F)    {
      a=cocor.indep.groups(r1.jk=corA[i,j], r2.hm=corB[i,j], n1=numA[i,j], n2=numB[i,j], alternative="two.sided", alpha=0.05, conf.level=0.95, null.value=0)
      resultCompare[i,j]=paste0(addStar(round(compareTemp$z[i,j],2),proCom[i,j])," (",round(a@zou2007$conf.int[1],2),"/",round(a@zou2007$conf.int[2],2),")")
    } 
  }
}

result=matrix(nrow=nVar*3,ncol=nVar,dimnames=list(rep("",nVar*3),colnames(corA)))


for (i in 1:nVar ){
  result[i*3-2,]=resultA[i,]
  rownames(result)[i*3-2]=rownames(resultA)[i]
  result[i*3-1,]=resultB[i,]
  result[i*3,]=resultCompare[i,]
}

result <- cbind( rep( c( paste0(nameGroupA," r(n)"),paste0(nameGroupB," r(n)"),"z value z(CI)"),nVar ),result )


#Experimental
# mapply(addStar(x,y),x=corMatrixA$r,y=proA)


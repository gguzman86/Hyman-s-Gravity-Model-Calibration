Production=matrix(c(400,350,250), byrow=TRUE)
Attraction=c(300,200,500)
Cij=matrix(c(5, 10, 18, 13, 5, 15, 20, 16, 6), nrow=3, ncol=3, byrow=TRUE)
Tij=matrix(c(148,79,173,99,82,169,53,39,158), nrow=3, ncol=3, byrow=TRUE)



##################Deterrence Function###########################
deterrence_function <- function(Cij, Beta){
  New_Cij=matrix(0,nrow(Cij),ncol(Cij))
  for (i in 1:length(Cij)){
    New_Cij[i]=exp(-(Beta)*Cij[i])
  }
  return(New_Cij)
}
####################END##########################################

##################Balancing Function###########################
Balancing<-function(Production, Attraction, Cij, Beta){
  
  Tolerance=0.01
  
  #Iteration 0
  Iter=deterrence_function(Cij,Beta)
  SumP=matrix(rowSums(Iter),byrow=TRUE)
  SumA=colSums(Iter)
  ai=Attraction/SumA
  bj=matrix(Production/SumP, byrow=TRUE)
  Ai=matrix(1,nrow(Iter),ncol=1,byrow=TRUE)
  Bj=matrix(1,nrow=1,ncol(Iter),byrow=TRUE)
  
  #Generate Bj
  Sum_Bj=matrix(0,nrow=1,ncol(Iter),byrow=TRUE)
  for (j in 1:nrow(Iter)){
    for (i in 1:nrow(Iter)){
      Sum_Bj[j]=sum(Ai[i]*Production[i]*Iter[i,j])+Sum_Bj[j]
    }
  }
  Bj=1/Sum_Bj
  
  #Generate Ai
  Sum_Ai=matrix(1,nrow(Iter),ncol=1,byrow=TRUE)
  Ai=1/Sum_Ai 
  
  #Loop
  z=3
  repeat{
    
    #Odds Iteration(1,3,5,...)
    
    if(z%%2==1){
      for(i in 1:nrow(Iter)){
        for(j in 1:ncol(Iter)){
          Iter[i,j]=Iter[i,j]*ai[j]
        }
      }
      SumP=matrix(rowSums(Iter),byrow=TRUE)
      SumA=colSums(Iter)
      ai=Attraction/SumA
      bj=matrix(Production/SumP, byrow=TRUE)
      
      #Generate Ai
      Sum_Ai=matrix(0,nrow(Iter),ncol=1,byrow=TRUE)
      for (i in 1:nrow(Iter)){
        for (j in 1:nrow(Iter)){
          Sum_Ai[i]=sum(Bj[j]*Attraction[j]*Iter[i,j])+Sum_Ai[i]
        }
      } 
      Ai=1/Sum_Ai
      
      #Generate Bj
      Sum_Bj=matrix(0,nrow=1,ncol(Iter),byrow=TRUE)
      for (j in 1:nrow(Iter)){
        for (i in 1:nrow(Iter)){
          Sum_Bj[j]=sum(Ai[i]*Production[i]*Iter[i,j])+Sum_Bj[j]
        }
      }
      Bj=1/Sum_Bj     
    }
    #Even Iteration(2,4,6...)
    else{
      for (i in 1:nrow(Iter)){
        for(j in 1:ncol(Iter)){
          Iter[i,j]=Iter[i,j]*bj[i]
        }
      }
      SumP=matrix(rowSums(Iter),byrow=TRUE)
      SumA=colSums(Iter)
      ai=Attraction/SumA
      bj=matrix(Production/SumP, byrow=TRUE)
      
      #Generate Bj
      i=1
      j=1
      Sum_Bj=matrix(0,nrow=1,ncol(Iter),byrow=TRUE)
      for (j in 1:nrow(Iter)){
        for (i in 1:nrow(Iter)){
          Sum_Bj[j]=sum(Ai[i]*Production[i]*Iter[i,j])+Sum_Bj[j]
        }
      }
      Bj=1/Sum_Bj
      
      #Generate Ai
      Sum_Ai=matrix(0,nrow(Iter),ncol=1,byrow=TRUE)
      for (i in 1:nrow(Iter)){
        for (j in 1:nrow(Iter)){
          Sum_Ai[i]=sum(Bj[j]*Attraction[j]*Iter[i,j])+Sum_Ai[i]
        }
      } 
      Ai=1/Sum_Ai   
      
    }
    
    #Generate Tolerance Check
    max_ai=abs(max(ai)-1)
    max_bj=abs(max(bj)-1)
    Max_T=max(max_ai,max_bj)
    
    #Tolerance Check
    if(Max_T<Tolerance){
      break
    }
    z=z+1
  }
  Results<-list("Iter"=Iter,"Ai"=Ai,"Bj"=Bj)
  return(Results)
}
####################END##########################################

################Hyman's Function######################################
Hyman<-function(Production, Attraction, Cij, Tij){
  Total_Cost=sum(Cij*Tij)
  Max_Iter=1000
  Tol=0.05
  Beta_Cost_Matrix=matrix(0,nrow=1,ncol=2,byrow=TRUE)
  
  #Iteration 1
  m=1
  z=1
  Beta=1.5/Total_Cost
  Balance<-Balancing(Production, Attraction, Cij, Beta)
  OD1 <-Balance$Iter
  A1 <- Balance$Ai
  B1 <- Balance$Bj
  New_Cost1=sum(OD1*Cij)
  Beta_Cost_Matrix=rbind(Beta_Cost_Matrix, c(Beta,New_Cost1) )
  Cost_Check=(abs(New_Cost1-Total_Cost))/Total_Cost
  if(Cost_Check<Tol){
    break
  }
  
  #Iteration 2
  m=2
  z=2
  Beta=(New_Cost1*Beta)/Total_Cost
  Balance<-Balancing(Production, Attraction, Cij, Beta)
  OD2 <-Balance$Iter
  A2 <- Balance$Ai
  B2 <- Balance$Bj 
  New_Cost2=sum(OD2*Cij)
  Beta_Cost_Matrix=rbind(Beta_Cost_Matrix, c(Beta,New_Cost2))
  Cost_Check=(abs(New_Cost2-Total_Cost))/Total_Cost
  if(Cost_Check<Tol){
    break
  }
  
  #m Iteration
  m=3
  z=3
  repeat{
    Beta=(((Total_Cost-Beta_Cost_Matrix[m-1,2])*Beta_Cost_Matrix[m,1])-((Total_Cost-Beta_Cost_Matrix[m,2])*Beta_Cost_Matrix[m-1,1]))/(Beta_Cost_Matrix[m,2]-Beta_Cost_Matrix[m-1,2])
    Balance<-Balancing(Production, Attraction, Cij, Beta)
    OD <-Balance$Iter
    A <- Balance$Ai
    B <- Balance$Bj 
    New_Cost=sum(OD*Cij)
    Beta_Cost_Matrix=rbind(Beta_Cost_Matrix, c(Beta,New_Cost))
    Cost_Check=(abs(New_Cost-Total_Cost))/Total_Cost
    if((Cost_Check<Tol)|z>Max_Iter){
      break
    }
    Beta_Cost_Matrix=Beta_Cost_Matrix[-1,]
    
    plot(z,Beta,xlim=c(0,Max_Iter),ylim=c(-1,1))
    Sys.sleep(0.1)
    z=z+1
  }
  Results<-list("Beta"=Beta, "OD"=OD,"A"=A,"B"=B)
  return(Results)
}
####################END##########################################

################Fernandez calibration######################################

Fernandez_Calibration <- function(Production, Attraction, Cij, Tij){
  Yij=matrix(0,nrow(Tij),ncol(Tij))
  for (j in 1:nrow(Yij)){
    for (i in 1:nrow(Yij)){
      Yij[i,j]=log(Tij[i,j])-log(Production[i])-log(Attraction[j])
    }
  } 
  
  Sij=1/Cij
  
  Max_Iter=100
  
  #Step 1
  Beta=0.5
  ro=0.5
  n=1
  
  repeat{
  Beta0=Beta
  ro0=ro
  
  #Step 2
  Balance<-Balancing(Production, Attraction, Cij, Beta)
  OD0 <-Balance$Iter
  A0 <- Balance$Ai
  B0 <- Balance$Bj
  
  
  #Step 3
  Yij0=matrix(0,nrow(Tij),ncol(Tij))
  for (j in 1:nrow(Yij)){
    for (i in 1:nrow(Yij)){
      Yij0[i,j]=Yij[i,j]-log(A0[i])-log(B0[j])
    }
  }   
  
  #Step 4
  r_Yij0=do.call(rbind,lapply(Yij0, unlist))
  r_Cij=do.call(rbind,lapply(Cij, unlist))
  r_Sij=do.call(rbind,lapply(Sij, unlist))
  
  reg <- lm(r_Yij0 ~ r_Cij + log(r_Sij))
  theta0 = summary(reg)$coefficients[1]
  Beta = -1*(summary(reg)$coefficients[2])
  ro = summary(reg)$coefficients[3]
  
  #Step 5
  
  if(Beta0-Beta<0.00000001 && ro0-ro<0.00000001){
    break
  }
  
  plot(n,Beta,xlim=c(0,Max_Iter),ylim=c(-1,1))
  Sys.sleep(0.1)
  n=n+1
  
  }
  Results<-list("Beta"=Beta, "OD"=OD0, "A"=A0,"B"=B0, "ro"=ro)
  return(Results) 
}
 

####################END##########################################

Fernandez_Calibrated_Beta<-Fernandez_Calibration(Production, Attraction, Cij, Tij)
Fernandez_New_Beta=Fernandez_Calibrated_Beta$Beta
Fernandez_New_Beta
Fernandez_New_ro=Fernandez_Calibrated_Beta$ro
Fernandez_New_ro
New_Ai=Fernandez_Calibrated_Beta$A
New_Ai
New_Bj=Fernandez_Calibrated_Beta$B
New_Bj

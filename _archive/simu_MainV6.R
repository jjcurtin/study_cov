tic = proc.time()[3] #to calculate elapsed time for simulations
########################################################
#General setup
S=100000                   #S=150000 #number of simulations per set of experiment characteristics
Ns=c(40,80,120,160)        #N
Bxs = c(0, .5)             #Population parameter for X
Kcs = c(4,8,16)            #number of covariates available
rCYs = c(.15, .30, .45)    #Population parameter for bivariate correlation between Y and each Cov with non-zero effect
nSigs = c(.25, .50, .75)   #proportion of non-zero Cov effects in set
rMin = -.6                 #minimum r between Covs for truncated normal distributio
rMax = .6                  #maximum r between Covs for truncated normal distributio
rSD = .1                   #SD of normal distribution of rs among Covs (MEAN =0)
ver = 6                    #for output file names to track versions
set.seed(01271969)         #for reproducable results

library(bestglm) #for bestglm() for BMX models with AIC selection
library(Matrix)  #for nearPD() to test SIGMA is positive definite
library(truncnorm) #for rtruncnorm for distribution of covariate covariances
library(MASS) #for mvnorm
########################################################
for (Bx in Bxs){
  for (N in Ns){
    for (Kc in Kcs){
      for (rCY in rCYs){
        for (nSig in nSigs){
          
          #report simulation characteristics
          print(sprintf(' Simulating: Bx= %0.1f; N= %0.0f; Kc= %0.0f; rCY= %0.2f; nSig= %0.2f', Bx, N, Kc, rCY, nSig))

          
          ########################################################
          #For simulation outcome storage
          
          #parameter estimates, etc for X from No covariate (NC) model
          bxNC = vector('numeric', S)
          SExNC = vector('numeric', S)
          txNC = vector('numeric', S) 
          pxNC = vector('numeric', S)
          
          #parameter estimates, etc for X from All covariate (AC) model
          bxAC = vector('numeric', S)
          SExAC = vector('numeric', S)
          txAC = vector('numeric', S) 
          pxAC = vector('numeric', S)
          
          #parameter estimates, etc for X from pCovX (PCX)  model
          bxPCX = vector('numeric', S)
          SExPCX = vector('numeric', S)
          txPCX = vector('numeric', S) 
          pxPCX = vector('numeric', S)    
          
          #parameter estimates, etc for X from pCovO (PCO)  model
          bxPCO = vector('numeric', S)
          SExPCO = vector('numeric', S)
          txPCO = vector('numeric', S) 
          pxPCO = vector('numeric', S)   
          
#           #parameter estimates, etc for X from Best Model w/X (BMX) model
#           bxBMX = vector('numeric', S)
#           SExBMX = vector('numeric', S)
#           txBMX = vector('numeric', S) 
#           pxBMX = vector('numeric', S)          
          
          #parameter estimates, etc for X from pX (PX) model
          bxPX = vector('numeric', S)
          SExPX = vector('numeric', S)
          txPX = vector('numeric', S) 
          pxPX = vector('numeric', S)          
          
          #Covariate parameter estimates from full model (useful as check)
          bc = matrix(NA,S,Kc, dimnames = list(1:S, paste('bc', 1:Kc, sep='')) )
          
          #Other indices to use for checks
          dx= vector('numeric', S) #effect size d for X effect
          VarY= vector('numeric', S)   #variance of Y
          VarYxl= vector('numeric', S) #variance of Y in low level of X
          VarYxh= vector('numeric', S) #variance of Y in high level of X
          rC1Y = vector('numeric', S) #correlation between C1 (always non-zero) and Y
          rC1Yxl = vector('numeric', S) #correlation between C1 (always non-zero) and Y in low level of X
          rC1Yxh = vector('numeric', S) #correlation between C1 (always non-zero) and Y in high level of X
          
          #Covariate ts from full model (useful for calculating effect size)
          tc= matrix(NA,S,Kc,dimnames = list(1:S, paste('tc_', 1:Kc, sep='')) )
          
          #Was covariate |t| > 1 in X + C model and therefore included in pCovX model
          icPCX= matrix(FALSE,S,Kc,dimnames = list(1:S, paste('icPCX_', 1:Kc, sep='')) )

          #Was covariate |t| > 1 in C model and therefore included in pCovO model
          icPCO= matrix(FALSE,S,Kc,dimnames = list(1:S, paste('icPCO_', 1:Kc, sep='')) )
          
#           #Was covariate included in BMX model
#           icBMX= matrix(FALSE,S,Kc,dimnames = list(1:S, paste('icBMX_', 1:Kc, sep='')) )

          #Did the covariate decrease p-value for X and therefore include in the pX model
          icPX= matrix(FALSE,S,Kc, dimnames = list(1:S, paste('icPX_', 1:Kc, sep='')) )
          
          cSumPCX = vector('numeric', S)  #binary variable to indicate which covariates included in pCovX model
          cSumPCO = vector('numeric', S)  #binary variable to indicate which covariates included in pCovO model
#          cSumBMX = vector('numeric', S)  #binary variable to indicate which covariates included in BMX model
          cSumPX = vector('numeric', S)   #binary variable to indicate which covariates included in pX model
          
          ########################################################
          #Set up X
          X = rep(c(-.5,.5),N/2)  #Set up X
          
          ########################################################
          #Initial setup for covariates and Y
          
          #Set up Bc based on Kc, rCY, and nSig
          Bc = c(rep(rCY,(Kc*nSig)),rep(0,(Kc*(1-nSig))))
          BNeg = rbinom(Kc,1,.5)  #generate binary numbers to indicate negative Bc
          Bc[which(BNeg==1)] = Bc[which(BNeg==1)] * -1  #set random (p=.5) Bcs to negative
          
          #Set up parameters and matrix for covariates themselves
          CovY = matrix(NA,nrow=N, ncol=Kc+1,dimnames = list(1:N,  c(paste('C', 1:Kc, sep=''),'Y')) )  #matrix of raw covariates and Y.  Values assigned in simulation loop
          MU = rep(0,Kc+1) #vector of means for Covs and Y (all 0)
          nCovars = (Kc * (Kc-1))/2  #number of covariances among Covs to set later in simulation.  rCY sets relationship between Covs and Y
          
          #######################################################
          #The simulations          
          for (s in (1:S)){
            
            if (s %% 1000 == 0) {print(s)}  #used to confirm simulation is progressing
            
            
            #Set up SIGMA which is var-covar matrix for Covs and Y
            SIGMA = list(mat = NA, converged = NA)
            SIGMA$mat = matrix(rep(0,(Kc+1)*(Kc+1)), nrow=Kc+1, ncol=Kc+1)  #variance-covariance matrix for Covs and Y
            SIGMA$converged = FALSE  #set to false before entering loop to establish postiive definite CovY matrix
            
            while(!SIGMA$converged) {
              Covars = c(rtruncnorm(nCovars,rMin, rMax,0,rSD), Bc) #Covariances among Covs from truncated (to be safe) normal distribution and Bc based on rCY
              SIGMA$mat[upper.tri(SIGMA$mat)]= Covars  #set upper triangle 
              SIGMA$mat = SIGMA$mat + t(SIGMA$mat)  #add lower triangle
              diag(SIGMA$mat) = 1  #set variances of Covs and Y to 1 in population (i.e., correlation matrix). 
              
              #check if SIGMA positive definite and modify if needed
              SIGMA = nearPD(SIGMA$mat, corr = TRUE, keepDiag = TRUE, ensureSymmetry = TRUE)  #SIGMA$converged = TRUE if successful
            }
            
            
            CovY = mvrnorm(N, MU, SIGMA$mat)  #make multivariate normal raw data for Covs and Y based on MU and SIGMA
            CovY = matrix(CovY,nrow=N, ncol=Kc+1,dimnames = list(1:N,  c(paste('C', 1:Kc, sep=''),'Y')) )  #add dimnames to CovY 
            
            #Calc final Y by adding X effect into Y from CovY
            Y = CovY[,'Y'] + Bx*X 
            
            #Full covariate model########################
            lmAC = 'Y ~ X'  #lm formula for AC model
            for (i in 1:Kc) lmAC = paste(lmAC, ' + CovY[,', i, ']', sep='') 
          
            mAC = lm(lmAC)
            tAC = summary(mAC)
            bxAC[s]=tAC$coefficients[2,1]  #parameter estimate for X
            SExAC[s]=tAC$coefficients[2,2] #Standard error
            txAC[s]=tAC$coefficients[2,3]  #t
            pxAC[s]=tAC$coefficients[2,4]  #p-value
            
            bc[s,]=tAC$coefficients[3:(Kc+2),1]    #parameter estimates for covs  
            tc[s,]=tAC$coefficients[3:(Kc+2),3]    #ts for covs
          
            
            #no covariate model##########################
            mNC = lm(Y ~ X)
            tNC = summary(mNC)
            bxNC[s]=tNC$coefficients[2,1]  #parameter estimate
            SExNC[s]=tNC$coefficients[2,2] #Standard error
            txNC[s]=tNC$coefficients[2,3]  #t
            pxNC[s]=tNC$coefficients[2,4]  #p-value
            
            #compute some checks (d effect size for X, var Y, and rC1Y)
            dx[s] = bxNC[s] / ((var(Y[X==-.5])+var(Y[X==.5]))/2)^.5  #cohen's d
            VarY[s] = var(Y)   #variance Y
            VarYxl[s] = var(Y[X==-.5])  #variance Y in low X level
            VarYxh[s] = var(Y[X== .5]) #variance  Y in low y level
            rC1Y[s] = cor(Y,CovY[,1])
            rC1Yxl[s] = cor(Y[X==-.5],CovY[,1][X==-.5])
            rC1Yxh[s] = cor(Y[X== .5],CovY[,1][X== .5])
            
          
            #X + C models to select pCX, pCO, and pX models###########################
            lmPCX = 'Y ~ X'  #lm formula for pCX model (modified below)
            lmPCO = 'Y ~ X'  #lm formula for pCO model (modified below)
            lmPX  = 'Y ~ X'  #lm formula for pX model (modified below)
            
            #loop through all covariates one at a time to develop pCX, pCO and pX models
            for (i in 1:Kc){
              m = lm(Y ~ X + CovY[,i])
              t = summary(m)

              #pCX
              if (abs(t$coefficients[3,3]) >1){  #Cov t-value > 1
                cSumPCX[s] = cSumPCX[s] + 10^(i-1)  #record selection in binary summary variable (1=C1, 10=C2, 100=C3, etc)
                icPCX[s,i] = TRUE  #Cov selected
                lmPCX = paste(lmPCX, ' + CovY[,', i, ']', sep='')
              }
              
              #pX
              if (t$coefficients[2,4] < tNC$coefficients[2,4]){  #pX is smaller than in No Cov model
                cSumPX[s] = cSumPX[s] + + 10^(i-1)
                icPX[s,i] = TRUE
                lmPX = paste(lmPX, ' + CovY[,', i, ']', sep='')
              }
              
              #pCO
              #Need model with only C for pCO so have to rerun new models
              m = lm(Y ~ CovY[,i])
              t = summary(m)
              if (abs(t$coefficients[2,3]) >1){  #Cov t-value > 1
                cSumPCO[s] = cSumPCO[s] + 10^(i-1)  #record selection in binary summary variable (1=C1, 10=C2, 100=C3, etc)
                icPCO[s,i] = TRUE  #Cov selected
                lmPCO = paste(lmPCO, ' + CovY[,', i, ']', sep='')
              }
            }
                     
            #run final pCX model###########################
            mPCX = lm(lmPCX)
            tPCX = summary(mPCX)
            bxPCX[s]=tPCX$coefficients[2,1]
            SExPCX[s]=tPCX$coefficients[2,2]
            txPCX[s]=tPCX$coefficients[2,3]
            pxPCX[s]=tPCX$coefficients[2,4] 
            
            #run final pCO model###########################
            mPCO = lm(lmPCO)
            tPCO = summary(mPCO)
            bxPCO[s]=tPCO$coefficients[2,1]
            SExPCO[s]=tPCO$coefficients[2,2]
            txPCO[s]=tPCO$coefficients[2,3]
            pxPCO[s]=tPCO$coefficients[2,4]                
          
            
            #run final pX model###########################
            mPX = lm(lmPX)
            tPX = summary(mPX)
            bxPX[s]=tPX$coefficients[2,1]
            SExPX[s]=tPX$coefficients[2,2]
            txPX[s]=tPX$coefficients[2,3]
            pxPX[s]=tPX$coefficients[2,4]
            
#             #run Best covariate model########################
#             mRes = lm(Y ~ X) #force in X by using residuals as Y
#             XY = as.data.frame(cbind(CovY[,1:Kc],mRes$residuals))
#             bestAIC <- bestglm(XY, IC="AIC", TopModels=1)
#             iCs = bestAIC$Subsets[which(bestAIC$Subsets[,Kc+3] == min(bestAIC$Subsets[,Kc+3])),2:(Kc+1)]  #row with < AIC, col 1 is intercept
#             iCs = as.matrix(iCs)          
#             
#             lmBMX = 'Y ~ X'  #lm formula for BMX model
#             for (i in 1:Kc){
#               if (iCs[i]){
#                 cSumBMX[s] = cSumBMX[s] + 10^(i-1)  #record selection in binary summary variable (1=C1, 10=C2, 100=C3, etc)
#                 icBMX[s,i] = TRUE  #Cov selected
#                 lmBMX = paste(lmBMX, ' + CovY[,', i, ']', sep='') 
#               }
#             }
#             
#             mBMX  = lm(lmBMX)
#             tBMX = summary(mBMX)
#             bxBMX[s]=tBMX$coefficients[2,1]
#             SExBMX[s]=tBMX$coefficients[2,2]
#             txBMX[s]=tBMX$coefficients[2,3]
#             pxBMX[s]=tBMX$coefficients[2,4]
            
          }#for s
          
#          d = data.frame(rCY, Kc, N, Bx, dx, VarY, VarYxl, VarYxh, rC1Y, rC1Yxl, rC1Yxh, bxNC, SExNC, txNC, pxNC, bxAC, SExAC, txAC, pxAC, bxPCX, SExPCX, txPCX, pxPCX, bxPCO, SExPCO, txPCO, pxPCO, bxBMX, SExBMX, txBMX, pxBMX, bxPX, SExPX, txPX, pxPX, bc, tc, icPCX, cSumPCX, icPCO, cSumPCO, icBMX, cSumBMX, icPX, cSumPX)
          d = data.frame(rCY, Kc, N, Bx, dx, VarY, VarYxl, VarYxh, rC1Y, rC1Yxl, rC1Yxh, bxNC, SExNC, txNC, pxNC, bxAC, SExAC, txAC, pxAC, bxPCX, SExPCX, txPCX, pxPCX, bxPCO, SExPCO, txPCO, pxPCO, bxPX, SExPX, txPX, pxPX, bc, tc, icPCX, cSumPCX, icPCO, cSumPCO,icPX, cSumPX)
          
          setwd('P:\\StudyData\\ANCOVA1\\Analysis\\Data')
          
          if(nSig==.25){nSigLabel = 'Q'}
          if(nSig==.50){nSigLabel = 'H'}
          if(nSig==.75){nSigLabel = 'T'}
          CovModel = paste('rCY', rCY*100, nSigLabel, Kc, sep='')

          filename = paste('v', ver, '_Bx', Bx*10, '_', CovModel, '_', N, '.rds', sep='')
          saveRDS(d,file=filename)
          
        }#for nSig
      }#for rCY
    } #for Kc
  }  #for N
} #for Bx

###################################
#report elapsed time for simulation
toc = proc.time()[3] - tic
sprintf('Elapsed time: %f hours', toc/3600)  
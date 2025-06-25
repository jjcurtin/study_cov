tic = proc.time()[3] #to calculate elapsed time for simulations
########################################################
#General setup
S=100000 #number of simulations per set of experiment characteristics (Used to be 150,000)
Ns = 120#Ns=c(40,80,120)            #N
Bxs = c(0, .5)             #Population parameter for X
Kcs=4# Kcs = c(4,8,16)            #number of covariates available
dR2s = .25# dR2s = c(.10, .25, .40)    #Population parameter for Covs
nSigs = .5# nSigs = c(.25, .50, .75)   #proportion of non-zero Cov effects in set
ver = 1                    #for output file names to track versions
set.seed(01271969)         #for reproducable results

########################################################
for (Bx in Bxs){
  for (N in Ns){
    for (Kc in Kcs){
      for (dR2 in dR2s){
        for (nSig in nSigs){
          
          #report simulation characteristics
          print(sprintf(' Simulating: Bx= %0.1f; N= %0.0f; Kc= %0.0f; dR2= %0.2f; nSig= %0.2f', Bx, N, Kc, dR2, nSig))

          
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
          
          #Covariate ts from full model (useful for calculating effect size)
          tc= matrix(NA,S,Kc,dimnames = list(1:S, paste('tc_', 1:Kc, sep='')) )
          
          #Was covariate |t| > 1 in X + C model and therefore included in pCovX model
          icPCX= matrix(FALSE,S,Kc,dimnames = list(1:S, paste('icPCX_', 1:Kc, sep='')) )

          #Was covariate |t| > 1 in C model and therefore included in pCovO model
          icPCO= matrix(FALSE,S,Kc,dimnames = list(1:S, paste('icPCO_', 1:Kc, sep='')) )
          
          #Did the covariate decrease p-value for X and therefore include in the pX model
          icPX= matrix(FALSE,S,Kc, dimnames = list(1:S, paste('icPX_', 1:Kc, sep='')) )
          
          cSumPCX = vector('numeric', S)  #binary variable to indicate which covariates included in pCovX model
          cSumPCO = vector('numeric', S)  #binary variable to indicate which covariates included in pCovO model
          cSumPX = vector('numeric', S)  #binary variable to indicate which covariates included in pX model
          
          ########################################################
          #Set up X
          X = rep(c(-.5,.5),N/2)  #Set up X
          
          ########################################################
          #Initial setup for covariates
          
          #Set up Bc based on Kc, dR2, and nSig
          Bc = c(rep(1,(Kc*nSig)),rep(0,(Kc*(1-nSig))))
          BCMag = (dR2/(Kc*nSig))^.5  #Calc Bc necessary for dR2 with nSig proportion of the Covs- others will be 0
          Bc = Bc*BCMag  #Population parameters for covariates

          #Set up parameters and matrix for covariates themselves
          Cov = matrix(NA,nrow=N, ncol=Kc,dimnames = list(1:N, paste('C', 1:Kc, sep='')) )  #matrix of covariates.  Values assigned in simulation loop
          MU = rep(0,Kc)
          SIGMA = matrix(rep(0,Kc*Kc), nrow=Kc, ncol=Kc)
          diag(SIGMA) = 1  #set variances of Covs to 1 in population.   for now covariances are 0 in population by default
          #NOTES on making a symetric matrix for when I need to have nonzero covariances
          #http://davetang.org/muse/2014/01/22/making-symmetric-matrices-in-r/
          
          
          ########################################################
          #Noise for Y.  Noise decreases as we "find" more covariates
          #This should maintain sd(Y) at 1 + the sd induced by manipulation of X
          #Therefore Bx = d and Bc = Rcy when calculated within each level of X
          YNoise = (1-sum(Bc^2))^.5
          
          #######################################################
          #The simulations          
          for (s in (1:S)){
            
            #Set up covariate values for this S 
            Cov = mvrnorm(N, MU, SIGMA)
            
            #Set up Y as function of X + Covs + Noise
            Y = Bx*X + Cov%*%Bc + rnorm(N,0,YNoise)  #X + Noise
          
            
            #Full covariate model########################
            lmAC = 'Y ~ X'  #lm formula for AC model
            for (i in 1:Kc) lmAC = paste(lmAC, ' + Cov[,', i, ']', sep='') 
          
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
            
            #computer d effect size for X and var Y checks
            dx[s] = bxNC[s] / ((var(Y[X==-.5])+var(Y[X==.5]))/2)^.5  #cohen's d
            VarY[s] = var(Y)   #variance Y
            VarYxl[s] = var(Y[X==-.5])  #variance Y in low X level
            VarYxh[s] = var(Y[X== .5]) #variance  Y in low y level
            
          
            #X + C models to select pCovX, PCovO, and pX models###########################
            lmPCX = 'Y ~ X'  #lm formula for pCovS model (modified below)
            lmPCO = 'Y ~ X'  #lm formula for pCov1 model (modified below)
            lmPX  = 'Y ~ X'  #lm formula for pX model (modified below)
            
            #loop through all covariates one at a time to develop pCovX, pCovO and pX models
            for (i in 1:Kc){
              m = lm(Y ~ X + Cov[,i])
              t = summary(m)

              #pCovX
              if (abs(t$coefficients[3,3]) >1){  #Cov t-value > 1
                cSumPCX[s] = cSumPCX[s] + 10^(i-1)  #record selection in binary summary variable (1=C1, 10=C2, 100=C3, etc)
                icPCX[s,i] = TRUE  #Cov selected
                lmPCX = paste(lmPCX, ' + Cov[,', i, ']', sep='')
              }
              
              #pX
              if (t$coefficients[2,4] < tNC$coefficients[2,4]){  #pX is smaller than in No Cov model
                cSumPX[s] = cSumPX[s] + + 10^(i-1)
                icPX[s,i] = TRUE
                lmPX = paste(lmPX, ' + Cov[,', i, ']', sep='')
              }
              
              #pCovO
              #Need model with only C for pCovO so have to rerun new models
              m = lm(Y ~ Cov[,i])
              t = summary(m)
              if (abs(t$coefficients[2,3]) >1){  #Cov t-value > 1
                cSumPCO[s] = cSumPCO[s] + 10^(i-1)  #record selection in binary summary variable (1=C1, 10=C2, 100=C3, etc)
                icPCO[s,i] = TRUE  #Cov selected
                lmPCO = paste(lmPCO, ' + Cov[,', i, ']', sep='')
              }
            }
            
            
            
            #run final pCovX model###########################
            mPCX = lm(lmPCX)
            tPCX = summary(mPCX)
            bxPCX[s]=tPCX$coefficients[2,1]
            SExPCX[s]=tPCX$coefficients[2,2]
            txPCX[s]=tPCX$coefficients[2,3]
            pxPCX[s]=tPCX$coefficients[2,4] 
            
            #run final pCovO model###########################
            mPCO = lm(lmPCO)
            tPCO = summary(mPCO)
            bxPCO[s]=tPCO$coefficients[2,1]
            SExPCO[s]=tPCO$coefficients[2,2]
            txPCO[s]=tPCO$coefficients[2,3]
            pxPCO[s]=tPCO$coefficients[2,4]                
          
            
            #fun final pX model###########################
            mPX = lm(lmPX)
            tPX = summary(mPX)
            bxPX[s]=tPX$coefficients[2,1]
            SExPX[s]=tPX$coefficients[2,2]
            txPX[s]=tPX$coefficients[2,3]
            pxPX[s]=tPX$coefficients[2,4]
          }#for s
          
          d = data.frame(dR2, Kc, N, Bx, dx, VarY, VarYxl, VarYxh, bxNC, SExNC, txNC, pxNC, bxAC, SExAC, txAC, pxAC, bxPCX, SExPCX, txPCX, pxPCX, bxPCO, SExPCO, txPCO, pxPCO, bxPX, SExPX, txPX, pxPX, bc, tc, icPCX, cSumPCX, icPCO, cSumPCO, icPX, cSumPX)
          setwd('P:\\StudyData\\ANCOVA1\\Analysis\\Data')
          
          if(nSig==.25){nSigLabel = 'Q'}
          if(nSig==.50){nSigLabel = 'H'}
          if(nSig==.75){nSigLabel = 'T'}
          CovModel = paste('dR', dR2*100, nSigLabel, Kc, sep='')

          filename = paste('v', ver, '_Bx', Bx*10, '_', CovModel, '_', N, '.rds', sep='')
          saveRDS(d,file=filename)
          
        }#for nSig
      }#for dR2
    } #for Kc
  }  #for N
} #for Bx

###################################
#report elapsed time for simulation
toc = proc.time()[3] - tic
sprintf('Elapsed time: %f hours', toc/3600)  
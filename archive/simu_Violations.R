#CHECK OPTIONS FOR GENERATING POSTIVE SKEW VIA F DISTRIBUTION
# n=1000
# df1=10
# df2=100
# f = rf(n,df1,df2)
# mean(f)
# sd(f)
# skew(f)
# kurtosi(f)
# plot(density(f))


########################################################
#THE SIMULATION
set.seed(01271969)  
S=100000 #number of simulations per set of experiment characteristics
P = rep(NA,S)  #p-values
D = rep(NA,S)  #decision
n=50
df1=2 #(For F)  skew of 2ish
df2=1000 #(For F)

X = rep(c(-.5,.5),n)

for (i in 1:S){
  #print(sprintf('Simulation: %0.0f', i))
  Y = c(rf(n,df1,df2),rf(n,df1,df2))  #skewed F distribution (mean D =  0.04827)
  #Y = c(rnorm(n,0,1),rnorm(n,0,2))  #heteroscedacity distribution 2X SD in one group (mean D = 0.05147)
  #Y = c(rnorm(n,0,1),rnorm(n,0,1))  #normal distribution (sanity check)
  
  m = lm(Y ~ X)
  t=summary(m)
  P[i] = t$coefficients[2,4]
  D[i] = P[i] < .05
}
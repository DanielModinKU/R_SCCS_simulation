#time grid 
time.grid <- seq(0,20,1)
NT = length(time.grid)
#sample size 
sample.size = 10 

#exposure duration
exp.duration = 2

#pre exposure 
X = rlnorm(n=sample.size,mean=2,sd=0.5)

#odds ratio exposure
oddsratio.exp = 1 

#oddsratioo for comoorbs
oddsrato.X = 1.05

#baseline risk 
odds.Y = 1/1000

#data matrix 
E = matrix(numeric(sample.size*NT),ncol=NT)
Y = matrix(numeric(sample.size*NT),ncol=NT)

#prob exposure
exp.prob = 0.1

#outcome model 
lin.pred = log(oddsrato.X)*X + log(odds.Y)

for (j in time.grid) {
  if (j >= exp.duration) {
    print (j) 
    E[,j+1] <- rbinom(sample.size,size=1,exp.prob)
    has.exposure = rowSums(E)
    E[has.exposure == 1, j+1] <- 0 #reset exspoure to 0 for previously exposed
    
    current.exp = rowSums(E[,(j+1-exp.duration):(j+1)])
    
    current.lin.pred = lin.pred + log(oddsratio.exp)*current.exp
    
    expit <- function(x) { exp(x)/(1+exp(x)) }
    
    Y[,j+1] <- rbinom(sample.size,size=1,expit(current.lin.pred))
    
  }
}

hist(X)





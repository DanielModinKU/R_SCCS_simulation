
################ DEFINE FUNCTIONS #######################

#expit function
expit <- function(x) { exp(x)/(1+exp(x)) }

#simulate data function 
sim <- function(time.grid.lower,  time.grid.upper,  time.step ,  sample.size,  exposure.duration,  frailty.mean,  frailty.sd ,  oddsratio.exposure,  oddsratio.X,  odds.baseline,  exposure.probability ) {
  
  #create time grid 
  time.grid <- seq(time.grid.lower,time.grid.upper,  time.step)
  time.grid.length = length(time.grid)
  
  #lav pre-study frailty 
  X = rlnorm(n= sample.size,  mean = frailty.mean,  sd = frailty.sd)
  
  #data matrix initialization
  E = matrix( numeric(sample.size*time.grid.length),  ncol = time.grid.length)
  Y = matrix( numeric(sample.size*time.grid.length),  ncol = time.grid.length)
  
  #linear prediction based on frailty and baseline risk 
  lin.pred = log(oddsratio.X)*X + log(odds.baseline)
  
  #simulate data
  for (j in time.grid) { 
    if (j>= exposure.duration) {
     print(j)
      
     #draw exposure 
     E[,j+1] <- rbinom(sample.size,  size=1, exposure.probability)
     
     #reset exposure to 0 for already exposed
     has.exposure = rowSums(E)
     E[has.exposure == 1, j+1] <- 0 #reset exp to 0 for prev exposed
     
     #update linear predictor ift exposed or not
     current.exposure = rowSums(E[,(j+1-exposure.duration):(j+1)])
     current.lin.pred = lin.pred + log(oddsratio.exposure)*current.exposure
    
     #update outcomes 
     Y[,j+1] <- rbinom(sample.size,size=1,expit(current.lin.pred))
    }
  }
  
  #plot pt frailty distribution
  hist(X)
  
  #return data in a list with E and Y 
  return.list = list("E"=E, "Y"=Y)
  return(return.list) 
}


################ RUN SIMULATIONS  #######################

#run sim 1
sim1 = sim( 0,   #time.grid.lower
        40,      #time.grid.upper
         1,      #time.step
         1000,   #sample.size
         2,      #exposure.duration
         2,      #frailty.mean 
         0.5,    #frailty.sd
         1,      #oddsratio.exposure
         1.05,   #oddsratio.X (odds associeret med comorbs/frailty)
         0.001,  #odds.baseline aka baseline risk 
         0.1     #exposure.probability 
)











################ DEFINE FUNCTIONS #######################

#expit function
expit <- function(x) { exp(x)/(1+exp(x)) }

#simulate data function 
sim <- function(time.grid.lower,  time.grid.upper,  time.step ,  sample.size,  exposure.duration,  frailty.mean,  frailty.sd ,  oddsratio.exposure,  oddsratio.X,  odds.baseline,  exposure.probability, baseline.probability.observation, oddsratio.observation ) {
  
  #create time grid 
  time.grid <- seq(time.grid.lower,time.grid.upper,  time.step)
  time.grid.length = length(time.grid)
  
  #lav pre-study frailty 
  X = rlnorm(n= sample.size,  mean = frailty.mean,  sd = frailty.sd)
  
  #data matrix initialization
  E = matrix( numeric(sample.size*time.grid.length),  ncol = time.grid.length)
  Y = matrix( numeric(sample.size*time.grid.length),  ncol = time.grid.length)
  E.observed = matrix( numeric(sample.size*time.grid.length),  ncol = time.grid.length)
  
  #outcome linear prediction based on frailty and baseline risk 
  b_0 = log(odds.baseline)
  b_1 = log(oddsratio.X)
  lin.pred = b_0 + b_1*X     
  
  #exposure coefficient
  b_exposure = log(oddsratio.exposure)
  
  
  #exposure observation probability
  b_baseline = log(baseline.probability.observation/(1-baseline.probability.observation))
  b_observation = log(oddsratio.observation)
  lin.pred.observation =b_baseline + b_observation*X
  observation.probability = expit(lin.pred.observation)

  
  #simulate data
  for (j in time.grid) { 
    if (j>= exposure.duration) {
     print(j)
      
     #draw exposure 
     E[,j+1] <- rbinom(sample.size,  size=1, exposure.probability)
     
     #reset exposure to 0 for already exposed
     has.exposure = rowSums(E[, 1:j])                                        #OBS: her i oprindelig kode fra møde stod rowSums(E), som gav fejl der når kode køres får alle exposures slettet (fordi exposure i E[,j+1] medtages)
     E[has.exposure == 1, j+1] <- 0 #reset exp to 0 for prev exposed
     
     #update linear predictor ift exposed or not
     current.exposure = rowSums(E[,(j+1-exposure.duration):(j+1)])
     current.lin.pred = lin.pred + b_exposure*current.exposure
     
     #update outcome probability
     outcome.probability = expit(current.lin.pred)
    
     #update outcomes 
     Y[,j+1] <- rbinom(sample.size,size=1, outcome.probability)
     
     #draw observed exposure
     has.current.exposure = E[,j] #logically, 1 if has and 0 if does not have
     E.observed[has.current.exposure == 1, j] <- rbinom( sum(has.current.exposure) , size=1, observation.probability[has.current.exposure==1] )  #draw observed exposure with probability modified by frailty / comorbidity status
    }
  }
  
  #plot pt frailty distribution
  hist(X)
  
  #return data in a list with E and Y 
  return.list = list("E"=E, "E.observed" = E.observed, "Y"=Y)
  return(return.list) 
}


################ RUN SIMULATIONS  #######################

#run sim 1
sim1 = sim( 0,      #time.grid.lower
            40,     #time.grid.upper
            1,      #time.step
            1000,   #sample.size
            2,      #exposure.duration
            2,      #frailty.mean 
            0.5,    #frailty.sd
            1,      #oddsratio.exposure
            1.05,   #oddsratio.X (odds ratio associeret med comorbs/frailty)
            0.001,  #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
            0.01,   #exposure.probability 
            0.5,    #baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
            1       #oddsratio for exposure observation depending on comorbidity/frailty
)

E = sim1$E
Y = sim1$Y
E.observed = sim1$E.observed

print(sum(Y)) #check hvor mange outcomes
print(sum(E)) #check hvor mange exposure 
print(sum(E.observed)) #check hvor mange observed exposures 











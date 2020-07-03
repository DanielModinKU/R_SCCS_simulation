
################ DEFINE FUNCTIONS #######################

#expit function
expit <- function(x) { 
  d = exp(x)/(1+exp(x))
  return(d)
}

#my match function skrevet fordi base R match() opfærte sig mærkeligt med apply over matrix 
my.match <- function(table, x, nomatches) {
  z = match (x, table, nomatch = nomatches )
  return(z)
}

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
  
  #lav cohort datasæt###############################
  
  #udhent exposure indices (return index / time period (for example day of exposure if time unit is days)) if pt was exposed, return index -10 if not ( removed later, set to outside study period ) ) 
  #OBS se funktion my.match, skrevet fordi der var nogle issues med rækkefølgen af argumenter som blev passet til match() i base R ved brug af apply funktionen
  E.time <- apply( E, 1, my.match, 1, -10) #vektor med time unit for exposure for hver pt (vil være dag for exposure hvis time unit er i dage ift start af studie periode eks) (-10 hvis pt ej exposed)
  E.observed.time <- apply( E.observed, 1, my.match, 1, -10) #samme som E.times, men blot for observed exposures
  Y.first.time <- apply( Y , 1, my.match, 1, -10) #vektor med time unit for outcome for hver pt (vil være dag for outcome  hvis time unit er i dage ift start af studie periode eks) (-10 hvis pt ej outcome)
  
  #konstruer matrix med 1 row = 1 pt og 3 kolonner ( Exposure time, observed exposure time og outcome time)
  data.full.cohort <- cbind(E.time, E.observed.time,Y.first.time)
  
  #add unikt pt ID
  ID <- 1:nrow(data.full.cohort)
  
  #flet på data
  data.full.cohort <- cbind(ID,data.full.cohort)
  
  #lav til dataframe 
  df <- data.frame(data.full.cohort)
  
  #return full cohort dataframe
  return(df) 
}


################ RUN SIMULATIONS  #######################

#run sim 1
sim1 = sim( 0,      #time.grid.lower
            100,     #time.grid.upper
            1,      #time.step
            1000,   #sample.size
            2,      #exposure.duration
            2,      #frailty.mean 
            0.5,    #frailty.sd
            1,      #oddsratio.exposure
            1.05,   #oddsratio.X (odds ratio associeret med comorbs/frailty)
            0.0005,  #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
            0.01,   #exposure.probability 
            0.5,    #baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
            1       #oddsratio for exposure observation depending on comorbidity/frailty
)


#case data med kun ptt der har outcome 
data.cases <- sim1[ sim1['Y.first.time']>0 , ]

data.cases

nrow(data.cases) 







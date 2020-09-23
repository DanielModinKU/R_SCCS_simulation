library(SCCS) #pakke fra farrington og co. til at fitte SCCS modeller

################ quick and dirty documentation ####################################### quick and dirty documentation  #######################
################ quick and dirty documentation ####################################### quick and dirty documentation  #######################

# INTRO #
# This file contains the function "sim" and helper functions used to generate the cohort datasets 

# USEAGE, INPUTS AND CALLING #
# Below an example call along with variable descriptions of the "sim" function can bee seen
# sim1 = sim( 0,      #time.grid.lower
#             130,     #time.grid.upper
#             1,      #time.step
#             25000,   #sample.size
#             14,      #exposure.duration
#             2,      #frailty.mean 
#             0.5,    #frailty.sd
#             1.3,      #oddsratio.exposure
#             1.15,   #oddsratio.X (odds ratio associeret med comorbs/frailty)
#             0.00015,  #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
#             0.005,   #exposure.probability 
#             0.15,    #baseline.probability.observation baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
#             1.5,       #oddsratio.observation oddsratio for exposure observation depending on comorbidity/frailty
#             #cases_only = T #optional, set true if only cases from the cohort with outcome to be returned in a df
#             exposed_cases_only = F, #optional, set true if one only wants exposed cases with outcome to be returned
#             #print_X = T #optional, to plot the comorbidity distribution 
#           )

# OUTPUT #
# The function call returns a dataframe and has options:
# 1: return entire cohort data (can be specified during call, se example call above)
# 2: return cases only (see example call)
# 3: return cases with observed exposure only (see example call)

################ quick and dirty documentation ####################################### quick and dirty documentation  #######################
################ quick and dirty documentation ####################################### quick and dirty documentation  #######################





################ DEFINE FUNCTIONS #######################

#expit function
expit <- function(x) { 
  d = exp(x)/(1+exp(x))
  return(d)
}

#my match function skrevet fordi base R match() opfærte sig mærkeligt med apply over matrix 
my.match <- function(table, match.element, nomatches) {
  z = match (match.element, table, nomatch = nomatches )
  return(z)
}

#funktion til udhent exposure indices (return index / time period (for example day of exposure if time unit is days)) if pt was exposed, return index -10 if not ( removed later, set to outside study period ) )
get_exposure_times = function(exposure.matrix) { #takes a matrix as input
  exposure.times = apply(exposure.matrix, 1, my.match, match.element = 1, nomatches = -10)
  return(exposure.times) #returns a vector with exposure time or -10 if not exposed with length = nrows in matrix and each element corresponding to exposure time or -10 for each row
}

#funktion til extract af outcome times fra outcome matrix. returnerer -10 for outcome dato hvis pt ej har haft outcome
get_outcome_times = function(outcome.matrix) {
  outcome.times = apply(outcome.matrix, 1, my.match, match.element = 1, nomatches = -10)
  return(outcome.times)
}

create_dataframe = function(...,start, end, arg_cases_only=F,arg_exposed_cases_only=F) { #take variable number of column arguemnts and return dataframe with unique ID 
  
  data = cbind(...)
  #add unique ID
  ID = 1:nrow(data)
  data = cbind(ID,data)
  data = data.frame(data)
  data['start'] = start #set start time point for SCCS analysis
  data['end'] = end #set end time point for SCCS analysis 
  
  if (arg_cases_only) {
    data = data[ data['Y.first.time'] > 0, ]
  }
  if (arg_exposed_cases_only) {
    data = data[ data['Y.first.time'] > 0 & data['E.observed.time'] > 0, ]
  }
  
  return(data)
  
}

#simulate data function 
sim <- function(time.grid.lower,  time.grid.upper,  time.step ,  sample.size,  exposure.duration,  frailty.mean,  frailty.sd ,  oddsratio.exposure,  oddsratio.X,  odds.baseline,  exposure.probability, baseline.probability.observation, oddsratio.observation, cases_only = F, exposed_cases_only =F, print_X =F ) {
  
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
    if (j>= exposure.duration) { #OBS dette  når man definerer observationsperioe i SCCS kaldet! 
      #print(j)
      
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
  if (print_X==T) {
    hist(X)
  }
  
  
  #lav cohort datasæt###############################
  
  #udhent exposure indices (return index / time period (for example day of exposure if time unit is days)) if pt was exposed, return index -10 if not ( removed later, set to outside study period ) ) 
  E.time <- get_exposure_times(E) #vektor med time unit for exposure for hver pt (vil være dag for exposure hvis time unit er i dage ift start af studie periode eks) (-10 hvis pt ej exposed)
  E.observed.time <- get_exposure_times(E.observed) #samme som E.times, men blot for observed exposures
  Y.first.time <- get_outcome_times(Y) #vektor med time unit for outcome for hver pt (vil være dag for outcome  hvis time unit er i dage ift start af studie periode eks) (-10 hvis pt ej outcome)
  
  #convert data til dataframe med unikt pt ID 
  df <- create_dataframe(E.time, 
                         E.observed.time,
                         Y.first.time, 
                         arg_cases_only = cases_only, 
                         arg_exposed_cases_only = exposed_cases_only,
                         start = time.grid.lower, #add timeline start of cohort data (impoprtant for sccs calll)
                         end = time.grid.upper) #add end 
  
  #return data
  return(df) 
}




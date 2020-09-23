library(SCCS) #pakke fra farrington og co. til at fitte SCCS modeller
library(parallel) #multi core processing
library(R6)

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

create_dataframe = function(...,start, end, arg_cases_only=F,arg_exposed_cases_only=F,arg_observed_exposed_cases_only=F) { #take variable number of column arguemnts and return dataframe with unique ID 
  
  data = cbind(...)
  #add unique ID
  ID = 1:nrow(data)
  data = cbind(ID,data)
  data = data.frame(data)
  data['start'] = start
  data['end'] = end 
  
  if (arg_cases_only) {
    data = data[ data['Y.first.time'] > 0, ]
  }
  else if (arg_exposed_cases_only) {
    data = data[ data['Y.first.time'] > 0 & data['E.time'] > 0, ]
  }
  else if (arg_observed_exposed_cases_only) {
    data = data[ data['Y.first.time'] > 0 & data['E.observed.time'] > 0, ]
  }
  
  return(data)
  
}

#simulate data function 
sim <- function(time.grid.lower,  time.grid.upper,  time.step ,  sample.size,  exposure.duration,  frailty.mean,  frailty.sd ,  oddsratio.exposure,  oddsratio.X,  odds.baseline,  exposure.probability, baseline.probability.observation, oddsratio.observation, cases_only = F, exposed_cases_only =F, observed_exposed_cases_only = F, print_X =F ) {
  
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
                         arg_observed_exposed_cases_only = observed_exposed_cases_only,
                         start = time.grid.lower, #add timeline start of cohort data (impoprtant for sccs calll)
                         end = time.grid.upper) #add end 
  
  #return data
  return(df) 
}


################ RUN SIMULATIONS  #######################
################ RUN SIMULATIONS  #######################


################ DEFINE TEMPORARY TEST FUNCTIONS  #######################

#top level function for temp testing
master_sim = function(parallel_dummy=F) {
  sim_data = run_sim()
  estimates = run_sccs(sim_data, risk_duration = 14)
  return(estimates)
}

#for temporary testing
run_sim = function(parallel_dummy=NULL) {
  #run sim to generate data 
  sim1 = sim( 0,      #time.grid.lower
              130,     #time.grid.upper
              1,      #time.step
              450000,   #sample.size
              14,      #exposure.duration
              2,      #frailty.mean 
              0.6,    #frailty.sd
              6,      #oddsratio.exposure
              1.15,   #oddsratio.X (odds ratio associeret med comorbs/frailty)
              0.00001,  #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
              0.0009,   #exposure.probability 
              0.2,    #baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
              1,       #oddsratio for exposure observation depending on comorbidity/frailty
              #cases_only = T #optional, set true if only cases from the cohort with outcome to be returned in a df
              #exposed_cases_only = T, #optional, set true if one only wants exposed cases with outcome to be returned
              observed_exposed_cases_only = T,
              print_X = F #optional, to plot the comorbidity distribution 
  )
  return(sim1)
}

#for temp testing, run sccs model
run_sccs = function(data, risk_duration = 14) {
  est1 = standardsccs(event~E.observed.time,
                      indiv=ID,
                      astart=start+risk_duration, ## OBS forklaring på "+risk_duration": se definition af for loop i kode til at genererre data: if-sætning gør at man i perioden 0-risk_duration intervallet ikke kan få exposure eller outcome, derfor startes risk period her
                      aend=end,
                      aevent=Y.first.time,
                      adrug=E.observed.time,
                      aedrug=E.observed.time+risk_duration,
                      data = data
  )
  return(est1)
}

################# simi1 ###########################

#antal cpu cores 
num_cores = detectCores()-2  #køre antal cores - 1 så computeren stadig kan bruges lidt ..
#create worker cluster
work_cluster = makeCluster(num_cores, setup_timeout = 0.5) #der er en fejl (se github: https://github.com/rstudio/rstudio/issues/6692 ) i makeCluster ift Mac OS, kræver setup_timeout < 1 for at virke, temporary fix
#export functions to workers (they cannot look up the function definition)
clusterExport(work_cluster, c('get_outcome_times','run_sccs','master_sim','run_sim','sim','standardsccs','expit','my.match','get_exposure_times','create_dataframe'))
#number of sims to run in parallel
nsims = 110


t1 = Sys.time()
######### RUN SIMS IN PARALLEL ##############
sims = parLapply(work_cluster, 1:nsims, function(x) {master_sim()})
#system.time( parLapply(work_cluster, 1:nsims, function(x) {master_sim()}))
######### RUN SIMS IN PARALLEL ##############
t2 = Sys.time()
print(paste('time to complete ',t2-t1))



#close workers
stopCluster(work_cluster) 

#extract p values of sims to a vector 
ps = sapply(sims, function(ele) {ele[[7]][5]})
#extract estimated coefficients
coefs = sapply(sims, function(ele) {ele[[7]][1]})
hist(coefs)
m = mean(coefs)
print(m)
print(exp(m))


#OR på 1.0 for exposure i sims call, så vi forventer at ca 5% af sims vil give signifikant resultat (p<0.05) pr. tillfældighed
v = ps < 0.05
proportion_sig = sum(v)/length(ps)
print(paste('Type I error rate: ',proportion_sig*100,'%'))



################# simi2 ###########################


#for temporary testing
run_sim = function(parallel_dummy=NULL) {
  #run sim to generate data 
  sim1 = sim( 0,      #time.grid.lower
              130,     #time.grid.upper
              1,      #time.step
              450000,   #sample.size
              14,      #exposure.duration
              2,      #frailty.mean 
              0.6,    #frailty.sd
              6,      #oddsratio.exposure
              1.15,   #oddsratio.X (odds ratio associeret med comorbs/frailty)
              0.00001,  #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
              0.0009,   #exposure.probability 
              0.2,    #baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
              1.15,       #oddsratio for exposure observation depending on comorbidity/frailty
              #cases_only = T #optional, set true if only cases from the cohort with outcome to be returned in a df
              #exposed_cases_only = T, #optional, set true if one only wants exposed cases with outcome to be returned
              observed_exposed_cases_only = T,
              print_X = F #optional, to plot the comorbidity distribution 
  )
  return(sim1)
}

## OBS function updated above, therefore must be exported again to workers ##

#antal cpu cores 
num_cores = detectCores()-2  #køre antal cores - 1 så computeren stadig kan bruges lidt ..
#create worker cluster
work_cluster = makeCluster(num_cores, setup_timeout = 0.5) #der er en fejl (se github: https://github.com/rstudio/rstudio/issues/6692 ) i makeCluster ift Mac OS, kræver setup_timeout < 1 for at virke, temporary fix
#export functions to workers (they cannot look up the function definition)
clusterExport(work_cluster, c('get_outcome_times','run_sccs','master_sim','run_sim','sim','standardsccs','expit','my.match','get_exposure_times','create_dataframe'))
#number of sims to run in parallel
nsims = 100



t1 = Sys.time()
######### RUN SIMS IN PARALLEL ##############
sims_2 = parLapply(work_cluster, 1:nsims, function(x) {master_sim()})
######### RUN SIMS IN PARALLEL ##############
t2 = Sys.time()
print(paste('time to complete ',t2-t1))

#close workers
stopCluster(work_cluster) 

#extract p values of sims to a vector 
ps_2 = sapply(sims_2, function(ele) {ele[[7]][5]})
#extract estimated coefficients
coefs_2 = sapply(sims_2, function(ele) {ele[[7]][1]})
hist(coefs_2)
m_2 = mean(coefs_2)
print(m_2)
print(exp(m_2))

#OR på 1.0 for exposure i sims call, så vi forventer at ca 5% af sims vil give signifikant resultat (p<0.05) pr. tillfældighed
v_2 = ps_2 < 0.05
proportion_sig_2 = sum(v_2)/length(ps_2)
print(paste('Type I error rate: ',proportion_sig_2*100,'%'))



################## SIM 3 ####################


#for temporary testing
run_sim = function(parallel_dummy=NULL) {
  #run sim to generate data 
  sim1 = sim( 0,      #time.grid.lower
              130,     #time.grid.upper
              1,      #time.step
              250000,   #sample.size
              14,      #exposure.duration
              2,      #frailty.mean 
              0.6,    #frailty.sd
              6,      #oddsratio.exposure
              1.15,   #oddsratio.X (odds ratio associeret med comorbs/frailty)
              0.00001,  #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
              0.0009,   #exposure.probability 
              0.2,    #baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
              2,       #oddsratio for exposure observation depending on comorbidity/frailty
              #cases_only = T #optional, set true if only cases from the cohort with outcome to be returned in a df
              #exposed_cases_only = T, #optional, set true if one only wants exposed cases with outcome to be returned
              observed_exposed_cases_only = T,
              print_X = F #optional, to plot the comorbidity distribution 
  )
  return(sim1)
}

## OBS function updated above, therefore must be exported again to workers ##

#antal cpu cores 
num_cores = detectCores()-2  #køre antal cores - 1 så computeren stadig kan bruges lidt ..
#create worker cluster
work_cluster = makeCluster(num_cores, setup_timeout = 0.5) #der er en fejl (se github: https://github.com/rstudio/rstudio/issues/6692 ) i makeCluster ift Mac OS, kræver setup_timeout < 1 for at virke, temporary fix
#export functions to workers (they cannot look up the function definition)
clusterExport(work_cluster, c('get_outcome_times','run_sccs','master_sim','run_sim','sim','standardsccs','expit','my.match','get_exposure_times','create_dataframe'))
#number of sims to run in parallel
nsims = 100



t1 = Sys.time()
######### RUN SIMS IN PARALLEL ##############
sims_3 = parLapply(work_cluster, 1:nsims, function(x) {master_sim()})
######### RUN SIMS IN PARALLEL ##############
t2 = Sys.time()
print(paste('time to complete ',t2-t1))

#close workers
stopCluster(work_cluster) 

#extract p values of sims to a vector 
ps_3 = sapply(sims_3, function(ele) {ele[[7]][5]})
#extract estimated coefficients
coefs_3 = sapply(sims_3, function(ele) {ele[[7]][1]})
hist(coefs_3)
m_3 = mean(coefs_3)
print(m_3)
print(exp(m_3))

#OR på 1.0 for exposure i sims call, så vi forventer at ca 5% af sims vil give signifikant resultat (p<0.05) pr. tillfældighed
v_3 = ps_3 < 0.05
proportion_sig_3 = sum(v_3)/length(ps_3)
print(paste('Type I error rate: ',proportion_sig_3*100,'%'))



################## SIM 4 ####################


#for temporary testing
run_sim = function(parallel_dummy=NULL) {
  #run sim to generate data 
  sim1 = sim( 0,      #time.grid.lower
              130,     #time.grid.upper
              1,      #time.step
              250000,   #sample.size
              14,      #exposure.duration
              2,      #frailty.mean 
              0.6,    #frailty.sd
              6,      #oddsratio.exposure
              1.15,   #oddsratio.X (odds ratio associeret med comorbs/frailty)
              0.00001,  #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
              0.0009,   #exposure.probability 
              0.2,    #baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
              5,       #oddsratio for exposure observation depending on comorbidity/frailty
              #cases_only = T #optional, set true if only cases from the cohort with outcome to be returned in a df
              #exposed_cases_only = T, #optional, set true if one only wants exposed cases with outcome to be returned
              observed_exposed_cases_only = T,
              print_X = F #optional, to plot the comorbidity distribution 
  )
  return(sim1)
}

## OBS function updated above, therefore must be exported again to workers ##

#antal cpu cores 
num_cores = detectCores()-2  #køre antal cores - 1 så computeren stadig kan bruges lidt ..
#create worker cluster
work_cluster = makeCluster(num_cores, setup_timeout = 0.5) #der er en fejl (se github: https://github.com/rstudio/rstudio/issues/6692 ) i makeCluster ift Mac OS, kræver setup_timeout < 1 for at virke, temporary fix
#export functions to workers (they cannot look up the function definition)
clusterExport(work_cluster, c('get_outcome_times','run_sccs','master_sim','run_sim','sim','standardsccs','expit','my.match','get_exposure_times','create_dataframe'))
#number of sims to run in parallel
nsims = 100



t1 = Sys.time()
######### RUN SIMS IN PARALLEL ##############
sims_4 = parLapply(work_cluster, 1:nsims, function(x) {master_sim()})
######### RUN SIMS IN PARALLEL ##############
t2 = Sys.time()
print(paste('time to complete ',t2-t1))

#close workers
stopCluster(work_cluster) 

#extract p values of sims to a vector 
ps_4 = sapply(sims_4, function(ele) {ele[[7]][5]})
#extract estimated coefficients
coefs_4 = sapply(sims_4, function(ele) {ele[[7]][1]})
hist(coefs_4)
m_4 = mean(coefs_4)
print(m_4)
print(exp(m_4))

#OR på 1.0 for exposure i sims call, så vi forventer at ca 5% af sims vil give signifikant resultat (p<0.05) pr. tillfældighed
v_4 = ps_4 < 0.05
proportion_sig_4 = sum(v_4)/length(ps_4)
print(paste('Type I error rate: ',proportion_sig_4*100,'%'))


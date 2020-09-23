
source('functions.R') #see 'functions.R' file for documentation on how to call 
library(parallel)
library(SCCS)

### write sim function wrapper ###
sim_wrapper = function(parallel_dummy=F,
                    time.grid.lower=0,
                    time.grid.upper=130,
                    time.step = 1,
                    sample.size = 20000,
                    exposure.duration = 14,
                    frailty.mean = 2,
                    frailty.sd = 0.5,
                    oddsratio.exposure = 1.3,
                    oddsratio.X = 1.15,
                    odds.baseline = 0.00015,
                    exposure.probability = 0.005,
                    baseline.probability.observation = 0.15,
                    oddsratio.observation = 1.5,
                    cases_only = F,
                    exposed.cases.only = F,
                    print_X = F) {
  
                    sim1 = sim( time.grid.lower,                    #time.grid.lower
                           time.grid.upper,                  #time.grid.upper
                           time.step,                        #time.step
                           sample.size,                      #sample.size
                           exposure.duration,                #exposure.duration
                           frailty.mean,                     #frailty.mean
                           frailty.sd,                       #frailty.sd
                           oddsratio.exposure,               #oddsratio.exposure
                           oddsratio.X,                      #oddsratio.X (odds ratio associeret med comorbs/frailty)
                           odds.baseline,                    #odds.baseline aka baseline risk (når odds er lille er odds ca lig sandsynlighed)
                           exposure.probability,             #exposure.probability
                           baseline.probability.observation, #baseline.probability.observation baseline probability of observing exposure given 0 effect of comorbidty (intercept term in log reg)
                           oddsratio.observation,            #oddsratio.observation oddsratio for exposure observation depending on comorbidity/frailty
                           #cases_only                       #optional, set true if only cases from the cohort with outcome to be returned in a df
                           exposed_cases_only = exposed.cases.only,           #optional, set true if one only wants exposed cases with outcome to be returned
                           #print_X                          #optional, to plot the comorbidity distribution 
                   )
                  return(sim1)
}

## write SCCS wrapper ##
sccs_wrapper = function(parallel_dummy=F,data,exposure_duration = 14,exposure_time_var='E.observed.time') {  #how long should the risk interval be? 
    formula = as.formula(paste('event','~',exposure_time_var,sep='')) #til at lave formula string 
    est = standardsccs(formula=formula,                          #fit model 
                      indiv=ID,                                  #unique ID of paitents in dataset
                      astart=start+exposure_duration,            ## OBS forklaring på "+risk_duration": se definition af for loop i kode til at genererre data: if-sætning gør at man i perioden 0-risk_duration intervallet ikke kan få exposure eller outcome, derfor startes risk period her
                      aend=end,                                  #slut af observationos periode
                      aevent=Y.first.time,                       #tidspunkt for event
                      adrug=E.observed.time,                     #tidspunkt for exposure start / start risikointerrval
                      aedrug=E.observed.time+exposure_duration,  #slut af exposure periode / slut af risikointerval 
                      data = data                                #dataset som model skal fittes på 
    )
    return(est)
} 



################# SET UP CLUSTERS  ###########################
################# SET UP CLUSTERS  ###########################

#antal cpu cores 
num_cores = detectCores()-1  #køre antal cores - 1 så computeren stadig kan bruges lidt ..
#create worker cluster
work_cluster = makeCluster(num_cores, setup_timeout = 0.5) #der er en fejl (se github: https://github.com/rstudio/rstudio/issues/6692 ) i makeCluster ift Mac OS, kræver setup_timeout < 1 for at virke, temporary fix
#export functions to workers (they cannot look up the function definition)
clusterExport(work_cluster, c(ls(),'standardsccs') ) #ls() list names of all obhects including functions in global space, so now we export these to workers 
#number of sims to run in parallel
nsims = 150


######### RUN SIMS TO GENERATE OR'S VARYING COMORB - OBS EXPOSURE FOR 2 DIFFERENT EXPOSURE - OUTCOME LEVELS ##############
######### RUN SIMS TO GENERATE OR'S VARYING COMORB - OBS EXPOSURE FOR 2 DIFFERENT EXPOSURE - OUTCOME LEVELS ##############

OR = 1 #bruges til at udregne coverage osv 

sims1 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 1, oddsratio.observation = 1, exposed.cases.only = T))})
coefs_sims1 = sapply(sims1, function(ele) {ele[[7]][2]})
lci_sim1 = sapply(sims1, function(ele) {ele[[8]][3]})
uci_sim1 = sapply(sims1, function(ele) {ele[[8]][4]})
coverage_sims1 = mean( OR <= uci_sim1 & OR >= lci_sim1 ) 
mean_absolute_bias1 = (1/nsims)*sum(abs(OR-coefs_sims1))
mean_bias1 = (1/nsims)*sum(OR-coefs_sims1)

sims2 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 1, oddsratio.observation = 1.5, exposed.cases.only = T))})
coefs_sims2 = sapply(sims2, function(ele) {ele[[7]][2]})
lci_sim2 = sapply(sims2, function(ele) {ele[[8]][3]})
uci_sim2 = sapply(sims2, function(ele) {ele[[8]][4]})
coverage_sims2 = mean( OR <= uci_sim2 & OR >= lci_sim2 ) 
mean_absolute_bias2 = (1/nsims)*sum(abs(OR-coefs_sims2))
mean_bias2 = (1/nsims)*sum(OR-coefs_sims2)

sims3 = parLapply(work_cluster, 1:nsims,function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 1, oddsratio.observation = 2, exposed.cases.only = T))})
coefs_sims3 = sapply(sims3, function(ele) {ele[[7]][2]})
lci_sim3 = sapply(sims3, function(ele) {ele[[8]][3]})
uci_sim3 = sapply(sims3, function(ele) {ele[[8]][4]})
coverage_sims3 = mean( OR <= uci_sim3 & OR >= lci_sim3 ) 
mean_absolute_bias3 = (1/nsims)*sum(abs(OR-coefs_sims3))
mean_bias3 = (1/nsims)*sum(OR-coefs_sims3)

sims4 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 1, oddsratio.observation = 4, exposed.cases.only = T))})
coefs_sims4 = sapply(sims4, function(ele) {ele[[7]][2]})
lci_sim4 = sapply(sims4, function(ele) {ele[[8]][3]})
uci_sim4 = sapply(sims4, function(ele) {ele[[8]][4]})
coverage_sims4 = mean( OR <= uci_sim4 & OR >= lci_sim4 ) 
mean_absolute_bias4 = (1/nsims)*sum(abs(OR-coefs_sims4))
mean_bias4 = (1/nsims)*sum(OR-coefs_sims4)

sims5 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 1, oddsratio.observation = 8, exposed.cases.only = T))})
coefs_sims5 = sapply(sims5, function(ele) {ele[[7]][2]})
lci_sim5 = sapply(sims5, function(ele) {ele[[8]][3]})
uci_sim5 = sapply(sims5, function(ele) {ele[[8]][4]})
coverage_sims5 = mean( OR <= uci_sim5 & OR >= lci_sim5 ) 
mean_absolute_bias5 = (1/nsims)*sum(abs(OR-coefs_sims5))
mean_bias5 = (1/nsims)*sum(OR-coefs_sims5)

#stopCluster(work_cluster) 

#PLOT FIGURE OF REESULTS
boxplot(coefs_sims1,coefs_sims2,coefs_sims3,coefs_sims4,coefs_sims5)


#########  CHANGE EXP - OUTCOME OR TO 8 ####################################################################################
#########  CHANGE EXP - OUTCOME OR TO 8 ####################################################################################

OR.2 = 8 #bruges til at udregne coverage osv 

sims11 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 8, oddsratio.observation = 1, exposed.cases.only = T))})
coefs_sims11 = sapply(sims11, function(ele) {ele[[7]][2]})
lci_sim11 = sapply(sims11, function(ele) {ele[[8]][3]})
uci_sim11 = sapply(sims11, function(ele) {ele[[8]][4]})
coverage_sims11 = mean( OR.2 <= uci_sim11 & OR.2 >= lci_sim11 ) 
mean_absolute_bias11 = (1/nsims)*sum(abs(OR.2-coefs_sims11))
mean_bias11 = (1/nsims)*sum(OR.2-coefs_sims11)


sims22 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 8, oddsratio.observation = 1.5, exposed.cases.only = T))})
coefs_sims22 = sapply(sims22, function(ele) {ele[[7]][2]})
lci_sim22 = sapply(sims22, function(ele) {ele[[8]][3]})
uci_sim22 = sapply(sims22, function(ele) {ele[[8]][4]})
coverage_sims22 = mean( OR.2 <= uci_sim22 & OR.2 >= lci_sim22 ) 
mean_absolute_bias22 = (1/nsims)*sum(abs(OR.2-coefs_sims22))
mean_bias22 = (1/nsims)*sum(OR.2-coefs_sims22)

sims33 = parLapply(work_cluster, 1:nsims,function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 8, oddsratio.observation = 2, exposed.cases.only = T))})
coefs_sims33 = sapply(sims33, function(ele) {ele[[7]][2]})
lci_sim33 = sapply(sims33, function(ele) {ele[[8]][3]})
uci_sim33 = sapply(sims33, function(ele) {ele[[8]][4]})
coverage_sims33 = mean( OR.2 <= uci_sim33 & OR.2 >= lci_sim33 ) 
mean_absolute_bias33 = (1/nsims)*sum(abs(OR.2-coefs_sims33))
mean_bias33 = (1/nsims)*sum(OR.2-coefs_sims33)

sims44 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 8, oddsratio.observation = 4, exposed.cases.only = T))})
coefs_sims44 = sapply(sims44, function(ele) {ele[[7]][2]})
lci_sim44 = sapply(sims44, function(ele) {ele[[8]][3]})
uci_sim44 = sapply(sims44, function(ele) {ele[[8]][4]})
coverage_sims44 = mean( OR.2 <= uci_sim44 & OR.2 >= lci_sim44 ) 
mean_absolute_bias44 = (1/nsims)*sum(abs(OR.2-coefs_sims44))
mean_bias44 = (1/nsims)*sum(OR.2-coefs_sims44)

sims55 = parLapply(work_cluster, 1:nsims, function(x) {sccs_wrapper(data=sim_wrapper(parallel_dummy=x, oddsratio.exposure = 8, oddsratio.observation = 8, exposed.cases.only = T))})
coefs_sims55 = sapply(sims55, function(ele) {ele[[7]][2]})
lci_sim55 = sapply(sims55, function(ele) {ele[[8]][3]})
uci_sim55 = sapply(sims55, function(ele) {ele[[8]][4]})
coverage_sims55 = mean( OR.2 <= uci_sim55 & OR.2 >= lci_sim55 ) 
mean_absolute_bias55 = (1/nsims)*sum(abs(OR.2-coefs_sims55))
mean_bias55 = (1/nsims)*sum(OR.2-coefs_sims55)

#PLOT FIGURE OF REESULTS
boxplot(coefs_sims11,coefs_sims22,coefs_sims33,coefs_sims44,coefs_sims55)


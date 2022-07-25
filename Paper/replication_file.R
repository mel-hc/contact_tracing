# Replication File *********************************************************#
#                                                                           #
#                                                                           #
# This file replicates figures and results for the pre-print in this folder.# 
# Code will be updated as paper is edited.                                  #
# 1) Define functions to run model & make plots                             #
# 2) Make and store plots                                                   #
# 3) Run code that informs ##s in paper                                     #
#***************************************************************************#

#### HOUSEKEEPING ####

# set working directory
setwd("/Users/melaniechitwood/Documents/GitHub/contact_tracing")
source("App/contact_tracing_v5.R")

# libraries
library(tidyverse)
library(scales)
library(viridis)

# Optimize for vaccine coverage at various levels of community detection (i), 
# contact tracing probabilities (j), and secondary attack rates (k). 

# Note that the model caclulates R0 as the produce of SAR, HiMSM_contacts, and 
# duration. To make it easier to test multiple R0s, HiMSM_contacts is set to 
# 1/21 and duration is set to 21, effectively making SAR = R0. 

result <- data.frame()

  for(j in seq(0,1, by = 0.05)){
    for(k in seq(1.2, 2, by = 0.2)){
      
  obj_fun = function(pars){ 
        z = make_params(vax = pars[1], 
                        contact_trace_prob = j,
                        SAR = k, 
                        duration = 21,  
                        HiMSM_prob.det_base = 0.1, 
                        HiMSM_prob.det_comm = 0.2,
                        HiMSM_prob.det_traced = 0.6, 
                        adh = 0.9, # adherence in detected
                        adh2 = 0.5, # adherence  in undetected
                        HiMSM_contacts =  1/21 #just making it easier to iter R0
                        )
        
        a = data.frame()
        a = bind_rows(a, calc_R(z[[1]], z[[2]], z[[3]], z[[4]]) ) %>% 
            dplyr::filter(Scenario == "Contact tracing") 
        
        return(abs(a$R - 1)) # = 0 when R = 1
      }
      
      # optimization routine
      res <- optim(par = c(0.2), #initial value
                   fn = obj_fun, # function
                   gr = NULL, 
                   method = "L-BFGS-B",
                   lower = 0, upper = 1, 
                   control = list(factr = 1e-7, 
                                  pgtol = 1e-7))
      # store results
      pre <- cbind(R0 = k*(1/21)*21, 
                   con.trac = j,
                   vax = (res[["par"]][1]), 
                   Rt_pm1 = (res[["value"]]), 
                   message = res[["convergence"]])
      
      result <- rbind(result, pre)
      }
    }


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
library(RColorBrewer)
library(gridExtra)

#### GET_R_PAPER: UMBRELLA FUNCTION ####
# 1) set variable over which to vary outcomes
# 2) call make_params to set up transition matrix
# 3) call calc_R to estimate outcomes
# 4) return data frame of R, det_frac, and R by symptom type

# think about a few scenarios that result in an R0 of 1.1

get_R_paper = function(duration = 21, # duration of infectiousness
                       HiMSM_prob.det = 0.1, # baseline detection
                       HiMSM_prob.det_traced = 0.85,
                       adh = 0.9, # adherence to isolation (% redux in contacts)
                       adh2 = 0.9) # adherence to isolation (% redux in contacts
{
  # parameters over which to vary
  param_vary = data.frame(
                expand.grid(
                  SAR = c(0.10, 0.14), 
                  HiMSM_contacts = c(7/14, 9/14),
                  HiMSM_prob.det_comm = seq(0.1, 0.7, length.out = 4),
                  vax = seq(0, 0.5, length.out = 6), 
                  contact_trace_prob = seq(0.2, 0.8, length.out = 4)
      ))
  
    # create data frame to store output 
    a = data.frame()
    
    # run model over each iteration of data frame
    for(i in 1:nrow(param_vary)){
      
      # make parameter vectors
      z = make_params(
        duration = duration,  
        HiMSM_prob.det = HiMSM_prob.det, 
        HiMSM_prob.det_traced = HiMSM_prob.det_traced, 
        adh = adh,
        adh2 = adh2, 
# varied over different iterations
        SAR = param_vary$SAR[i], 
        HiMSM_contacts = param_vary$HiMSM_contacts[i],
        HiMSM_prob.det_comm = param_vary$HiMSM_prob.det_comm[i],
        vax = param_vary$vax[i], 
        contact_trace_prob = param_vary$contact_trace_prob[i])
        
  a = bind_rows(a, calc_R(z[[1]], z[[2]], z[[3]], z[[4]]) %>% 
                      mutate(
                        # store variable values
                        SAR = param_vary$SAR[i], 
                        HiMSM_contacts = param_vary$HiMSM_contacts[i],
                        HiMSM_prob.det_comm = param_vary$HiMSM_prob.det_comm[i],
                        vax = param_vary$vax[i], 
                        contact_trace_prob = param_vary$contact_trace_prob[i]))
      
    }
    
    # Post-processing
    # calculate relative R
    # make variable labels
    #a = a %>% dplyr::group_by(SAR, HiMSM_contacts, HiMSM_prob.det_comm, vax, 
    #                          contact_trace_prob) %>%
    #  # find maximum R
    #  # then take ratio compared to this
    #  # recall eigenvalues scale linearly
    #  # so the percent reduction does not depend on the base value of R(t)
    #  dplyr::mutate(maxR = max(R), 
    #                perc_red = 100*(1-R/max(R))) %>% 
    #  ungroup() %>%
    #  # label variables
    #  dplyr::mutate(program = factor(Scenario, 
    #                                 levels = c("No contact tracing", 
    #                                            "Contact tracing")),
    #    ctract = paste(contact_trace_prob*100, "% contacts traced"))
    return(a)
    
  }

all_output <- get_R_paper(adh = 1, adh2 = 1) %>% 
              mutate(R_scenario = case_when(
                SAR == 0.1 & HiMSM_contacts == (7/14) ~  "A. R0 = 1.05", 
                SAR == 0.1 & HiMSM_contacts == (9/14) ~  "B. R0 = 1.35", 
                SAR == 0.14 & HiMSM_contacts == (7/14) ~ "C. R0 = 1.47", 
                SAR == 0.14 & HiMSM_contacts == (9/14) ~ "D. R0 = 1.89"),
                HiMSM_contacts = HiMSM_contacts*14) 
# old
ggplot(filter(all_output, 
              Scenario != "No contact tracing", 
              HiMSM_prob.det_comm == 0.3),
       aes(x = vax, y = R)) +
  geom_line(aes(color = as.factor(contact_trace_prob)#,
                #alpha = as.factor(contact_trace_prob)
  )) + 
  facet_wrap(.~R_scenario) + 
  geom_hline(yintercept = 1, color = "forest green", size = 0.5) + 
  theme_bw() + 
  labs(y = "Estiamted R Effective", 
       x = "Vaccine Coverage (assume 100% efficacy)", 
       #color = "Community Detection Rate (High Contact)", 
       color = "Probability that Contacts are Traced")






# optimization attempt (it does not do what we want it to)
result <- data.frame()

for(i in seq(0.1, 0.8, by = 0.1)){
  for(j in seq(0,1, by = 0.2)){
    for(k in seq(1.2, 1.8, by = 0.1)){
      
  obj_fun = function(pars){ 
        # make parameter vectors
        z = make_params(vax = pars[1], 
                        HiMSM_prob.det_comm = i,
                        contact_trace_prob = j,
                        SAR = k, 
                        duration = 21,  
                        HiMSM_prob.det = 0.1, 
                        HiMSM_prob.det_traced = 0.85, 
                        adh = 0.9,
                        adh2 = 0.9, 
                        HiMSM_contacts =  1/21#just making it easier to iter R0
                        )
        
        a = data.frame()
        a = bind_rows(a, calc_R(z[[1]], z[[2]], z[[3]], z[[4]]) ) %>% 
          dplyr::mutate(program = ifelse(Scenario == "Contact tracing and increased testing", 
                                         "Test all", "No contact tracing")) %>%
          dplyr::filter(program == "Test all") 
        
        return(abs(a$R - 1))
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
                   comm.det = i,
                   vax = (res[["par"]][1]), 
                   Rt_pm1 = (res[["value"]]), 
                   message = res[["convergence"]])
      
      result <- rbind(result, pre)
      }
    }
  }



# placeholder figures
ggplot(filter(result, comm.det == "0.3"), 
       aes(y = vax, x = R0)) + 
  geom_line(aes(color = as.factor(con.trac)), size = 1) +
  theme_bw()

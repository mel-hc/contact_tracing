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
                        HiMSM_prob.det_traced = 0.5, 
                        adh = 0.9, # adherence in detected
                        HiMSM_contacts =  1/21, 
                        rel_trans = 0.5
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


# figure, by Jiye Kwon
pal = c("#2D1160FF",  "#51127CFF", "#932B80FF", "#B63679FF", "#F1605DFF")

ggplot(filter(result, vax > 0,),
       aes(x = con.trac, y = vax)) + 
  geom_line(aes(color = as.factor(R0)), size = 1,  show.legend = F) +
  scale_color_manual(values = rev(pal))+
  geom_label(
    aes(label = paste0("R0 = ", R0), 
        color = as.factor(R0)),
    data = filter(result,vax > 0,
                  con.trac == 0.5),
    size = 3, show.legend = F) + 
  theme_bw() + 
  labs(y = "Vaccine Coverage", 
       x = "Fraction of Cases Contact Traced")


# sensitivity analyses

sens1 <- data.frame()
for(i in c(0.1, 0.2, 0.4))
for(j in seq(0,1, by = 0.05)){
  for(k in seq(1.2, 2, by = 0.2)){
    
    obj_fun = function(pars){ 
      z = make_params(vax = pars[1], 
                      contact_trace_prob = j,
                      SAR = k, 
                      duration = 21,  
                      HiMSM_prob.det_base = 0.1, 
                      HiMSM_prob.det_comm = i,
                      HiMSM_prob.det_traced = 0.5, 
                      adh = 0.9, # adherence in detected
                      HiMSM_contacts =  1/21, 
                      rel_trans = 0.5
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
                 comm.det = i,
                 vax = (res[["par"]][1]), 
                 Rt_pm1 = (res[["value"]]), 
                 message = res[["convergence"]])
    
    sens1 <- rbind(sens1, pre)
  }
}

ggplot(filter(sens1, vax > 0,),
       aes(x = con.trac, y = vax)) + 
  geom_line(aes(color = as.factor(R0)), size = 1,  show.legend = F) +
  scale_color_manual(values = rev(pal))+
  facet_wrap(.~comm.det, 
             labeller = labeller(comm.det = c(
               "0.1" = "10% community detection", 
               "0.2" = "20% community detection", 
               "0.4" = "40% community detection"
             ))) +
  geom_label(
    aes(label = paste0("R0 = ", R0), 
        color = as.factor(R0)),
    data = filter(sens1,vax > 0,
                  con.trac == 0.25),
    size = 3, show.legend = F) + 
  theme_bw() +  
  labs(y = "Vaccine Coverage", 
       x = "Fraction of Cases Contact Traced")


sens2 <- data.frame()
for(i in c(0.5, 0.7, 0.9)){
for(j in seq(0,1, by = 0.05)){
  for(k in seq(1.2, 2, by = 0.2)){
    
    obj_fun = function(pars){ 
      z = make_params(vax = pars[1], 
                      contact_trace_prob = j,
                      SAR = k, 
                      duration = 21,  
                      HiMSM_prob.det_base = 0.1, 
                      HiMSM_prob.det_comm = 0.2,
                      HiMSM_prob.det_traced = i, 
                      adh = 0.9, # adherence in detected
                      HiMSM_contacts =  1/21, 
                      rel_trans = 0.5
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
                 prob.det.trac = i,
                 vax = (res[["par"]][1]), 
                 Rt_pm1 = (res[["value"]]), 
                 message = res[["convergence"]])
    
    sens2 <- rbind(sens2, pre)
  }
 }
}
pal2 =  c("#2D1160FF",  "#932B80FF",  "#F1605DFF")
ggplot(filter(sens2, R0 == "1.6", vax > 0),
       aes(x = con.trac, y = vax)) + 
  geom_line(aes(color = as.factor(prob.det.trac)), 
            size = 0.75,  show.legend = T) +
  scale_color_manual(values = rev(pal2))+
  geom_label(aes(x = 0.9, y = 0.35), 
            label = "R0 = 1.6") +
  theme_bw() + 
  scale_y_continuous(limits = c(0,0.35)) +
  theme(legend.position = "bottom") +
  labs(y = "Vaccine Coverage", 
       x = "Fraction of Cases Contact Traced", 
       color = "Probablity a Traced Contact is Detected")

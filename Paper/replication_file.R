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
library(forcats)
library(scales)
library(viridis)

# Optimize for vaccine coverage at various levels of community detection (i), 
# contact tracing probabilities (j), secondary attack rates (k), and vaccine
# efficacy (l). 

# Note that the model calculates R0 as the product of SAR, HiMSM_contacts, and 
# duration. To make it easier to test multiple R0s, HiMSM_contacts is set to 
# 1/21 and duration is set to 21, effectively making SAR = R0. 

result <- data.frame()
for(i in c(0.1, 0.2, 0.4)){ # community detection 
 for(j in seq(0,1, by = 0.05)){ # contact tracing
  for(k in seq(1.4, 2, by = 0.2)){ # R0
    for(l in c(0.5, 0.75, 0.85, 1)){ # vaccine efficacy
    
    obj_fun = function(pars){ 
      z = make_params(vax = pars[1], 
                      contact_trace_prob = j,
                      SAR = k, 
                      vax.eff = l, 
                      duration = 21,  
                      HiMSM_prob.det_base = 0.1, 
                      HiMSM_prob.det_comm = i,
                      HiMSM_prob.det_traced = 0.5, 
                      adh = 0.8, # adherence in detected AND traced == 90% redux
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
                 vax.eff = l,
                 vax = (res[["par"]][1]), 
                 Rt_pm1 = (res[["value"]]),  # check that this is < 1e-14
                 message = res[["convergence"]])
    
    result <- rbind(result, pre)
      }
    }
  }
}

pal = c("#2D1160FF",  "#51127CFF", "#932B80FF", "#B63679FF", "#F1605DFF")

# Main text figure
ggplot(filter(result, vax > 0, vax.eff == 1),
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
    data = filter(result,vax > 0, vax.eff == 1,
                  con.trac == 0.25),
    size = 3, show.legend = F) + 
  theme_bw() +  
  labs(y = "Vaccine Coverage", 
       x = "Fraction of Cases Contact Traced")


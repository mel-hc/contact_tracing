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

# Optimize for vaccine coverage at various levels of community detection (i), 
# contact tracing probabilities (j), and secondary attack rates (k). 

# Note that the model caclulates R0 as the produce of SAR, HiMSM_contacts, and 
# duration. To make it easier to test multiple R0s, HiMSM_contacts is set to 
# 1/21 and duration is set to 21, effectively making SAR = R0. 

result <- data.frame()

for(i in seq(0.1, 0.4, by = 0.1)){
  for(j in seq(0,1, by = 0.1)){
    for(k in seq(1.2, 2, by = 0.1)){
      
  obj_fun = function(pars){ 
        z = make_params(vax = pars[1], 
                        HiMSM_prob.det_comm = i,
                        contact_trace_prob = j,
                        SAR = k, 
                        duration = 21,  
                        HiMSM_prob.det_base = 0.1, 
                        HiMSM_prob.det_traced = 0.8, 
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
                   comm.det = i,
                   vax = (res[["par"]][1]), 
                   Rt_pm1 = (res[["value"]]), 
                   message = res[["convergence"]])
      
      result <- rbind(result, pre)
      }
    }
  }

# placeholder figures
p1 <- ggplot(filter(result, 
                    (con.trac == "0" | con.trac == "0.3" | con.trac == "0.6" | 
                      con.trac == "0.9"),
                    (comm.det == "0.2" | comm.det == "0.4"),
                    vax > 0),
             aes(y = R0, x = vax)) + 
  geom_line(aes(color = as.factor(con.trac)), size = 1, show.legend = F) +
  facet_wrap(.~comm.det, 
             labeller = 
               labeller(comm.det = 
                          c(#"0.1" = "10% Community Detection", 
                            "0.2" = "20% Community Detection",
                            #"0.3" = "30% Community Detection",
                            "0.4" = "40% Community Detection"))) + 
  geom_label(aes(label = con.trac, 
                 color = as.factor(con.trac)),
             data = filter(result, (con.trac == "0" | con.trac == "0.3" | 
                                      con.trac == "0.6" | con.trac == "0.9"), 
                           (comm.det == "0.2" | comm.det == "0.4"), R0 == 1.6),
             size = 2, show.legend = F) +
  theme_bw() + 
  labs(title = "Frontier where effective reproduction number (Rt) equals 1", 
       y = "Basic Reproduction Number (R0)",
       x = "Vaccine Coverage", 
       subtitle = "By fraction of cases contact traced", 
       caption = "Assumes 80% of traced individual's secondary infections are detected")

p2 <- ggplot(filter(result, 
                    (R0 == "1.2" | R0 == "1.4" | R0 == "1.6" | R0 == "1.8" |
                    R0 == "2"), (comm.det == "0.2" | comm.det == "0.4"),
                    vax > 0),
       aes(y = con.trac, x = vax)) + 
  geom_line(aes(color = as.factor(R0)), size = 1,  show.legend = F) +
  facet_wrap(.~comm.det, 
             labeller = 
             labeller(comm.det = 
                        c(#"0.1" = "10% Community Detection", 
                          "0.2" = "20% Community Detection",
                          #"0.3" = "30% Community Detection",
                          "0.4" = "40% Community Detection"))) + 
  geom_label(
    aes(label = R0, 
        color = as.factor(R0)),
    data = filter(result, 
           (R0 == "1.2" | R0 == "1.4" | R0 == "1.6" | R0 == "1.8" |
              R0 == "2"), (comm.det == "0.2" | comm.det == "0.4"),
           con.trac == 0.5),
    size = 2, show.legend = F) + 
  theme_bw() + 
  labs(title = "Frontier where effective reproduction number (Rt) equals 1", 
       subtitle = "By basic reproduction number (R0)",
       x = "Vaccine Coverage", 
       y = "Fraction of Cases Contact Traced", 
       caption = "Assumes 80% of traced individual's secondary infections are detected")

p3 <- ggplot(filter(result, 
                    (con.trac == "0.3" | con.trac == "0.6"),
                    vax > 0),
             aes(y = R0, x = vax)) + 
  geom_line(aes(color = as.factor(comm.det)), size = 1,  show.legend = F) +
  facet_wrap(.~con.trac, 
             labeller = 
               labeller(con.trac = 
                          c(#"0.1" = "10% Community Detection", 
                            "0.3" = "30% Contact Tracing",
                            #"0.3" = "30% Community Detection",
                            "0.6" = "60% Contact Tracing"))) + 
  geom_label(
    aes(label = comm.det, 
        color = as.factor(comm.det)),
    data = filter(result, 
                  (con.trac == "0.3" | con.trac == "0.6"), R0 == 1.6),
    size = 2, show.legend = F) + 
  theme_bw() + 
  labs(title = "Frontier where effective reproduction number (Rt) equals 1", 
       subtitle = "By community detection rate",
       x = "Vaccine Coverage", 
       y = "Basic Reproduciton number", 
       caption = "Assumes 80% of traced individual's secondary infections are detected")


pdf("test.plots.pdf", width = 8, height = 6)
p2
p1
p3
dev.off()




#************************************ Replication File *********************#
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
get_R_paper = function(SAR = 0.2,
                       HiMSM_contacts = 10/14, # 10 contacts in 14 days
                       LoMSM_contacts = 2/14,  # 2 contacts in 14 days
                       duration = 21, # Estimated time from lesion appearance 
                                      # to fully healed (no scabs)
                       LoMSM_prob = 0.6, # very made up
                       # the above parameters work out to an R0=1.5 for MSM pop.
                       comparator = "Contact tracing only",
                       baseline_HiMSM_prob.det = 0.1, 
                       baseline_LoMSM_prob.det = 0.05, 
                       test_uptake = 0.5, 
                       rel_trans = 0.5)
{
  # parameters over which to vary
  param_vary = data.frame(
                expand.grid(
                  HiMSM_prob.det = seq(0.2, 0.4, length.out = 3),
                  LoMSM_prob.det = seq(0.1, 0.3, length.out = 3),                    
                  adh = c(0.75, 0.8, 0.95),
                  adh2 = c(0.75, 0.8, 0.95),
                  vax = c(0.1, 0.3, 0.5), # FIX LATER, should increase over time
                  contact_trace_prob = seq(0.1, 0.7, length.out = 4))) 
  
    # create data frame to store output 
    a = data.frame()
    
    # run model over each iteration of data frame
    for(i in 1:nrow(param_vary)){
      
      # make parameter vectors
      z = make_params(
        SAR = SAR, 
        HiMSM_contacts = HiMSM_contacts, 
        LoMSM_contacts = LoMSM_contacts,
        duration = duration,  
        LoMSM_prob = LoMSM_prob,
        comparator = comparator, 
        baseline_HiMSM_prob.det = baseline_HiMSM_prob.det, 
        baseline_LoMSM_prob.det = baseline_LoMSM_prob.det,  
        test_uptake = test_uptake, 
        rel_trans = rel_trans,
# varied over different iterations
        HiMSM_prob.det = param_vary$HiMSM_prob.det[i], 
        LoMSM_prob.det = param_vary$LoMSM_prob.det[i], 
        adh = param_vary$adh[i],
        adh2 = param_vary$adh2[i], 
        vax = param_vary$vax[i], 
        contact_trace_prob = param_vary$contact_trace_prob[i])
        
  a = bind_rows(a, calc_R(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]],z[[6]]) %>% 
                      mutate(
                        # store variable values
                        HiMSM_prob.det = param_vary$HiMSM_prob.det[i], 
                        LoMSM_prob.det = param_vary$LoMSM_prob.det[i], 
                        adh = param_vary$adh[i],
                        adh2 = param_vary$adh2[i], 
                        vax = param_vary$vax[i], 
                        contact_trace_prob = param_vary$contact_trace_prob[i]))
      
    }
    
    # Post-processing
    # calculate relative R
    # make variable labels
    a = a %>%
      dplyr::group_by(HiMSM_prob.det, contact_trace_prob, adh) %>%
      # find maximum R
      # then take ratio compared to this
      # recall eigenvalues scale linearly
      # so the percent reduction does not depend on the base value of R(t)
      dplyr::mutate(maxR = max(R), perc_red = 100*(1-R/max(R)))  %>% ungroup() %>%
      # label variables
      dplyr::mutate(program = ifelse(Scenario == "Contact tracing\n(Test all)", 
                                     "Test all", "No contact tracing"),
                    program = ifelse(Scenario == "Contact tracing\n(Test symptomatic)", 
                                     "Test symptomatic", program),
                    program = factor(program, levels = c("No contact tracing", 
                                                         "Test symptomatic", 
                                                         "Test all")),
        var = paste(adh*100, "% reduction in transmission due to detection", sep = ""))
    return(a)
    
  }

## JIYE TO MAKE ADDITIONAL CHANGES HERE

#### MAKE HEATMAP ####
make_heatmap = function(R, title, perc_asymp = 0.4, save = T, show_legend = T) {
  
  # make color palette
  pal = brewer.pal(9, "Greens")[1:9]
  
  # make the plot
  plot = ggplot(R %>% filter(Scenario != "No contact \ntracing" & A_prob == perc_asymp), 
         aes(x = S_prob.det, y = contact_trace_prob, fill = perc_red)) + geom_tile() +
    geom_text(aes(label = round(perc_red))) + 
    theme(axis.line = element_blank()) + 
    scale_fill_gradient(name = expression("Percentage reduction in "*R[t]), low = pal[1], high = pal[8]) + 
    labs(x = "Fraction of symptomatic cases detected in community", y = "Fraction of contacts successfully traced") + 
    facet_grid(var~program)
  
  if(save){
    # store output as png
    png(paste0(title, ".png"),width=11, height=9, units = "in", res = 300)
      print(plot)
    dev.off()
  }
  
  # return plot 
  return(plot)
  
}

#### MAKE LINEGRAPH ####
make_linegraph = function(R, title, perc_asymp = .4, save = T){
  
  # make color palette
  pal = brewer.pal(9, "Blues")[c(2,3,5,7,9)]
  
  # select lines to show
  ctrace_keep = unique(R$contact_trace_prob)[c(1,3,5,7,9)]
  
  # make plot
  plot = ggplot(R %>% filter(Scenario != "No contact \ntracing" & A_prob == .4 & 
                         (contact_trace_prob%in%ctrace_keep)), 
         aes(x = S_prob.det, y = perc_red, group = factor(contact_trace_prob), 
             col = factor(contact_trace_prob))) + geom_line() +
    theme(axis.line = element_blank()) + 
    scale_color_manual(name = "Fraction of contacts\nsuccessfully traced", values= pal) + 
    labs(x = "Fraction of symptomatic cases detected in community", y = expression("Percentage reduction in "*R[t])) + 
    facet_grid(var~program) +  guides(colour = guide_legend(reverse=T))
  
  if(save){
    # store output as png
    png(paste0("base_case.png"),width=9, height=6.5, units = "in", res = 300)
      print(plot)
    dev.off()
  }
  
  # return plot
  return(plot)
  
}

#### FIGURES ####

# base case
  # run model
  f = get_R_paper()
  
  # base case heatmap
  h1 = make_heatmap(f, "heatmap_base_case")
  l1 = make_linegraph(f, "linegraph_base_case")
  
  # asymptomatic sensitivity analysis
  h2 = make_heatmap(f, "heatmap_sens_asymp", perc_asymp = .2)
  l2 = make_linegraph(f, "linegraph_sens_asymp", perc_asymp = .2)

# tracing + testing scale-up
  # run model
  g = get_R_paper(comparator = "Scale-up")
  h3 = make_heatmap(g, "heatmap_sens_test")
  l3 = make_linegraph(g, "linegraph_sens_test")
  
# varying multiple on transmission
  i = get_R_paper(rel_trans = .75)
  h4 = make_heatmap(i, "heatmap_rel_trans")

#### PAPER NUMBERS ####

  # filter to be > 50% testing of symptomatics + > 50% contact tracing
  # set at base case for asymptomatics (40%)
  f2 = f %>% filter(Scenario!="No contact \ntracing" & S_prob.det > 0.5 & contact_trace_prob > 0.5 & A_prob==.4)
  
  # increase from testing
  print(f2 %>% 
    group_by(contact_trace_prob, S_prob.det, adh, A_prob) %>% 
    dplyr::summarize(out = perc_red[2]/perc_red[1]) %>% ungroup() %>% 
      summarize(median(out), quantile(out, .25), quantile(out, .75)))
  
  # increase from testing by isolation and quarantine efficacy
  print(f2 %>% 
          group_by(contact_trace_prob, S_prob.det, adh, A_prob) %>% 
          dplyr::summarize(out = perc_red[2]/perc_red[1]) %>% group_by(adh) %>% 
          dplyr::summarize(median(out)))
  
  # benefits by isolation and quarantine efficacy
  print(f2 %>% 
          group_by(Scenario, adh) %>% 
          dplyr::summarize(median(perc_red)))
  
  # compare to lower asymptomatic probability
  print(f %>% filter(Scenario!="No contact \ntracing" & S_prob.det > 0.5 & contact_trace_prob > 0.5) %>% 
    group_by(contact_trace_prob, S_prob.det, adh, Scenario) %>% 
    dplyr::summarize(out = perc_red[1]/perc_red[2]) %>% ungroup() %>% 
    summarize(median(out), quantile(out, .25), quantile(out, .75)))
  

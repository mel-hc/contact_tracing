### CONTACT TRACING MODEL### 

#### GET_R: UMBRELLA FUNCTION ####
# 1) set variable over which to vary outcomes
# 2) call make_params to set up transition matrix
# 3) call calc_R to estimate outcomes
# 4) return data frame of R, det_frac, and R by symptom type

get_R = function(P_RR, # relative infectiousness presymptomatic to symptomatic
                 P_dur, # duration of infectiousness, presymptomatic
                 HiMSM_RR, # infectiousness of symptomatic group (set = 1)
                 HiMSM_dur, # duration of infectiousness, symptomatic
                 LoMSM_RR, # relative infectiousness asymptomatic to symptomatic
                 LoMSM_dur, # duration of infectiousness, asymptomatic 
                 HiMSM_prob.det, # detection probability, symptomatic
                 LoMSM_prob.det, # detection probability, asymptomatic
                 LoMSM_prob, # fraction of infections that are asymptomatic 
                 contact_trace_prob, # probability of contact tracing 
                 comparator, # character string, assigns a scenario
                 baseline_HiMSM_prob.det, # probability detection w/o tracing (sym)
                 baseline_LoMSM_prob.det, # probability detection w/o tracing (asym)
                 test_uptake, # probability tested if traced 
                 adh, # adherence to isolation measures (% redux in contacts)
                 adh2, # adherence to isolation measures (% redux in contacts)
                 rel_trans, # relative number of secondary infections (detected 
                            # compared to undetected)
                 xaxis){
  
  # set x-axis
  var = c(eval(parse(text = xaxis)), seq(.1, .9, length.out = 9), .95)
  
  # create data frame to store output 
  a = data.frame()

  # run model over each iteration of data frame
  for(i in 1:length(var)){
    
    # set relevant variable to outcome
    assign(xaxis, var[i])
    
    # make parameter vectors
    z = make_params(P_RR, 
                    P_dur, 
                    HiMSM_RR, 
                    HiMSM_dur, 
                    LoMSM_RR, 
                    LoMSM_dur, 
                    HiMSM_prob.det,
                    LoMSM_prob.det,
                    LoMSM_prob, 
                    contact_trace_prob,
                    comparator, 
                    baseline_HiMSM_prob.det, 
                    baseline_LoMSM_prob.det, 
                    test_uptake, 
                    adh, 
                    adh2, 
                    rel_trans)
    
    a = bind_rows(a, 
                  calc_R(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]],z[[6]]) %>% 
                    mutate(
                      # store variable values
                      var = var[i], 
                      # note if these were selected inputs
                      point = ifelse(i == 1, "point", "omit")))
    
  }
  
  return(a)
  
}

#### make_params: Make parameter set associated with each scenario ####
# 1) Take in model inputs (see above)
# 2) Return set of parameters for each scenario
  # Scenarios: 
  # Counterfactual: No contact tracing
  # Base Case: Symptomatic testing, contact tracing
  # 
make_params = function(P_RR, # relative infectiousness presym to symp
                       P_dur, # duration of infectiousness, presymp
                       HiMSM_RR, # infectiousness of symptomatic group (set = 1)
                       HiMSM_dur, # duration of infectiousness, symptomatic
                       LoMSM_RR, # relative infectiousness asym to sym
                       LoMSM_dur, # duration of infectiousness, asymptomatic 
                       HiMSM_prob.det, # detection probability, symptomatic
                       LoMSM_prob.det, # detection probability, asymptomatic
                       LoMSM_prob, # fraction of infections that are asymptomatic 
                       contact_trace_prob, # probability of contact tracing 
                       comparator, # character string, assigns a scenario
                       baseline_HiMSM_prob.det, # prob detection w/o tracing (sym)
                       baseline_LoMSM_prob.det,  # prob detection w/o tracing (asym)
                       test_uptake, # probability tested if traced 
                       adh, # adherence to isolation (% redux in contacts)
                       adh2, # adherence to isolation (% redux in contacts)
                       rel_trans # Redux in number of secondary infections, 
                                 # detected vs. undetected
                       ){
  # BASE CASE 
  # Here, 'U' refers to those who are not detected and 'D' refers to those 
  # who are detected. In most cases, new param name is for consistency. 
  params = data.frame(P_RR, 
                      P_dur, 
                      HiU_RR = HiMSM_RR, 
                      HiU_dur = HiMSM_dur, 
                      HiD_RR = HiMSM_RR*rel_trans, 
                      HiD_dur = HiMSM_dur,
                      LoU_RR = LoMSM_RR, 
                      LoU_dur = LoMSM_dur,
                      LoD_RR = LoMSM_RR*rel_trans, 
                      LoD_dur = LoMSM_dur,
                      HiMSM_prob.det, 
                      LoMSM_prob.det, 
                      LoMSM_prob,
                      symp_contact_trace_prob = contact_trace_prob, 
                      asymp_contact_trace_prob = contact_trace_prob)
  
  # COUNTERFACTUAL: base case with no contact tracing
  params_cf = params %>% 
    mutate(symp_contact_trace_prob = 0, 
           asymp_contact_trace_prob = 0,
           HiMSM_prob.det = ifelse(comparator=="Contact tracing only", 
                               HiMSM_prob.det, 
                               baseline_HiMSM_prob.det), 
           LoMSM_prob.det = ifelse(comparator=="Contact tracing only", 
                               LoMSM_prob.det, 
                               baseline_LoMSM_prob.det))
  
  # CONTACT TRACING GEN 1
  # 1) Increase testing of symptomatic contacts (test_uptake)
  # 2) Decrease transmission according to adherence (adh)
  params_ctrace_1 = params %>% 
    mutate(HiMSM_prob.det = test_uptake,
           P_RR       = P_RR*(1-adh),
           HiU_RR      = HiU_RR*(1-adh),
           HiD_RR      = HiD_RR*(1-adh),
           LoU_RR      = LoU_RR*(1-adh), 
           LoD_RR      = LoD_RR*(1-adh))
  
  # CONTACT TRACING GEN 2
  # 1) Increase testing of symptomatic contacts (test_uptake)
  # 2) Decrease transmission according to adherence (adh2)
  params_ctrace_2plus = params %>% 
    mutate(HiMSM_prob.det = test_uptake,
           P_RR       = P_RR*(1-adh2),
           HiU_RR      = HiU_RR*(1-adh2),
           HiD_RR      = HiD_RR*(1-adh2),
           LoU_RR      = LoU_RR*(1-adh2),
           LoD_RR      = LoD_RR*(1-adh2))
  
  # TEST ASYMPTOMATICS GEN 1
  # 1) Increase asymptomatic testing
  # 2) Allow  transmission reduction for detected presymptomatics 
  params_test_ctrace_1 = params_ctrace_1 %>% 
    mutate(LoMSM_prob.det = test_uptake, 
           P_RR = P_RR*.5) # drawing from uniform dist, reduce by half???
  
  # TEST ASYMPTOMATICS GEN 2
  # Same as gen 1, but adapting params_ctrace_2plus
  params_test_ctrace_2_plus = params_ctrace_2plus %>% 
    mutate(LoMSM_prob.det = test_uptake, 
           P_RR = P_RR*.5)
  
  # RETURN OUTPUT
  return(list(params_cf, params, params_ctrace_1, params_ctrace_2plus, 
              params_test_ctrace_1, params_test_ctrace_2_plus))
}

#### calc_R: run over each scenario ####
# 1) No contact tracing (params_cf)
# 2) Test symptomatic (params and params_ctrace_XX)
# 3) Test all (params and params_test_ctrace_XX)
calc_R = function(params_cf,# NO CONTACT TRACING (counterfactual)
                 
                  params, # BASE CASE (first generation for CT scenarios)
                  
                  params_ctrace_1, #CT and test symptomatic
                  params_ctrace_2plus, #CT and test symptomatic
                  
                  params_test_ctrace_1, #CT and test all
                  params_test_ctrace_2_plus #CT and test all 
                  ) {
  
  # run different methods
  out = bind_rows(dom_eigen(params_cf, params_cf, params_cf) %>% 
                    mutate(Scenario = "No contact \ntracing"),
                  dom_eigen(params, params_ctrace_1, params_ctrace_2plus) %>% 
                    mutate(Scenario = "Contact tracing\n(Test symptomatic)"),
                  dom_eigen(params, params_test_ctrace_1, 
                            params_test_ctrace_2_plus) %>% 
                    mutate(Scenario = "Contact tracing\n(Test all)"))
  return(out) 
}


#### dom_eigen: Make transition matrix and calculate R + steady state ####
# 1) Call transition probs for each case
# 2) Make transition matrix
# 3) Take eigenvalues
# 4) Take eigenvectors of transpose
# 5) Return data from of relevant quantities
# Note that we don't care where R(t) starts in this function
# because the relative R(t) will be the same, no matter what the baseline
# i.e. eigenvalues scale linearly when multiplied by a constant
# so we adjust to user-input R(t) later on in make_plots
dom_eigen = function(params, params_ctrace_1, params_ctrace_2plus){
  
  # call transition probs
  x = get_trans_probs(params, first_gen = T, ctrace = F)
  y = get_trans_probs(params_ctrace_1, first_gen = F, ctrace = T) 
  z = get_trans_probs(params_ctrace_2plus,  first_gen = F, ctrace = T)
  
  # put in matrix
  mat = matrix(unlist(c(rep(x, 9), 
                        rep(y, 3),
                        rep(z, 3))), ncol = 15, byrow = T)
  
  # estimate eigenvectors/values
  vec = Re(eigen(t(mat))$vectors[,1])
  
  # return output
  return(data.frame(
    # R(t)
    R = max(Re(eigen(mat)$values)),
    
    # detection fraction
    det_frac = sum(vec[3], vec[6], vec[9], vec[13:15])/sum(vec[2:3], 
                   vec[5:6], vec[8:9], vec[10:15]),
    
    # transmission by symptom status
    presymp = sum(vec[1:3],vec[10], vec[13])/sum(vec),
    symp = sum(vec[4:6], vec[11], vec[14])/sum(vec),
    asymp = sum(vec[7:9], vec[12], vec[15])/sum(vec)))
}


#### get_trans_probs: Pull together transition probabilities ####
# by symptom status, detection status, whether originated from contact tracing 
# + gen, whether traced

# This is our next generation matrix -- a way of representing our compartments
# as a matrix. In each compartment, what is happening in the next generation?

# This maps onto an estimate of Rt -- how many infections per infected in the 
# next generation. The matrix lets us aggregate of multiple groups

# see wiki pages on next generation matrices and eigenvalues 

# Here, we calculate the transition probabilities -- e.g. what should go in 
# the next generation matrix. It should be 6x6 (double check)
get_trans_probs = function(params, first_gen = F, ctrace = F) {
    
    # pre-symptomatic
    psymp_D_T_1  =
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$P_RR*params$P_dur)*params$symp_contact_trace_prob*first_gen
    psymp_D_T_2  = 
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$P_RR*params$P_dur)*params$symp_contact_trace_prob*(1-first_gen)
    psymp_D_NT_noctrace = 
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$P_RR*params$P_dur)*(1-params$symp_contact_trace_prob)*(1-ctrace)
    psymp_D_NT_ctrace = 
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$P_RR*params$P_dur)*(1-params$symp_contact_trace_prob)*ctrace
    psymp_U = 
      (1-params$LoMSM_prob)*(1-params$HiMSM_prob.det)*(params$P_RR*params$P_dur)
  
    # symptomatic
    symp_D_T_1  = 
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$HiD_RR*params$HiD_dur)*params$symp_contact_trace_prob*first_gen
    symp_D_T_2  = 
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$HiD_RR*params$HiD_dur)*params$symp_contact_trace_prob*(1-first_gen)
    symp_D_NT_noctrace = 
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$HiD_RR*params$HiD_dur)*(1-params$symp_contact_trace_prob)*(1-ctrace)
    symp_D_NT_ctrace = 
      (1-params$LoMSM_prob)*(params$HiMSM_prob.det)*
      (params$HiD_RR*params$HiD_dur)*(1-params$symp_contact_trace_prob)*(ctrace)
    symp_U = 
      (1-params$LoMSM_prob)*(1-params$HiMSM_prob.det)*(params$HiU_RR*params$HiU_dur)
    
    # asymptomatic
    asymp_D_T_1  = 
      (params$LoMSM_prob)*(params$LoMSM_prob.det)*
      (params$LoD_RR*params$LoD_dur)*params$asymp_contact_trace_prob*first_gen
    asymp_D_T_2  = 
      (params$LoMSM_prob)*(params$LoMSM_prob.det)*
      (params$LoD_RR*params$LoD_dur)*params$asymp_contact_trace_prob*(1-first_gen)
    asymp_D_NT_noctrace = 
      (params$LoMSM_prob)*(params$LoMSM_prob.det)*
      (params$LoD_RR*params$LoD_dur)*(1-params$asymp_contact_trace_prob)*(1-ctrace)
    asymp_D_NT_ctrace = 
      (params$LoMSM_prob)*(params$LoMSM_prob.det)*
      (params$LoD_RR*params$LoD_dur)*(1-params$asymp_contact_trace_prob)*ctrace
    asymp_U = 
      (params$LoMSM_prob)*(1-params$LoMSM_prob.det)*(params$LoU_RR*params$LoU_dur)
  
  # return values
  trans = c(
    # not traced
    psymp_U, psymp_D_NT_noctrace, psymp_D_NT_ctrace, 
    symp_U, symp_D_NT_noctrace, symp_D_NT_ctrace,
    asymp_U, asymp_D_NT_noctrace,asymp_D_NT_ctrace,
    
    # traced first gen
    psymp_D_T_1, symp_D_T_1, asymp_D_T_1,
    
    # traced second gen
    psymp_D_T_2,  symp_D_T_2, asymp_D_T_2)
    
}

#### make_plots: Make plots ####
make_plots = function(R_plot, xaxis = "test", R0 = 2, Rt = 1) {
  
  # CLEANING
  # modify axis label
  xaxis = paste(" \n", xaxis, collapse = "")
  
  # set theme
  t = theme(panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 13),
            axis.text = element_text(size = 11)) 
  
  # estimate convergence  
  R_plot = R_plot %>% dplyr::group_by(var) %>% 
    dplyr::mutate(
           # R as a ratio
           maxR = R[1], ratio = R/maxR, 
           # percent reduction
           Rt_new = Rt*ratio, 
           
           # containment margin
           
           c_margin = (R0 - 1/(ratio*Rt))/(R0 - Rt),
           c_margin = ifelse(c_margin < 0, 0, c_margin),
           
           # percent from each group
           Symptomatic = symp*ratio*Rt,
           Asymptomatic = asymp*ratio*Rt,
           Presymptomatic = presymp*ratio*Rt, 
           
           # set as factor for cleaning
           Scenario2 = factor(Scenario, levels = 
                                c("No contact \ntracing", "Testing scale-up", 
           "Contact tracing\n(Test symptomatic)", "Contact tracing\n(Test all)")))

  # TOP ROW
  # process data
  R1 = R_plot %>%
    gather(chk, value, Rt_new, c_margin) %>% 
    mutate(var2 = ifelse(chk == "Rt_new", "R(t) with contact tracing", 
                         "Fraction of current physical distancing needed for \n R(t)<1 with contact tracing"),
           var3 = factor(var2, levels = c("R(t) with contact tracing",
                                          "Fraction of current physical distancing needed for \n R(t)<1 with contact tracing")),
           txt =  paste(Scenario2, ": \n x=", round(var,2), "\n y=",  
                        round(value,2), sep = "")) %>% 
    filter(Scenario != "No contact \ntracing") 
  
  # make plot
  ymax = max(1, R1$value)
  a = ggplot(R1, aes(x = var, y = value)) +
    geom_line(lwd = 1, aes(group = Scenario2, col = Scenario2, text = txt)) + 
    theme_minimal(base_size = 20) +
    geom_point(data = R1 %>% filter(point=="point"), 
               aes(x = var, y = value, group = Scenario2, text = txt), 
               size = 2, col = "black") + 
    facet_wrap(.~var3, ncol = 2, scales = "free_y") + t +
    scale_color_brewer(name = "", palette = "Set1") + 
    labs(x = xaxis, y = "", title = "") + ylim(0, ymax)
  
  # BOTTOM LEFT
  # process data
  R2 = R_plot %>% filter(!Scenario%in%c("No contact \ntracing", "Testing scale-up")) %>% 
    mutate(temp = "Fraction of confirmed cases who \nare known contacts",
           txt = paste(Scenario2, ": \n x=", round(var, 2), "\n y=",  
                       round(det_frac, 2), sep = ""))
  
  # make plot
  b = ggplot(R2, 
             aes(x = var, y = det_frac, group = Scenario2, col = Scenario2,
                                                           text = txt)) + 
    geom_line(lwd = 1) + 
    geom_point(data = R2 %>% filter(point=="point"), 
               aes(x = var, y = det_frac, group = Scenario2, text = txt), 
               size = 2, col = "black") + 
    theme_minimal(base_size = 20) + 
    scale_color_brewer(name = "", palette = "Set1") + 
    labs(x = xaxis, y = "", title = "") + t + facet_grid(.~temp) + 
    ylim(0,100) + theme(legend.position='none') + ylim(0,1)
  
  # BOTTOM RIGHT
  # process data
  R3 = R_plot %>% filter(point=="point") %>%
    gather(chk, value, Asymptomatic, Presymptomatic, Symptomatic) %>%
    mutate(txt = paste(chk, ": ", round(value, 2), sep = ""),
           temp = "R(t) by symptom status")
  
  # make plot
  c = ggplot(R3,
             aes(x = Scenario2, y = value, fill = chk, text = txt)) + 
    geom_bar(stat = "identity") + 
    geom_text(position = position_stack(vjust = .5), aes(label = round(value, 2))) + 
    geom_text(aes(y = ratio*Rt*1.1, label = round(ratio*Rt,2))) + 
    theme_minimal(base_size = 20) + 
    scale_fill_brewer(name = "", palette = "Set1") + 
    labs(x = xaxis, y = "", title = "") + t + facet_grid(.~temp)
  
  # return output
  return(list(ggplotly(a, tooltip = c("text")) %>%
                layout(margin = list(b = 50, t = 80)) %>% config(displayModeBar = F),
              ggplotly(b, tooltip = c("text")) %>%
                layout(margin = list(b = 50, t = 80)) %>% config(displayModeBar = F),
              ggplotly(c, tooltip = c("text")) %>%
                layout(margin = list(b = 50, t = 80)) %>% config(displayModeBar = F)
              ))
}


### CONTACT TRACING MODEL### 

#### GET_R: UMBRELLA FUNCTION ####
# 1) set variable over which to vary outcomes
# 2) call make_params to set up transition matrix
# 3) call calc_R to estimate outcomes
# 4) return data frame of R, det_frac, and R by contact group type

get_R = function(SAR, # secondary attack rate
                 HiMSM_contacts, # avg. (daily) contacts, High Contact group 
                 LoMSM_contacts, # avg. (daily) contacts, Low Contact group 
                 duration, # duration of infectiousness
                 HiMSM_prob.det, # detection probability, High Contact
                 LoMSM_prob.det, # detection probability, Low Contact
                 LoMSM_prob, # fraction of infections that are Low Contact. Moot?
                 contact_trace_prob, # probability of contact tracing 
                 comparator, # character string, assigns a scenario
                 baseline_HiMSM_prob.det, # prob detection w/o tracing (Hi)
                 baseline_LoMSM_prob.det, # prob detection w/o tracing (Lo)
                 test_uptake, # probability tested (if traced? without traced?)
                 adh, # adherence / % redux in transmission
                 adh2, # adherence / % redux in transmission
                 vax, # fraction of MSM population with vax
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
    z = make_params(SAR,
                    HiMSM_contacts, 
                    LoMSM_contacts,  
                    duration, 
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
                    vax,
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
  # Base Case: High Contact testing, contact tracing
  # 
make_params = function(SAR, # secondary attack rate
                       HiMSM_contacts, # avg. contacts, High Contact group 
                       LoMSM_contacts, # avg. contacts, Low Contact group 
                       duration, # duration of infectiousness
                       HiMSM_prob.det, # detection probability, High Contact
                       LoMSM_prob.det, # detection probability, Low Contact
                       LoMSM_prob, # fraction of infections that are Low Contact 
                       contact_trace_prob, # probability of contact tracing 
                       comparator, # character string, assigns a scenario
                       baseline_HiMSM_prob.det, # prob detection w/o tracing (sym)
                       baseline_LoMSM_prob.det,  # prob detection w/o tracing (asym)
                       test_uptake, # probability tested if traced 
                       adh, # adherence to isolation (% redux in contacts)
                       adh2, # adherence to isolation (% redux in contacts)
                       vax, # fraction of MSM pop. vaccinated
                       rel_trans # Redux in number of secondary infections, 
                                 # detected vs. undetected
                       ){
  # BASE CASE 
  # Here, 'U' refers to those who are not detected and 'D' refers to those 
  # who are detected. In most cases, new param name is for consistency. 
  params = data.frame(HiMSM_U_RR = HiMSM_contacts*SAR, 
                      HiMSM_D_RR = HiMSM_contacts*SAR*rel_trans, 
                      LoMSM_U_RR = LoMSM_contacts*SAR, 
                      LoMSM_D_RR = LoMSM_contacts*SAR*rel_trans, 
                      duration, 
                      HiMSM_prob.det, 
                      LoMSM_prob.det, 
                      LoMSM_prob,
                      HiMSM_contact_trace_prob = contact_trace_prob, 
                      LoMSM_contact_trace_prob = contact_trace_prob,
                      vax)
  
  # COUNTERFACTUAL: base case with no contact tracing
  params_cf = params %>% 
    mutate(HiMSM_contact_trace_prob = 0, 
           LoMSM_contact_trace_prob = 0,
           HiMSM_prob.det = ifelse(comparator=="Contact tracing only", 
                               HiMSM_prob.det, 
                               baseline_HiMSM_prob.det), 
           LoMSM_prob.det = ifelse(comparator=="Contact tracing only", 
                               LoMSM_prob.det, 
                               baseline_LoMSM_prob.det))

### !!!! ADD NEW SCENARIOS FOR MONKEY POX !!! #
# Name changes are *encouraged*, please update names in calc_R fx, etc. 
  #line 159, line 167, line 179, line 204
############################################################################### 
  # CONTACT TRACING GEN 1
  # 1) Increase testing of High Contact contacts (test_uptake)
  # 2) Decrease transmission according to adherence (adh)
  params_ctrace_1 = params %>% 
    mutate(HiMSM_prob.det = test_uptake,
           HiMSM_U_RR      = HiMSM_U_RR*(1-adh),
           HiMSM_D_RR      = HiMSM_D_RR*(1-adh),
           LoMSM_U_RR      = LoMSM_U_RR*(1-adh), 
           LoMSM_D_RR      = LoMSM_D_RR*(1-adh))
  
  # CONTACT TRACING GEN 2
  # 1) Increase testing of High Contact contacts (test_uptake)
  # 2) Decrease transmission according to adherence (adh2)
  params_ctrace_2plus = params %>% 
    mutate(HiMSM_prob.det = test_uptake,
           HiMSM_U_RR      = HiMSM_U_RR*(1-adh2),
           HiMSM_D_RR      = HiMSM_D_RR*(1-adh2),
           LoMSM_U_RR      = LoMSM_U_RR*(1-adh2),
           LoMSM_D_RR      = LoMSM_D_RR*(1-adh2))
  
  # !!!! ADD NEW SCENARIOS FOR MONKEY POX !!! #
  # CONTACT TRACTING AND PASSIVE TESTING GEN 1
  # 1) Increase Low Contact testing
  # 2) Allow  transmission reduction for detected presymptomatics 
  params_test_ctrace_1 = params_ctrace_1 %>% 
    mutate(LoMSM_prob.det = test_uptake, 
           P_RR = P_RR*.5) # drawing from uniform dist, reduce by half???
  
  # CONTACT TRACTING AND PASSIVE TESTING GEN 2
  # Same as gen 1, but adapting params_ctrace_2plus
  params_test_ctrace_2_plus = params_ctrace_2plus %>% 
    mutate(LoMSM_prob.det = test_uptake, 
           P_RR = P_RR*.5)
  
  # RETURN OUTPUT
  return(list(params_cf, params, params_ctrace_1, params_ctrace_2plus, 
              params_test_ctrace_1, params_test_ctrace_2_plus))
}

#### calc_R: run over each scenario #### NAME CHANGES
# 1) No contact tracing (params_cf)
# 2) Test High Contact (params and params_ctrace_XX)
# 3) Test all (params and params_test_ctrace_XX)
calc_R = function(params_cf,# NO CONTACT TRACING (counterfactual)
                 
                  params, # BASE CASE (first generation for CT scenarios)
           
                  params_ctrace_1, #CT and test High Contact
                  params_ctrace_2plus, #CT and test High Contact
                  
                  params_test_ctrace_1, #CT and test all
                  params_test_ctrace_2_plus #CT and test all 
                  ) {
  
  # run different methods NAME CHANGES HEREA
  out = bind_rows(dom_eigen(params_cf, params_cf, params_cf) %>% 
                    mutate(Scenario = "No contact \ntracing"),
                  dom_eigen(params, params_ctrace_1, params_ctrace_2plus) %>% 
                    mutate(Scenario = "Contact tracing\n(Test High Contact)"),
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
# because eigenvalues scale linearly when multiplied by a constant
# We adjust to user-input R(t) later on in make_plots

dom_eigen = function(params, params_ctrace_1, params_ctrace_2plus){
 # !!! Have Alyssa double check updates here
  
# call transition probs -- NAME CHANGES HERE
  x = get_trans_probs(params, first_gen = T, ctrace = F)
  y = get_trans_probs(params_ctrace_1, first_gen = F, ctrace = T) 
  z = get_trans_probs(params_ctrace_2plus,  first_gen = F, ctrace = T)
  
  # put in matrix
  mat = matrix(unlist(c(rep(x, 6), 
                        rep(y, 2),
                        rep(z, 2))), ncol = 10, byrow = T)
  
  # estimate eigenvectors/values
  vec = Re(eigen(t(mat))$vectors[,1])
  
  # return output
  return(data.frame(
    # R(t)
    R = max(Re(eigen(mat)$values)),
    # detection fraction
    det_frac = sum(vec[3], vec[6], vec[9:10])/
                sum(vec[2:3], vec[5:6], vec[7:10]),
    # transmission by contact group status
    HiMSM = sum(vec[1:3], vec[7], vec[9])/sum(vec),
    LoMSM = sum(vec[4:6], vec[8], vec[10])/sum(vec)))
}

#### get_trans_probs: Pull together transition probabilities ####
# by symptom status, detection status, whether originated from contact tracing 
# + gen, whether traced

# This is our next generation matrix -- a way of representing our compartments
# as a matrix. This maps onto an estimate of Rt: how many infections per 
# infected in the next generation? 

# Here, we calculate the transition probabilities -- e.g. what should go in 
# the next generation matrix. It should be 6x6 (double check)
get_trans_probs = function(params, first_gen = F, ctrace = F) {

    # High Contact group
    HiMSM_D_T_1  = # first gen detected through contact tracing
      params$HiMSM_contact_trace_prob*first_gen* # contact tracing prob
      (params$HiMSM_prob.det)* # detection rate
      (1-params$LoMSM_prob)* #prob in the HiMSM pop. 
      (params$HiMSM_D_RR*params$duration*(1-params$vax)) # Rt (RR * duration)

    HiMSM_D_T_2  = # second gen detected through contact tracing
      params$HiMSM_contact_trace_prob*(1-first_gen)*
      (1-params$LoMSM_prob)* 
      (params$HiMSM_prob.det)*
      (params$HiMSM_D_RR*params$duration*(1-params$vax)) 

    HiMSM_D_NT_noctrace = # gen 1 not traced, gen 2 not traced
      (1-params$HiMSM_contact_trace_prob)* # prob. gen 1 not contact traced
      (1-params$LoMSM_prob)*
      (params$HiMSM_prob.det)*
      (1-ctrace)* # prob. gen 2 not contact traced
      (params$HiMSM_D_RR*params$duration*(1-params$vax))* 

    HiMSM_D_NT_ctrace = # gen 1 not traced, gen 2 traced
      (1-params$HiMSM_contact_trace_prob)*
      (1-params$LoMSM_prob)*
      (params$HiMSM_prob.det)*
      (ctrace)* # prob. gen 2 contact traced
      (params$HiMSM_D_RR*params$duration*(1-params$vax))
    
    HiMSM_U = # undetected -- second generation cannot be traced
      (1-params$HiMSM_prob.det)*
      (1-params$LoMSM_prob)*
      (params$HiMSM_U_RR*params$duration*(1-params$vax))
    
    # Low Contact group
    LoMSM_D_T_1  = 
      params$LoMSM_contact_trace_prob*
      (params$LoMSM_prob)*
      (params$LoMSM_prob.det)*
      (params$LoMSM_D_RR*params$duration*(1-params$vax))*
      first_gen
    
    LoMSM_D_T_2  = 
      params$LoMSM_contact_trace_prob*
      (params$LoMSM_prob)*
      (params$LoMSM_prob.det)*
      (params$LoMSM_D_RR*params$duration*(1-params$vax))*
      (1-first_gen)
    
    LoMSM_D_NT_noctrace = 
      (1-params$LoMSM_contact_trace_prob)*
      (params$LoMSM_prob)*
      (params$LoMSM_prob.det)*
      (1-ctrace)*
      (params$LoMSM_D_RR*params$duration*(1-params$vax))

    LoMSM_D_NT_ctrace = 
      (1-params$LoMSM_contact_trace_prob)*
      (params$LoMSM_prob)*
      (params$LoMSM_prob.det)*
      ctrace*
      (params$LoMSM_D_RR*params$duration*(1-params$vax))

    LoMSM_U = 
      (1-params$LoMSM_prob.det)*
      (params$LoMSM_prob)*
      (params$LoMSM_U_RR*params$duration*(1-params$vax))
    
# return values
  trans = c(
    # not traced
    HiMSM_U, HiMSM_D_NT_noctrace, HiMSM_D_NT_ctrace,
    LoMSM_U, LoMSM_D_NT_noctrace, LoMSM_D_NT_ctrace,
    # traced first gen
    HiMSM_D_T_1, LoMSM_D_T_1,
    # traced second gen
    HiMSM_D_T_2, LoMSM_D_T_2)
}

#### make_plots: Make plots ####
#make_plots = function(R_plot, xaxis = "test", R0 = 2, Rt = 1) {
#  
#  # CLEANING
#  # modify axis label
#  xaxis = paste(" \n", xaxis, collapse = "")
#  
#  # set theme
#  t = theme(panel.background = element_blank(),
#            panel.grid.major = element_blank(),
#            panel.grid.minor = element_blank(),
#            axis.title = element_text(size = 13),
#            axis.text = element_text(size = 11)) 
#  
#  # estimate convergence  
#  R_plot = R_plot %>% dplyr::group_by(var) %>% 
#    dplyr::mutate(
#           # R as a ratio
#           maxR = R[1], ratio = R/maxR, 
#           # percent reduction
#           Rt_new = Rt*ratio, 
#           
#           # containment margin
#           
#           c_margin = (R0 - 1/(ratio*Rt))/(R0 - Rt),
#           c_margin = ifelse(c_margin < 0, 0, c_margin),
#           
#           # percent from each group
#           HiMSM_pct = symp*ratio*Rt,
#           LoMSM_pct = LoMSM*ratio*Rt,
#           Presymptomatic = presymp*ratio*Rt, 
#           
#           # set as factor for cleaning
#           Scenario2 = factor(Scenario, levels = 
#                                c("No contact \ntracing", "Testing scale-up", 
#           "Contact tracing\n(Test High Contact)", "Contact tracing\n(Test all)")))
#
#  # TOP ROW
#  # process data
#  R1 = R_plot %>%
#    gather(chk, value, Rt_new, c_margin) %>% 
#    mutate(var2 = ifelse(chk == "Rt_new", "R(t) with contact tracing", 
#                         "Fraction of current physical distancing needed for \n R(t)<1 with contact tracing"),
#           var3 = factor(var2, levels = c("R(t) with contact tracing",
#                                          "Fraction of current physical distancing needed for \n R(t)<1 with contact tracing")),
#           txt =  paste(Scenario2, ": \n x=", round(var,2), "\n y=",  
#                        round(value,2), sep = "")) %>% 
#    filter(Scenario != "No contact \ntracing") 
#  
#  # make plot
#  ymax = max(1, R1$value)
#  a = ggplot(R1, aes(x = var, y = value)) +
#    geom_line(lwd = 1, aes(group = Scenario2, col = Scenario2, text = txt)) + 
#    theme_minimal(base_size = 20) +
#    geom_point(data = R1 %>% filter(point=="point"), 
#               aes(x = var, y = value, group = Scenario2, text = txt), 
#               size = 2, col = "black") + 
#    facet_wrap(.~var3, ncol = 2, scales = "free_y") + t +
#    scale_color_brewer(name = "", palette = "Set1") + 
#    labs(x = xaxis, y = "", title = "") + ylim(0, ymax)
#  
#  # BOTTOM LEFT
#  # process data
#  R2 = R_plot %>% filter(!Scenario%in%c("No contact \ntracing", "Testing scale-up")) %>% 
#    mutate(temp = "Fraction of confirmed cases who \nare known contacts",
#           txt = paste(Scenario2, ": \n x=", round(var, 2), "\n y=",  
#                       round(det_frac, 2), sep = ""))
#  
#  # make plot
#  b = ggplot(R2, 
#             aes(x = var, y = det_frac, group = Scenario2, col = Scenario2,
#                                                           text = txt)) + 
#    geom_line(lwd = 1) + 
#    geom_point(data = R2 %>% filter(point=="point"), 
#               aes(x = var, y = det_frac, group = Scenario2, text = txt), 
#               size = 2, col = "black") + 
#    theme_minimal(base_size = 20) + 
#    scale_color_brewer(name = "", palette = "Set1") + 
#    labs(x = xaxis, y = "", title = "") + t + facet_grid(.~temp) + 
#    ylim(0,100) + theme(legend.position='none') + ylim(0,1)
#  
#  # BOTTOM RIGHT
#  # process data
#  R3 = R_plot %>% filter(point=="point") %>%
#    gather(chk, value, LoMSM_pct, Presymptomatic, HiMSM_pct) %>%
#    mutate(txt = paste(chk, ": ", round(value, 2), sep = ""),
#           temp = "R(t) by symptom status")
#  
#  # make plot
#  c = ggplot(R3,
#             aes(x = Scenario2, y = value, fill = chk, text = txt)) + 
#    geom_bar(stat = "identity") + 
#    geom_text(position = position_stack(vjust = .5), aes(label = round(value, 2))) + 
#    geom_text(aes(y = ratio*Rt*1.1, label = round(ratio*Rt,2))) + 
#    theme_minimal(base_size = 20) + 
#    scale_fill_brewer(name = "", palette = "Set1") + 
#    labs(x = xaxis, y = "", title = "") + t + facet_grid(.~temp)
#  
#  # return output
#  return(list(ggplotly(a, tooltip = c("text")) %>%
#                layout(margin = list(b = 50, t = 80)) %>% config(displayModeBar = F),
#              ggplotly(b, tooltip = c("text")) %>%
#                layout(margin = list(b = 50, t = 80)) %>% config(displayModeBar = F),
#              ggplotly(c, tooltip = c("text")) %>%
#                layout(margin = list(b = 50, t = 80)) %>% config(displayModeBar = F)
#              ))
#}

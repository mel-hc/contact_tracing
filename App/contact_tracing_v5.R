### CONTACT TRACING MODEL### 

#### GET_R: UMBRELLA FUNCTION ####
# 1) set variable over which to vary outcomes
# 2) call make_params to set up transition matrix
# 3) call calc_R to estimate outcomes
# 4) return data frame of R, det_frac, and R by contact group type

get_R = function(SAR, # secondary attack rate
                 HiMSM_contacts, # avg. (daily) contacts, High Contact group p 
                 duration, # duration of infectiousness
                 HiMSM_prob.det, # detection probability, High Contact
                 contact_trace_prob, # probability of contact tracing 
                 comparator, # character string, assigns a scenario
                 test_uptake, # fraction of contacts tested
                 adh, # adherence / % redux in transmission
                 adh2, # adherence / % redux in transmission
                 vax,     # compared to undetected)
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
                    duration, 
                    HiMSM_prob.det,
                    contact_trace_prob,
                    comparator, 
                    test_uptake, 
                    adh, 
                    adh2, 
                    vax)

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

make_params = function(SAR, # secondary attack rate
                       HiMSM_contacts, # avg. contacts, High Contact group 
                       duration, # duration of infectiousness
                       HiMSM_prob.det, # detection probability, High Contact
                       contact_trace_prob, # probability of contact tracing 
                       comparator, # character string, assigns a scenario
                       test_uptake, # probability tested if traced 
                       adh, # adherence to isolation (% redux in contacts)
                       adh2, # adherence to isolation (% redux in contacts)
                       vax # fraction of MSM pop. vaccinated
                       ){
  # BASE CASE 
  # Here, 'U' refers to those who are not detected and 'D' refers to those 
  # who are detected. In most cases, new param name is for consistency. 
  params = data.frame(HiMSM_U_RR = HiMSM_contacts*SAR, 
                      HiMSM_D_RR = HiMSM_contacts*SAR,
                      duration, 
                      HiMSM_prob.det, 
                      contact_trace_prob,
                      vax) 
  # COUNTERFACTUAL: base case with no contact tracing
  params_cf = params %>% 
    mutate(contact_trace_prob = 0)

  # CONTACT TRACING GEN 1
  # 1) Increase detection in High Contact group (test_uptake)
  # 2) Decrease transmission according to adherence (adh)
  params_ctrace_1 = params %>% 
    mutate(HiMSM_prob.det = test_uptake,
           HiMSM_U_RR      = HiMSM_U_RR*(1-(adh*0.75)),
           HiMSM_D_RR      = HiMSM_D_RR*(1-adh))
  
  # CONTACT TRACING GEN 2
  # 1) Increase detection in High Contact group (test_uptake)
  # 2) Decrease transmission according to adherence (adh2)
  params_ctrace_2plus = params %>% 
    mutate(HiMSM_prob.det = test_uptake,
           HiMSM_U_RR      = HiMSM_U_RR*(1-(adh*0.75)),
           HiMSM_D_RR      = HiMSM_D_RR*(1-adh))
  
  # RETURN OUTPUT
  return(list(params_cf, params, params_ctrace_1, params_ctrace_2plus))
}

#### calc_R: run over each scenario #### NAME CHANGES
# 1) Testing, no contact tracing (params_cf)
# 2) Contact tracing and increased testing (params and params_ctrace_XX)

calc_R = function(params_cf,# NO CONTACT TRACING (counterfactual)
                  params, # BASE CASE (first generation for CT scenarios)
                  params_ctrace_1, #CT and test 
                  params_ctrace_2plus, #CT and test 
                  ) {
  
  # run different methods 
  out = bind_rows(#dom_eigen(params_cf, params_cf, params_cf) %>% 
                  #  mutate(Scenario = "No contact tracing"),
                  dom_eigen(params, params_ctrace_1, params_ctrace_2plus) %>% 
                    mutate(Scenario = "Contact tracing and increased testing"))
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
  mat = matrix(unlist(c(rep(x, 3), 
                        rep(y, 1), 
                        rep(z, 1))), 
               ncol = 5, byrow = T)
  
  # estimate eigenvectors/values
  vec = Re(eigen(t(mat))$vectors[,1])
  
  # return output
  return(data.frame(
    # R(t)
    R = max(Re(eigen(mat)$values)))
)
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
      params$contact_trace_prob*first_gen* # contact tracing prob
      (params$HiMSM_prob.det)* # detection rate
      (params$HiMSM_D_RR*params$duration*(1-params$vax)) # Rt (RR * duration)

    HiMSM_D_T_2  = # second gen detected through contact tracing
      params$contact_trace_prob*(1-first_gen)*
      (params$HiMSM_prob.det)*
      (params$HiMSM_D_RR*params$duration*(1-params$vax)) 

    HiMSM_D_NT_noctrace = # gen 1 not traced, gen 2 not traced
      (1-params$contact_trace_prob)* # prob. gen 1 not contact traced
      (params$HiMSM_prob.det)*
      (1-ctrace)* # prob. gen 2 not contact traced
      (params$HiMSM_D_RR*params$duration*(1-params$vax))

    HiMSM_D_NT_ctrace = # gen 1 not traced, gen 2 traced
      (1-params$contact_trace_prob)*
      (params$HiMSM_prob.det)*
      (ctrace)* # prob. gen 2 contact traced
      (params$HiMSM_D_RR*params$duration*(1-params$vax))
    
    HiMSM_U = # undetected -- second generation cannot be traced
      (1-params$HiMSM_prob.det)*
      (params$HiMSM_U_RR*params$duration*(1-params$vax))
  
  
# return values
  trans = c(
    HiMSM_U, HiMSM_D_NT_noctrace, HiMSM_D_NT_ctrace, # not traced
    HiMSM_D_T_1, # traced first gen
    HiMSM_D_T_2) # traced second gen
}

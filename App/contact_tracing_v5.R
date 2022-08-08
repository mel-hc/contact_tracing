### CONTACT TRACING MODEL### 

#### GET_R: UMBRELLA FUNCTION ####
# 1) set variable over which to vary outcomes
# 2) call make_params to set up transition matrix
# 3) call calc_R to estimate outcomes
# 4) return data frame of R, det_frac, and R by contact group type
get_R = function(SAR, # secondary attack rate
                 HiMSM_contacts, # avg. contacts, High Contact group 
                 duration, # duration of infectiousness
                 HiMSM_prob.det_base, # baseline detection
                 HiMSM_prob.det_comm, # community detection 
                 HiMSM_prob.det_traced, # traced detection
                 contact_trace_prob, # probability of contact tracing 
                 adh, # adherence to isolation (% redux in contacts)
                 vax, # fraction of MSM pop. vaccinated
                 vax.eff, # the efficacy of the vaccine against infection
                 rel_trans, # relative transmissibility of detected cases
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
                    HiMSM_prob.det_base, 
                    HiMSM_prob.det_comm, 
                    HiMSM_prob.det_traced,
                    contact_trace_prob,
                    adh, 
                    adh2, 
                    vax, 
                    vax.eff,
                    rel_trans)

    a = bind_rows(a, 
                  calc_R(z[[1]], z[[2]], z[[3]], z[[4]]) %>% 
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
                       HiMSM_prob.det_base, # baseline detection
                       HiMSM_prob.det_comm, # community detection 
                       HiMSM_prob.det_traced, # traced detection
                       contact_trace_prob, # probability of contact tracing 
                       adh, # adherence to isolation (% redux in contacts)
                       vax, # fraction of MSM pop. vaccinated
                       vax.eff, # the efficacy of the vaccine against infection
                       rel_trans # relative transmissibility of detected cases
                       ){
  # BASE CASE 
  # Here, 'U' refers to those who are not detected and 'D' refers to those 
  # who are detected. 
  params = data.frame(HiMSM_U_RR = HiMSM_contacts*SAR, 
                    HiMSM_D_RR = HiMSM_contacts*SAR*(1-rel_trans),
                    HiMSM_prob.det = HiMSM_prob.det_comm,
                    duration, 
                    contact_trace_prob,
                    vax, 
                    vax.eff) 
  
  # COUNTERFACTUAL: base case with no contact tracing
  params_cf = params %>% 
    mutate(contact_trace_prob = 0, 
           HiMSM_prob.det = HiMSM_prob.det_base)

  # CONTACT TRACING GEN 1
  # 1) Increase detection in High Contact group (HiMSM_prob.det)
  # 2) Decrease transmission according to adherence (adh)
  params_ctrace_1 = params %>% 
    mutate(HiMSM_prob.det = HiMSM_prob.det_traced,
           HiMSM_U_RR      = HiMSM_U_RR*(1-adh),
           HiMSM_D_RR      = HiMSM_D_RR*(1-adh))
  
  # CONTACT TRACING GEN 2
  # 1) Increase detection in High Contact group (incr_HiMSM_prob.det)
  # 2) Decrease transmission according to adherence (adh2)
  params_ctrace_2plus = params_ctrace_1 %>%
    mutate(HiMSM_prob.det  = HiMSM_prob.det_traced)
  
  # RETURN OUTPUT
  return(list(params_cf, params, params_ctrace_1, params_ctrace_2plus))
}

#### calc_R: run over each scenario #### NAME CHANGES
# 1) Testing, no contact tracing (params_cf)
# 2) Contact tracing and increased testing (params and params_ctrace_XX)
calc_R = function(params_cf,# NO CONTACT TRACING (counterfactual)
                  params, # BASE CASE (first generation for CT scenarios)
                  params_ctrace_1, #CT and test 
                  params_ctrace_2plus #CT and test 
                  ) {
  
  # run different methods 
  out = bind_rows(dom_eigen(params_cf, params_cf, params_cf) %>% 
                    mutate(Scenario = "No contact tracing"),
                  dom_eigen(params, params_ctrace_1, params_ctrace_2plus) %>% 
                    mutate(Scenario = "Contact tracing"))
  return(out)  
}

#### dom_eigen: Make transition matrix and calculate R + steady state ####
# 1) Call transition probs for each case
# 2) Make transition matrix
# 3) Take eigenvalues
# 4) Take eigenvectors of transpose
# 5) Return data from of relevant quantities

dom_eigen = function(params, params_ctrace_1, params_ctrace_2plus){
# call transition probs 
  x = get_trans_probs(params, first_gen = T, ctrace = F)
  y = get_trans_probs(params_ctrace_1, first_gen = F, ctrace = T) 
  z = get_trans_probs(params_ctrace_2plus, first_gen = F, ctrace = T)
  
  # put in matrix
  mat = matrix(unlist(c(rep(x, 3), rep(y, 1), rep(z, 1))), 
               ncol = 5, byrow = T)
  # return output
  return(data.frame(
    # R(t)
    R = max(Re(eigen(mat)$values)))
)
}

#### get_trans_probs: Pull together transition probabilities ####
# by symptom status, detection status, whether originated from contact tracing 
# + gen, whether traced

# Probabilities of transitioning INTO each state

get_trans_probs = function(params, first_gen = F, ctrace = F) {

    HiMSM_D_T_1  = 
      first_gen* 
      params$contact_trace_prob*# contact tracing prob
      params$HiMSM_prob.det* # detection rate
      (params$HiMSM_D_RR*params$duration*(1-params$vax*params$vax.eff)) # Rt (RR * duration)

    HiMSM_D_T_2  = 
      (1-first_gen)*
      params$contact_trace_prob*
      params$HiMSM_prob.det*
      (params$HiMSM_D_RR*params$duration*(1-params$vax*params$vax.eff)) 

    HiMSM_D_NT_noctrace = 
      (1-params$contact_trace_prob)* # prob. gen 1 not contact traced
      params$HiMSM_prob.det*
      (1-ctrace)* #  contact trace indicator
      (params$HiMSM_D_RR*params$duration*(1-params$vax*params$vax.eff))

    HiMSM_D_NT_ctrace = 
      (1-params$contact_trace_prob)*
      (params$HiMSM_prob.det)*
      (ctrace)* #  contact trace indicator
      (params$HiMSM_D_RR*params$duration*(1-params$vax*params$vax.eff))
    
    HiMSM_U = 
      (1-params$HiMSM_prob.det)*
      (params$HiMSM_U_RR*params$duration*(1-(params$vax*params$vax.eff)))

# return values
  trans = c(
    HiMSM_U, HiMSM_D_NT_noctrace, HiMSM_D_NT_ctrace, # not traced
    HiMSM_D_T_1, # traced first gen
    HiMSM_D_T_2) # traced second gen
}

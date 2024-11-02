# This R script fit several univariable generalized linear models
# and outputs the results -Estimate, 95% CI, p-
# in a merged table for any outcome and risk factor

# The user needs to provide:
# 1. The dataset
# 2. A vector of the outcomes
# 3. A vector of the risk factors
# 4. The family (gaussian or binomial)
# 5. Whether they want a Firth penalty

library(dplyr)
library(rstatix)

# ---------------------------------------------------------
#           Logistic Regression
# ---------------------------------------------------------


fit_logistic_univ = function(data, outcome, rf, Firth){
  
  
  if (Firth){
    
    myfit = brglm::brglm(unlist(data[,outcome]) ~ rf, data = data, family = 'binomial') %>%
      
      summary
    
  } else {
    
    myfit = glm(unlist(data[,outcome]) ~ rf, data = data, family = 'binomial') %>%
      
      summary
    
  }
  
  
  model_results = myfit$coefficients %>%
    
    as.data.frame() %>%
    
    mutate(
      
      RF = row.names(.),
      
      low = round(exp(Estimate - 1.96*`Std. Error`), 2),
      
      upper = round(exp(Estimate + 1.96*`Std. Error`), 2)
      
    ) %>%
    
    filter(RF != '(Intercept)') %>%
    
    transmute(
      
      Category = gsub('rf', '', RF),
      
      OR = round(exp(Estimate), 2),
      
      `95% CI` = paste0('(', low, ', ', upper, ')'),
      
      p = round(`Pr(>|z|)`, 3)
      
      
    )
  
  rownames(model_results) = NULL
  
  return(model_results)
  
  
}

# ---------------------------------------------------------
#           Linear Regression
# ---------------------------------------------------------


fit_linear_univ = function(data, outcome, rf){
  
  
  myfit = glm(unlist(data[,outcome]) ~ rf, data = data) %>%
    
    summary
  
  
  model_results = myfit$coefficients %>%
    
    as.data.frame() %>%
    
    mutate(
      
      RF = row.names(.),
      
      low = round(Estimate - 1.96*`Std. Error`, 2),
      
      upper = round(Estimate + 1.96*`Std. Error`, 2)
      
    ) %>%
    
    filter(RF != '(Intercept)') %>%
    
    transmute(
      
      Category = gsub('rf', '', RF),
      
      beta = round(Estimate, 2),
      
      `95% CI` = paste0('(', low, ', ', upper, ')'),
      
      p = round(`Pr(>|t|)`, 3)
      
      
    )
  
  rownames(model_results) = NULL
  
  return(model_results)
  
}

# ---------------------------------------------------------
#           GLM Model Fit
# ---------------------------------------------------------

glm_univ = function(data, outcome, riskfactor, family = 'gaussian', Firth = F){
  
  mydata = data %>%
    
    select(any_of(outcome), any_of(riskfactor)) %>%
    
    gather(., key = 'Outcome', value = 'Outcome_Value', any_of(outcome)) 
  
  
  mydata1 = split(mydata, f = mydata$Outcome)
  
  tmp = lapply(mydata1, function(x){
    
    y = x %>%
      
      select(Outcome_Value, everything()) %>%
      
      select(-Outcome)
    
    y1 = y %>% select(-Outcome_Value)
    
    result_within_outcome = lapply(y1, function(z){
      
      
      # fit_logistic_univ(y, 'Outcome_Value', z, Firth)
      if (family == 'gaussian'){
        
        tryCatch(
          
          fit_linear_univ(y, 'Outcome_Value', z),
          error = function(e){NULL})
        
      } else {
        
        tryCatch(
          
          fit_logistic_univ(y, 'Outcome_Value', z, Firth),
          error = function(e){NULL})
        
      }
      
    })
    
    tryCatch(
      
      data.table::rbindlist(result_within_outcome, idcol = 'Risk_Factor'),
      
      error = function(e){NULL}
      
    )
    
  }) 
  
  results = tryCatch(
    
    data.table::rbindlist(tmp, idcol = 'Outcome'),
    
    error = function(e){NULL}
    
  )
  
  
  return(results)
  
}

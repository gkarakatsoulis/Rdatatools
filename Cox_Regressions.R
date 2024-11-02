# This R script conducts univariable Cox regressions
# for multiple risk factors separately
# and reports the results in a merged table
# It also allows the inclusion of control covariates


# The user needs to provide:
# 1. The dataset
# 2. The time variable
# 3. The event variable
# 4. A vector with the names of the risk factors
# 5. In case of controlling, the vector of controlling covariates


library(survival)
library(survminer)


# ---------------------------------------------------------
#           Unadjusted Univariable Testing
# ---------------------------------------------------------

fit_cox_univ = function(data, time, event, rf){
  
  myfit = coxph(
    
    Surv(time = data[, time], event = data[,event]) ~ rf
    
  )
  
  tmp = summary(myfit)$coefficients %>%
    
    as.data.frame() %>%
    
    mutate(Category = rownames(.))
  
  
  tmp$Category = gsub('rf', '', tmp$Category)
  
  tmp$Overall_p = round(summary(myfit)$sctest[3], 3)
  
  
  tmp$L = round(exp(tmp$coef - 1.96 * tmp$`se(coef)`), 2)
  tmp$U = round(exp(tmp$coef + 1.96 * tmp$`se(coef)`), 2)
  
  tmp$p = round(tmp$`Pr(>|z|)`, 3)
  
  tmp$Estimate = round(exp(tmp$coef), 2)
  
  tmp = tmp %>%
    
    transmute(
      
      Category = Category,
      HR = Estimate,
      `95% CI` = paste0('(', L, ', ', U, ')'),
      p = ifelse(p < 0.001, '<0.001', as.character(p)),
      Overall_p = ifelse(Overall_p < 0.001, '<0.001', as.character(Overall_p))
      
    )
  
  rownames(tmp) = NULL
  
  tmp$Overall_p = ifelse((duplicated(tmp$Category)), NA,
                         tmp$Overall_p)
  
  
  return(tmp)
  
}


f_cox_univ = function(data, time, event, riskfactor){
  
  
  mydata = data %>%
    
    dplyr::select(all_of(time), all_of(event), any_of(riskfactor)) %>%
    
    as.data.frame()
  
  mydata1 = mydata %>%
    
    select(-all_of(time), -all_of(event))
  
  tmp = lapply(mydata1, function(x){
    
    tryCatch(
      
      fit_cox_univ(mydata, time, event, x),
      error = function(e){NULL})
    
  })
  
  return(
    
    tryCatch(
      
      data.table::rbindlist(tmp, idcol = 'Risk_Factor'),
      
      error = function(e){NULL}
      
    )
    
    
  )
  
}


# ---------------------------------------------------------
#           Ajusted Univariable Testing
# ---------------------------------------------------------

fit_cox_adjusted = function(data, outcome, rf, adjusted){
  
  data = data %>%
    
    select(all_of(outcome), all_of(adjusted))
  
  data = cbind(data, rf)
  
  colnames(data)[-1] = paste0(colnames(data)[-1], '_rf_')
  
  myfit = coxph(outcome ~ ., data = data)
  
  tmp = summary(myfit)$coefficients %>%
    
    as.data.frame() %>%
    
    mutate(Category = rownames(.))
  
  tmp$Risk_Factor = gsub('_rf_.*', '', tmp$Category)
  
  tmp$Category = gsub('.*_rf_', '', tmp$Category)
  
  tmp$Overall_p = round(summary(myfit)$sctest[3], 3)
  
  
  tmp$L = round(exp(tmp$coef - 1.96 * tmp$`se(coef)`), 2)
  tmp$U = round(exp(tmp$coef + 1.96 * tmp$`se(coef)`), 2)
  
  tmp$p = round(tmp$`Pr(>|z|)`, 3)
  
  tmp$Estimate = round(exp(tmp$coef), 2)
  
  tmp = tmp %>%
    
    transmute(
      
      Risk_Factor = Risk_Factor,
      Category = Category,
      HR = Estimate,
      `95% CI` = paste0('(', L, ', ', U, ')'),
      p = ifelse(p < 0.001, '<0.001', as.character(p))
      
    )
  
  rownames(tmp) = NULL
  
  return(tmp)
  
}


f_cox_adjusted = function(data, time, event, riskfactor, adjusted){
  
  
  mydata = data %>%
    
    dplyr::select(all_of(time), all_of(event), any_of(riskfactor), any_of(adjusted)) %>%
    
    as.data.frame()
  
  mydata$outcome = Surv(time = mydata[, time], event = mydata[,event])
  
  mydata = mydata %>%
    
    select(-all_of(time), - all_of(event))
  
  mydata1 = mydata %>%
    
    select(-outcome)
  
  tmp = lapply(mydata1, function(x){
    
    
    tryCatch(
      
      fit_cox_adjusted(mydata, 'outcome', x, adjusted),
      
      error = function(e){NULL}
      
    )
    
    
  })
  
  return(
    
    tryCatch(
      
      data.table::rbindlist(tmp, idcol = 'Var') %>%
        
        mutate(
          
          Risk_Factor = case_when(
            
            Risk_Factor == 'rf' ~ Var,
            .default = Risk_Factor
          )
          
        ),
      
      error = function(e){NULL}
      
    )
    
    
  )
  
}



# ---------------------------------------------------------
#           Multivariable Testing
# ---------------------------------------------------------

fit_cox_multiv = function(data, time, event, riskfactor, step, trace){
  
  
  mydata1 = data %>%
    
    dplyr::select(all_of(riskfactor))
  
  myfit = coxph(
    
    Surv(time = data[, time], event = data[,event]) ~ ., data = mydata1
  )
  
  
  if (!is.null(step)){
    
    myfit = myfit %>%
      
      step(trace = trace)
    
  } 
  
  mymatrix = summary(myfit)$coefficients %>%
    
    as.data.frame() %>%
    
    mutate(RF = rownames(.))
  
  
  mymatrix$Category = gsub('.*_CAT_', '', mymatrix$RF)
  
  mymatrix$RF = gsub('_CAT_.*', '', mymatrix$RF)
  
  mymatrix$L = round(exp(mymatrix$coef - 1.96 * mymatrix$`se(coef)`), 2)
  mymatrix$U = round(exp(mymatrix$coef + 1.96 * mymatrix$`se(coef)`), 2)
  
  mymatrix$p = round(mymatrix$`Pr(>|z|)`, 3)
  
  mymatrix$Estimate = round(exp(mymatrix$coef), 2)
  
  mymatrix = mymatrix %>%
    
    transmute(
      RF = RF,
      Category = Category,
      HR = Estimate,
      `95% CI` = paste0('(', L, ', ', U, ')'),
      p = ifelse(p < 0.001, '<0.001', as.character(p))
      
    )
  
  rownames(mymatrix) = NULL
  
  mymatrix$RF[which(duplicated(mymatrix$RF))] = NA
  
  return(mymatrix)
  
  
}

f_cox_multiv = function(data, time, event, riskfactor, step = NULL, trace = F){
  
  
  mydata = data %>%
    
    select(all_of(time), all_of(event), any_of(riskfactor)) %>%
    
    na.omit() %>%
    
    as.data.frame()
  
  # Apply a pattern to distinguish the risk factors from the categories
  
  colnames(mydata)[-c(1,2)] = paste0(colnames(mydata)[-c(1,2)], '_CAT_')
  
  riskfactor = paste0(riskfactor, '_CAT_')
  
    
  return(
    
    tryCatch(
      
      fit_cox_multiv(mydata, time, event, riskfactor, step, trace),
      
      error = function(e){ NULL}
      
    )
    
  )    
      
  
}





# Test for trend
fit_cox_ordinal = function(data, time, event, rf){
  

  myfit = coxph(
    
    Surv(time = data[, time], event = data[,event]) ~ rf
    
  )
  
  # tmp = surv_pvalue(myfit, data = ddd, test.for.trend = T) %>%
  # 
  #   as.data.frame() %>%
  # 
  #   transmute(variable = rf,
  #             Overall_p = round(pval, 3))

  return(myfit)
  
}

f_cox_ordinal = function(data, time, event, riskfactor){
  
  
  mydata = data %>%
    
    dplyr::select(all_of(time), all_of(event), all_of(riskfactor)) %>%
    
    as.data.frame()
  
  mydata1 = mydata %>%
    
    select(-all_of(time), -all_of(event))
  
  tmp = lapply(mydata1, function(x){
    
    fit_cox_ordinal(mydata, time, event, x)
    
    # tryCatch(
    #   
    #   fit_cox_ordinal(mydata, time, event, x),
    #   error = function(e){NULL})
    
  })
  
  return(
    
    tryCatch(
      
      data.table::rbindlist(tmp, idcol = T),
      
      error = function(e){NULL}
      
    )
    
    
  )
  
}




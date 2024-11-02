# This R script applies univariable Mixture Cure models
# as proposed by Kuk and Chen (1992),
# for a given time-to-event outcome and multiple risk factors.

# For the implementation, the intsurv package is used.

# The output is a merged table that reports:
# Hazard ratios, 95% CI and p-value for the latency
# and
# Odds ratios, 95% CI and p-value for the incidence component.

# The user needs to provide:
# 1. The dataset
# 2. The time and the censoring variable
# 3. A vector with the risk factors
# 4. Whether Firth's penalty to be applied
# 5. The number of bootstraps for the standard error estimation

library(intsurv)

set.seed(99)

fit_cure_univ = function(data, time, event, rf, firth, n_boots){
    
  mydata2 = cbind(data, rf) %>%
    
    na.omit()
  
  myfit = cox_cure(~ rf,
                   ~ rf,
                   time = mydata2[, time],
                   event = mydata2[, event],
                   data = mydata2,
                   firth = firth,
                   bootstrap = n_boots)
    
  surv_matrix = summary(myfit)$surv_coef_mat %>%
      
      as.data.frame() %>%
      
      mutate(Type = 'Survival', Category = rownames(.))
    
  cure_matrix = summary(myfit)$cure_coef_mat %>%
      
      as.data.frame() %>%
      
      mutate(Type = 'Cure', Category = rownames(.))
    
    
  mymatrix = rbind(surv_matrix, cure_matrix)
    
  rownames(mymatrix) = NULL
    
  mymatrix = mymatrix[-grep('Intercept', mymatrix$Category), ]
  
  mymatrix$Category = gsub('.*rf', '', mymatrix$Category)
  
  mymatrix$CI = paste0(
    
    '(',
    round(exp(mymatrix$coef - 1.96 * mymatrix$`se(coef)`), 2),
    '-',
    round(exp(mymatrix$coef + 1.96 * mymatrix$`se(coef)`), 2),
    ')'
    
  )
  
  
  mymatrix$p = round(mymatrix$`Pr(>|z|)`, 3)
  
  mymatrix$Estimate = round(exp(mymatrix$coef), 2)
  
  mymatrix = mymatrix %>%
    
    select(
      
      Category, Type, Estimate, CI, p
      
    )
  
  mymatrix = split(mymatrix, mymatrix$Type)

  mymatrix = lapply(mymatrix, function(x){

    colnames(x)[-c(1:2)] = paste0(x$Type, '_', colnames(x)[-c(1:2)])

    x = x %>% select(-Type)

    return(x)


  })

  mymatrix = merge(mymatrix$Cure, mymatrix$Survival) %>%

    as.data.frame()

  
  return(mymatrix)
  
  
  
}


f_cure_updated <- function(data, time, event, riskfactor, firth = T, n_boots = 0){
  
  
  mydata = data %>%
    
    dplyr::select(any_of(time), any_of(event)) %>%
    
    as.data.frame()
  
  
  mydata1 = data %>%
    
    select(any_of(riskfactor))

    
  tmp = lapply(mydata1, function(x){
  
    tryCatch(
      
      fit_cure_univ(mydata, time, event, x, firth, n_boots),
      error = function(e){NULL})
    
  })
  
  
  return(
    
    tryCatch(
      
      data.table::rbindlist(tmp, idcol = T) %>%
        
        rename('Risk_Factor' = '.id'),
      
      error = function(e){NULL}
      
    )
    
    
  )
    
  
}





f_cure <- function(data, time, event, riskfactor){
  
  mymatrix = matrix(ncol = 7)
  
  colnames(mymatrix) = c('coef',  'exp(coef)', 'se(coef)',  'z', 'Pr(>|z|)', 'Type', 'RF')
  
  for (i in 1:length(riskfactor)){
    
    mydata <- data %>%
      dplyr::select(any_of(time), any_of(event), any_of(riskfactor[i])) %>%
      na.omit() %>%
      as.data.frame()
    
    myfit <- cox_cure(~ mydata[,riskfactor[i]],
                      ~ mydata[,riskfactor[i]],
                      time = mydata[,time],
                      event = mydata[,event],
                      data = mydata,
                      firth = T,
                      bootstrap = 50)
    
    surv_matrix = summary(myfit)$surv_coef_mat %>%
      
      as.data.frame() %>%
      
      mutate(Type = 'Survival',
             RF = riskfactor[i])
    
    cure_matrix = summary(myfit)$cure_coef_mat %>%
      
      as.data.frame() %>%
      
      mutate(Type = 'Cure',
             RF = riskfactor[i])
    
    
    
    mymatrix = rbind(mymatrix, surv_matrix, cure_matrix)
    
    tmp_matrix <- mymatrix
    
    tmp_matrix$Rows = rownames(mymatrix)
    
    tmp_matrix <<- tmp_matrix
    
  }
  
  mymatrix = mymatrix[-1, ]
  
  mymatrix$Category = rownames(mymatrix)
  
  mymatrix = mymatrix[-grep('Intercept', mymatrix$Category),]
  
  mymatrix$CI = paste0(
    
    '(',
    round(exp(mymatrix$coef - 1.96 * mymatrix$`se(coef)`), 2),
    '-',
    round(exp(mymatrix$coef + 1.96 * mymatrix$`se(coef)`), 2),
    ')'
    
  )
  
  # mymatrix$L = round(exp(mymatrix$coef - 1.96 * mymatrix$`se(coef)`), 2)
  # mymatrix$U = round(exp(mymatrix$coef + 1.96 * mymatrix$`se(coef)`), 2)
  
  mymatrix$p = round(mymatrix$`Pr(>|z|)`, 3)
  
  mymatrix$Estimate = round(exp(mymatrix$coef), 2)
  
  mymatrix = mymatrix %>%
    
    select(
      
      RF, Category, Type, Estimate, CI, p
      
    )
    
  
  return(mymatrix)
  
}

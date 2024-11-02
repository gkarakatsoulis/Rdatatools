# This R script automates the process for comparing differences
# between independent groups

# Step 1: Check normality assumption
# Step 2: Check Variance homogeneity
# Step 3: Choose the optimal test (anova, welch, kruskal-wallis)
# Step 4: Conduct the optimal test and get the results

library(dplyr)
library(rstatix)

# ---------------------------------------------------------
#           Check Normality assumption
# ---------------------------------------------------------


normality_function = function(data, y, x){
  
  mydata = data %>%
    
    select(any_of(y), any_of(x)) %>%
    
    gather(., key = 'Outcome', value = 'Outcome_Value', any_of(y)) %>%
    
    gather(., key = 'RiskFactor', value = 'RiskFactor_Value', any_of(x))
  
  
  groups_for_normality  = mydata %>%
    
    group_by(Outcome, RiskFactor, RiskFactor_Value) %>%
    
    count() %>%
    
    na.omit() %>%
    
    filter(n >=3 & n < 300) %>%
    
    transmute(
      
      group = paste(Outcome, RiskFactor, RiskFactor_Value)
      
    )
  
  
  
  if (nrow(groups_for_normality) > 0){
    
    non_normal = mydata %>%
      
      mutate(
        
        group = paste(Outcome, RiskFactor, RiskFactor_Value)
        
      ) %>%
      
      filter(group %in% groups_for_normality$group) %>%
      
      group_by(Outcome, RiskFactor, RiskFactor_Value) %>%
      
      summarise(
        
        pvalue = shapiro.test(Outcome_Value)$p.value
        
      ) %>%
      
      ungroup() %>%
      
      filter(pvalue < 0.05) %>%
      
      transmute(
        
        group = paste(Outcome, RiskFactor)
        
      ) %>%
      
      unlist() %>%
      
      unique()
    
    
  } else {
    
    non_normal = NA
    
  }
  
  
  mydata$Comparison = paste(mydata$Outcome, mydata$RiskFactor)
  
  
  mydata$normal = ifelse(mydata$Comparison %in% non_normal, F, T)
  
  total_results = mydata %>%
    
    select(Comparison, normal) %>%
    
    filter(!duplicated(Comparison))
  
  
  return(total_results)
  
  
}

# ---------------------------------------------------------
#           Check Variance Homogeneity
# ---------------------------------------------------------



levene_function = function(data, y, x){
  
  mydata = data %>%
    
    select(any_of(y), any_of(x)) %>%
    
    gather(., key = 'Outcome', value = 'Outcome_Value', any_of(y)) %>%
    
    gather(., key = 'RiskFactor', value = 'RiskFactor_Value', any_of(x)) %>%
    
    mutate_at('RiskFactor_Value', as.character)
  
  
  total_results = mydata %>%
    
    group_by(Outcome, RiskFactor) %>%
    
    mutate_if(is.character, as.factor) %>%
    
    levene_test(Outcome_Value ~ RiskFactor_Value) %>%
    
    transmute(
      
      Comparison = paste(Outcome, RiskFactor),
      
      Vareq = ifelse(p < 0.05, F, T)
      
    )
  
  
  return(total_results)
  
  
}

# ---------------------------------------------------------
#           Choose the optimal test (par vs nonpar)
# ---------------------------------------------------------


t_test_optimal = function(data, y, x){
  
  norm = normality_function(data, y, x)
  
  variance = levene_function(data, y, x)
  
  
  total_result = merge(norm, variance, by = 'Comparison')
  
  total_result$test = ifelse(!total_result$normal, 'nonpar',
                             ifelse(total_result$normal & total_result$Vareq,
                                    'anova', 'welch'))
  
  
  total_result = total_result %>%
    
    select(Comparison, test)
  
  return(total_result)
  
}



# ---------------------------------------------------------
#           Non-parametric test
# ---------------------------------------------------------


kruskal_function = function(data, y, x, long = F){
  
 if (long) {
   
   mydata = data %>%
     
     mutate_if(is.factor, as.numeric)
   
 } else {
   
   mydata = data %>%
     
     select(any_of(y), any_of(x)) %>%
     
     mutate_if(is.factor, as.numeric) %>%
     
     gather(., key = 'Outcome', value = 'Outcome_Value', any_of(y)) %>%
     
     gather(., key = 'RiskFactor', value = 'RiskFactor_Value', any_of(x))
   
 }
    
  
  total_comparisons = mydata %>%
    
    group_by(Outcome, RiskFactor) %>%
    
    kruskal_test(Outcome_Value ~ RiskFactor_Value) %>%
    
    transmute(Outcome, RiskFactor, p, Comparison = paste(Outcome, RiskFactor))
    
  comparisons_for_pairwise = c()
  
  if (any(total_comparisons$p < 0.05)){
    
    comparisons_for_pairwise = total_comparisons$Comparison[total_comparisons$p < 0.05]
    
  }
  
  mydata$Comparison = paste(mydata$Outcome, mydata$RiskFactor)
  
  
  if (length(comparisons_for_pairwise) >= 1){
    
    pairwise_comparisons = mydata %>%
      
      filter(Comparison %in% comparisons_for_pairwise) %>%
      
      group_by(Outcome, RiskFactor) %>%
      
      
      pairwise_wilcox_test(Outcome_Value ~ RiskFactor_Value,
                           
                           p.adjust.method = 'bonferroni') %>%
      
      
      select(Outcome, RiskFactor, group1, group2, p, p.adj)
    
    group1 = apply(pairwise_comparisons[,c('group1', 'group2')], 1, function(x){
      
      return(sort(x)[1])
      
    })
    
    group2 = apply(pairwise_comparisons[,c('group1', 'group2')], 1, function(x){
      
      return(sort(x)[2])
      
    })
    
    
    pairwise_comparisons$group1 = group1
    
    pairwise_comparisons$group2 = group2
    
    
  }
  
  
    descriptives = mydata %>%
      
      group_by(Outcome, RiskFactor, RiskFactor_Value) %>%
      
      na.omit() %>%
      
      summarise(
        
        result = paste0(
          
          round(mean(Outcome_Value, na.rm = T), 2),
          ' (',
          round(sd(Outcome_Value, na.rm = T), 2),
          ')'
          
        )
        
      ) %>%
      
      ungroup() %>%
      
      spread(., key = 'Outcome', value = 'result') 
    
    if (length(comparisons_for_pairwise) >= 1){
      
      total_results = list(
        
        # 'Descriptives' = descriptives,
        'Total Comparisons' = total_comparisons,
        'Post-Hoc comparisons' = pairwise_comparisons
        
      )
      
      
    } else {
      
      
      total_results = list(
        
        # 'Descriptives' = descriptives,
        'Total Comparisons' = total_comparisons
        
      )
      
  
      
    }
    
    
    return(total_results)
    
    
}


# ---------------------------------------------------------
#           Parametric test
# ---------------------------------------------------------


anova_function = function(data, y, x, long = F, vareq = T){
  
  if (long){
    
    mydata = data %>%
      
      mutate_if(is.factor, as.numeric)
    
  } else {
    
    mydata = data %>%
      
      mutate_if(is.factor, as.numeric) %>%
      
      select(any_of(y), any_of(x)) %>%
      
      gather(., key = 'Outcome', value = 'Outcome_Value', any_of(y)) %>%
      
      gather(., key = 'RiskFactor', value = 'RiskFactor_Value', any_of(x))
    
  }
  
  # ANOVA
  if (vareq){
    
    total_comparisons = mydata %>%
      
      filter(!is.na(RiskFactor_Value)) %>%
      
      group_by(Outcome, RiskFactor) %>%
      
      anova_test(Outcome_Value ~ RiskFactor_Value) %>%
      
      ungroup() %>%
      
      as.data.frame() %>%
      
      transmute(Outcome, RiskFactor, p, Comparison = paste(Outcome, RiskFactor))
    
    
  } else { # Welch 
    
    total_comparisons = mydata %>%
      
      filter(!is.na(RiskFactor_Value)) %>%
      
      group_by(Outcome, RiskFactor) %>%
      
      welch_anova_test(Outcome_Value ~ RiskFactor_Value) %>%
      
      mutate(Comparison = paste(Outcome, RiskFactor)) %>%
      
      select(Outcome, RiskFactor, p, Comparison)
    
    
  }
  

  comparisons_for_pairwise = c()
  
  if (any(total_comparisons$p < 0.05)){
    
    comparisons_for_pairwise = total_comparisons$Comparison[total_comparisons$p < 0.05]
    
  }
  
  
  mydata$Comparison = paste(mydata$Outcome, mydata$RiskFactor)
  
  
  if (length(comparisons_for_pairwise) >= 1){
    
    # Students t-test
    if (vareq){
      
      pairwise_comparisons = mydata %>%
        
        filter(Comparison %in% comparisons_for_pairwise) %>%
        
        group_by(Outcome, RiskFactor) %>%
        
        
        pairwise_t_test(Outcome_Value ~ RiskFactor_Value,
                        
                        p.adjust.method = 'bonferroni') %>%
        
        select(Outcome, RiskFactor, group1, group2, p, p.adj)
      
    } else { # Welch t-test
      
      pairwise_comparisons = mydata %>%
        
        filter(Comparison %in% comparisons_for_pairwise) %>%
        
        group_by(Outcome, RiskFactor) %>%
        
        pairwise_t_test(
        
        Outcome_Value ~ RiskFactor_Value,
        
        p.adjust.method = 'bonferroni') %>%
        
        select(Outcome, RiskFactor, group1, group2, p, p.adj)
      
    }
    
    
    
    group1 = apply(pairwise_comparisons[,c('group1', 'group2')], 1, function(x){
      
      return(sort(x)[1])
      
    })
    
    group2 = apply(pairwise_comparisons[,c('group1', 'group2')], 1, function(x){
      
      return(sort(x)[2])
      
    })
    
    
    pairwise_comparisons$group1 = group1
    
    pairwise_comparisons$group2 = group2
    
    
    
    
  }
  
  
  
  if (length(comparisons_for_pairwise) >= 1){
    
    total_results = list(
      
      # 'Descriptives' = descriptives,
      'Total Comparisons' = total_comparisons,
      'Post-Hoc comparisons' = pairwise_comparisons
      
    )
    
    
  } else {
    
    
    total_results = list(
      
      # 'Descriptives' = descriptives,
      'Total Comparisons' = total_comparisons
      
    )
    
    
    
  }
  
  return(total_results)
  
  
}



# ---------------------------------------------------------
#           Between-group differences
# ---------------------------------------------------------

between_group_diff_function = function(data, y, x){
  
 mydata = data %>%
    
    select(any_of(y), any_of(x)) %>%
    
    gather(., key = 'Outcome', value = 'Outcome_Value', any_of(y)) %>%
    
    gather(., key = 'RiskFactor', value = 'RiskFactor_Value', any_of(x)) %>%
    
    mutate(
      
      Comparison = paste(Outcome, RiskFactor)
      
    )
  
  mytest = t_test_optimal(data, y, x)
  
  nonpar_opt = mytest %>%
    
    filter(test == 'nonpar')
  
  data_nonpar = mydata %>%
    
    filter(Comparison %in% nonpar_opt$Comparison)
  
  
  anova_test_opt = mytest %>%
    
    filter(test == 'anova')
  
  data_anova = mydata %>%
    
    filter(Comparison %in% anova_test_opt$Comparison)
  
  welch_opt = mytest %>%
    
    filter(test == 'welch')
  
  data_welch = mydata %>%
    
    filter(Comparison %in% welch_opt$Comparison)
  
  total_results = list()
  
  
  total_results[['Descriptives']] = mydata %>%
    
    group_by(Outcome, RiskFactor, RiskFactor_Value) %>%
    
    na.omit() %>%
    
    summarise(
      
      'M (SD)' = paste0(
        
        round(mean(Outcome_Value, na.rm = T),2),
        ' (',
        round(sd(Outcome_Value, na.rm = T),2),
        ')'
        
      )
      
    ) %>%
    
    ungroup()
  
  total_results[['Total Comparisons']] = data.frame(matrix(ncol = 4))
  
  colnames(total_results[['Total Comparisons']]) = c('Outcome', 'RiskFactor', 'p', 'Comparison')
  
  
  total_results[['Post-Hoc comparisons']] = data.frame(matrix(ncol = 6))
  
  colnames(total_results[['Post-Hoc comparisons']]) = c('Outcome',
                                                     'RiskFactor',
                                                     'group1',
                                                     'group2',
                                                     'p', 'p.adj')
  
  
  
  if (nrow(data_nonpar) > 0){
    
    nonpar_result = kruskal_function(data_nonpar, y, x, long = T)
    
    total_results[['Total Comparisons']] = rbind(total_results[['Total Comparisons']],
                                                 nonpar_result[['Total Comparisons']])
    
    
    total_results[['Post-Hoc comparisons']] = rbind(total_results[['Post-Hoc comparisons']],
                                                 nonpar_result[['Post-Hoc comparisons']])
    
    
  }
  
  
  if (nrow(data_anova) > 0){
    
    anova_result = anova_function(data_anova, y, x, long = T, vareq = T)
    
    total_results[['Total Comparisons']] = rbind(total_results[['Total Comparisons']],
                                                 anova_result[['Total Comparisons']])
    
    total_results[['Post-Hoc comparisons']] = rbind(total_results[['Post-Hoc comparisons']],
                                                 nonpar_result[['Post-Hoc comparisons']])
  }
  
  if (nrow(data_welch) > 0){
    
    welch_result = anova_function(data_welch, y, x, long = T, vareq = F)
    total_results[['Total Comparisons']] = rbind(total_results[['Total Comparisons']],
                                                 welch_result[['Total Comparisons']])
    
    total_results[['Post-Hoc comparisons']] = rbind(total_results[['Post-Hoc comparisons']],
                                                 welch_result[['Post-Hoc comparisons']])
    
  }
  
  total_results[['Total Comparisons']] = total_results[['Total Comparisons']][-1, ]
  
  total_results[['Total Comparisons']] = merge(
    
    total_results[['Descriptives']],
    total_results[['Total Comparisons']]
    
  ) %>%
    
    mutate(p = ifelse(duplicated(Comparison), NA, p)) %>%
    
    select(-Comparison)
  
  total_results[['Descriptives']] = NULL
  
  total_results[['Post-Hoc comparisons']] = total_results[['Post-Hoc comparisons']][-1, ]
  
  
  return(total_results)
  
}


# This R script investigates independence between qualitative variables

# To compare independence between two categorical variables
# the dominant tests are the chi-squared and the Fisher's exact test.
# The former can be used only if some further assumptions hold.

# This script takes two sets of variables: 
# one acting as the outcomes and
# one acting as the risk factors that we want to associate with the outcomes

# It first checks whether the chi-squared test is appropriate.
# If yes, it is applied, otherwise the Fisher's exact test is applied.

# The output is a list with:
# (1) A table with the p-value and the test conducted for 
# for each combination between an outcome and a risk factor.
# (2) A merged contingency table for each outcome with the risk factors.


# The user needs to provide:
# 1. The dataset
# 2. A vector of the outcome
# 3. A vector of the risk factors


library(dplyr)
library(rstatix)


chisq_function = function(data, y, x){
  
  total_results = list()
  
  mydata = data %>%
    
    select(any_of(y), any_of(x)) %>%
    
    gather(., key = 'Outcome', value = 'Outcome_Value', any_of(y)) %>%
    
    gather(., key = 'RiskFactor', value = 'RiskFactor_Value', any_of(x)) %>%
    
    mutate(
      
      Comparison = paste(Outcome, RiskFactor)
      
    )
    
  
  total_comparisons = mydata %>%
    
    group_by(Outcome, RiskFactor) %>%
    
    summarise(
      
      prop.expected = any(
        
        as.vector(
          
          
          chisq.test(
            
            x = .data[['Outcome_Value']],
            y=  .data[['RiskFactor_Value']]
          )$expected
          
        ) < 5
        
        
      ),
      
      p = chisq.test(x = .data[['Outcome_Value']], y=  .data[['RiskFactor_Value']])$p.value,
      
      test = 'Chi-square'
      
    ) %>%
    
    ungroup()
  
  
  total_comparisons$Comparison = paste(
    
    total_comparisons$Outcome, total_comparisons$RiskFactor
    
  )
  
  
  if (any(total_comparisons$prop.expected)){
    
    fisher_comparisons = total_comparisons$Comparison[which(
      
      total_comparisons$prop.expected
      
    )]
    
    total_comparisons_chisq  = total_comparisons %>%
      
      filter(!prop.expected)
    
    
    total_comparisons_fisher = mydata %>%
      
      filter(Comparison %in% fisher_comparisons) %>%
      
      group_by(Outcome, RiskFactor) %>%
      
      summarise(
        
        p = fisher.test(
          
          x = .data[['Outcome_Value']],
          y =  .data[['RiskFactor_Value']]
          
          )$p.value,
        
        test = 'Fisher',
        
        prop.expected = T
        
      ) %>%
      
      ungroup() %>%
      
      mutate(
        
        Comparison = paste(Outcome, RiskFactor)
        
      )
    
    
    total_comparisons = rbind(
      
      total_comparisons_chisq,
      total_comparisons_fisher
      
    )
    
  }
  
  total_comparisons = total_comparisons %>%
    
    select(-Comparison, - prop.expected)
    
  
  descriptives = mydata %>%
    
    freq_table(Outcome, RiskFactor, Outcome_Value, RiskFactor_Value)
  
  descriptives = split(descriptives, f = descriptives$Outcome) %>%
    
    lapply(., function(x){
      
      return(x %>%
               
               mutate(
                 
                 res = paste0(
                   n,
                   ' (',
                   prop,
                   '%)'
                 )
                 
                 
               ) %>%
               
               select(RiskFactor, Outcome_Value, RiskFactor_Value, res) %>%
               
               spread(., key = 'Outcome_Value', value = 'res', fill = '0 (0%)')
      )
      
      
    })
  
  
  
  
  total_results = c(
    
    list(total_comparisons),
    
    descriptives
    
  )
  
  names(total_results)[1] <- 'Total Comparisons'
  
  return(total_results)
  
  
}






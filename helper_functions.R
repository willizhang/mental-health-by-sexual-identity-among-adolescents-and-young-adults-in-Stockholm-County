# title: "Helper functions"
# author: Willi Zhang
# email: willi.zhang@ki.se


### Function to Categorize Age ###

categorize_age_group <- function( age ) {
  factor(
    
    case_when(
      age <= 29 ~ "<=29",
      age >= 30 & age <= 44 ~ "30–44",
      age >= 45 & age <= 59 ~ "45–59",
      age >= 60             ~ "\u226560" # \u2265 is the Unicode for ≥
      ),
    
    levels = c( "<=29", "30–44", "45–59", "\u226560" ) )
  }


##### Individual Surveys #####

# overall prevalence
cal_prev_indi_overall_cc <- function( variables_list, design, year ) {
  
  results <- list()
  
  for ( var in variables_list ) {
    prop <- svyciprop( formula = as.formula( paste0( "~ I(", var$variable, " == '", var$condition, "')" ) ),
                       design = subset( design, !is.na( get( var$variable ) ) ),
                       vartype = "ci",
                       method = "beta" )
    
    ci <- attr( prop, "ci" )
    
    df <- data.frame( point_estimate = as.vector( prop ),
                      lower_ci = ci[ 1 ],
                      upper_ci = ci[ 2 ],
                      row.names = NULL )
    
    colnames( df ) <- c( paste0( var$name, "_", year, "_point_estimate" ),
                         paste0( var$name, "_", year, "_lower_ci" ),
                         paste0( var$name, "_", year, "_upper_ci" ) )
    
    results[[ var$name ]] <- df
    }
  
  return( results )
  }


# prevalence by sexual identity
cal_prev_indi_cc <- function( variables_list, design, year, group_var, subgroup = FALSE ) {
  
  results <- list()
  
  for ( var in variables_list ) {
    svyby_result <- svyby(
      formula = as.formula( paste0( "~ I(", var$variable, " == '", var$condition, "')" ) ),
      by = as.formula( paste0( "~", group_var ) ),
      design = subset( design, !is.na( get( var$variable ) ) ),
      FUN = svyciprop,
      vartype = "ci",
      method = "beta"
      ) %>% remove_rownames()
    
    if ( subgroup ) {
      colnames( svyby_result )[ -c( 1, 2 ) ] <- c(
        paste0( var$name, "_", year, "_point_estimate" ),
        paste0( var$name, "_", year, "_lower_ci" ),
        paste0( var$name, "_", year, "_upper_ci" )
      )
    } else {
      colnames( svyby_result )[ -1 ] <- c(
        paste0( var$name, "_", year, "_point_estimate" ),
        paste0( var$name, "_", year, "_lower_ci" ),
        paste0( var$name, "_", year, "_upper_ci" )
      )
    }
    
    results[[ var$name ]] <- svyby_result
  }
  
  return( results )
  }


# calculate prevalence ratio
cal_pr_indi_cc <- function( design, exposure, outcome_list, year, covariates = NULL ) {
  
  model_list <- list()
  
  for ( outcome in outcome_list ) {
    
    formula_str <- paste0( outcome$variable, 
                           " ~ ",
                           exposure,
                           if ( !is.null( covariates ) ) paste0( " + ", covariates ) else "" )
    
    mod <- svyglm( formula = as.formula( formula_str ),
                   design = design, 
                   family = quasipoisson( link = "log" ) ) # Poisson regression
    
    model_list[[ paste0( outcome$name, "_", year ) ]] <- mod
    }
  return( model_list ) 
}


# calculate prevalence difference
cal_pd_indi_cc <- function( design, exposure, outcome_list, year, covariates = NULL ) {
  
  model_list <- list()
  
  for ( outcome in outcome_list ) {
    
    formula_str <- paste0( outcome$variable, 
                           " ~ ",
                           exposure,
                           if ( !is.null( covariates ) ) paste0( " + ", covariates ) else "" )
    
    mod <- svyglm( formula = as.formula( formula_str ),
                   design = design, 
                   family = gaussian( link = "identity" ) )
    
    model_list[[ paste0( outcome$name, "_", year ) ]] <- mod
  }
  return( model_list ) 
}



##### Pooled Sample #####

# annual prevalence by sexual identity
cal_prev_poo_cc <- function( year_start, year_end, data, year_var, variables_list, group_var, subgroup = FALSE ) {
  
  results <- list()
  
  for ( yr in seq( year_start, year_end ) ) {
    
    data_year <- data %>% filter( .data[[ year_var ]] == yr )
    
    design <- svydesign( ids = ~ 1, data = data_year )
    
    for ( var in variables_list ) {
      
      svyby_result <- svyby(
        formula = as.formula( paste0( "~ I(", var$variable, " == '", var$condition, "')" ) ),
        by = as.formula( paste0( "~", group_var ) ),
        design = subset( design, !is.na( get( var$variable ) ) ),
        FUN = svyciprop,
        vartype = "ci",
        method = "beta" 
      ) %>% 
        remove_rownames()
      
      if ( subgroup ) {
        colnames( svyby_result )[ -c( 1, 2 ) ] <- c( "point_estimate", "lower_ci", "upper_ci" )

      } else {
        colnames( svyby_result )[ -1 ] <- c( "point_estimate", "lower_ci", "upper_ci" )
      }
      
      svyby_result <- svyby_result %>%
        mutate(
          outcome = var$name,
          year = yr
        )
      
      results[[ paste0( var$name, "_", yr ) ]] <- svyby_result
    }
  }
  
  return( results )
  }

# cumulative prevalence by sexual identity
cal_prev_poo_cc_cum <- function( year_start, year_end, data, year_var, variables_list, group_var, subgroup = FALSE ) {
  
  results <- list()
  
  for ( yr in seq( year_start, year_end ) ) {
    
    # cumulative follow-up up to current year
    data_year <- data %>% filter( .data[[ year_var ]] <= yr )
    
    for ( var in variables_list ) {
      
      # ever had outcome up to year yr
      data_cum <- data_year %>%
        group_by( lopnr ) %>%
        mutate(
          cum_outcome = if_else(
            any( .data[[ var$variable ]] == var$condition ),
            "Yes",
            "No"
            )
          ) %>%
        slice_tail( n = 1 ) %>%   # retain one row per individual
        ungroup()
      
      design <- svydesign( ids = ~1, data = data_cum )
      
      svyby_result <- svyby(
        formula = ~I( cum_outcome == "Yes" ),
        by = as.formula( paste0( "~", group_var ) ),
        design = design,
        FUN = svyciprop,
        vartype = "ci",
        method = "beta"
        ) %>%
        remove_rownames()
      
        if ( subgroup ) {
          colnames( svyby_result )[ -c( 1, 2 ) ] <- c( "point_estimate", "lower_ci", "upper_ci" )
          
        } else {
          colnames( svyby_result )[ -1 ] <- c( "point_estimate", "lower_ci", "upper_ci" )
        }
        
        svyby_result <- svyby_result %>%
          mutate(
            outcome = var$name,
            cum_year = yr
          )
        
        results[[ paste0( var$name, "_cum_", yr ) ]] <- svyby_result
    }
  }
  
  return( results )
}


# calculate prevalence ratio
cal_pr_poo_cc <- function( year_start, year_end, data, year_var, exposure, outcome_list, covariates = NULL ) {
  
  model_list <- list()
  
  for ( yr in seq( year_start, year_end ) ) {
    
    data_year <- data %>% filter( .data[[ year_var ]] == yr )
    
    design <- svydesign( ids = ~1, data = data_year )
    
    for ( outcome in outcome_list ) {
      
      formula_str <- paste0( outcome$variable, 
                             " ~ ", 
                             exposure, 
                             if ( !is.null( covariates ) ) paste0( " + ", covariates ) else "" )
      
      mod <- svyglm( formula = as.formula( formula_str ),
                     design = design, 
                     family = quasipoisson( link = "log" ) ) # Poisson regression
  
      model_list[[ paste0( outcome$name, "_", yr ) ]] <- mod
      
    }
  }
  
  return( model_list )
}


# calculate prevalence difference
cal_pd_poo_cc <- function( year_start, year_end, data, year_var, exposure, outcome_list, covariates = NULL ) {
  
  model_list <- list()
  
  for ( yr in seq( year_start, year_end ) ) {
    
    data_year <- data %>% filter( .data[[ year_var ]] == yr )
    
    design <- svydesign( ids = ~1, data = data_year )
    
    for ( outcome in outcome_list ) {
      
      formula_str <- paste0( outcome$variable, 
                             " ~ ", 
                             exposure, 
                             if ( !is.null( covariates ) ) paste0( " + ", covariates ) else "" )
      
      mod <- svyglm( formula = as.formula( formula_str ),
                     design = design, 
                     family = gaussian( link = "identity" ) )
      
      model_list[[ paste0( outcome$name, "_", yr ) ]] <- mod
      
    }
  }
  
  return( model_list )
}


##### Individual Surveys and Pooled Sample #####

# extract prevalence ratio
extract_pr_cc <- function( model_list, exposure_var ) {
  
  results <- lapply( names( model_list ), function( name ) {
    
    mod <- model_list[[ name ]]
    
    outcome <- sub( "_[0-9]{4}$", "", name )
    year <- as.numeric( sub( ".*_", "", name ) )
    
    coef_est <- coef( mod )
    ci_est <- confint( mod, ddf = degf( mod$survey.design ) )
    
    exposure_terms <- grep( paste0( "^", exposure_var ),
                            names( coef_est ),
                            value = TRUE )
    
    data.frame(
      outcome = outcome,
      year = year,
      exposure_level = sub( paste0( exposure_var ), "", exposure_terms ),
      estimate = exp( coef_est[ exposure_terms ] ),
      lower_ci = exp( ci_est[ exposure_terms, 1 ] ),
      upper_ci = exp( ci_est[ exposure_terms, 2 ] ),
      row.names = NULL
    )
  } )
  
  bind_rows( results )
}

# extract prevalence ratio by sex
extract_pr_by_sex_cc <- function( model_list, exposure_var ) {
  
  results <- lapply( names( model_list ), function( name ) {
    
    mod <- model_list[[ name ]]
    
    outcome <- sub( "_[0-9]{4}$", "", name )
    year <- as.numeric( sub( ".*_", "", name ) )
    
    coef_est <- coef( mod )
    ci_est <- confint( mod, ddf = degf( mod$survey.design ) )
    
    exposure_terms <- grep( paste0( "^", exposure_var ),
                            names( coef_est ),
                            value = TRUE )
    
    exposure_terms <- exposure_terms[ !grepl( ":", exposure_terms ) ]
    
    data.frame(
      outcome = outcome,
      year = year,
      exposure_level = sub( paste0( exposure_var ), "", exposure_terms ),
      estimate = exp( coef_est[ exposure_terms ] ),
      lower_ci = exp( ci_est[ exposure_terms, 1 ] ),
      upper_ci = exp( ci_est[ exposure_terms, 2 ] ),
      row.names = NULL
    )
  } )
  
  bind_rows( results )
}

# extract prevalence difference
extract_pd_cc <- function( model_list, exposure_var ) {
  
  results <- lapply( names( model_list ), function( name ) {
    
    mod <- model_list[[ name ]]
    
    outcome <- sub( "_[0-9]{4}$", "", name )
    year <- as.numeric( sub( ".*_", "", name ) )
    
    coef_est <- coef( mod )
    ci_est <- confint( mod, ddf = degf( mod$survey.design ) )
    
    exposure_terms <- grep( paste0( "^", exposure_var ),
                            names( coef_est ),
                            value = TRUE )
    
    data.frame(
      outcome = outcome,
      year = year,
      exposure_level = sub( paste0( exposure_var ), "", exposure_terms ),
      estimate = coef_est[ exposure_terms ],
      lower_ci = ci_est[ exposure_terms, 1 ],
      upper_ci = ci_est[ exposure_terms, 2 ],
      row.names = NULL
    )
  } )
  
  bind_rows( results )
}

# extract prevalence difference by sex
extract_pd_by_sex_cc <- function( model_list, exposure_var ) {
  
  results <- lapply( names( model_list ), function( name ) {
    
    mod <- model_list[[ name ]]
    
    outcome <- sub( "_[0-9]{4}$", "", name )
    year <- as.numeric( sub( ".*_", "", name ) )
    
    coef_est <- coef( mod )
    ci_est <- confint( mod, ddf = degf( mod$survey.design ) )
    
    exposure_terms <- grep( paste0( "^", exposure_var ),
                            names( coef_est ),
                            value = TRUE )
    
    exposure_terms <- exposure_terms[ !grepl( ":", exposure_terms ) ]
    
    data.frame(
      outcome = outcome,
      year = year,
      exposure_level = sub( paste0( exposure_var ), "", exposure_terms ),
      estimate = coef_est[ exposure_terms ],
      lower_ci = ci_est[ exposure_terms, 1 ],
      upper_ci = ci_est[ exposure_terms, 2 ],
      row.names = NULL
    )
  } )
  
  bind_rows( results )
}

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



##### Pooled Cohort #####

# prevalence by sexual identity
cal_prev_poo_cc <- function( year_start, year_end, data, variables_list, group_var, subgroup = FALSE ) {
  
  results <- list()
  
  for ( yr in seq( year_start, year_end ) ) {
    
    data_year <- data %>% filter( year == yr )
    
    design <- svydesign( ids = ~1, data = data_year )
    
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
      
      group_vars <- strsplit( group_var, "\\s*\\+\\s*")[[1]]
      
      sample_size <- data_year %>%
        drop_na( all_of( c( var$variable, group_vars ) ) ) %>%
        group_by( across( all_of( group_vars ) ) ) %>%
        summarise( n = n(), .groups = "drop" )
      
      svyby_result <- svyby_result %>%
        left_join( sample_size, by = group_vars )
      
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


# calculate prevalence ratio
cal_pr_poo_cc <- function( year_start, year_end, data, exposure, outcome_list, covariates = NULL ) {
  
  model_list <- list()
  
  for ( yr in seq( year_start, year_end ) ) {
    
    data_year <- data %>% filter( year == yr )
    
    design <- svydesign( ids = ~1, data = data_year )
    
    for ( outcome in outcome_list ) {
      
      formula_str <- paste0( outcome$variable, 
                             " ~ ", 
                             exposure, 
                             if ( !is.null( covariates ) ) paste0( " + ", covariates ) else "" )
      
      mod <- svyglm( formula = as.formula( formula_str ),
                     design = subset( design, !is.na( get( outcome$variable ) ) & !is.na( get( exposure ) ) ), 
                     family = quasipoisson( link = "log" ) ) # Poisson regression
  
      model_list[[ paste0( outcome$name, "_", yr ) ]] <- mod
      
    }
  }
  
  return( model_list )
}


# calculate prevalence difference
cal_rd_poo_cc <- function( year_start, year_end, data, exposure, outcome_list, covariates = NULL ) {
  
  model_list <- list()
  
  for ( yr in seq( year_start, year_end ) ) {
    
    data_year <- data %>% filter( year == yr )
    
    design <- svydesign( ids = ~1, data = data_year )
    
    for ( outcome in outcome_list ) {
      
      formula_str <- paste0( outcome$variable, 
                             " ~ ", 
                             exposure, 
                             if ( !is.null( covariates ) ) paste0( " + ", covariates ) else "" )
      
      mod <- svyglm( formula = as.formula( formula_str ),
                     design = subset( design, !is.na( get( outcome$variable ) ) & !is.na( get( exposure ) ) ), 
                     family = gaussian( link = "identity" ) )
      
      model_list[[ paste0( outcome$name, "_", yr ) ]] <- mod
      
    }
  }
  
  return( model_list )
}



##### Individual Surveys and Pooled Cohort #####

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

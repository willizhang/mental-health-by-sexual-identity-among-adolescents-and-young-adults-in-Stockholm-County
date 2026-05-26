# characteristics of participants from individual SPHC-B surveys, summarized across all five surveys
# run each survey analysis before running the following scripts

d_2002_young <- d_2002 %>%
  select( lopnr, age, sex, country_of_birth, education, sexual_identity_2010 ) %>%
  filter( age <= 29 ) %>%
  rename( sexual_identity = sexual_identity_2010 ) %>%
  mutate( survey_year = "SPHC-B 2002" )

d_2006_young <- d_2006 %>%
  select( lopnr, age, sex, country_of_birth, education, sexual_identity_2010 ) %>%
  filter( age <= 29 ) %>%
  rename( sexual_identity = sexual_identity_2010 ) %>%
  mutate( survey_year = "SPHC-B 2006" )

d_2010_young <- d_2010 %>%
  select( lopnr, age, sex, country_of_birth, education, sexual_identity_2010 ) %>%
  filter( age <= 29 ) %>%
  rename( sexual_identity = sexual_identity_2010 ) %>%
  mutate( survey_year = "SPHC-B 2010" )

d_2014_young <- d_2014 %>%
  select( lopnr, age, sex, country_of_birth, education, sexual_identity_2014 ) %>%
  filter( age <= 29 ) %>%
  rename( sexual_identity = sexual_identity_2014 ) %>%
  mutate( survey_year = "SPHC-B 2014" )

d_2021_young <- d_2021 %>%
  select( lopnr, age, sex, country_of_birth, education, sexual_identity_2021 ) %>%
  filter( age <= 29 ) %>%
  rename( sexual_identity = sexual_identity_2021 ) %>%
  mutate( survey_year = "SPHC-B 2021" )

d_pooled_cha <- bind_rows(
  d_2002_young,
  d_2006_young,
  d_2010_young,
  d_2014_young,
  d_2021_young ) %>%
  mutate( survey_year = as.factor( survey_year ) )

summary( d_pooled_cha )

# check for overlapping participants across SPHC-B 2002 to 2014
# unable to check for SPHC-B 2021, because of its independent lopnr number system

sum( duplicated( d_pooled_cha$lopnr[ !is.na( d_pooled_cha$lopnr ) ] ) ) # 62 duplicates

length( intersect( d_2002_young$lopnr, d_2006_young$lopnr ) ) # 60 individuals participated in both SPHC-B 2002 and 2006
length( intersect( d_2002_young$lopnr, d_2010_young$lopnr ) ) # 0
length( intersect( d_2002_young$lopnr, d_2014_young$lopnr ) ) # 0
length( intersect( d_2006_young$lopnr, d_2010_young$lopnr ) ) # 0
length( intersect( d_2006_young$lopnr, d_2014_young$lopnr ) ) # 2
length( intersect( d_2010_young$lopnr, d_2014_young$lopnr ) ) # 0
length( Reduce( intersect, list( d_2002_young$lopnr, d_2006_young$lopnr, d_2014_young$lopnr ) ) ) # 0

# make characteristics table by sexual identity
explanatory =  c( "age", "sex", "country_of_birth", "education" )
dependent = "sexual_identity"

d_poo_table <- d_pooled_cha %>% 
  summary_factorlist( dependent,
                      explanatory, 
                      na_include = TRUE,
                      na_include_dependent = TRUE, 
                      total_col = TRUE,
                      add_col_totals = TRUE,
                      column = FALSE )

writexl::write_xlsx( d_poo_table, "results_output/characteristics_table_poo_indi.xlsx" )

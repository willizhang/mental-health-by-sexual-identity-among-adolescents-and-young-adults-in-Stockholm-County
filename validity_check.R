##### validity check of constructing study sample #####

# exclude move-in/out/death during the year
event_during_year_ids <- eligible_dates %>%
  filter( date >= followup_start, 
          date <= followup_end,
          event_type != "death" ) %>%
  pull( lopnr ) %>%
  unique()

# exclude those whose last record before the year was a move-out or death
exit_before_year_ids <- eligible_dates %>%
  filter( date < followup_start ) %>%
  arrange( lopnr, date ) %>%
  group_by( lopnr ) %>%
  slice_tail( n = 1 ) %>%
  ungroup() %>%
  filter( event_type %in% c( "move-out", "death" ) ) %>%
  pull( lopnr )

# exclude those whose first move-in was after the year
enter_after_year_ids <- eligible_dates %>%
  arrange( lopnr, date ) %>%
  group_by( lopnr ) %>%
  slice( 1 ) %>%
  ungroup() %>%
  filter( date > followup_end, event_type == "move-in" ) %>%
  pull( lopnr )

excluded_year_ids_check <- Reduce( union, list( 
  event_during_year_ids, exit_before_year_ids, enter_after_year_ids ) 
)

identical(
  sort( excluded_year_ids ),
  sort( excluded_year_ids_check )
  ) # true



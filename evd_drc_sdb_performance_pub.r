#..........................................................................................
### +++++++ FACTORS ASSOCIATED WITH EBOLA SAFE AND DIGNIFIED BURIAL PERFORMANCE +++++++ ###
#..........................................................................................

#..........................................................................................
## ----------------- R CODE TO PREPARE DATA AND FIT STATISTICAL MODELS ----------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (July 2020)
                                          # francesco.checchi@lshtm.ac.uk 


#..........................................................................................
### Preparatory steps
#..........................................................................................

  #...................................      
  ## Install or load required R packages

    # List of required packages
    x1 <- c("broom.mixed", "data.table", "ggplot2", "ggpubr", "lme4", "lubridate", 
      "RColorBrewer", "readxl", "scales")
    
    # Install any packages not yet installed
    x2 <- x1 %in% row.names(installed.packages())
    if (any(x2 == FALSE)) { install.packages(x1[! x2]) }

    # Load all packages    
    lapply(x1, library, character.only = TRUE)
    
  #...................................      
  ## Starting setup

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
    # Set font
    windowsFonts(Arial=windowsFont("Arial"))

    # Set working directory to where this file is stored
    current_path = rstudioapi::getActiveDocumentContext()$path 
    setwd(dirname(current_path ))
    print( getwd() )
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    show_col(palette_cb)

#.........................................................................................
### Reading in required files
#.........................................................................................
    
  #...................................      
  ## Variable dictionary
  dict <- read_excel("evd_drc_sdb_performance_datasets_pub.xlsx", sheet = "dictionary")
    # remove tibble nonsense
    dict <- as.data.frame(dict)
    # dataset names
    dnames <- unique(dict[, "tab"])
      
  #...................................      
  ## Read in all the datasets
    # for each dataset...
    for (i in dnames) {
      # read in
      assign(i, read_excel("evd_drc_sdb_performance_datasets_pub.xlsx", sheet = i) )
        # remove tibble nonsense
        assign(i, as.data.frame(get(i)) )
      # only keep needed columns
      x1 <- subset(dict, tab == i)[, "use"]
      x1 <- which(! x1 %in% c("no") )
      x2 <- get(i)[, x1]
      assign(i, x2)
    }
    

#.........................................................................................                            
### Generating time series of health zones and time
#.........................................................................................
    
  #...................................       
  ## Determine start and end of analysis period (based on extent of SDB dataset)  
    # Start
    year_start <- min(sdb_dataset[, "epi_year"])
    month_start <- min(sdb_dataset[sdb_dataset$year == year_start, "month"]) 
    x1 <- sdb_dataset[sdb_dataset$epi_year == year_start, "epi_week"]
    week_start <- min(x1)
  
    # End
    year_end <- max(sdb_dataset[, "epi_year"])
    month_end <- max(sdb_dataset[sdb_dataset$year == year_end, "month"]) 
    x1 <- sdb_dataset[sdb_dataset$epi_year == year_end, "epi_week"]
    week_end <- max(x1)    

  #...................................    
  ## Time series of health zones and epidemiological weeks
    
    # Create time series
    hz <- unique(admin_units[, "hz"])
    weeks <- c(c((year_start*100 + week_start) : (year_start*100 + 52) ) , 
               c((year_end*100 + 1) : (year_end*100 + week_end)) )
    tsw <- expand.grid(hz, weeks)
    colnames(tsw) <- c("hz", "weeks")
    tsw[, "epi_year"] <- floor(tsw[, "weeks"] / 100)
    tsw[, "epi_week"] <- tsw[, "weeks"] %% 100
    tsw <- tsw[order(tsw[, "hz"], tsw[, "weeks"]), ]
    tsw[, "hz"] <- as.character(tsw[, "hz"])

    # Remove health areas from administrative database
    admin_units <- unique(admin_units[, c("province", "territory", "hz")])


#.........................................................................................                            
### Making additional preparations / creating additional variables in explanatory datasets
#.........................................................................................

  #...................................   
  ## Merge explanatory datasets and population denominators onto time series
    # Merge admin units
    tsw <- merge(tsw, admin_units, by="hz", all.x = TRUE)
 
    # Merge explanatory datasets and population denominator
    tsw <- merge(tsw, evd_cases_weeks, by=c("hz", "epi_year", "epi_week"), all.x = TRUE)
    tsw <- merge(tsw, attacks_evd_weeks, by=c("hz", "epi_year", "epi_week"), all.x = TRUE)
    tsw <- merge(tsw, insecurity_acled_weeks, by=c("territory", "epi_year", "epi_week"), all.x = TRUE)
    tsw <- merge(tsw, etc_presence[, ! colnames(etc_presence) %in% c("month", "year")], by=c("hz", "epi_year", "epi_week"), all.x = TRUE)
    tsw <- merge(tsw, health_facilities, by="hz", all.x = TRUE)
    tsw <- merge(tsw, cell_coverage, by="hz", all.x = TRUE)
    tsw <- merge(tsw, radio_coverage, by="hz", all.x = TRUE)
    tsw <- merge(tsw, road_coverage, by="hz", all.x = TRUE)
    tsw <- merge(tsw, mines, by="hz", all.x = TRUE)
    tsw <- merge(tsw, population, by="hz", all.x = TRUE)
    
        
  #...................................   
  ## Manage data on EVD cases
    # Set all NA case totals to zero
    x1 <- c("confirmed", "probable", "suspected", "alive", "dead", "unknown")
    tsw[, x1][is.na(tsw[, x1])] <- 0
    tsm[, x1][is.na(tsm[, x1])] <- 0
    
    # Calculate running cumulative number of confirmed cases (as a cumulative incidence rate per 100000 population)
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    hz <- sort(hz)
    x1 <- c()
    for (i in hz) {
      # select health zone time series
      x2 <- subset(tsw, hz == i)
      # compute running sum
      x3 <- c()
      for (j in 1:nrow(x2) ) { x3[j] <- sum(x2[max(0, j-6):max(0, j-1), "confirmed"], na.rm = TRUE) }
      # add results
      x1 <- c(x1, x3)
    }
    # add results
    tsw[, "confirmed_6weeks_rate"] <- x1 * 100000 / tsw[, "pop"]

    # Stage in the epidemic for each health zone: pre-epidemic, during epidemic, epidemic over
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    hz <- sort(hz)
    x1 <- c()
    for (i in hz) {
      # compute cumulative sum of confirmed cases for health zone time series
      x2 <- cumsum(subset(tsw, hz == i)[, "confirmed"])
      # first label all time points as pre-epidemic
      x3 <- rep("pre-epidemic", length(x2) )
      # then detect first time point with cases and label time points from then on as epidemic
      if (sum(x2) > 0) {x3[min(which(x2 > 0)):length(x3)] <- "epidemic"}
      # then detect last time point with an increases in cases and label subsequent time points as post-epidemic
      if (sum(x2) > 0 & which.max(x2) < length(x2) ) {x3[(which.max(x2) + 1):length(x2)] <- "post-epidemic"}
      # add results
      x1 <- c(x1, x3)
    }
      # add results
      tsw[, "epidemic_stage"] <- x1
   
      
  #...................................   
  ## Manage data on attacks against EVD response
    # Set all NA case totals to zero
    x1 <- c("against_evd", "suspend_activities")
    tsw[, x1][is.na(tsw[, x1])] <- 0
    tsm[, x1][is.na(tsm[, x1])] <- 0
    
    # Whether any attack occurred recently
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    hz <- sort(hz)
    x1 <- c()
    for (i in hz) {
      # select health zone time series
      x2 <- subset(tsw, hz == i)
      # compute running sum
      x3 <- c()
      for (j in 1:nrow(x2) ) { x3[j] <- ifelse(sum(x2[max(0, j-6):max(0, j-1), "against_evd"], na.rm = TRUE) > 0, "yes", "no") }
      # add results
      x1 <- c(x1, x3)
    }
    # add results
    tsw[, "against_evd_6weeks"] <- x1

  # Whether any suspension occurred recently
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    hz <- sort(hz)
    x1 <- c()
    for (i in hz) {
      # select health zone time series
      x2 <- subset(tsw, hz == i)
      # compute running sum
      x3 <- c()
      for (j in 1:nrow(x2) ) { x3[j] <- ifelse(sum(x2[max(0, j-6):max(0, j-1), "suspend_activities"], na.rm = TRUE) > 0, "yes", "no") }
      # add results
      x1 <- c(x1, x3)
    }
    # add results
    tsw[, "suspend_activities_6weeks"] <- x1

        
  #...................................   
  ## Manage data on insecurity events
    # Set all NA case totals to zero
    x1 <- c("n_events", "fatalities")
    tsw[, x1][is.na(tsw[, x1])] <- 0
    tsm[, x1][is.na(tsm[, x1])] <- 0
    
    # Calculate running cumulative number of events (as a rate per 100,000 population)
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    hz <- sort(hz)
    x1 <- c()
    for (i in hz) {
      # select health zone time series
      x2 <- subset(tsw, hz == i)
      # compute running sum
      x3 <- c()
      for (j in 1:nrow(x2) ) { x3[j] <- sum(x2[max(0, j-2):max(0, j-1), "n_events"], na.rm = TRUE) }
      # add results
      x1 <- c(x1, x3)
    }
    # add results
    tsw[, "n_events_2weeks_rate"] <- x1 * 100000 / tsw[, "pop"] 

    # Calculate running cumulative number of fatalities (as a rate per 100,000 population)
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    hz <- sort(hz)
    x1 <- c()
    for (i in hz) {
      # select health zone time series
      x2 <- subset(tsw, hz == i)
      # compute running sum
      x3 <- c()
      for (j in 1:nrow(x2) ) { x3[j] <- sum(x2[max(0, j-2):max(0, j-1), "fatalities"], na.rm = TRUE) }
      # add results
      x1 <- c(x1, x3)
    }
    # add results
    tsw[, "fatalities_2weeks_rate"] <- x1 * 100000 / tsw[, "pop"] 


  #...................................   
  ## Manage date on the presence of an ETC/transit centre
    # Whether an ETC or TC was open within the HZ at the given time
    tsw <- tsw[order(tsw[, "hz"], tsw[, "epi_year"], tsw[, "epi_week"]), ]
    hz <- sort(hz)
    x1 <- c()
    for (i in hz) {
      # select ETC open column for health zone time series
      x2 <- subset(tsw, hz == i)[, "open"]
      # base is that an ETC/TC was not open
      x3 <- rep("no", length(x2))
      # if an ETC/TC was ever open, set value from that time point to the end of the time series as 'yes'
      if ("yes" %in% x2) {x3[match("yes", x2):length(x3)] <- "yes"}
      # add results
      x1 <- c(x1, x3)
    }
    # add results
    tsw[, "etc_open"] <- x1
     
        
  #...................................   
  ## Convert other explanatory variables to population rates per 100,000      
    # Number of health facilities
    tsw[, "n_hf_rate"] <- tsw[, "n_hf"] * 100000 / tsw[, "pop"]  

    # Road coverage
    tsw[, "road_length_rate"] <- tsw[, "road_length"] * 100000 / tsw[, "pop"]  

    # Number of mines
    tsw[, "n_mines_rate"] <- tsw[, "n_mines"] * 100000 / tsw[, "pop"]  

        
#.........................................................................................                            
### Merging SDB and explanatory datasets, creating categorical variables
#.........................................................................................
    
  #...................................   
  ## Merge SDB with explanatory time series dataset
  sdb <- merge(sdb_dataset[, ! colnames(sdb_dataset) %in% c("month", "year")], 
                 tsw, by = c("hz", "epi_year", "epi_week"), all.x = TRUE )

  
  #...................................   
  ## Explore distributions and create categorical exposure variables     
    # Function to plot histograms
    f_hist <- function(f_var, f_data) {
      # start plot
      pl <- ggplot(f_data)
      
      # if the variable has >= 20 unique values...
      if (length(unique(na.omit(f_data[, f_var]))) >= 20) {
        pl <- pl + geom_histogram(aes(x = f_data[, f_var]) ) + theme_bw() + 
              # xlim(NA, quantile(f_data[, f_var], probs = 0.95, na.rm=TRUE)) +
              xlab(f_var) + expand_limits(x = 0)
      }
 
      # otherwise...
      if (length(unique(na.omit(f_data[, f_var]))) < 20) {
        pl <- pl +geom_histogram(aes(x = f_data[, f_var]), stat="count") + theme_bw()+ xlab(f_var)
      }
        
      print(pl)  
    }
  
    # Sub-coordination hub
    table(sdb$hub)
    sdb[, "hub_cat"] <- sdb[, "hub"]
    sdb[, "hub_cat"] <- ifelse(sdb[, "hub"] %in% c("Biakato", "Mambasa", "Mangina", "Mwenga", "Tchomia"), "other", sdb[, "hub"])
    table(sdb$hub_cat)
    
    # Age of deceased
    f_hist("age", sdb)
    sdb[, "age_cat"] <- cut(sdb[, "age"], breaks = c(0, 0.999, 4.999, 17.999, 59.999, 200),
                              labels = c("0y", "1-4y", "5-17y", "18-59y", ">=60y"), include.lowest = TRUE)
    table(sdb$age_cat)

    # Number of networks
    f_hist("n_networks", sdb)
    sdb[, "n_networks_cat"] <- cut(sdb[, "n_networks"], breaks = c(0, 3.999, 100),
                                 labels = c("< 4 networks", "4 networks"), include.lowest = TRUE)
    table(sdb$n_networks_cat)
    
    # Number of radio frequencies
    f_hist("n_frequencies", sdb)
    sdb[, "n_frequencies_cat"] <- cut(sdb[, "n_frequencies"], breaks = c(0, 9.999, 19.999, 100),
                                 labels = c("< 10 frequencies", "10-19 frequencies", ">= 20 frequencies"), include.lowest = TRUE)
    table(sdb$n_frequencies_cat)
  
    # Rate of confirmed cases over the time window
    f_hist("confirmed_6weeks_rate", sdb)
    sdb[, "confirmed_6weeks_rate_cat1"] <- cut(sdb[, "confirmed_6weeks_rate"], breaks = c(0, 0.001, 19.999, 39.999, 100),
                                 labels = c("0", "0.1 to 19.9", "20.0 to 39.9", ">= 40.0"), include.lowest = TRUE)
    sdb[, "confirmed_6weeks_rate_cat2"] <- cut(sdb[, "confirmed_6weeks_rate"], breaks = c(0, 0.001, 100),
                                 labels = c("0", "> 0"), include.lowest = TRUE)
    table(sdb$confirmed_6weeks_rate_cat1)
    table(sdb$confirmed_6weeks_rate_cat2)

    # Rate of insecurity events over the time window
    f_hist("n_events_2weeks_rate", sdb)
    sdb[, "n_events_2weeks_rate_cat"] <- cut(sdb[, "n_events_2weeks_rate"], breaks = c(0, 0.001, 4.999, 100),
                                 labels = c("0", "0.1 to 4.9", ">= 5.0"), include.lowest = TRUE)
    table(sdb$n_events_2weeks_rate_cat)

    # Rate of insecurity deaths over the time window
    f_hist("fatalities_2weeks_rate", sdb)
    sdb[, "fatalities_2weeks_rate_cat1"] <- cut(sdb[, "fatalities_2weeks_rate"], breaks = c(0, 0.001, 9.999, 100),
                                 labels = c("0", "0.1 to 9.9", ">= 10.0"), include.lowest = TRUE)
    sdb[, "fatalities_2weeks_rate_cat2"] <- cut(sdb[, "fatalities_2weeks_rate"], breaks = c(0, 0.001, 100),
                                 labels = c("0", "> 0"), include.lowest = TRUE)
    table(sdb$fatalities_2weeks_rate_cat1)
    table(sdb$fatalities_2weeks_rate_cat2)

    # Number of health facilities per population
    f_hist("n_hf_rate", sdb)
    sdb[, "n_hf_rate_cat"] <- cut(sdb[, "n_hf_rate"], breaks = c(0, 24.999, 49.999, 100),
                                 labels = c("< 25.0", "25.0 to 49.9", ">= 50.0"), include.lowest = TRUE)
    table(sdb$n_hf_rate_cat)

    # Road length per population
    f_hist("road_length_rate", sdb)
    sdb[, "road_length_rate_cat"] <- cut(sdb[, "road_length_rate"], breaks = c(0, 199.99, 399.99, 1000),
                                 labels = c("< 200 Km", "200-399 Km", ">= 400 Km"), include.lowest = TRUE)
    table(sdb$road_length_rate_cat)
    
    # Number of mines per population
    f_hist("n_mines_rate", sdb)
    sdb[, "n_mines_rate_cat"] <- cut(sdb[, "n_mines_rate"], breaks = c(0, 0.001, 1000),
                                 labels = c("no mining", "some mining"), include.lowest = TRUE)
    table(sdb$n_mines_rate_cat)
  
        
  #...................................   
  ## Create outcome categories (failure = 1; success = 0)
  table(sdb$outcome_lshtm)
  sdb[, "outcome_cat"] <- NA
  sdb[, "outcome_cat"] <- ifelse(sdb[, "outcome_lshtm"] %in% c("sdb not needed", "success") , 0, sdb[, "outcome_cat"] )        
  sdb[, "outcome_cat"] <- ifelse(sdb[, "outcome_lshtm"] %in% c("failure") , 1, sdb[, "outcome_cat"] )        
  table(sdb$outcome_cat)  
        

#.........................................................................................      
### Doing descriptive analysis
#.........................................................................................    

  #...................................   
  ## Profile of decedents in SDB database (gender, age, location of death, hub), by EVD status
    # Variables to tabulate by EVD status
    x1 <- c("gender", "age_cat", "origin_cat", "hub_cat")
    
    # Tabulate totals and percentages
    x2 <- c()
    x3 <-c()
    for (i in x1) {
      x2 <- rbind(x2, table(sdb[, i], sdb[, "evd_status"], useNA = "ifany") )
      x3 <- rbind(x3, prop.table(table(sdb[, i], sdb[, "evd_status"], useNA = "ifany"), margin = 2) * 100 )
    }
    x2 <- cbind(x2, x3)
    x2 <- rbind(x2, c(table(sdb[, "evd_status"], useNA = "ifany"), 
                      prop.table(table(sdb[, "evd_status"], useNA = "ifany")) * 100 ) )
    x2[, 4:6] <- round(x2[, 4:6], 1)
    
    # Create table
    out <- c()
    for (i in 1:3) {
      out <- cbind(out, paste(x2[, i], " (", x2[, i+3], ")" , sep="") )
    }
    out <-cbind(rownames(x2), out)
    out <- as.data.frame(out)
    colnames(out) <- c("category", colnames(x2[, 1:3]) )
    out[, "category"] <- ifelse(is.na(out[, "category"]), "unknown / missing", out[, "category"] )
    out[nrow(out), "category"] <- "totals"
    
    # Save table
    write.csv(out, "table_char_by_evd_status.csv", row.names = FALSE)
    
    
  #...................................   
  ## Trends by week in SDBs responded to, by province and team type, superimposed onto trends in EVD confirmed cases    
    
    # Aggregate SDB data
      # select SDBs that were actually responded to and were not false alerts
      x1 <- subset(sdb, responded == "yes")
      # aggregate by province, team type and week
      x1[, "n_sdb"] <- 1
      x1 <- aggregate(x1[ , "n_sdb"], by = x1[, c("province", "team_type", "epi_year", "epi_week")], FUN = sum )
      colnames(x1) <- c("province", "team_type", "epi_year", "epi_week", "n_sdb")
      # relabel and factorise team type
      x1[, "team_type"] <- ifelse(x1[, "team_type"] == "CEHRB", "Red Cross-supported CEHRB", x1[, "team_type"])
      x1[, "team_type"] <- ifelse(x1[, "team_type"] == "IFRC", "Red Cross mobile team", x1[, "team_type"])
      x1[, "team_type"] <- ifelse(x1[, "team_type"] == "Civil protection", "Civil Protection (mobile or CEHRB)", x1[, "team_type"])
      x1[, "team_type"] <- factor(x1[, "team_type"], 
        levels = c("Red Cross mobile team", "Red Cross-supported CEHRB", "Civil Protection (mobile or CEHRB)"))
      
    # Prepare EVD case data
     # aggregate by province
      x2 <- merge(evd_cases_weeks, admin_units, by="hz", all.x = TRUE)
      x2 <- aggregate(x2[ , "confirmed"], by = x2[, c("province", "epi_year", "epi_week")], FUN = sum )
      colnames(x2) <- c("province", "epi_year", "epi_week", "confirmed")
        
    # Preparatory steps for plotting  
     # create a date variable for the x axis
     x1[, "date"] <- as.Date(paste(x1[, "epi_year"], x1[, "epi_week"], 7), format = "%Y %U %u")
     x2[, "date"] <- as.Date(paste(x2[, "epi_year"], x2[, "epi_week"], 7), format = "%Y %U %u")      

     # Eliminate South Kivu province (only 1 observation)
     x1 <- subset(x1, province != "South Kivu")
     x2 <- subset(x2, province != "South Kivu")
     
    # Draw plot
    plot <- ggplot(x1) +
      geom_bar(mapping = aes(fill = team_type, x = date, y = n_sdb), position="stack", 
        stat="identity", alpha = 0.5) +
      geom_step(data = x2, mapping = aes(x = date, y = confirmed ), colour = palette_cb[6], size = 1) +
      scale_fill_manual(values = palette_cb[c(7, 2, 4)]) +
      scale_y_continuous("SDB alerts responded to", sec.axis = sec_axis(~ ., name="confirmed EVD cases")) +
      theme_bw() +
      facet_wrap(~province, nrow=2) +
      theme(legend.position="bottom", legend.direction="horizontal") +
      scale_x_date("", expand=c(0,0) , minor_breaks=NULL, date_breaks="2 months",
        date_labels = "%b-%Y", limits = c(as.Date("2018-07-01"), max(x1$date, na.rm=TRUE))  ) +
      labs(fill = "Responder:  ") +
      theme(legend.title = element_text(color="grey20", size=11),
        strip.text.x = element_text(color="grey20", size=11),
        legend.text = element_text(color="grey20", size=11),
        axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=11),               
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_line(color = "grey20"),
        axis.text.y = element_text(color = "grey20", size=11),
        axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
        axis.line.y.right = element_line(color = "grey20"),
        axis.ticks.y.right = element_line(color = "grey20"),
        axis.text.y.right = element_text(color = "grey20", size=11),
        axis.title.y.right = element_text(color="grey20", margin = margin(l = 10), size=11 )
      )
    plot
    ggsave("figure_1.png", units = "cm", height = 25, width = 22, dpi = "print", device = "png")
   
  #...................................   
  ## Trends in SDB percent success over time, by province, superimposed onto trends in SDB delay (% with delay <24h)
    
    # Prepare data
      # identify successes and failures, acceptable delays
      x1 <- sdb
      x1[, "n_sdb"] <- ifelse(is.na(x1$outcome_cat), 0, 1)
      x1[, "n_success"] <- ifelse(x1$outcome_cat == 0, 1, 0)
      x1[, "n_success"] <- ifelse(is.na(x1$n_success), 0, x1$n_success)
      x1[, "n_delay"] <- ifelse(is.na(x1$delay_sdb), 0, 1)
      x1[, "n_delay_ok"] <- ifelse(x1$delay_sdb == "<= 24h", 1, 0)
      x1[, "n_delay_ok"] <- ifelse(is.na(x1$delay_sdb), 0, x1$n_delay_ok)
      
      # aggregate by province and week
      x1 <- aggregate(x1[ , c("n_success", "n_sdb", "n_delay_ok", "n_delay")],
        by = x1[, c("province", "epi_year", "epi_week")], FUN = sum )
      colnames(x1) <- c("province", "epi_year", "epi_week", "n_success", "n_sdb", "n_delay_ok", "n_delay")
      
      # calculate percent indicators
      x1[, "percent_success"] <- (x1$n_success / x1$n_sdb) * 100
      x1[, "percent_delay_ok"] <- (x1$n_delay_ok / x1$n_delay) * 100

    # Preparatory steps for plotting  
     # create a date variable for the x axis
     x1[, "date"] <- as.Date(paste(x1[, "epi_year"], x1[, "epi_week"], 7), format = "%Y %U %u")
     
     # eliminate South Kivu province (only 1 observation)
     x1 <- subset(x1, province != "South Kivu")
     
     # reshape database so that there is one line for each time point and indicator
     x2 <- x1[, c("province", "date", "percent_success")]
     colnames(x2) <- c("province", "date", "percent")
     x3 <- x1[, c("province", "date", "percent_delay_ok")]
     colnames(x3) <- c("province", "date", "percent")     
     x1 <- rbind(x2, x3)
     x1[, "variable"] <- c(rep("success", nrow(x2)), rep("delay_ok", nrow(x3)))
     
    # Draw plot
    plot <- ggplot(x1) +
      # geom_point(aes(x = date, y = percent, colour = variable), alpha = 0.8, size = 3) +
      geom_step(mapping = aes(x = date, y = percent, colour = variable, linetype = variable ), 
        size = 2, alpha = 0.6) +
      scale_colour_manual(name = "Outcome:  ", labels = c("delay < 24 h  ", "successful SDB"), values = palette_cb[c(6, 7)]) +
      scale_linetype_manual(name = "Outcome:  ", labels = c("delay < 24 h  ", "successful SDB"), values=c("solid", "solid")) +
      scale_y_continuous("percentage of all SDB responses with non-missing data") +
      theme_bw() +
      facet_wrap(~province, nrow=2) +
      theme(legend.position="bottom", legend.direction="horizontal") +
      scale_x_date("", expand=c(0,0) , minor_breaks=NULL, date_breaks="2 months",
                   date_labels = "%b-%Y", limits = c(as.Date("2018-07-01"), max(x1$date, na.rm=TRUE))  ) +
      geom_hline(yintercept = 80, colour = palette_cb[4], linetype = "21", alpha = 0.7, size = 1) +
      annotate(geom = "text", x = as.Date("2018-08-01"), y = 85, colour = palette_cb[4], label = "target") +
      theme(legend.title = element_text(color="grey20", size=11),
        strip.text.x = element_text(color="grey20", size=11),
        legend.text = element_text(color="grey20", size=11),
        axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=11),               
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_line(color = "grey20"),
        axis.text.y = element_text(color = "grey20", size=11),
        axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 )
     )
    plot
    ggsave("figure_3.png", units = "cm", height = 20, width = 20, dpi = "print", device = "png")

      
#.........................................................................................      
### Doing univariate analysis (GLM, logistic, random effect = sub-coordination hub)
#.........................................................................................
      
  #...................................   
  ## Remove unnecessary variables
  x1 <- c("territory", "weeks", "months", "confirmed", "probable", "suspected", "alive", "dead", "unknown", "against_evd",
    "suspend_activities", "n_events", "fatalities", "open", "n_hf", "road_length", "n_mines", "pop", "hz", "epi_year",
    "epi_week", "month", "year", "hub", "age", "outcome_lshtm", "province", "n_networks", "n_frequencies",
    "confirmed_6weeks_rate", "n_events_2weeks_rate", "fatalities_2weeks_rate", "n_hf_rate", "road_length_rate",
    "n_mines_rate", "confirmed_2months_rate", "n_events_1month_rate", "fatalities_1month_rate",
    subset(dict, use == "descriptive")$variable )
  sdb <- sdb[, ! colnames(sdb) %in% x1]

  #...................................   
  ## Identify exposure variables
    
    # Data frame to denote which variables are exposures to include ("yes")
    exposures <- data.frame(colnames(sdb), rep("yes", length(colnames(sdb))) )
    colnames(exposures) <- c("variable", "include_exposure")  
    exposures[exposures[, "variable"] %in% c("date", "outcome_cat", "hub_cat"), "include_exposure"] <- "no"
    exposures[exposures[, "variable"] %in% c("fatalities_2weeks_rate_cat1", "confirmed_6weeks_rate_cat1"), "include_exposure"] <- "no"
    exposures <- subset(exposures, include_exposure == "yes")
    
    # Specify reference categories for categorical variables
    sdb <- as.data.frame(lapply(sdb, as.factor))
    exposures[, "ref"] <- c("IFRC", "M", "no", "ETC", "pre-epidemic", "no", "no", "yes", "18-59y", "< 4 networks",
      "< 10 frequencies", "0", "0", "0", ">= 50.0", ">= 400 Km", "no mining")
    for (i in 1:nrow(exposures) ) {
      x1 <- exposures[i, "variable"]
      sdb[, x1] <- relevel(sdb[, x1], exposures[i, "ref"])
    }
    
    # Specify distal, intermediate and proximate exposure variables
    exposures[, "level"] <- c("proximate", "distal", "proximate", "proximate", "intermediate", "intermediate",
      "intermediate", "intermediate", "distal", "intermediate", "intermediate", "intermediate",
      "intermediate", "intermediate", "distal", "distal", "distal")    
  
  #...................................   
  ## Function to fit the random effect logistic model and display clean results
  f_model <- function(f_vars, f_data) {
    # write the model formula
    form <- as.formula( paste("outcome_cat", "~", paste(f_vars, collapse= " + "), "+ (1|hub_cat)", sep="")  )
    # fit GLM with random effect
    fit <- glmer(form, data = f_data, family="binomial" )
    # compute ORs in linear form
    f_out <- tidy(fit, conf.int=TRUE, exponentiate=TRUE, effects="fixed")
    return(f_out)
  }
      
  #...................................   
  ## Fit univariate models and screen exposures in/out
    # Define p-value threshold for screening in
    p_level <- 0.20
    
    # Variables to keep in or out
    exposures[, "univariate_pass"] <- "no"
    
    # Run models and store output
    out_univariate <- data.frame()
    
    for (i in exposures$variable) {
      print(paste("now running univariate model for variable ", i))
      # fit model
      out <- f_model(i, sdb)
      # if at least one p-value among all categories is less than threshold, keep variable
      if (min(out[-1, "p.value"]) < p_level) {exposures[exposures[, "variable"] == i, "univariate_pass"] <- "yes"}
      # store output
      out_univariate <-rbind(out_univariate, out)
    }
    
    out_univariate
    
    # Save output
    write.csv(out_univariate, "out_univariate.csv", row.names = FALSE)

        
#.........................................................................................      
### Doing multivariate analysis (GLM, logistic, random effect = sub-coordination hub)
#.........................................................................................
    
  #...................................   
  ## Distal model

    # Select exposure variables
    x1 <- subset(exposures, level == "distal" & univariate_pass == "yes")[, "variable"]
    
    # Fit model
    out <- f_model(x1, sdb)
    print(out, n = Inf)

    # Remove variables and refit iteratively, looking for most parsimonious model
    x1 <- x1[x1 != "n_mines_rate_cat"]
    out <- f_model(x1, sdb)
    print(out, n = Inf)

  #...................................   
  ## Distal + intermediate model
    
    # Select exposure variables
    x1 <- c(x1, subset(exposures, level == "intermediate" & univariate_pass == "yes")[, "variable"] )
    
    # Fit model
    out <- f_model(x1, sdb)
    print(out, n = Inf)
      
    # Remove variables and refit iteratively, looking for most parsimonious model
    x1 <- x1[x1 != "n_events_2weeks_rate_cat"]
    out <- f_model(x1, sdb)
    print(out, n = Inf)
      
  #...................................   
  ## Distal + intermediate + proximate model
    
    # Select exposure variables
    x1 <- c(x1, subset(exposures, level == "proximate" & univariate_pass == "yes")[, "variable"] )
    
    # Fit model
    out <- f_model(x1, sdb)
    print(out, n = Inf)
    
    # Remove variables and refit iteratively, looking for most parsimonious model
    x1 <- x1[x1 != "road_length_rate_cat"]
    out <- f_model(x1, sdb)
    print(out, n = Inf)  
    
    x1 <- x1[x1 != "road_length_rate_cat"]
    out <- f_model(x1, sdb)
    print(out, n = Inf) 
    
    x1 <- x1[x1 != "age_cat"]
    out <- f_model(x1, sdb)
    print(out, n = Inf)     
  
  #...................................   
  ## Save final model
  write.csv(out, "out_multivariate.csv", row.names = FALSE)
    
        
#.........................................................................................
### ENDS
#.........................................................................................


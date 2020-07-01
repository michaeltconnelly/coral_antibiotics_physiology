volume <- 0.00065 * 0.40

temperature <- 24.3
salinity <- 35

exp_wells <- c("A2", "A4", "A5", "A6", "B2", "B5", "C1", "C2", "C4", "D1", "D3", "D5")
ctl_wells <- c("A1", "A3", "B1", "B3", "B4", "B6", "C3", "C5", "C6", "D2", "D4", "D6")



# Import Data -------------------------------------------------------------

path <- "~/Documents/Github repositories/respR Example Analyses/Mike Connelly/pocillopora_368_connelly_032619_b_Oxygen.xlsx"

## import
coral <- parse_multiplate(path)

## Convert well names to column index. 
data_cols <- which(names(coral) %in% exp_wells)
ctrl_cols <- which(names(coral) %in% ctl_wells)


# Background --------------------------------------------------------------

bg <- calc_rate.bg(coral, xcol = 2, ycol = ctrl_cols)

# Inspect and save data ---------------------------------------------------

## empty list to save results to
coral_insp <- list()

for (i in 1:length(data_cols)){
  coral_insp[[i]] <- inspect(coral,
                             time = 2,
                             oxygen = data_cols[i])
}


# Calculate rates ---------------------------------------------------------

coral_rates_cr <- list()

for(i in 1: length(data_cols)){
  
  coral_rates_cr[[i]] <- calc_rate(coral_insp[[i]],
                                   from = 500,
                                   to = 1500,
                                   by = "time")
}

# Correct for background --------------------------------------------------

coral_rates_cr_adj <- list()

for(i in 1:length(data_cols)){
  
  coral_rates_cr_adj[[i]] <- adjust_rate(coral_rates_cr[[i]], by = bg)
  
}


# Convert -----------------------------------------------------------------

coral_rates_cr_adj_conv <- list()

for(i in 1:length(data_cols)){
  
  coral_rates_cr_adj_conv[[i]] <- convert_rate(coral_rates_cr_adj[[i]], 
                                               o2.unit = "%",
                                               time.unit = "min",
                                               volume = volume,
                                               output.unit = "mg/h",
                                               t = temperature, # set at start
                                               S = salinity)
}


# extract all
final_rates <- sapply(coral_rates_cr_adj_conv, function(x) x$output)



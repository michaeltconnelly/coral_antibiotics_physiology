## Load required packages
library(respR)
library(readxl)


## Create function to import these types of file
## (This will be in the next version, out in a week or so)
parse_multiplate <- function(path) {
  out <- readxl::read_xlsx(path, skip = 10) # skip top 10 rows
  return(out)
}


# Experiment parameters ---------------------------------------------------

# These could all be in a single file of course

## Volume - for conversion later this should be in L
## Corrected for the approx. 40% you said the corals take up
## I'll show you later how to calculate this more precisely
volume <- 0.00065 * 0.40

temperature <- 24.3
salinity <- 35

## these are not quite what you told me, but i worked it out.
exp_wells <- c("A2", "A4", "A5", "A6", "B2", "B5", "C1", "C2", "C4", "D1", "D3", "D5")
ctl_wells <- c("A1", "A3", "B1", "B3", "B4", "B6", "C3", "C5", "C6", "D2", "D4", "D6")



# Import Data -------------------------------------------------------------

## If you have lots of files you only need to change this path and the code
## should work as long as everything else is the same

## Set file path
path <- "~/Documents/Github repositories/respR Example Analyses/Mike Connelly/pocillopora_368_connelly_032619_b_Oxygen.xlsx"

## import
coral <- parse_multiplate(path)

## Convert well names to column index. This just converts the data and ctrl
## columns to an index (i.e. column numbers) which makes selecting the correct
## columns easier
data_cols <- which(names(coral) %in% exp_wells)
ctrl_cols <- which(names(coral) %in% ctl_wells)

## quick check there are no duplicates - should be FALSE
any(data_cols %in% ctrl_cols)



# Background --------------------------------------------------------------

# The calc_rate.bg function is different to the others in that it will accept
# several columns of oxygen, and return an average rate for all of them (as well
# as individually). Using the average may or may bnot be what you want depending
# on your experiment, but I'm going to use it here.

# save bg
bg <- calc_rate.bg(coral, xcol = 2, ycol = ctrl_cols)

# view
print(bg)
plot(bg)

# Usually bg rates are negative, that is the oxygen decreases because
# microorganisms are using it. Here, it's increasing....  I don't know exactly
# what's going on in your experiment. I presume you don't have anything in there
# that will produce oxygen, so it suggests the wells might not be sealed
# correctly. If you are using water that is not 100% saturated and there is a
# leak it will slowly absorb oxygen and % will increase. Either that or the
# probes are drifting for some reason.

## Anyway, you can see the rates are pretty consistent, so presumably the same
## thing is happening in the experimental wells, so we can still use it to
## correct them. But i would run some trials to try and investigate what is
## going on. If this is going on while your corals are respiring you will be
## underestimating their uptake rate.

## Anyway, we have saved the background rate to the `bg` object (it contains
## both the 12 individual bg rates and a mean) and will use that later to
## correct.


# Inspect and save data ---------------------------------------------------

## this will inspect all data columns
inspect(coral, time = 2, oxygen = data_cols)

## However it doesn't let you look at them individually (depending on respR
## version you are using - that will change in next version), so I'm going to
## inspect and save them as separate objects and save them to a list

## empty list to save results to
coral_insp <- list()

for (i in 1:length(data_cols)){
  coral_insp[[i]] <- inspect(coral,
                             time = 2,
                             oxygen = data_cols[i])
}

## You will see it has plotted each one - too fast for you to see, but you can
## scroll back through plots to view each

## Or use these to look at individual ones
## (double square brackets are for referencing list() elements)
print(coral_insp[[1]])
plot(coral_insp[[1]])

## They all look ok (though see comments at very end). Don't worry about the
## unevenly spaced time warning - it's because you are using decimal minutes.



# Calculate rates ---------------------------------------------------------

## There's a number of different ways you can do this. Depends on your
## experiemnt and what you are looking for. I'll show you a couple.

## First, using calc_rate to use the same time period from each. I picked a time
## after things have appeared to have stabilised a little. 

## calc_rate

coral_rates_cr <- list()

for(i in 1: length(data_cols)){

  coral_rates_cr[[i]] <- calc_rate(coral_insp[[i]],
                                   from = 500,
                                   to = 1500,
                                   by = "time")
  }

## Again, it will go too fast to see all the plots, but you can scroll back through
## them, or use these on the list the calc_rate objects are saved to.

print(coral_rates_cr[[1]])
plot(coral_rates_cr[[1]])

## You can view or extract all rates using this

sapply(coral_rates_cr, function(x) x$rate)

## Seem fairly consistent!



## auto_rate

## Secondly using auto_rate which finds the most linear rate in each (see
## documentation and vignette).

coral_rates_ar <- list()

for(i in 1:length(data_cols)){

  coral_rates_ar[[i]] <- auto_rate(coral_insp[[i]])

}

## View all
## Need to add the [1] because auto_rate may identify several linear
## sections, and [1] is the top ranked one (i.e. most linear)
sapply(coral_rates_ar, function(x) x$rate[1])


## Again, you can scroll back through the plots to look at each, or view
## particular ones with print() and plot().

## Which of these methods you use is up to you. It depends on the experiment and
## what you are looking at. For yours, i would be careful with auto_rate. In
## most of your data, you have a clear section up to 500s that looks different.
## My guess is the corals are retracted or less active after the disturbance of
## setting up and starting the experiment, so you don't really want to use these
## sections.

## So there's a good chance auto_rate may select linear regions within these
## regions that you don't want to use. Going back through the auto_rate results,
## you'll see this has happened for at least one -

plot(coral_rates_ar[[10]])

## Not really representative of the experiment!

## However, you can go to the next ranked result using (pos =) to see if
## it fits better. This one is more consistent with the others.
plot(coral_rates_ar[[10]], pos = 2)
print(coral_rates_ar[[10]], pos = 2)

## Another thing you could do is subset out only the regions of data you are
## interested in and run auto_rate on those. 

coral_rates_ar_sub <- list()

for(i in 1:length(data_cols)){

  subs <- subset_data(coral_insp[[i]], from = 500, to = 1772, by = "time")
  coral_rates_ar_sub[[i]] <- auto_rate(subs)

}

## these look a lot better, and rates more consistent
sapply(coral_rates_ar_sub, function(x) x$rate[1])



# Correct for background --------------------------------------------------

## I'll use the calc_rate results and correct them. Remember there are 12
## results saved in a list.

## adjust first one - not saved
adjust_rate(coral_rates_cr[[1]], by = bg)

## This just shows you what's happening. Because your bg rate is positive the
## specimen rate actually gets corrected upwards. i.e. get bigger (more
## negative). By default this uses the mean bg rate of all 12 controls, but you
## can enter a different one or any value.
adjust_rate(coral_rates_cr[[1]], by = bg$bgrate[1])
adjust_rate(coral_rates_cr[[1]], by = 0.0005)

## Can adjust them all like this

coral_rates_cr_adj <- list()

for(i in 1:length(data_cols)){
  
  coral_rates_cr_adj[[i]] <- adjust_rate(coral_rates_cr[[i]], by = bg)

}

## new adjusted rates - make sure you extract the $corrected rate
sapply(coral_rates_cr_adj, function(x) x$corrected)


# Convert -----------------------------------------------------------------

## Ok, last thing is to convert the rates to units. 

## One to show you what's happening
convert_rate(coral_rates_cr_adj[[1]],
             o2.unit = "%",
             time.unit = "min",
             volume = volume,
             output.unit = "mg/h",
             t = temperature, # set at start
             S = salinity) # set at start


## All
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


# view one
print(coral_rates_cr_adj_conv[[4]])

# extract all
final_rates <- sapply(coral_rates_cr_adj_conv, function(x) x$output)

## these are your final rates in mg/h
final_rates



# Comments ----------------------------------------------------------------

## I know this is just a pilot experiment, but one thing with your data that
## stands out is a lot of the curves are getting steeper as the experiment goes
## on. That is they haven't levelled out to a consistent slope. You can see this
## in the bottom panel of the inspect() plots - these show the rate over the
## experiment and in nearly all it is constantly rising. It should level off at
## some point if you are looking for a standard or routine metabolic rate.

## This suggests to me your experiments are not long enough to reach a steady
## state. In other words the corals are using up the available O2 too quickly.
## At 25 mins or so, these experiments are really brief! I know you are in
## tropical temperatures, but my own experiments on inverts generally last
## several hours.  The only thing you can do if this is the case (given your
## respirometry system) is try to reduce the size of your coral fragments. That
## means they will take up less proportional volume. This is also good thing -
## generally you want a much lower specimen:chamber volume ratio than what you
## have at 40%. It's slightly different for mobile inverts, but i generally go
## for the specimen taking up around 10% of the volume, or at most 20% or so. It
## just gives the animal more time to become accustomed to the chamber, handling
## etc. and so revert eventually to routine behaviour. Alternatively you could
## have them in the chambers for a longer time but open within a tank of
## water, and only seal them at the last minute. 

## Basically you need to do a lot of trials. What you want is an inspect() plot
## where the rolling rate eventually flattens out.  

## For converting mass and density of the corals to volume, and calculating the
## effective volume see my other R package
## https://github.com/nicholascarey/respfun


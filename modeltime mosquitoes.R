##########################################
### Modeltime for mosquito control ###
##########################################

# install packages
library(tidymodels)
library(modeltime)
library(timetk)   
library(lubridate)
library(tidyverse)
library(parsnip)
library(ggplot2)

library(prophet)
library(randomForest)

#################################################################################
#
#   NOTES: 
#   - run "FS trap data merge and summarize v20" for exporting FieldSeeker data
#   - see modeltime tutorial for detailed information:
#   - https://www.business-science.io/code-tools/2020/06/29/introducing-modeltime.html?utm_content=bufferd20d1&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer#
#
#################################################################################

#### Load required packages ###

library(ggplot2)
library(reshape2)
library(lubridate)
library(zoo)
library(dplyr)
library(tidyr)


#### Pull coordinate data from local hard drive ###

dat.coord <- read.csv("H:/Sample_scripts_(Becky)/TrapLocation_3.csv")  #temp home location


#### Import data from ArcGIS ###

# install.packages("arcgisbinding", repos = "https://r.esri.com", type = "win.binary")  #if misbehaving
library(arcgisbinding)
arc.check_product()  #installation check
arc.check_portal()   #portal check


## Trap-related layers

TrapData <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/15")  #pulls from server
dat.tr <- arc.select(object = TrapData, fields = names(TrapData@fields))  #puts all fields into data.frame
dat.tr <- dat.tr[-grep("fail", dat.tr$COMMENTS, ignore.case = TRUE),]     #removes records marked as failed
colnames(dat.tr)

SpeciesAbundance <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/10")
dat.sa <- arc.select(object = SpeciesAbundance, fields = names(SpeciesAbundance@fields))
dat.sa$SPECIES <- as.factor(dat.sa$SPECIES)
colnames(dat.sa)


TrapLocation <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/3")
dat.loc <- arc.select(object = TrapLocation, fields = names(TrapLocation@fields))
colnames(dat.loc)


## Pool-related layers (also download trap layers above)

Pool <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/9")
dat.pool <- arc.select(object = Pool, fields = names(Pool@fields))

PoolDetail <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/8")
dat.pdt <- arc.select(object = PoolDetail, fields = names(PoolDetail@fields))

PointLocation <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/11")
dat.ptloc <- arc.select(object = PointLocation, fields = names(PointLocation@fields))


## Treatment- and inventory-related layers

MosquitoInspection <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/14")
dat.inspect <- arc.select(object = MosquitoInspection, fields = names(MosquitoInspection@fields))

Treatment <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/18")
dat.treat <- arc.select(object = Treatment, fields = names(Treatment@fields))

ProductInventoryTransactions <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/43")
dat.transact <- arc.select(object = ProductInventoryTransactions, fields = names(ProductInventoryTransactions@fields))

ProductInventory <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/42")
dat.prodinv <- arc.select(object = ProductInventory, fields = names(ProductInventory@fields))


## LRC layers

LandingCountLocation <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/27")
dat.lrcloc <- arc.select(object = LandingCountLocation, fields = names(LandingCountLocation@fields))

LandingCount <- arc.open(path = "https://services9.arcgis.com/pxQng7oePvh2pa6G/arcgis/rest/services/CollierFieldSeeker3/FeatureServer/16")
dat.lrc <- arc.select(object = LandingCount, fields = names(LandingCount@fields))

# ## TEMPORARY: Test data subset
# dat.tr <- subset(dat.tr, as.Date(as.POSIXct(ENDDATETIME)) == as.Date("2020-11-17"))
# dat.sa <- subset(dat.sa, as.Date(as.POSIXct(created_date)) == as.Date("2020-11-18"))



#### Functions ###

full_merge <- function(date.start = "2018-01-01", date.end = Sys.time()){  
  d.start <- as.POSIXct(date.start)
  d.end <- as.POSIXct(date.end)
  retrievals <- subset(dat.tr,
                       TRAPACTIVITYTYPE == "R" & ENDDATETIME > d.start & ENDDATETIME < d.end)  #all retrievals in date range
  species.sub <- subset(dat.sa,
                        select = c("MALES", "FEMALES", "SPECIES", "TRAPDATA_ID", "GlobalID"))
  colnames(species.sub)[colnames(species.sub) == "GlobalID"] <- "GlobalID_Species"
  merge1 <- merge(retrievals, species.sub,
                  by.x = "GlobalID",
                  by.y = "TRAPDATA_ID",
                  all.x = TRUE)
  merge1 <- subset(merge1,
                   select = c("LOC_ID", "ENDDATETIME", "COMMENTS", "LOCATIONNAME", "FEMALES",
                              "MALES", "SPECIES", "TRAPTYPE", "GlobalID", "GlobalID_Species"))
  colnames(merge1)[names(merge1) == "GlobalID"] <- "GlobalID_tr"  #retain TrapData GlobalID
  merge2 <- merge(merge1, dat.loc,
                  by.x = "LOC_ID",
                  by.y = "GlobalID",
                  all.x = TRUE)
  merge2 <- subset(merge2,
                   select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "NAME", "SPECIES",
                              "FEMALES", "MALES", "COMMENTS.x", "LOCATIONNAME", "GlobalID_tr",
                              "HABITAT", "ACTIVE", "LOC_ID", "GlobalID_Species"))  #select final columns and reorder for output
  colnames(merge2)[names(merge2) == "COMMENTS.x"] <- "COMMENTS_TrapData"
  merge2$FEMALES[is.na(merge2$FEMALES) == TRUE] <- 0
  merge2$ENDDATETIME <- as.POSIXct(merge2$ENDDATETIME)  #formatting important for trap_totals
  return(merge2)
}  #merges TrapLocation, TrapData and SpeciesAbundance tables, long format, no added zeroes

trap_totals <- function(trap.start = "2018-01-01", trap.end = Sys.time()){  #optional date filter; defaults to full FS dataset
  merged.traps <- full_merge(trap.start, trap.end)
  output <- merged.traps %>%
    group_by(ENDDATETIME, TRAPTYPE, ZONE, NAME, COMMENTS_TrapData, 
             LOCATIONNAME, GlobalID_tr, HABITAT, ACTIVE, LOC_ID) %>%  #keeps all columns from full_merge except SPECIES and FEMALES
    summarise(sum_females = sum(FEMALES),
              sum_males = sum(MALES))
  output <- as.data.frame(output)
  return(output)
}  #sum females of all species collected by site and date

species_counts <- function(sp.ref, trap.start = "2018-01-01", trap.end = Sys.time()){
  merged.traps <- full_merge(trap.start, trap.end)
  totals <- trap_totals(trap.start, trap.end)
  species.sub <- subset(merged.traps, 
                        SPECIES == sp.ref,
                        select = c("GlobalID_tr", "SPECIES", "FEMALES", "MALES"))
  output <- merge(totals, species.sub,  
                  by = "GlobalID_tr",
                  all.x = TRUE,
                  all.y = TRUE)    #using totals as "master" creates zero values for traps without target species
  output <- output[, c("ENDDATETIME", "TRAPTYPE", "ZONE", "NAME", "SPECIES",
                       "FEMALES", "MALES", "sum_females", "sum_males", "COMMENTS_TrapData",
                       "LOCATIONNAME", "GlobalID_tr", "HABITAT", "ACTIVE", "LOC_ID")]
  sel <- grepl("males", names(output), ignore.case = TRUE)  #lists "female" and "male" columns and totals
  output[sel] <- lapply(output[sel], function(x) replace(x, x %in% NA, 0))  #in sel columns, replace any NAs with 0
  output$SPECIES[which(is.na(output$SPECIES) == TRUE)] <- sp.ref  #fills in empty SPECIES values
  return(output)
}  #return selected species females and males along with fem/male totals (all species) by collection date and site

full_merge0 <- function(date.start = "2018-01-01", date.end = Sys.time()){  
  d.start <- as.POSIXct(date.start)
  d.end <- as.POSIXct(date.end)
  retrievals <- subset(dat.tr,
                       TRAPACTIVITYTYPE == "R" & ENDDATETIME > d.start & ENDDATETIME < d.end)  #all retrievals in date range
  retrievals <- subset(retrievals,
                       TRAPTYPE != "OV")  #removes ovicup collections
  species.sub <- subset(dat.sa,
                        SPECIES != "Eggs - Aedes aegypti",  #removes ovicup values
                        select = c("FEMALES", "SPECIES", "TRAPDATA_ID"))
  species.wide <- dcast(species.sub,
                        TRAPDATA_ID ~ SPECIES,
                        fun.aggregate = sum,    #shouldn't be necessary once data corrected?
                        value.var = "FEMALES")  #long to wide, adds zeroes
  species.wide$Total_females <- rowSums(species.wide[,-1])  #calculate total trap catch
  species.long <- melt(species.wide,
                       id.vars = "TRAPDATA_ID",
                       variable.name = "SPECIES",
                       value.name = "FEMALES")  #wide to long
  merge1 <- merge(retrievals, species.sub,
                  by.x = "GlobalID",
                  by.y = "TRAPDATA_ID",
                  all.x = TRUE)
  merge1 <- subset(merge1,
                   select = c("LOC_ID", "ENDDATETIME", "COMMENTS", "LOCATIONNAME", "FEMALES",
                              "SPECIES", "TRAPTYPE", "GlobalID"))
  colnames(merge1)[names(merge1) == "GlobalID"] <- "GlobalID_tr"  #retain TrapData GlobalID
  merge2 <- merge(merge1, dat.loc,
                  by.x = "LOC_ID",
                  by.y = "GlobalID",
                  all.x = TRUE)
  merge2 <- subset(merge2,
                   select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "NAME", "SPECIES",
                              "FEMALES", "COMMENTS.x", "LOCATIONNAME", "GlobalID_tr",
                              "HABITAT", "ACTIVE", "LOC_ID"))  #select final columns and reorder for output
  colnames(merge2)[names(merge2) == "COMMENTS.x"] <- "COMMENTS_TrapData"
  merge2$FEMALES[is.na(merge2$FEMALES) == TRUE] <- 0  #replace NAs with zeroes
  merge2$ENDDATETIME <- as.POSIXct(merge2$ENDDATETIME)  #formatting important for trap_totals
  return(merge2)
}  #merges trap-related tables (long format) and adds zeroes, keeps all species, no ovicups

wide_merge <- function(date.start = "2018-01-01", date.end = Sys.time()){  
  d.start <- as.POSIXct(date.start)
  d.end <- as.POSIXct(date.end)
  retrievals <- subset(dat.tr,
                       TRAPACTIVITYTYPE == "R" & ENDDATETIME > d.start & ENDDATETIME < d.end)  #all retrievals in date range
  retrievals <- subset(retrievals,
                       TRAPTYPE != "OV")  #removes ovicup collections
  colnames(retrievals)[names(retrievals) == "GlobalID"] <- "GlobalID_tr"  #to distinguish from dat.loc$GlobalID
  merge1 <- merge(retrievals, dat.loc,
                  by.x = "LOC_ID",
                  by.y = "GlobalID",
                  all.x = TRUE)
  merge1 <- subset(merge1,
                   select = c("ENDDATETIME", "TRAPTYPE", "ZONE.y", "NAME", "COMMENTS.x", "LOCATIONNAME", "GlobalID_tr",
                              "HABITAT", "ACTIVE", "LOC_ID"))  #select final columns and reorder for output
  colnames(merge1)[names(merge1) == "COMMENTS.x"] <- "COMMENTS_TrapData"
  colnames(merge1)[names(merge1) == "ZONE.y"] <- "ZONE"
  merge1$ENDDATETIME <- as.POSIXct(merge1$ENDDATETIME)  #formatting important for trap_totals
  species.sub <- subset(dat.sa,
                        SPECIES != "Eggs - Aedes aegypti",  #removes ovicup values
                        select = c("FEMALES", "SPECIES", "TRAPDATA_ID"))
  species.wide <- dcast(species.sub,
                        TRAPDATA_ID ~ SPECIES,
                        fun.aggregate = sum,    #shouldn't be necessary once data corrected?
                        value.var = "FEMALES")  #long to wide, adds zeroes
  species.wide$Total_females <- rowSums(species.wide[,-1])  #calculate total trap catch
  merge2 <- merge(merge1, species.wide,
                  by.x = "GlobalID_tr",
                  by.y = "TRAPDATA_ID",
                  all.x = TRUE)
  sel <- 10:length(colnames(merge2))  #all count columns
  merge2[sel] <- lapply(merge2[sel], function(x) replace(x, x %in% NA, 0))  #in sel columns, replace any NAs with 0
  return(merge2)
}  #merges tables (wide format) and adds zeroes, keeps all species, no ovicups

graph_prep <- function(df){
  x <- df  #data frame
  x$Date <- as.Date(x$ENDDATETIME)
  x$DateRef <- x$Date
  year(x$DateRef) <- 2000  #this was a leap year; use this to overplot multiple years
  x$Year <- year(x$Date)   #useful for color coding overplots
  x$Year <- factor(x$Year)
  # TO DO: add epi week column
  x$NAME <- factor(x$NAME)  #get rid of extraneous factor values following subsetting
  x$ZONE <- factor(x$ZONE)
  #  x$SPECIES <- factor(x$SPECIES)
  x$TRAPTYPE <- factor(x$TRAPTYPE)
  return(x)
}  #adds useful columns for graphing

genus_counts <- function(trap.start = "2018-01-01", trap.end = Sys.time()){
  merged.traps <- full_merge(trap.start, trap.end)
  merged.traps <- subset(merged.traps, select = c("GlobalID_tr", "SPECIES", "FEMALES", "MALES"))
  
  merged.traps <- cbind(merged.traps, colsplit(merged.traps$SPECIES, " ", c("Genus", "Species")))
  genus.abbrev <- c("Ae.", "Aed.", "An.", "Cq.", "Cs.", "Cx.", "Man.", "Ps.", "Ur.", "Wy.", "Aedes", "Aedeomyia", "Anopheles", "Coquillettidia", "Culiseta", "Culex", "Mansonia", "Psorophora", "Uranotaenia", "Wyeomyia", "Unknown")
  genus.full <- c(rep(c("Aedes", "Aedeomyia", "Anopheles", "Coquillettidia", "Culiseta", "Culex", "Mansonia", "Psorophora", "Uranotaenia", "Wyeomyia"), times = 2), "Unknown")
  map <- setNames(genus.full, genus.abbrev)
  merged.traps$Genus[] <- map[merged.traps$Genus]  #corrects genus names; NAs are empty traps
  
  genus.sums <- merged.traps %>%
    group_by(GlobalID_tr, Genus) %>%
    summarise(FEMALES = sum(FEMALES),
              MALES = sum(MALES))  #so that this will work with general plot code below
  
  totals <- trap_totals(trap.start, trap.end)
  
  output <- merge(totals, genus.sums,  
                  by = "GlobalID_tr",
                  all.x = TRUE)    #using totals as "master" creates zero values for traps without target species
  output <- output[, c("ENDDATETIME", "TRAPTYPE", "ZONE", "NAME", "Genus",
                       "FEMALES", "MALES", "sum_females", "sum_males", "COMMENTS_TrapData",
                       "LOCATIONNAME", "GlobalID_tr", "HABITAT", "ACTIVE", "LOC_ID")]
  sel <- grepl("males", names(output), ignore.case = TRUE)  #lists "female" and "male" columns and totals
  output[sel] <- lapply(output[sel], function(x) replace(x, x %in% NA, 0))  #in sel columns, replace any NAs with 0
  
  return(output)
}  #return counts by genus, long format, sum_females/males are total trap counts (for calculating proportions)

pool_merge <- function(trap.start = "2018-01-01", trap.end = Sys.time()){
  sub.pdt <- subset(dat.pdt, 
                    select = c("POOL_ID", "SPECIES", "FEMALES"))  #avoids duplication of colnames
  merge1 <- merge(sub.pdt, dat.pool,
                  by.x = "POOL_ID",
                  by.y = "GlobalID",  #PoolDetail is sub to Pool and TrapData
                  all.x = TRUE)       #data.frame of PoolDetail and Pool data
  merge1 <- subset(merge1,
                   select = c("SPECIES", "FEMALES", "TRAPDATA_ID", "DATETESTED", "TESTTECH", 
                              "COMMENTS", "TESTMETHOD", "DISEASETESTED", "DISEASEPOS", "POOL_ID"))
  
  merge1b <- trap_totals(trap.start = trap.start, trap.end = trap.end)  #TrapData and TrapLocation, no species
  
  merge2 <- merge(merge1, merge1b,
                  by.x = "TRAPDATA_ID",
                  by.y = "GlobalID_tr")  #merge all together; all.x would have extraneous traps, all.y isn't subset by date
  return(merge2)
}  #return pool results with trap collection data, long format

join_field <- function(target, tfield, joinx, jfield){  #e.g. (dat.ptloc, "GlobalID", dat.inspect, "POINTLOCID")
  target <- target
  tfield <- tfield
  joinx <- joinx
  jfield <- jfield
  merge1 <- merge(target, joinx,
                  by.x = tfield,
                  by.y = jfield
  )
  return(merge1)
}  #mimics ArcGIS Join Field function

add_trapcoord <- function(df){
  temp.coords <- subset(dat.coord, select = c(GlobalID, Name, x, y))
  df$LOC_ID <- gsub("[{]", "", df$LOC_ID)
  df$LOC_ID <- gsub("[}]", "", df$LOC_ID)  #curly brackets interfere with merge below
  merge(df, temp.coords,
        by.x = "LOC_ID", by.y = "GlobalID",
        all.x = TRUE)
}  #adds coordinates to trap data based on LOC_ID


# #### Testing: once duplicate species entries for trap events removed, use to check species_counts function ###
# vec.species <- levels(dat.sa$SPECIES)
# output <- rep(NA, times = length(vec.species))
# for(i in 1:length(vec.species)){
#   x <- vec.species[i]
#   temp <- species_counts(x)
#   output[i] <- length(temp[,1])
# }
# (temp <- data.frame(vec.species, output))  #test for species_counts changes! (temporary troubleshooting)
# 
# temp <- rbind(species_counts("Ae. aegypti"),
#               species_counts("Ae. albopictus"),
#               species_counts("Cx. nigripalpus"))
# temp$Date <- as.Date(as.POSIXct(temp$ENDDATETIME))
# temp$DateLoc <- paste(temp$Date, temp$LOCATIONNAME)
# temp2 <- temp %>%
#   group_by(DateLoc, SPECIES) %>%
#   summarise(N = n())
# temp2 <- as.data.frame(temp2)
# temp3 <- temp2[which(temp2$N > 1), ]









#####################
#
# Collate the data
#
#####################


# Cx. nigripalpus from CDC traps

sp.cdc <- "Cx. nigripalpus"  #set target species for section
dat.cdc <- species_counts(sp.cdc)

dat.cdc$NAME <- gsub("CDC AM Del Webb 2", "CDC AM Del Webb", as.character(dat.cdc$NAME))  #fixes site crossover issue; coordinates not affected

CDC.list <- c("CDC AM Ave Entrance",
              "CDC AM Del Webb",
              "CDC AM Davey Lndscp Yrd",
              "CDC AM Dog Park",
              "CDC AM Dog Park 2", # whats up with this site?
              "CDC AM NORTH PRK PASTURE",
              "CDC AM Maple Ridge",
              "CDC AM Maple Ridge 11", # whats up with this site?
              "CDC AM North Park E",
              "CDC Oil Well Grade Rd",
              "CDC AM Trap 7", # whats up with this site?
              "CDC AM Trap 8", # whats up with this site?
              "CDC AM Water Park",
              "CDC CS Collier Seminole 1",
              "CDC CS Collier Seminole 2", 
              "CDC E Bayshore Dr",
              "CDC E Bayview",
              "CDC E Kings Lake",
              "CDC E Sports Park",
              "CDC E Wendy",
              "CDC GGC Whippoorwill",
              "CDC GGE 15th St NW",
              "CDC GGE 29th Ave SW",
              "CDC GGE 43rd and 40th",
              "CDC GGE 47th and 12th",
              "CDC GGE Richard St",
              "CDC GGE Shady Hollow",
              "CDC GGE Weber",
              "CDC GGE 66th Ave NE",
              "CDC GGE North Desoto",
              "CDC GGE South Desoto ",
              "CDC HC 1", 
              "CDC HC 2", 
              "CDC HC Fiddlers Creek",
              "CDC HC KOA",
              "CDC HC Lipman House",
              "CDC HC Lipman Office",
              "CDC HC Trevisio Bay Dr",
              "CDC IM Airport Ditch",
              "CDC IM Anez",
              "CDC IM Carson Rd",
              "CDC IM Hector Ramirez", # whats up with this site?
              "CDC IM Indian Camp Rd",
              "CDC IM Little League",
              "CDC IM Farm Village",
              "CDC IM LCEC",
              "CDC IM Water Plant",
              "CDC KW Keewaydin 1",
              "CDC KW Keewaydin S Ctrl",
              "CDC KW Windstar",
              "CDC M Stevens Landing",
              "CDC N Nicholas Blvd",
              "CDC N Veterans Park",
              "CDC NA Keewaydin 2", # note the wrong zone label?
              "CDC PY Canal Rd",
              "CDC PY Newman Dr",
              "CDC PY Picayune Trailhead",
              "CDC T Crayton")

# create the dataset
dat.cdc <- dat.cdc[dat.cdc$NAME %in% CDC.list,]



# calculate average female abundance per day
data <- dat.cdc %>%
  mutate(day = floor_date(ENDDATETIME, "day")) %>%
  group_by(day) %>%
  summarize(avg = mean(FEMALES))




######## TESTING ROLLING AVERAGES (this works!)

library(zoo)
data <- data %>%
  arrange(day) %>% 
  mutate(avg3 = zoo::rollmean(avg, k = 3, fill = NA))  %>%
  select(day,
         avg,
         avg3)
head(data1)




# visualize time series data set
data %>%
  plot_time_series(day, avg3, .interactive = FALSE)


# create a csv to evaluate the raw data
# write.csv(data, "Z:/Research Dept/Rob Straser/cx.nig.avgs")




#####################
#
# Train and test data
#
#####################

# split data to train and test model
# setting 'assess = 3 months" tells function to use last 3 months to data as a testing set
# setting 'cumulative = TRUE' tells the sampling to use all prior data as training set
splits <- data %>%
  time_series_split(assess = "12 months", cumulative = TRUE)

# visualize the train/test split
splits %>%
  tk_time_series_cv_plan() %>%
  plot_time_series_cv_plan(day, avg, .interactive = FALSE)






#####################
#
# Modeling data
#
#####################

#### Automatic models ####

# Auto ARIMA
model_fit_arima <- arima_reg() %>%
  set_engine("auto_arima") %>%
  fit(avg ~ day, training(splits))
model_fit_arima

# Prophet
model_fit_prophet <- prophet_reg(seasonality_yearly = TRUE) %>%
  set_engine("prophet") %>%
  fit(avg ~ day, training(splits))
model_fit_prophet





#### Machine learning models ####

# Preprocessing recipe

# Create a preprocessing recipe using recipe() and adding time series steps. 
# The process uses the "date" column to create 45 new features that I'd like to model. 
# These include time-series signature features and fourier series.

recipe_spec <- recipe(avg ~ day, training(splits)) %>%
  step_timeseries_signature(day) %>%
  step_rm(contains("am.pm"), contains("hour"), contains("minute"),
          contains("second"), contains("xts")) %>%
  step_fourier(day, period = 365, K = 5) %>%
  step_dummy(all_nominal())

recipe_spec %>% prep() %>% juice()



# Elastic net
model_spec_glmnet <- linear_reg(penalty = 0.01, mixture = 0.5) %>%
  set_engine("glmnet")

workflow_fit_glmnet <- workflow() %>%
  add_model(model_spec_glmnet) %>%
  add_recipe(recipe_spec %>% step_rm(day)) %>%
  fit(training(splits))



# Random forest
model_spec_rf <- rand_forest(trees = 500, min_n = 50) %>%
  set_engine("randomForest")

workflow_fit_rf <- workflow() %>%
  add_model(model_spec_rf) %>%
  add_recipe(recipe_spec %>% step_rm(day)) %>%
  fit(training(splits))




# Hybrid ML models
# hybrid models that combine both automated algorithms with machine learning

# Prophet boost (prophet & XGBoost)
model_spec_prophet_boost <- prophet_boost(seasonality_yearly = TRUE) %>%
  set_engine("prophet_xgboost") 

workflow_fit_prophet_boost <- workflow() %>%
  add_model(model_spec_prophet_boost) %>%
  add_recipe(recipe_spec) %>%
  fit(training(splits))
workflow_fit_prophet_boost






#####################
#
# Model comparison
#
#####################


# modeltime table organizes models with IDs and creates
# generic descriptions to help keep track of models.

model_table <- modeltime_table(
  model_fit_arima, 
  model_fit_prophet,
  workflow_fit_glmnet,
  workflow_fit_rf,
  workflow_fit_prophet_boost
) 

model_table


# calibration is used to quantify error and estimate confidence intervals.
# perform model calibration on the testing set
# wwo new columns are generated (".type" and ".calibration_data"), 
# the most important of which is the ".calibration_data". This includes the actual values, 
# fitted values, and residuals for the testing set.


calibration_table <- model_table %>%
  modeltime_calibrate(testing(splits))
calibration_table



# visualize forecast models
calibration_table %>%
  modeltime_forecast(actual_data = data) %>%
  plot_modeltime_forecast(.interactive = FALSE)


# test accuracy of testing set
# MAE - Mean absolute error, mae()
# MAPE - Mean absolute percentage error, mape()
# MASE - Mean absolute scaled error, mase()
# SMAPE - Symmetric mean absolute percentage error, smape()
# RMSE - Root mean squared error, rmse()
# RSQ - R-squared, rsq()

calibration_table %>%
  modeltime_accuracy() %>%
  table_modeltime_accuracy(.interactive = FALSE)



#####################
#
# Refit and forecasting
#
#####################

# Refit and Forecast Forward
calibration_table %>%
  modeltime_refit(data) %>%
  modeltime_forecast(h = "12 months", actual_data = data) %>%
  plot_modeltime_forecast(.interactive = FALSE)




























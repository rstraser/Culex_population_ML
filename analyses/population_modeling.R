##################################################
### Modeltime for mosquito population modeling ###
##################################################



### Load required packages ###
library(tidymodels)
library(modeltime)
library(timetk)   
library(lubridate)
library(tidyverse)
library(parsnip)
library(randomForest)
library(zoo)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(gridExtra)
library(leaflet)
library(forecast)
library(prophet)






#####################
#
# Pull and collate data
#
#####################


#### Pull coordinate and population data from local drive ####

dat.coord <- read.csv("coordinated_PROTECTED")
fs.coord <- read.csv("fieldseeker_PROTECTED")
geopro.coord <- read.csv("geopro_PROTECTED")
paper.coord <- read.csv("paper_PROTECTED")


# subset columns to join
sub.fs.coord <- subset(fs.coord,
                       select = c("Name", "x", "y"))
sub.geopro.coord <- subset(geopro.coord,
                           select = c("Name", "x", "y"))
sub.paper.coord <- subset(paper.coord,
                          select = c("Name", "x", "y"))

# join all coordinate data
dat.coord <- rbind(sub.paper.coord, sub.geopro.coord, sub.fs.coord)




#### Import data from ArcGIS ####

# install.packages("arcgisbinding", repos = "https://r.esri.com", type = "win.binary")  #run if error
library(arcgisbinding)
arc.check_product()  #installation check
arc.check_portal()   #portal check


# Trap-related layers
TrapData <- arc.open(path = "server_PROTECTED")  #pulls from server
dat.tr <- arc.select(object = TrapData, fields = names(TrapData@fields))  #puts all fields into data.frame
dat.tr <- dat.tr[-grep("fail", dat.tr$COMMENTS, ignore.case = TRUE),]     #removes records marked as failed
colnames(dat.tr)

SpeciesAbundance <- arc.open(path = "server_PROTECTED")
dat.sa <- arc.select(object = SpeciesAbundance, fields = names(SpeciesAbundance@fields))
dat.sa$SPECIES <- as.factor(dat.sa$SPECIES)
colnames(dat.sa)

TrapLocation <- arc.open(path = "server_PROTECTED")
dat.loc <- arc.select(object = TrapLocation, fields = names(TrapLocation@fields))
colnames(dat.loc)














#####################
#
# Helpful functions (Heinig)
#
#####################


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












#####################
#
# Clean and structure data
#
#####################

# Cx. nigripalpus from CDC traps
sp.cdc <- "Cx. nigripalpus"  #set target species for selection
dat.cdc <- species_counts(sp.cdc)


# join in historical data (GeoPro and Paper records)
paper <- read.csv("Historical data/2011-2014 FS compatible paper trap data.csv")
geopro <- read.csv("Historical data/2015-2019 GeoPro and paper consolidated data.csv")

# subset variables of interest
sub.paper <- subset(paper,
                    select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "SPECIES",
                               "FEMALES", "MALES", "LOCATIONNAME"))
sub.geopro <- subset(geopro,
                     select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "SPECIES",
                                "TOTAL", "MALES", "LOCATIONNAME"))

# no delineation between FEMALES and TOTAL in GeoPro data
sub.geopro <- sub.geopro %>%
  rename(FEMALES = TOTAL)

# subset joined data
sub.dat.cdc <- subset(dat.cdc,
                      select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "SPECIES",
                                 "FEMALES", "MALES", "LOCATIONNAME"))

# reformat dates
sub.paper$ENDDATETIME <- as.Date(sub.paper$ENDDATETIME)
sub.geopro$ENDDATETIME <- as.Date(sub.geopro$ENDDATETIME)
sub.dat.cdc$ENDDATETIME <- as.Date(sub.dat.cdc$ENDDATETIME)

# join full data set
full.data <- rbind(sub.paper, sub.geopro, sub.dat.cdc)




# create uniform species naming
full.data$SPECIES <- recode_factor(full.data$SPECIES, "Culex nigripalpus" = "Cx. nigripalpus")
table(full.data$SPECIES)


# assess site location names annomalies
table(full.data$LOCATIONNAME)

# correct duplicate site names
dat.cdc$NAME <- gsub("CDC AM Davey Lndscp Yrd", "CDC AM Davey Landscape Yard", as.character(dat.cdc$NAME))
dat.cdc$NAME <- gsub("CDC AM WATER PLANT SWAMP", "CDC AM Water Plant Swamp", as.character(dat.cdc$NAME))
dat.cdc$NAME <- gsub("CDC GGE South Desoto ", "CDC GGE South Desoto", as.character(dat.cdc$NAME))


# create list of site locations of interest
CDC.list <- c("CDC AM Ave Entrance",
              "CDC AM Ave Maria",
              "CDC AM Bellera",
              "CDC AM Davey Landscape Yard",
              "CDC AM Del Webb",
              "CDC AM Del Webb 2",
              "CDC AM Del Webb North",
              "CDC AM Dog Park",
              "CDC AM Dog Park 2",
              "CDC AM Emerson Park",
              "CDC AM EMS 32",
              "CDC AM Maple Ridge",
              "CDC AM Maple Ridge 11",
              "CDC AM North Park",
              "CDC AM North Park E",
              "CDC AM North Park Pasture",
              "CDC AM Oil Well Grade Rd",
              "CDC AM Patrick",
              "CDC AM Trap 7",
              "CDC AM Trap 8",
              "CDC AM Water Park",
              "CDC AM Water Plant",
              "CDC AM Water Plant Swamp",
              "CDC CS Collier Seminole 1",
              "CDC CS Collier Seminole 2",
              "CDC E Bayshore Dr",
              "CDC E Bayview",
              "CDC E Becca Ave",
              "CDC E Cope Ln",
              "CDC E Florida Sports Park",
              "CDC E Kings Lake",
              "CDC E Lely",
              "CDC E Lely High School",
              "CDC E Sports Park",
              "CDC E Weeks",
              "CDC E Wendy",
              "CDC E Whitaker",
              "CDC GGC Golden Gate High School",
              "CDC GGC Hickory Wood Dr",
              "CDC GGC Lancewood Way",
              "CDC GGC Whippoorwill",
              "CDC GGE 15th St NW",
              "CDC GGE 20th St SE",
              "CDC GGE 29th Ave SW",
              "CDC GGE 43rd and 40th",
              "CDC GGE 47th and 12th",
              "CDC GGE 66th Ave NE",
              "CDC GGE Buekowski",
              "CDC GGE Corkscrew Fire Station",
              "CDC GGE Desoto",
              "CDC GGE Everglades North",
              "CDC GGE Everglades South",
              "CDC GGE Fire Station 12",
              "CDC GGE Landfill",
              "CDC GGE Leann",
              "CDC GGE Morris",
              "CDC GGE North Desoto",
              "CDC GGE Oil Well Rd",
              "CDC GGE Phillips",
              "CDC GGE Quarry",
              "CDC GGE Richard St",
              "CDC GGE Ryan ",
              "CDC GGE Shady Hollow",
              "CDC GGE South Desoto",
              "CDC GGE Stivers",
              "CDC GGE Weber",
              "CDC GGE Willis",
              "CDC HC 6L Farms",
              "CDC HC Brandy Ln",
              "CDC HC Fiddlers Creek",
              "CDC HC KOA",
              "CDC HC Lipman House",
              "CDC HC Lipman Office",
              "CDC HC Sabal Palm",
              "CDC HC Trevisio Bay Dr",
              "CDC IM Airport Ditch",
              "CDC IM Anez",
              "CDC IM Animal Control",
              "CDC IM Audobon",
              "CDC IM Career Source",
              "CDC IM Carson Rd",
              "CDC IM Farm Village",
              "CDC IM Hector Ramirez",
              "CDC IM Immokalee",
              "CDC IM Immokalee Airport",
              "CDC IM Indian Camp Rd",
              "CDC IM James Brown",
              "CDC IM Lake Trafford",
              "CDC IM Lake Trafford Fire Dept",
              "CDC IM LCEC",
              "CDC IM LCEC 11",
              "CDC IM Little League",
              "CDC IM Little League Field",
              "CDC IM Little League Park",
              "CDC IM Naples Children's Hospital",
              "CDC IM Racetrack",
              "CDC IM Village Oaks Elementary",
              "CDC IM Water Plant",
              "CDC IM Williams Farms Swamp",
              "CDC KW Keewaydin 1",
              "CDC KW Keewaydin S Ctrl",
              "CDC KW Windstar",
              "CDC M Capri Blvd",
              "CDC M Goodland",
              "CDC M Isles of Capri",
              "CDC M Stevens Landing",
              "CDC N Collier Reserve",
              "CDC N Little Hickory",
              "CDC N Nicholas Blvd",
              "CDC N Palm River",
              "CDC N Pine Ridge",
              "CDC N Quail Creek",
              "CDC N Rose",
              "CDC N Veterans Park",
              "CDC N Victoria Park",
              "CDC N Water Park",
              "CDC N Willoughby",
              "CDC PY Canal Rd",
              "CDC PY Newman Dr",
              "CDC PY Picayune Trailhead",
              "CDC T Crayton",
              "CDC T Crayton ",
              "CDC T Klein",
              "CDC T test trap",
              "CDC T Town",
              "CDC T Wilderness",
              "CDC T Windmill Academy",
              "CDC XO 66th St SW",
              "CDC XO 92nd Ave",
              "CDC XO Caxambas",
              "CDC XO Fiddlers Creek",
              "CDC XO Golf Course",
              "CDC XO Mediterra",
              "CDC XO North Lake Dr",
              "CDC XO Picayune",
              "CDC XO Rail Head",
              "CDC XO Seminole",
              "CDC XO Willow",
              "NJL AM Ave Maria",
              "NJL AM Del Webb",
              "NJL AM Dog Park",
              "NJL AM Emerson Park",
              "NJL AM EMS 32",
              "NJL AM North Park",
              "NJL AM Pump House",
              "NJL AM Water Plant Swamp",
              "NJL E Boat Place",
              "NJL E Boca Ciega",
              "NJL E County Barn",
              "NJL E East",
              "NJL E Fire Station 42",
              "NJL E Hughes",
              "NJL E Lely",
              "NJL E Palm Springs",
              "NJL E Royal Harbor",
              "NJL E Tower Rd",
              "NJL E Weeks",
              "NJL GGC Collier Tire",
              "NJL GGC Dogwood Way",
              "NJL GGC Fraternal Order of Eagles",
              "NJL GGC Logan",
              "NJL GGC Rose",
              "NJL GGC Wyndemere",
              "NJL GGE 10th Ave SE Desoto",
              "NJL GGE 6th Ave S",
              "NJL GGE Buekowski",
              "NJL GGE Church",
              "NJL GGE EMS 12",
              "NJL GGE Everglades South",
              "NJL GGE Fire Station 12",
              "NJL GGE Landfill",
              "NJL GGE May",
              "NJL GGE Olde Florida",
              "NJL GGE Radio and 951",
              "NJL GGE Ryan",
              "NJL GGE Salinas",
              "NJL GGE South Desoto",
              "NJL GGE Stivers",
              "NJL GGE Wilson",
              "NJL HC Hidden Oaks Ln",
              "NJL HC KOA",
              "NJL HC Lord's Way",
              "NJL HC Naples Shores",
              "NJL HC Sabal Palm",
              "NJL IM Anez",
              "NJL IM Animal Control",
              "NJL IM Bass & Bass",
              "NJL IM Ernesto",
              "NJL IM Farm Village",
              "NJL IM Health Department",
              "NJL IM Immokalee Airport",
              "NJL IM Immokalee Tire",
              "NJL IM James Brown",
              "NJL IM Lake Trafford",
              "NJL IM Lake Trafford Fire Dept",
              "NJL IM Naples Children's Hospital",
              "NJL IM NCH & Del Webb",
              "NJL IM Nunez",
              "NJL IM Racetrack",
              "NJL IM Raulerson Rd",
              "NJL IM Water Plant",
              "NJL M Goodland",
              "NJL M Island Club",
              "NJL M Isles of Capri",
              "NJL M Tiger Tail Beach",
              "NJL M Woodway",
              "NJL N Bonita Bay",
              "NJL N Little Hickory",
              "NJL N Livingston",
              "NJL N Palm River",
              "NJL N Pelican Bay",
              "NJL N Pine Ridge",
              "NJL N Rose",
              "NJL N Venetian Way",
              "NJL T Coach House",
              "NJL T Park Shore",
              "NJL T Parker Aero",
              "NJL T Reynolds",
              "NJL T Town",
              "NJL XO Naples South",
              "NJL XO Palm",
              "NJL XO Seminole",
              "NJL XO Shadow Ridge",
              "NJL XO Shepherd of the Glades",
              "NJL XO Sod Farm",
              "NJL XO Unknown",
              "NJL XO Willow")


# create dataset containing sites of interest
complete.dat <- full.data[full.data$LOCATIONNAME %in% CDC.list,]



# add climate data to "complete.dat"
clim <- read.csv("climate_PROTECTED", fileEncoding="UTF-8-BOM")
clim$ENDDATETIME <- as.Date(clim$ENDDATETIME)
comp.data <- left_join(complete.dat, clim, by=c("ENDDATETIME"))



# subset data to include only Zones of interest (eg. AM, IM and GEE)
sub.data <- subset(comp.data, ZONE == "Zone AM" |
                     ZONE == "Zone IM" |
                     ZONE == "Zone GGE")


# subset data to include only Cx. nigripalpus
sub.data2 <- subset(sub.data, SPECIES == "Cx. nigripalpus")

# subset data to include only post 2017
my.data <- subset(sub.data2, ENDDATETIME > "2017-01-01")












#####################
#
# Model with 'prophet'
#
#####################

# plot population daily averages

df_tidy_mean <- my.data %>%
  filter(!is.na(FEMALES)) %>%
  group_by(ENDDATETIME) %>%
  summarise(n = n(),
            mean = mean(FEMALES),
            median = median(FEMALES)) 

ggplot(df_tidy_mean, aes(x=ENDDATETIME, y=mean)) +
  geom_line(aes(x=ENDDATETIME, y=mean)) +
  theme_bw()


# plot data coverage per day

ggplot(df_tidy_mean, aes(x=ENDDATETIME, y=n)) +
  geom_point(aes(x=ENDDATETIME, y=n)) +
  geom_smooth(aes(x=ENDDATETIME, y=n))+
  theme_bw()




### test with daily averages first
# "prophet" must include 'ds' (date) and 'y' (abundance)
# rename date and abundance in dataset

df_tidy_mean$ds <- df_tidy_mean$ENDDATETIME
df_tidy_mean$y <- df_tidy_mean$mean

head(df_tidy_mean)

m <- prophet(df_tidy_mean)
future <- make_future_dataframe(m, periods = 365)
forecast <- predict(m, future)

plot(m, forecast) +
  xlab("Date") +
  ylab("Cx. nigripalpus abundance") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Prophet model for Cx. nig. seasonal trends")




# Basic predictions machine learning

my.data$ds <- my.data$ENDDATETIME
my.data$y <- my.data$FEMALES

head(my.data)

m <- prophet(my.data)
future <- make_future_dataframe(m, periods = 365)
forecast <- predict(m, future)

plot(m, forecast) +
  xlab("Date") +
  ylab("Cx. nigripalpus abundance") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Prophet model for Cx. nig. seasonal trends")






# add climate regressors

pdat <- data.frame(ds = my.data$ENDDATETIME,
                   y = my.data$FEMALES,
                   temp = my.data$Thermometer,
                   humid = my.data$Hygrometer,
                   precip = my.data$Rain.Gauge)

# remove the few dates that have missing climate data
pdat <- na.omit(pdat)

# Forecast weather regressors
pfdat <- data.frame(ds=max(my.data$ENDDATETIME) + 1:25)

pprecip <- pdat %>% 
  select(ds,y=precip) %>% 
  prophet() %>%
  predict(pfdat)

phumid <- pdat %>% 
  select(ds,y=humid) %>% 
  prophet() %>%
  predict(pfdat)

ptemp <- pdat %>% 
  select(ds,y=temp) %>% 
  prophet() %>%
  predict(pfdat)

fdat <-  data.frame(ds=pfdat$ds,
                    precip=pprecip$yhat,
                    humid=phumid$yhat,
                    temp=ptemp$yhat)



# Fit the model (Seasonality automatically determined)
fit6 <- prophet() %>% 
  add_regressor('precip') %>% 
  add_regressor('humid') %>% 
  add_regressor('temp') %>% 
  fit.prophet(pdat)

future.fit <- make_future_dataframe(fit6, periods = 365)
forecast.fit <- predict(m, future.fit)

# remove negative values from projected yhats
forecast.fit$yhat <- pmax(forecast.fit$yhat,0)
forecast.fit$yhat_lower <- pmax(forecast.fit$yhat_lower,0)

# plot
plot(m, forecast.fit) +
  xlab("Date") +
  ylab("Cx. nigripalpus abundance") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Prophet model for Cx. nig. seasonal trends")




# evaluate weather regressors influence on mosquito populations

fit1 <- lm(y~humid+temp+precip,data=pdat)
summary(fit1)

pdat$resid[!is.na(pdat$y)] <- resid(fit1)

# Plot the residuals
ggplot(pdat,aes(ds,resid)) + 
  geom_point() + geom_smooth() +
  ggtitle("Linear Regression Residuals",
          subtitle = paste0("RMSE: ",round(sqrt(mean(pdat$resid^2,na.rm=TRUE)),2)))

Acf(pdat$resid, main="ACF of OLS Residuals")











#####################
#
# Model with 'modeltime'
#
#####################


# summarize mosquito count averages (cumulative across AM, IM, GGE)!!!
data <- my.data %>%
  mutate(day = floor_date(ENDDATETIME, "day")) %>%
  group_by(day) %>%
  summarize(avg = mean(FEMALES),
            temp = Thermometer,
            humid = Hygrometer,
            precip = Rain.Gauge)


# calculate the moving averages (3 moving average)
data3 <- data %>%
  arrange(day) %>% 
  mutate(avg3 = zoo::rollmean(avg, k = 3, fill = NA))  %>%
  select(day,
         avg,
         avg3)



# visualize time series data set
data %>%
  plot_time_series(day, avg, .interactive = FALSE)


# plot as boxplots (visualize variation)
data %>%
  plot_time_series_boxplot(day, avg,
                           .period = "1 month",
                           .facet_ncol  = 2, 
                           .interactive = FALSE)

# detect anomalies
data %>%
  plot_anomaly_diagnostics(day, avg,)


# frequency and trends
data %>%
  plot_stl_diagnostics(day, avg, 
                       .frequency = "auto", .trend = "auto",
                       .interactive = FALSE)

















#####################
#
# Test and train models
#
#####################


# split data to train and test model
# setting 'assess = 9 months" tells function to use last 9 months to data as a testing set
# setting 'cumulative = TRUE' tells the sampling to use all prior data as training set
splits <- data %>%
  time_series_split(assess = "9 months", cumulative = TRUE)

# visualize the train/test split
splits %>%
  tk_time_series_cv_plan() %>%
  plot_time_series_cv_plan(day, avg, .interactive = FALSE)



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
  #workflow_fit_rf,
  workflow_fit_prophet_boost
) 

model_table


# calibration is used to quantify error and estimate confidence intervals.
# perform model calibration on the testing set
# new columns are generated (".type" and ".calibration_data"), 
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
# Refit model and run forecast
#
#####################

# Refit and Forecast Forward
calibration_table %>%
  modeltime_refit(data) %>%
  modeltime_forecast(h = "12 months", actual_data = data) %>%
  plot_modeltime_forecast(.interactive = FALSE)




























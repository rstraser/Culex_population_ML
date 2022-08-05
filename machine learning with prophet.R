##########################################
### Deep learning for mosquito control ###
##########################################

# install packages
install.packages('prophet')
library(prophet)
library(dplyr)
library(ggplot2)
library(dplyr)

#################################################################################
#
#   NOTES: 
#   - run "FS trap data merge and summarize v20" for exporting FieldSeeker data
# 
#################################################################################


#### Cx. nigripalpus in operations CDCs ####

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

dat.cdc <- dat.cdc[dat.cdc$NAME %in% CDC.list,]


# structue summary data set
dat.cdc <- graph_prep(dat.cdc)


## Plot nigripalpus time series by trap site

ggplot(dat.cdc, aes(Date, FEMALES, color = NAME)) +
  geom_line() +
  facet_grid(ZONE ~ .) +
  ylab(paste(sp.cdc, "females", sep = " ")) +
  xlab("Collection day") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_date(date_labels = "%b-%Y", date_breaks = "1 month", limits = as.Date(c("2019-01-01", "2022-04-28")))
















# plot the raw data
# Plot the residuals
ggplot(dat.cdc,aes(Date,FEMALES)) + 
  geom_point() + geom_smooth() +
  ylim(0, 10000)












# # testing 'prophet' package with Culex nigripalpus
# # subset Cx. quinquefasciatus
# cx <- subset(dat.sa, SPECIES == "Cx. nigripalpus")
# names(cx)
# 
# 
# 
# exploratory plot
plot(FEMALES ~ ENDDATETIME, dat.cdc, type = "l")



# Basic predictions machine learning

# "prophet" must include 'ds' (date) and 'y' (abundance)
# rename date and abunance in dataset

dat.cdc$ds <- dat.cdc$ENDDATETIME
dat.cdc$y <- dat.cdc$FEMALES



m <- prophet(dat.cdc)
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







# inspecting model components
prophet_plot_components(m, forecast)






# you can account for peaks in seasonality to improve model
# but probabaly best predicts known spikes (eg. saltmash mozzie emergance)

# this examples accounts for playoff finals in searching Lebron James
# playoff_brackets <- data_frame(
#   holiday = 'playoffs',
#   ds = seq(as.Date("2014/09/04"), by = "day", length.out = 5),
#  lower_window = 0,
#  upper_window = 45
# )
# playoff_finals <- data_frame(
#  holiday = 'playoff_finals',
#  ds = as.Date(c('2016-06-02', '2015-06-04', '2014-06-05')),
#  lower_window = 0,
#  upper_window = 20
# )





# removing outliers

# this methods of outlier detection is based on percentiles, based on 1 and 99 percentiles
#find Q1, Q3, and interquartile range for values in column A
Q1 <- quantile(dat.cdc$y, .01)
Q3 <- quantile(dat.cdc$y, .99)
IQR <- IQR(dat.cdc$y)

#only keep rows in dataframe that have values within 1.5*IQR of Q1 and Q3
dat.cdc.ed <- subset(dat.cdc, dat.cdc$y> (Q1 - 1.5*IQR) & dat.cdc$y< (Q3 + 1.5*IQR))

#view row and column count of new data frame
dim(dat.cdc.ed) 









### rerun the model

m <- prophet(dat.cdc.ed)
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




# inspecting model components
prophet_plot_components(m, forecast)











##############################


fcast_results <- 
  train_valid_xreg %>%
  mutate(prophet_fit = map2(
    train,
    xreg,
    function(x, y)
      prophet(
        df = x,
        holidays = y,
        growth = "auto",
        yearly.seasonality = "auto",
        weekly.seasonality = "auto"
      )
  ),
  pred_prophet = map2(prophet_fit, valid, predict)
  )


cv_slice = 1

slice_training <-
  fcast_results$train[[cv_slice]] %>%
  mutate(key = "Training") %>%
  rename(value = y)

slice_validation <-
  fcast_results$pred_prophet[[cv_slice]] %>%
  transmute(
    ds = ymd(ds),
    predicted = yhat,
    actual = fcast_results$valid[[cv_slice]]$y
  ) %>%
  gather(key, value, -ds) %>%
  mutate(key = str_to_title(key))


bind_rows(
  slice_training,
  slice_validation
) %>%
  ggplot(aes(ds, value, color = key)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.3) +
  my_plot_theme() +
  facet_zoom(x = ds %in% tail(slice_validation$ds, 90)) +
  labs(
    x = "Date (Day)",
    y = "Units",
    color = "Part"
  )










##########################################
### Modeltime for mosquito control ###
##########################################



#################################################################################
#
#   NOTES: 
#   - see modeltime tutorial for detailed information:
#   - https://www.business-science.io/code-tools/2020/06/29/introducing-modeltime.html?utm_content=bufferd20d1&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer#
#
#################################################################################

#### Load required packages ###
library(tidymodels)
library(modeltime)
library(timetk)   
library(lubridate)
library(tidyverse)
library(parsnip)
library(ggplot2)
library(prophet)
library(randomForest)
library(zoo)
library(ggplot2)
library(reshape2)
library(lubridate)
library(zoo)
library(dplyr)
library(tidyr)


#### Pull coordinate data from local hard drive ###

# original FieldSeeker coordinates Becky used
#dat.coord <- read.csv("Z:/Research Dept/Research/Surveillance data records/Field Seeker/CollierFieldSeeker3_20220508/TrapLocation_3.csv")  #temp home location

fs.coord <- read.csv("Z:/Research Dept/Rob Straser/Historical data/Historical trap coordinates/FS trap export 20220508.csv")  #temp home location
geopro.coord <- read.csv("Z:/Research Dept/Rob Straser/Historical data/Historical trap coordinates/GeoPro coords.csv")  #temp home location
paper.coord <- read.csv("Z:/Research Dept/Rob Straser/Historical data/Historical trap coordinates/site correction ref list.csv")  #temp home location


# subset columns to join
sub.fs.coord <- subset(fs.coord,
                      select = c("Name", "x", "y"))
sub.geopro.coord <- subset(geopro.coord,
                       select = c("Name", "x", "y"))
sub.paper.coord <- subset(paper.coord,
                       select = c("Name", "x", "y"))

# join all coordinate data
dat.coord <- rbind(sub.paper.coord, sub.geopro.coord, sub.fs.coord)
#write.csv(dat.coord, "Z:/Research Dept/Rob Straser/totalcoords")

































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



























#### Functions (from Heinig FieldSeeker analyses) ###

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

graph_prep1 <- function(df){
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

graph_prep <- function(df){
  x <- df  #data frame
  x$Date <- as.Date(x$ENDDATETIME)
  x$DateRef <- x$Date
  year(x$DateRef) <- 2000  #this was a leap year; use this to overplot multiple years
  x$Year <- year(x$Date)   #useful for color coding overplots
  x$Year <- factor(x$Year)
  # TO DO: add epi week column
  x$LOCATIONNAME <- factor(x$LOCATIONNAME)  #get rid of extraneous factor values following subsetting
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
#write.csv(dat.cdc, "Z:/Research Dept/Rob Straser/test.dat.cdc")




# join in historical data (GeoPro and Paper records)

# read in Geopro and Paper records
setwd("Z:/Research Dept/Rob Straser/")
paper <- read.csv("Historical data/2011-2014 FS compatible paper trap data.csv")
geopro <- read.csv("Historical data/2015-2019 GeoPro and paper consolidated data.csv")


sub.paper <- subset(paper,
                            select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "SPECIES",
                                       "FEMALES", "MALES", "LOCATIONNAME"))


sub.geopro <- subset(geopro,
                             select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "SPECIES",
                                       "TOTAL", "MALES", "LOCATIONNAME"))

# no delineation between FEMALES and TOTAL in GeoPro data
# assuming TOTALS is equal to FEMALES
sub.geopro <- sub.geopro %>%
                    rename(FEMALES = TOTAL)

sub.dat.cdc <- subset(dat.cdc,
                             select = c("ENDDATETIME", "TRAPTYPE", "ZONE", "SPECIES",
                                        "FEMALES", "MALES", "LOCATIONNAME"))


sub.paper$ENDDATETIME <- as.Date(sub.paper$ENDDATETIME)
sub.geopro$ENDDATETIME <- as.Date(sub.geopro$ENDDATETIME)
sub.dat.cdc$ENDDATETIME <- as.Date(sub.dat.cdc$ENDDATETIME)
#write.csv(sub.paper, "Z:/Research Dept/Rob Straser/paper")
#write.csv(sub.geopro, "Z:/Research Dept/Rob Straser/geopro")
#write.csv(sub.dat.cdc, "Z:/Research Dept/Rob Straser/datcdc")



full.data <- rbind(sub.paper, sub.geopro, sub.dat.cdc)
#write.csv(full.data, "Z:/Research Dept/Rob Straser/fulldata")








#################################
################################# 
################################# 
#################################








# align naming for culex nigripaplus
full.data$SPECIES <- recode_factor(full.data$SPECIES, "Culex nigripalpus" = "Cx. nigripalpus")
table(full.data$SPECIES)

# assess site location names
a <- table(full.data$LOCATIONNAME)
#write.csv(a, "Z:/Research Dept/Rob Straser/AAAAA")


# corrert duplicate site names
dat.cdc$NAME <- gsub("CDC AM Davey Lndscp Yrd", "CDC AM Davey Landscape Yard", as.character(dat.cdc$NAME))  #fixes site crossover issue; coordinates not affected
dat.cdc$NAME <- gsub("CDC AM WATER PLANT SWAMP", "CDC AM Water Plant Swamp", as.character(dat.cdc$NAME))  #fixes site crossover issue; coordinates not affected
dat.cdc$NAME <- gsub("CDC GGE South Desoto ", "CDC GGE South Desoto", as.character(dat.cdc$NAME))  #fixes site crossover issue; coordinates not affected


              
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





















######################
# THIS IS ONLY WORKS FOR FEILDSEEKER DATA
#####################
#
## create the dataset
#dat.cdc <- dat.cdc[dat.cdc$NAME %in% CDC.list,]
#
## join trap location coordinates to dataset
#dat.coord <- dat.coord %>% 
#                rename(NAME = Name)
#data.coord <- merge(dat.cdc, dat.coord, by='NAME')
#
#
## calculate average female abundance per day
#data <- dat.cdc %>%
#  mutate(day = floor_date(ENDDATETIME, "day")) %>%
# group_by(day) %>%
#  summarize(avg = mean(FEMALES))
#
#
## group by ZONE
#data1 <- dat.cdc %>%
#  mutate(day = floor_date(ENDDATETIME, "day")) %>%
#  group_by(ZONE, day) %>%
#  summarize(avg = mean(FEMALES))
#  
#########################




# create the dataset
complete.dat <- full.data[full.data$LOCATIONNAME %in% CDC.list,]

# join trap location coordinates to dataset
#dat.coord <- dat.coord %>% 
#  rename(LOCATIONNAME = Name)
#complete.dat <- merge(full.dat, dat.coord, by='LOCATIONNAME')


#   Notes:
#   "complete.dat' = includes cx nig counts from CDC and NJL traps from paper, geopro, feildseeker including coords

#write.csv(complete.dat, "Z:/Research Dept/Rob Straser/completedat")
























##################################
### Add climate data to "complete.dat"
##################################


clim <- read.csv("Z:/Research Dept/Rob Straser/Historical data/Parkside Elementary School weather station.csv", fileEncoding="UTF-8-BOM")  #temp home location

# join complete data (mosquito abundnace) and climate data
# Remember to rename "Timestamp" to "ENDDATETIME"

clim$ENDDATETIME <- as.Date(clim$ENDDATETIME)
comp.data <- left_join(complete.dat, clim, by=c("ENDDATETIME"))

#write.csv(comp.data, "Z:/Research Dept/Rob Straser/complete dataset with climate")




# subset data to include only Zones AM, IM and GEE
sub.data <- subset(comp.data, ZONE == "Zone AM" |
                 ZONE == "Zone IM" |
                 ZONE == "Zone GGE")


# subset data to include only Zones AM, IM and GEE
sub.data2 <- subset(sub.data, SPECIES == "Cx. nigripalpus")

# subset data to include only post 2017
my.data <- subset(sub.data2, ENDDATETIME > "2017-01-01")




















#####################################
#####################################
#####################################
###########
###########   Prophet package
###########
#####################################
#####################################
#####################################






# plot the daily averages

df_tidy_mean <- my.data %>%
  filter(!is.na(FEMALES)) %>%
  group_by(ENDDATETIME) %>%
  summarise(n = n(),
            mean = mean(FEMALES),
            median = median(FEMALES)) 



ggplot(df_tidy_mean, aes(x=ENDDATETIME, y=mean)) +
  geom_line(aes(x=ENDDATETIME, y=mean)) +
  theme_bw()



ggplot(df_tidy_mean, aes(x=ENDDATETIME, y=n)) +
  geom_point(aes(x=ENDDATETIME, y=n)) +
  geom_smooth(aes(x=ENDDATETIME, y=n))+
  theme_bw()




### test with daily averages first
# "prophet" must include 'ds' (date) and 'y' (abundance)
# rename date and abunance in dataset

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


names(m)









# Basic predictions machine learning

# "prophet" must include 'ds' (date) and 'y' (abundance)
# rename date and abunance in dataset

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

#write.csv(pdat, "Z:/Research Dept/Rob Straser/pdat")

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


plot(m, forecast.fit) +
  xlab("Date") +
  ylab("Cx. nigripalpus abundance") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Prophet model for Cx. nig. seasonal trends")






























#################
##### REMOVING OUTLIERS CREATES MAY DUPLICATE DAILY VALUES FOR SOME REASON - DONT RUN
#################

# removing outliers

# this methods of outlier detection is based on percentiles, based on 1 and 99 percentiles
#find Q1, Q3, and interquartile range for values in column A
Q1 <- quantile(pdat$y, .01)
Q3 <- quantile(pdat$y, .99)
IQR <- IQR(pdat$y)

#only keep rows in dataframe that have values within 1.5*IQR of Q1 and Q3
pdat2 <- subset(pdat, pdat$y> (Q1 - 1.5*IQR) & pdat$y< (Q3 + 1.5*IQR))

#view row and column count of new data frame
dim(pdat2) 
head(pdat2)







### rerun the model


# remove the few dates that have missing climate data
pdat2 <- na.omit(pdat2)

# Forecast weather regressors
pfdat2 <- data.frame(ds=max(pdat2$ds) + 1:25)

pprecip <- pdat2 %>% 
  select(ds,y=precip) %>% 
  prophet() %>%
  predict(pfdat2)



phumid <- pdat2 %>% 
  select(ds,y=humid) %>% 
  prophet() %>%
  predict(pfdat2)

ptemp <- pdat2 %>% 
  select(ds,y=temp) %>% 
  prophet() %>%
  predict(pfdat2)

fdat2 <-  data.frame(ds=pfdat2$ds,
                    precip=pprecip$yhat,
                    humid=phumid$yhat,
                    temp=ptemp$yhat)



# Fit the model (Seasonality automatically determined)
fit7 <- prophet() %>% 
  add_regressor('precip') %>% 
  add_regressor('humid') %>% 
  add_regressor('temp') %>% 
  fit.prophet(pdat2)





future.fit <- make_future_dataframe(fit7, periods = 365)
forecast.fit <- predict(m, future.fit)



plot(fit7, forecast.fit) +
  xlab("Date") +
  ylab("Cx. nigripalpus abundance") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ggtitle("Prophet model for Cx. nig. seasonal trends")









# inspecting model components
prophet_plot_components(fit6, forecast.fit)






# plot other variables (including climate data)

pdat


temp.plot <- ggplot(my.data, aes(x=ds, y=temp), color = "black") +geom_line()+ theme_bw()
temp.plot

humid.plot <- ggplot(pdat, aes(x=ds, y=humid), color = "black") +geom_line()+ theme_bw()
humid.plot

precip.plot <- ggplot(pdat, aes(x=ds, y=precip), color = "black") +geom_line()+ theme_bw()
precip.plot












############
## TESTING A NEW APPROACH
#############


library(tidyverse)
library(gridExtra)
library(lubridate)
library(leaflet)
library(randomForest)
library(forecast)
library(prophet)



# How well weather regressors are at predicting mosquitoes?

fit1 <- lm(y~humid+temp+precip,data=pdat)
summary(fit1)

pdat$resid[!is.na(pdat$y)] <- resid(fit1)


library(ggplot2)
# Plot the residuals
ggplot(pdat,aes(ds,resid)) + 
  geom_point() + geom_smooth() +
  ggtitle("Linear Regression Residuals",
          subtitle = paste0("RMSE: ",round(sqrt(mean(pdat$resid^2,na.rm=TRUE)),2)))


Acf(pdat$resid, main="ACF of OLS Residuals")
















#####################################
#####################################
#####################################
###########
###########   Second take on Prophet Model
###########
#####################################
#####################################
#####################################


# use the "pdat" data set - includes ds, y, temp, humid, precip


# check proportion of missing values in each variable
# proportion of NA values for these columns

df <- pdat
names(df)

library(ggcorrplot) #for plotting correlations
#calculating and plotting correlation 
corr <- cor(df %>% select(y, temp, humid, precip))
ggcorrplot(corr , hc.order = TRUE, method = "circle")




# creating train-test split with split having exactly 60 days
train <- df %>% slice(1:(dim(df)[1] - (60*24)))  # train having the rest
test <- df %>% slice((dim(df)[1] - (60*24)) + 1:dim(df)[1]) # test having last 60 days




### RANDOM FOREST
#partial dependence plot approach
library(randomForest) # for randomForest, partialPlot, and varImpPlot functions

set.seed(100) # for reproducibility

train.rf <- randomForest(y ~ ds + temp + precip + humid, data = df, importance = TRUE)

#plotting relative importance of variables
varImpPlot(train.rf)







# aggregating data at daily level
train_day <- train

test_day <- test

# creating prophet train and train_validation sets
prophet_train <- train_day[1:(nrow(train_day) - 60),]
prophet_train_val <- train_day[(nrow(train_day) - 60 + 1):nrow(train_day),] 

# prophet test set
prophet_test <- test_day





# library for prophet model
library(prophet)

#creating an improved prophet model
im <- (prophet(
  df = NULL,   # Dataframe containing the history
  growth = "linear",    # trend change/growth can't be logistic or flat for this TS
  #changepoints = c('2013-08-10','2014-01-10', '2014-08-02', '2015-01-13', '2015-07-13',  
  #'2016-01-22', '2016-08-04'), 
  n.changepoints = 25, #more than default, might overfit 
  changepoint.range = 0.80, # Proportion of history in which trend changepoints will be estimated
  yearly.seasonality = TRUE, # Default Fourier Order
  weekly.seasonality = FALSE, # no evidence for temp change for days of a week
  daily.seasonality = FALSE, # Daily seasonality locked as off as data is daily leveled
  holidays = NULL, # no evidence that holidays affect temp 
  seasonality.mode = "additive", # by observation
  seasonality.prior.scale = 10, # default
  holidays.prior.scale = 10, # default for regressors too
  changepoint.prior.scale = 0.05, # default
  mcmc.samples = 0, # default
  interval.width = 0.80, #default
  uncertainty.samples = 1000, #default
  fit = TRUE
))

# adding external regressors
im = add_regressor(im,'temp',standardize = FALSE) # added temp as an external regressor
im = add_regressor(im,'precip', standardize = FALSE) # added precip as an external regressor
im = add_regressor(im,'humid', standardize = FALSE) # added humid as an external regressor

# fitting a prophet model
im = fit.prophet(im, df = prophet_train)

#making future dataframes
future_im <- make_future_dataframe(im, periods = 60, 
                                   freq = "day", include_history = TRUE) #predictions for two months

future_im$temp = head(train_day$temp ,nrow(future_im))
future_im$precip = head(train_day$precip ,nrow(future_im))
future_im$humid = head(train_day$humid, nrow(future_im))




# prediction
fcst_im <- predict(im,future_im) # creating forecast for 60 days
tail(fcst_im[c('ds', 'yhat', 'yhat_lower', 'yhat_upper')]) #observing tail observations
dyplot.prophet(im,fcst_im , uncertainty = TRUE) # creating interactive plots for the forecast



# remove negative values from projected yhats
fcst_im$yhat <- pmax(fcst_im$yhat,0)
fcst_im$yhat_lower <- pmax(fcst_im$yhat_lower,0)




# plotting components of the forecast
prophet_plot_components(
  im,
  fcst_im,
  uncertainty = TRUE,
  plot_cap = TRUE,
  yearly_start = 0,
  render_plot = TRUE
)




# out of sample (validation set) assessment
RMSE_im = sqrt(mean((tail(train_day$y,60) - tail(fcst_im$yhat,60))^2))
MAE_im = mean(abs( tail(train_day$y,60) - tail(fcst_im$yhat,60) ))
MAPE_im = mean(abs(  (tail(train_day$y,60) - tail(fcst_im$yhat,60))  / tail(train_day$y,60) ))

tibble("RMSE"= c(round(RMSE_im,4)), "MAE" = round(MAE_im,4), "MAPE" = round(MAPE_im,4))









#making future dataframes
future_fm <- make_future_dataframe(im, periods = 60, 
                                   freq = "day", include_history = FALSE) #predictions for two months

future_fm$WSPM = prophet_test$WSPM 
future_fm$DEWP = prophet_test$DEWP
future_fm$O3 = prophet_test$O3
future_fm$PRES = prophet_test$PRES

# prediction
fcst_fm <- predict(fm,future_fm) # creating forecast for 60 days

#regressor coefficients 
regressor_coefficients(im)






#####################################
#####################################
#####################################
###########
###########   Modeltime package
###########
#####################################
#####################################
#####################################


# summarize mosquito count averages (cummulative across AM, IM, GGE)!!!
data <- my.data %>%
  mutate(day = floor_date(ENDDATETIME, "day")) %>%
  group_by(day) %>%
  summarize(avg = mean(FEMALES),
            temp = Thermometer,
            humid = Hygrometer,
            precip = Rain.Gauge)



#write.csv(data, "Z:/Research Dept/Rob Straser/cignigavg")



# calculate the moving averages (3 moving average)
data3 <- data %>%
  arrange(day) %>% 
  mutate(avg3 = zoo::rollmean(avg, k = 3, fill = NA))  %>%
  select(day,
         avg,
         avg3)









##########
## USE THE PDAT DATASET INSTEAD
###########


#data <- pdat %>% 
#          rename(day = ds,
#                 avg = y)






# visualize time series data set
data %>%
 plot_time_series(day, avg, .interactive = FALSE)






# can also plot as boxplots
data %>%
  plot_time_series_boxplot(day, avg,
                           .period = "1 month",
                   .facet_ncol  = 2, 
                   .interactive = FALSE)

# can detect anomalies
data %>%
  plot_anomaly_diagnostics(day, avg,)


# can detect frequency and trends
data %>%
  plot_stl_diagnostics(day, avg, 
                       .frequency = "auto", .trend = "auto",
                       .interactive = FALSE)



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
  time_series_split(assess = "9 months", cumulative = TRUE)

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
  #workflow_fit_rf,
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




























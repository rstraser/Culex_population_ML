#### Mapping


str(data.coord)
names(data.coord)





#Load the library
library("ggmap")

# Get the latest Install
# if(!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("dkahle/ggmap", ref = "tidyup", force=TRUE)

#Set your API Key
ggmap::register_google(key = "AIzaSyAT1lO6EJiCMOD7WIOQKGkofKTSFxZunKA")








# generate base map for Collier County
collier_basemap <- get_map(location=c(lon = -81.679235, lat = 26.188009), zoom=10, maptype = 'satellite', color = "bw", source = 'google')
ggmap(collier_basemap)





# structue summary data set
df.coord <- graph_prep(data.coord)
head(df.coord)


#produce the final map
ggmap(collier_basemap) + geom_point(data=df.coord, aes(x = x, y = y, colour= ZONE, size = FEMALES, alpha = 0.25)) +
  facet_wrap(~ Year)











#Le script qui synthètise les covariables effort, wsprop, suitability avec les présences

library(pacman)
p_load(httr, httpuv,brms, StanHeaders,  spOccupancy, jsonlite, dplyr, data.table, readr, sf, ggplot2, sp, terra, mapview, MCMCvis, raster, viridisLite, viridis, tidyr, tibble, tidyverse, boot)
set.seed(500)

your.path ="C:/Users/IF_240606/Desktop/Projet/" #replace


#### auth.code requirements####
if(file.exists(paste0(your.path,".Renviron"))){
  readRenviron(paste0(your.path,".Renviron"))
} else {
  print("No .Renviron existing!")
}

### set-up the authentication rights
client_key = Sys.getenv("client_key")
client_secret =  Sys.getenv("client_secret")
file.remove(".httr-oauth")
auth.code <- httr::oauth2.0_token(endpoint = httr::oauth_endpoint(authorize="https://auth.infoflora.ch/oauth2/authorize",
                                                                  acces="https://auth.infoflora.ch/oauth2/token"),
                                  app = httr::oauth_app(
                                    appname = "infoflora",
                                    key = client_key,
                                    secret = client_secret
                                  ),
                                  scope = "observation:read-all",
                                  use_oob = FALSE, cache = TRUE, client_credentials = TRUE) #Default to FALSE. Set to TRUE to use Client Credentials Grant instead of Authorization Code Grant.
auth.code$credentials


####object for specifying taxa & model dates####
taxon = c(1028530, 1019640, 1041760,1006340)#if what to loop
l=2
taxa = taxon[l]

taxon_names = c("Medicago","Festuca","Saxifraga","Baldellia")

suitatbility_strings = c(paste0(your.path,"map_eco_tiled/1028530_Medicago_lupulina.tif"),paste0(your.path,"map_eco_tiled/1019640_Festuca_valesiaca.tif"),paste0(your.path,"map_eco_tiled/1041760_saxifraga_hirculus.tif"), paste0(your.path,"map_eco_tiled/1041760_saxifraga_hirculus.tif"))

SWISS <- st_union(st_read(paste0(your.path,"swissBOUNDARIES3D_1_5_LV95_LN02.gpkg"), layer = "tlm_landesgebiet"))
SWISS <- st_as_sf(data.frame("inout"=1),geom=SWISS) 

SWISS <- st_union(st_read(paste0(your.path,"swissBOUNDARIES3D_1_5_LV95_LN02.gpkg"), layer = "tlm_landesgebiet"))
SWISS <- st_as_sf(data.frame("inout"=1),geom=SWISS) 

dates <-c("1000-01-01","1995-01-01","1996-01-01","1997-01-01","1998-01-01","1999-01-01","2000-01-01","2001-01-01","2002-01-01","2003-01-01","2004-01-01","2005-01-01","2006-01-01","2007-01-01","2008-01-01","2009-01-01","2010-01-01","2011-01-01","2012-01-01","2013-01-01","2014-01-01","2015-01-01","2016-01-01","2017-01-01","2018-01-01","2019-01-01","2020-01-01","2021-01-01","2022-01-01","2023-01-01","2024-01-01")

source(paste0(your.path,"WSProp_55.R")) #loading in the wsprop function

####CREATING THE 5X5 GRID#### 
lower_left <- data.frame(xmin=5.959163471, ymin = 45.816873303) #GRID 5X5 OPERATIONS
lower_left_sf <- st_as_sf(lower_left, coords = c('xmin', 'ymin'), crs = 4326) %>% 
  st_transform(crs = 2056)
grid <- st_make_grid(SWISS, cellsize = c(5000, 5000), what = "polygons",offset = st_bbox(lower_left_sf)[c("xmin", "ymin")])
GRID <- st_as_sf(data.frame(id_grid=1:length(grid)),geom=grid,crs=2056)
GRID <- st_join(GRID,SWISS,join=st_intersects) #output should be 3150x3, CHECK, if 3700 rerun swiss <- and grid


####Loading the effort DF for 30 years 5x5####
GRID_taxa_effort <- readRDS(paste0(your.path,"GRID_taxa_effort.30years.Rdata")) 

#if want to see code for generating these refer to script in folder ""

#### presence grid for each year ####
presence_taxa_30 <- list()
#get all obs for a specific taxa and intersect with GRID 
my.request <- httr::GET(paste0("https://obs.infoflora.ch/rest/v4/observations?",
                               "srid=2056&",
                               "limit=4000&",
                               "with_total=true&",
                               "taxa=",taxa,
                               "&taxa_matching=synonyms&",
                               "taxa_relations_references=3,5,6,7,8"),
                        config(token = auth.code))

#find total #of obs for taxa then loop and add offset
total_num_obs <- as.numeric(my.request$headers$`x-filtered-total`)
total_loops <- ceiling(total_num_obs/4000)

for(i in 1:total_loops){
  
  my.request <- tryCatch({
    httr::GET(paste0("https://obs.infoflora.ch/rest/v4/observations?",
                     "srid=2056&",
                     "offset=",as.character((i-1)*4000),
                     "&limit=4000&",
                     "with_total=true&",
                     "taxa=",taxa,
                     "&taxa_matching=synonyms&",
                     "taxa_relations_references=3,5,6,7,8"),
              config(token = auth.code))}, 
    
    error= function(e){return(NULL)})
  
  if(!is.null(my.request)){
    my.request <- httr::content(my.request, as="text",encoding = "UTF-8") # convert as text 
    my.request <- jsonlite::fromJSON(my.request) 
    my.request <- my.request$data #Table of all available taxa in CH (and their IDs)
  }
  
  #row bind the new 4000 set of obs
  if(i == 1){
    my_observations <- my.request
  }
  
  else{
    my_observations <- bind_rows(my_observations, my.request)
  }
}

#Checking validation state
validation.status<-unique(my.request$v_validation_status[
  my.request$v_validation_status>=200 &
    my.request$v_validation_status != 250 &
    my.request$v_validation_status != 260])

for(c in 1:(length(dates)-1)){
  
  #Filtering observations 
  my.sp.filter<-sp.filter(my_observations, 
                          with.geom = TRUE, 
                          exclude.prec = 2500, #parameters for xy radius max
                          presence.status = 2, #keep only presences 
                          introduction.status = c(0,1,2,3), #introduction status not 4
                          doubt.status = c(0,1), #eliminate error and doubtful observations
                          date.interval = c(dates[c],dates[c+1]) 
  )
  taxa_data <- my.sp.filter[!is.na(my.sp.filter$x),c("x","y","obs_id")]
  
  taxa_data <- taxa_data %>% 
    st_as_sf(coords = c('x', 'y')) %>% 
    st_set_crs(2056)
  
  placeholder <- st_intersection(GRID,taxa_data)
  id_presence <- (unique(placeholder$id_grid))
  
  grid_presence <- rep(0,nrow(GRID))
  
  for(i in 1:length(id_presence)){
    grid_presence[id_presence[i]] =1 #finding only presence in grid
  }
  presence_yr <- paste0("presence",dates[c+1]) #want it to be called by year 94
  presence_taxa_30[[presence_yr]] <- grid_presence   
}

####Synthesise all data into 1 data####
#this df will be use to create the list of object required for out stPGOcc model
taxa_all_info <- (wsprop55(taxa))

#generalizing sutability factor for 5x5 grid  
r_taxa <- terra::rast(suitatbility_strings[l])
r_taxa <- terra::project(r_taxa,"epsg:2056")
suita_taxa <- terra::zonal(r_taxa, vect(GRID) , fun='mean') #make sure rasterize function maintains the proper resolution size of r_medicago grid
suitability <- suita_taxa$layer
suitability[is.na(suitability)] <- 0

####Combining all our info: wsprop, presence, suita, effort####
taxa_all_info <- cbind(taxa_all_info, suitability, GRID_taxa_effort[,2:31], presence_taxa_30[1:30], st_drop_geometry(GRID$inout))

taxa_all_info <- subset(taxa_all_info, !is.na(st_drop_geometry.GRID.inout.)) #inside CH border

#object taxa_all_info is what we will use in scripts "stPGOcc_model_1849.R" & "stPGOcc.model_reduced.R" to create respective models

####extra: if you can to plot observations from a specific year over the grid####
c=1 #j is the index for the selection of dates, j=1 -> 1800-1994

my.sp.filter<-sp.filter(my_observations, 
                        with.geom = TRUE, 
                        exclude.prec = 2500, #parameters for xy radius max
                        presence.status = 2, #keep only presences 
                        introduction.status = c(0,1,2,3), #introduction status not 4
                        doubt.status = c(0,1), #eliminate error and doubtful observations
                        date.interval = c(dates[c],dates[c+1]) 
)
taxa_data <- my.sp.filter[!is.na(my.sp.filter$x),c("x","y","obs_id")]

taxa_data <- taxa_data %>% 
  st_as_sf(coords = c('x', 'y')) %>% 
  st_set_crs(2056)
subgrid <- subset(GRID, !is.na(inout))
grid_bbox <- st_bbox(subgrid)

plot2save <- ggplot() +
  geom_sf(data=subgrid, fill = NA, color = "gray") +
  geom_sf(data=taxa_data, color="blue", size=2)+
  geom_sf(data=ws, color="red", alpha=0.01)+
  theme_minimal() +
  labs(title = "Spatial Plot of Observations for Medicago Lupilina L. 2023",
       x = "Longitude", 
       y = "Latitude") +
  coord_sf(xlim = c(grid_bbox["xmin"], grid_bbox["xmax"]), 
           ylim = c(grid_bbox["ymin"], grid_bbox["ymax"])) 

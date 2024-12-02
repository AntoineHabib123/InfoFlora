
####To Load before running function####
library(pacman)
source("C:/Users/IF_240606/Desktop/Projet/functions.R")
p_load(httr, httpuv, StanHeaders,  eSDM, spOccupancy, jsonlite, dplyr, data.table, readr, sf, ggplot2, sp, terra, mapview, MCMCvis, raster, viridisLite, viridis, tidyr, tibble, tidyverse)

if(file.exists("C:/Users/IF_240606/Desktop/Projet/.Renviron")){
  readRenviron("C:/Users/IF_240606/Desktop/Projet/.Renviron")
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

#TO LOAD BEFORE RUNNING FUNCTION!
atlas27_ids <- httr::GET(paste0('https://obs.infoflora.ch/rest/v4/communities/0/public/atlases/27/names?offset=0&limit=10000&with_total=true',
                                "lmit=10000&",
                                "with_total=true&",
                                "id=0&",
                                "atlas_id=27"),
                         httr::config(token = auth.code)) %>% 
  httr::content(as="text",encoding = "UTF-8") %>% 
  jsonlite::fromJSON() %>% 
  .$data #Table of all available taxa in CH (and their IDs)

# get taxon id from species name
taxonid<-httr::GET(paste0("https://obs.infoflora.ch/rest/public/v4/taxonomy/taxa?",
                          "lang=en&",
                          "limit=10000&",
                          "references=3&",
                          "ranks=2263%2C696"),
                   httr::config(token = auth.code)) %>% 
  httr::content(as="text",encoding = "UTF-8") %>% 
  jsonlite::fromJSON() %>% 
  .$data #Table of all available taxa in CH (and their IDs)

SWISS <- st_read("C:/Users/IF_240606/Desktop/Projet/swissBOUNDARIES3D_1_5_LV95_LN02.gpkg", layer = "tlm_kantonsgebiet")

ws <- st_read("C:/Users/IF_240606/Desktop/Projet/welten_sutter_with_id.gpkg", layer = "welten_sutter")




wsprop55 <- function(taxa){ 
  
  #### Condition for method to apply ####
  stopifnot(is.numeric(taxa), taxa %in% taxonid$id, taxa %in% atlas27_ids$taxon_id) #checking if function inputs are compatible with method

  #### Extracting all taxa observations from DB ####
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
    my.request <- httr::GET(paste0("https://obs.infoflora.ch/rest/v4/observations?",
                                   "srid=2056&",
                                   "offset=",as.character((i-1)*4000),
                                   "&limit=4000&",
                                   "with_total=true&",
                                   "taxa=", taxa,
                                   "&taxa_matching=synonyms&",
                                   "taxa_relations_references=3,5,6,7,8"),
                            config(token = auth.code))
    
    my.request <- httr::content(my.request, as="text",encoding = "UTF-8") # convert as text 
    my.request <- jsonlite::fromJSON(my.request) 
    my.request <- my.request$data #Table of all available taxa in CH (and their IDs)
    
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
  
  #Filtering observations 
  my.sp.filter<-sp.filter(my_observations, 
                          with.geom = TRUE, 
                          exclude.prec = 2500, #parameters for xy radius max
                          presence.status = 2, #keep only presences 
                          introduction.status = c(0,1,2,3), #introduction status not 4
                          doubt.status = c(0,1), #eliminate error and doubtful observations
                          date.interval = c("1000-01-01","1994-12-31") #NOTE: 1994, last year of WS supplement 
  )
  taxa_obs <- my.sp.filter[!is.na(my.sp.filter$x),c("x","y","obs_id","date")]
  taxa_obs <- st_as_sf(taxa_obs, coords = c("x","y"), crs= 2056) #convert them to spatial object
  
  ####Handling WS Atlas####  
  ws_atlas<-httr::GET(paste0('https://obs.infoflora.ch/rest/v4/communities/0/public/atlases/27/names/',taxa,'/units?srid=4326&with-legend=true'))%>% 
    httr::content(as="text",encoding = "UTF-8") %>% 
    jsonlite::fromJSON() %>% 
    .$data #Table of all available taxa in CH (and their IDs)
  
  ws_atlas_grpd <- ws_atlas %>% group_by(year_min, legend) %>% #group by lengend and year to gather all of same frequency
    summarise(Ws_ids = list(gis_id))
  
  ws_freq_table <- subset(ws, ws$id %in% unique(ws_atlas$gis_id))#getting geometries for all sectors with WS INFO
  st_geometry(ws_freq_table) <- "geometry"
  
  ####creating 5X5 grid directly:####
  grid.extent<- st_sf(a = 1:2, geom = st_sfc(st_point(x = c(5.959163471,45.816873303)), #these coordinates are provided by infoflora API
                                             st_point(x = c(10.576523538944393,47.80822259322885))),
                      crs = 4326) %>%
    st_transform(st_crs(taxa_obs))
  grid <- st_coordinates(grid.extent) %>%
    vect %>%
    rast(res = 5000,crs = "epsg:2056") # create a 5000 m raster                        
  values(grid) <- ncell(grid):1 #create id_grid cells backwards to match that of sf_grid we create later
  #intersection ws and grid
  grid_ws_int <-  terra::extract(grid, vect(ws_freq_table) ,method="simple",touches=TRUE, exact=TRUE,xy=TRUE, ID=TRUE)
  grid_ws_int$ws_id <-  ws_freq_table$id[grid_ws_int$ID] #assigning ws_ids
  grid_ws_int <- st_as_sf(grid_ws_int, coords = c("x","y"), crs=2056) #matching squares and ws_ids
  grid_ws_int <- st_drop_geometry(grid_ws_int)
  ####Create DF with information for WS_prop calculation #####
  tt <- grid_ws_int %>% 
    group_by(ws_id) %>% 
    summarise(tsquares = sum(fraction)) #total squares calculation per sector
  
  grid_obs_int <- terra::extract(grid, vect(taxa_obs), touches=TRUE,xy=TRUE, ID=TRUE)  #finding the grid squares that have observations in them
  grid_obs_int <- grid_obs_int %>%   #grouping by each square
    group_by(lyr.1) %>% 
    summarise(N_obs_in_sqr = length(ID)) #number of observations in each square
  
  #create 1 df that has number of squares and the total proportion each square in each WS sector
  grid_ws_int <- merge(grid_ws_int, grid_obs_int, by.x="lyr.1", by.y="lyr.1", all.x=TRUE) #merging both intersection tables
  
  ####Calculation of area covered per WS ####
  #isolate all of the grids with obs in them with their proportions
  grid_obs_int<- subset(grid_ws_int, lyr.1 %in% grid_obs_int$lyr.1) #so all of these have presence of obs from DB
  sqrsws_obs <- grid_obs_int %>% 
    group_by(ws_id) %>% 
    summarise(propsqrs = sum(fraction)) #add the proportions of squares in ws and dividing by the number of squares in that ws
  sqrsws_obs <- merge(sqrsws_obs, tt, by.x="ws_id", by.y="ws_id", all.x=TRUE) #merge tt and new object to calculate squareswobs/totalsquares
  sqrsws_obs$wsprop <- sqrsws_obs$propsqrs/sqrsws_obs$tsquares
  
  ####Finding the median across both categories:frequent,rare####
  #want to calculate the median frequencies of the sectors for each classification of in the WS atlas:
  freq_ws <- subset(sqrsws_obs, ws_id %in% ws_atlas$gis_id[ws_atlas$legend==3]) #for all the sectors in atlas with legend == 3
  medf_wsprop <- median(freq_ws$wsprop)
  
  rare_ws <- subset(sqrsws_obs, ws_id %in% ws_atlas$gis_id[ws_atlas$legend!=3])#repeat for rare sectors
  medr_wsprop <- median(rare_ws$wsprop)
  ####Actually assigning a wsprop to all squares on grid####
  GRID_PROP_TAXA<- grid_ws_int  
  
  ws_int_obs <- grid_ws_int %>%  #calculating the number of observations in each WS with atlas information
    group_by(ws_id) %>% #doing this to find which sectors don't have any observations
    summarise(nobs_ws = sum(N_obs_in_sqr, na.rm = TRUE))
  
  ws_wo_obs <- ws_int_obs$ws_id[ws_int_obs$nobs_ws == 0] #all sectors without obs in them
  GRID_PROP_TAXA$wsprop <- 0
  GRID_PROP_TAXA2 <- merge(GRID_PROP_TAXA, ws_atlas[,c("gis_id","legend")], by.y="gis_id", by.x="ws_id", all.x=TRUE)
  #Method were only squares where ws_atlas info is present obtain a wsprop
  for(i in 1:nrow(GRID_PROP_TAXA2)){
    if(is.na(GRID_PROP_TAXA2$legend[i])){GRID_PROP_TAXA2$wsprop[i] <- 0} #CAS pas d'info wsatlas
    #else if(!is.na(GRID_PROP_TAXA2$N_obs_in_sqr[i])){GRID_PROP_TAXA2$wsprop[i] <- GRID_PROP_TAXA$fraction[i]}
    else{
      if(GRID_PROP_TAXA2$legend[i]==3){GRID_PROP_TAXA2$wsprop[i] <- medf_wsprop*GRID_PROP_TAXA2$fraction[i]}
      else{GRID_PROP_TAXA2$wsprop[i] <- medr_wsprop*GRID_PROP_TAXA2$fraction[i]}
    }
  }
  
  ####SUMMING TO GET UNIQUE SQUARES, and sythesising data####
  GRID_PROP_TAXA <- GRID_PROP_TAXA2 %>% #grid cells not unique in GRID_PROP_taxa2 since can touch multiple 
    group_by(lyr.1) %>% 
    summarise(wsprop = sum(wsprop, na.rm = TRUE))
  ####creating final 5x5 grid####
  lower_left <- data.frame(xmin=5.959163471, ymin = 45.816873303)
  lower_left_sf <- st_as_sf(lower_left, coords = c('xmin', 'ymin'), crs = 4326) %>% 
    st_transform(crs = 2056)
  grid2 <- st_make_grid(SWISS, cellsize = c(5000, 5000), what = "polygons",offset = st_bbox(lower_left_sf)[c("xmin", "ymin")])
  grid2 <- st_as_sf(data.frame(id_grid=1:length(grid2)),geom=grid2) %>% 
    st_transform(crs=2056)
  
  #aligning the id grid of the grid w that of raster
  centroidgrid2 <- st_centroid(grid2)
  ovrlcentro <- terra::extract(grid, vect(centroidgrid2))
  GRID_PROP_TAXA$lyr.1 <- ovrlcentro$ID[match(GRID_PROP_TAXA$lyr.1,ovrlcentro$lyr.1)] 
  mgrid <- merge(grid2, GRID_PROP_TAXA, by.x="id_grid", by.y="lyr.1", all.x=TRUE)

  ####wsprop for model calibration without injecting into the squares which don't have observations in associated weltens ####
  
  GRID_PROP_TAXA2$wsprop[GRID_PROP_TAXA2$ws_id%in%ws_wo_obs] <-NA
  
  GRID_PROP_TAXA_MODEL <- GRID_PROP_TAXA2 %>% 
    group_by(lyr.1, .drop=FALSE) %>% 
    summarise(wsprop = if(all(is.na(wsprop))) NA else sum(wsprop, na.rm = TRUE)) # Handle NA group explicitly, squares in weltens without observations)
  
  GRID_PROP_TAXA_MODEL$lyr.1 <- ovrlcentro$ID[match(GRID_PROP_TAXA_MODEL$lyr.1,ovrlcentro$lyr.1)]
  
  mgrid_w_predws <- merge(mgrid, GRID_PROP_TAXA_MODEL, by.x="id_grid", by.y="lyr.1", all.x=TRUE)
  names(mgrid_w_predws) <- c("id_grid", "wsprop_pred", "wsprop_model", "geom")
  
  #turning all sqaures outside swiss boundaries to 0 not na:
  mgrid_w_predws$wsprop_pred[is.na(mgrid_w_predws$wsprop_pred) & is.na(mgrid_w_predws$wsprop_model)] <- 0 #making only the grids which have to be NA, NA
  mgrid_w_predws$wsprop_model[(mgrid_w_predws$wsprop_pred)==0 & is.na(mgrid_w_predws$wsprop_model)] <- 0 

  return(mgrid_w_predws)
}




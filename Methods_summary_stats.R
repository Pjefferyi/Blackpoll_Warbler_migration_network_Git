# Thesis methods details 

# load necessary libraries 

# Load the helper functions script  
source("C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_Warbler_migration_network_Git/Geolocator_data_analysis_scripts/Geolocator_analysis_helper_functions.R")

# Load geolocator reference data 
ref_path <- "C:/Users/Jelan/OneDrive/Desktop/University/University of Guelph/Thesis/Blackpoll_data/Geolocator_reference_data_consolidated.csv"
ref_data <- read.csv(ref_path)

# location data (Only for geolocator used in the analysis )
geo.list <- c("V8757_010",
              "V8296_004",
              "V8296_005",
              "V8296_006",
              "V8757_055",
              "V8757_018",
              "V8757_021",
              "V8296_015",
              "V8296_017",
              "V8296_026",
              "V8296_025",
              "V8296_007",
              "V8296_008",
              "V8757_019",
              "V8757_096",
              "V8757_134",
              "V8757_029",
              "V8757_078",
              "V7638_001",
              "V7638_005",
              "V7638_009",
              "V7638_010",
              "V7638_011",
              "blpw09",
              "blpw12",
              "3254_001",
              "4068_014",
              "blpw14",
              "3254_003",
              "3254_008",
              "3254_011",
              "3254_057",
              "blpw15",
              "blpw25",
              "4105_008",
              "4105_009",
              "4105_016",
              "4105_017",
              "4210_002",
              "4210_004",
              "4210_006",
              "4210_010",
              "WRMA04173",
              "A",
              "B",
              "C",
              #"E",
              "D")

geo.all <- findLocData(geo.ids = geo.list, check_col_length = F)

# filter reference data to include only geolocators used in the analysis 
analysis_ref <- ref_data %>% dplyr::filter(geo.id %in% unique(geo.all$geo_id))

# mean Zenith angle measured using in-habitat calibration
IH_zentih <- mean(analysis_ref$In_habitat_median_zenith_angle, na.rm = T)

# mean zenith angle measured using Hill-Ekstrom Calibration 
He_zenith <- mean(analysis_ref$Hill_Ekstrom_median_angle, na.rm = T)

#mean difference between the Ih and HE zenith
zenith_diff <- mean(analysis_ref$In_habitat_median_zenith_angle - analysis_ref$Hill_Ekstrom_median_angle, na.rm = T)


hist(analysis_ref$In_habitat_median_zenith_angle - analysis_ref$Hill_Ekstrom_median_angle, breaks = 20)

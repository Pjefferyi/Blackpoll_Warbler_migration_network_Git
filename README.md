# Blackpoll warbler migration network project

This repository contains code  used as part of the analysis described in the following manuscript: 

Duali, J., DeLuca, W.V., Mackenzie, S.A. et al. Range-wide post- and pre-breeding migratory networks of a declining neotropical–nearctic migratory bird, the blackpoll warbler. *Sci Rep* 14, 30229 (2024). https://doi.org/10.1038/s41598-024-80838-9

Using this code requires light-level measurements from geolocators deployed on blackpoll warblers. This data is publicly available on Movebank, on pages with the following ID numbers:

- 126313959 (https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study126313959)
	
	Includes data collected as part of this study: DeLuca WV, Woodworth BK, Rimmer CC, Marra PP, Taylor PD, McFarland KP, Mackenzie SA, Norris DR. 2016. Data from: Transoceanic migration 	by a 12 g songbird. Movebank Data Repository. https://royalsocietypublishing-org.subzero.lib.uoguelph.ca/doi/10.1098/rsbl.2014.1045

- 959756713 (https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study959756713)

	Includes data collected as part of this study: DeLuca, WV, BK Woodworth, SA Mackenzie, AEM Newman, HA Cooke, LM Phillips, NE Freeman, AO Sutton, L Tauzer, C McIntyre, IJ Stenhouse, S 	Weidensaul, PD Taylor, DR Norris. 2019. A boreal songbirds's 20,00 km migration across North America and the Atlantic Ocean. https://esajournals-onlinelibrary-wiley-com.subzero.lib.uoguelph.ca/doi/10.1002/ecy.2651

- 4715484202 (https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study4715484202)
	
	Includes data from geolocators deployed after 2018.

Other sources of data used in this analysis include:

- Distribution polygons for the blackpoll warbler, used for plotting and for the spatial mask of the geolocator analysis. These were obtained via request to BirdLife International 	(https://datazone.birdlife.org/species/factsheet/blackpoll-warbler-setophaga-striata/distribution). This polygons can be substituted by range polygons from the eBird Status and 	Trends dataset, which are available upon request with the raster data describe below. 

- Raster of the blackpoll warbler's abundance across its breeding and nonbreeding range, which were obtained from the eBird Status and Trends database and are publicly available from ebird upon request (https://science.ebird.org/en/status-and-trends/species/bkpwar/downloads?week=1)

## Contact information 

If you would like to use the code in this repository but encounter some difficulties, you can reach out to me (the corresponding author of the manuscript) at jelanyduali@hotmail.com or jduali@uoguelph.ca


## Running this analysis 

You'll need to perform these steps in order

1. Download the geolocator data from eBird. Select the raw data (i.e., light level measurements)

2. Clone this repository in Rstudio 

4. Using the geolocator analysis scripts, estimate the locations for each birds from the geolocator data. You will need to modify the geolocator analysis script so they can read the raw light level data you have downloaded and ouput the location data in a desired location on your PC. All other scripts will need no modifications to run as long as the location data is generated, stored in "Network Construction", and named "All.locations.csv". 

5. Use the "Network_construction.R"" script to build the networks 

6. Use the "strenght of Strength_of_migratory_connectivity.R"" script to estimate the MC metric

7. Run the "blackpoll_warbler_project_figures.R"" script to calculate the network metrics and generate the figures

## Index of folders and files containing the code used for different parts of the analysis: 

### Geolocator analysis 

See this folder for related scripts: https://github-com.subzero.lib.uoguelph.ca/Pjefferyi/Blackpoll_Warbler_migration_network_Git/tree/main/Geolocator_data_analysis_scripts

The folder contains the scripts used to extract location data from each geolocator using the TwGeos and SGAT packages. There is one script per bird.

The scripts share their names with the ID of the birds/geolocators as shown in movebank and the manuscript.  

The script "Geolocator_analysis_helper_functions.R" is also located in the folder and contains functions used during the geolocator analysis as well as in other operations (notably to retriveve location data) 

### Network construction 

See this folder for related scripts (Network_construction.R): https://github-com.subzero.lib.uoguelph.ca/Pjefferyi/Blackpoll_Warbler_migration_network_Git/tree/main/Network_construction

The folder Contains the script that combines the location data derived from the geolocators to build the post- and pre-breeding migration network (Network_construction.R)

The data generated during the construction process is also stored in the Network Construction folder. This includes edge lists for the networks constructed for the post- and pre-breeding migrations and estimates of relative abundance used to weight each nodes.  \
The script includes code for a nonbreeidng movement network, but this network was not used in the analysis. 

### Strength of migratory connectivity metric 

The strength of migratory connectivity was calculated in this script (Strength_of_migratory_connectivity.R): https://github-com.subzero.lib.uoguelph.ca/Pjefferyi/Blackpoll_Warbler_migration_network_Git/blob/main/Strength_of_migratory_connectivity.R

### Figures, network metric calculation (betweenness centrality, time-adjusted weight, community analysis)

This script was used to generate the paper's figure and to calculated network metrics (blackpoll_warbler_project_figures.R): https://github-com.subzero.lib.uoguelph.ca/Pjefferyi/Blackpoll_Warbler_migration_network_Git/blob/main/blackpoll_warbler_project_figures.R 


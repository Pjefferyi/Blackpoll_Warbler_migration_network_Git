# Blackpoll warbler migration network project

This repository contains code that used as part of the analysis described in the following manuscript: 

Duali J., Deluca, W. V. Mackenzie, S. A., Tremblay, J. A., Drolet, B., Haché, S., Roberto-Charron, A., Ruiz, M. A., Boardman, R., Cooke, H. A., Rimmer, C. C., McFarland, K. P., Marra, P. P., Taylor, P. D., and Norris, D. R. Range-wide post- and pre-breeding migratory networks of a declining Neotropical-Nearctic migratory bird, the blackpoll warbler. In review.


Using this code requires light-level measurements from geolocators deployed on blackpoll warblers. This data is publicly available on Movebank, on pages with the following ID numbers:

- 126313959 (https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study126313959)
	
	Includes data collected as part of this study: DeLuca WV, Woodworth BK, Rimmer CC, Marra PP, Taylor PD, McFarland KP, Mackenzie SA, Norris DR. 2016. Data from: Transoceanic migration 	by a 12 g songbird. Movebank Data Repository. https://www.doi.org/10.5441/001/1.jb182ng4

- 959756713 (https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study959756713)

	Includes data collected as part of this study: DeLuca, WV, BK Woodworth, SA Mackenzie, AEM Newman, HA Cooke, LM Phillips, NE Freeman, AO Sutton, L Tauzer, C McIntyre, IJ Stenhouse, S 	Weidensaul, PD Taylor, DR Norris. 2019. A boreal songbirds's 20,00 km migration across North America and the Atlantic Ocean

- 4715484202 (https://www.movebank.org/cms/webapp?gwt_fragment=page=studies,path=study4715484202)
	
	Includes data from geolocators deployed after 2018.

Other sources of data used in this analysis include:

- Distribution polygons for the blackpoll warbler, used for plotting and for the spatial mask of the geolocator analysis. These were obtained via request to BirdLife International 	(https://datazone.birdlife.org/species/factsheet/blackpoll-warbler-setophaga-striata/distribution). This polygons can be substituted by range polygons from the eBird Status and 	Trends dataset, which are also available upon request with the raster data describe below. 

- Raster's of the blackpoll warbler's abundance across its breeding and nonbreeding range, which were obtained from the eBird Status and Trends database and are publicly available upon 	request (https://science.ebird.org/en/status-and-trends/species/bkpwar/downloads?week=1)


## Index of folders and files containing the code used for different parts of the analysis: 

### Geolocator analysis 

See this folder for related scripts: https://github.com/Pjefferyi/Blackpoll_Warbler_migration_network_Git/tree/main/Geolocator_data_analysis_scripts

The folder contains the scripts used to extract location data from each geolocator using the TwGeos and SGAT packages.

The script "Geolocator_analysis_helper_functions.R" is also located in the folder and contains functions used during the geolocator analysis as well as in other operations (notably to retriveve location data) 

### Network construction 

See this folder for related scripts: https://github.com/Pjefferyi/Blackpoll_Warbler_migration_network_Git/tree/main/Network_construction

The folder Contains the script that combines the location data derived from the geolocators to build the post- and pre-breeding migration network (Network_construction.R)

### Strenght of migratory connectivity metric 

The strength of migratory connectivity was calculated in this script: 

### Figures, network metric calculation (betweenness centrality, time-adjusted weight, Community analysis)

This script was used to generate the paper's figure and to calculated network metrics: 


#' The MASTER DICE DATA SET
#'
#' diceData is a dataset structured to organize global Influenza Like Illness (ILI) data on a multitude of spatial scales.  
#'
#' \describe{
#' 
#'The following is a brief description of each list entry:
#'   \item{$attr}{A dataframe that associates each population with a continent, country, region, etc.  Also includes latitude/longitude coordinates for the population-density-weighted centroid for that region. For example the United States has the following na}
#'   \item{$ili}{A dataframe that contains Google Flu Trends (GFT) \%ILI data for flu seasons starting in 2003-2014. Data is weekly and covers the United States at the national, regional, and state levels. Each column corresponds to a specific nation/region/state. Each row corresponds to a specific week.}
#'   \item{$sh}{A dataframe that contains Specific Humidity (SH) data for each region and week. Hourly/Daily SH data on a spatial grid has been averaged spatially and temporally to give a single value for each region-week. Each column corresponds to a specific nation/region/state. Each row corresponds to a specific week. SH is in units kg/kg.}
#'   \item{$school}{A dataframe with approximated school schedules for each region and week. Values range from 0 to 1. 0 indicates all schools are in-session for the entire week. 1 indicates that all schools are out-of-session for the entire week. Each column corresponds to a specific nation/region/state. Each row corresponds to a specific week.}
#'   \item{$pop}{A dataframe with census-reported populations by year. Each column corresponds to a specific nation/region/state. Each row corresponds to a specific year. Values are number of residents.}
#'   \item{$CDCili}{A dataframe containing Centers for Disease Control (CDC) \%ILI data. Currently contains data for flu seasons starting 2003-2015. Each column corresponds to a specific nation/region/state. Each row corresponds to a specific week. Note: when attempting to access the current flu season, get.DICE.data()/get.subset() will first attempt to get CDC data directly from the CDC server.}
#'   \item{$CDCbaseline}{A dataframe containing CDC onset-values for each region and year. CDC onset values provide a quantitative method for determining when a flu season begins and ends. Weeks with \%ILI > CDCbaseline are considered `flu weeks' that are part of the flu season.}
#' }
#' An important aspect of the datastructure organization is the idea of spatial `levels'. The levels are defined as: 0-Global, 1-Continent, 2-Country, 3-Region, 4-State, 5-County, 6-City in the United States. Levels smaller than `country' may be defined differently in each country.
#' 
#' In diceData$attr, the NAME_X fields indicate the association at level X. For example the United States has the following values: level=2, NAME_1='North.America', NAME_2='United.States', NAME_3='', NAME_4='', NAME_5='', NAME_6=''. An empty character '' indicates a NAME that is not applicable. Another example is San Diego: level=6, NAME_1='North.America', NAME_2='United.States', NAME_3='Region9', NAME_4='California', NAME_5='San.Diego.County', NAME_6='San.Diego'
#' 
#' @examples 
#' # Load dice
#' require(DICE)
#' # Load the dataset
#' data(DICE_dataset)
#' 
#' # to view the attributes for CDC Region 1
#' diceData$attr[diceData$attr$level==3 & diceData$attr$NAME_3=='Region1',]
#' 
#' # to view the first 20 weeks of GFT ILI data for all states in CDC Region 1
#' NamesInd = diceData$attr$NAME_4[diceData$attr$level==4 & diceData$attr$NAME_3=='Region1']
#' NamesInd
#' diceData$ili[1:20,NamesInd]
#' 
"diceData"
#' 
#' #' Gaussian Fit Parameters for Posterior Distributions
#' #' 
#' #' Results of Gaussian fits (mean and standard deviation) for each of the ten HHS regions using each of the five models
#' #' supported by \code{dice}. We only fit the basic reproduction number \eqn{R_0} and the percent clinical \eqn{pC}.  
#' #' In both cases we do a Gaussian fit to the posterior distribution.
#' #' This is done for each season, for each region and for each model.
#' #' These are needed if the MCMC procedure is done using a prior.  
#' #' The mean and standard deviation are packed into a named list with two elements:
#' 
#' #' \describe{
#' #'   \item{mean}{The mean value}
#' #'   \item{sigma}{The standard deviation}
#' #' }
#' #' In each case there are 550 rows and 5 columns and each entry in the fourth and fifth column defines the \emph{mean} or \emph{sigma} for \eqn{R0} and \eqn{pC} given a specific 
#' #' set of: 
#' #' \describe{
#' #'   \item{Region}{Starting year for the prior flu season}
#' #'   \item{year.start}{Starting year for the prior flu season}
#' #'   \item{model}{model number 1-5}
#' #' }
#' 
#' "postFit"

#' The GADM_2.8 level 1 spatial polygon shape file for the USA
#'
#' The file was downloaded from GADM using the getData command:
#' port1 = suppressMessages(getData('GADM', country = 'USA', level = 1))
#' and then converted to 'rda' format so it can be loaded with the package when LazyData is set to yes
#' The port1 object has 13 entries for each polygon
#' \describe{
#'  \item{OBJECTID}{Integer, polygon number}
#'  \item{ID_0}{Integer, country ID number}
#'  \item{ISO}{String, Three letter ISO3 country code}
#'  \item{NAME_0}{String, Country Name}
#'  \item{ID_1}{Integer, level 1 of country-state for the USA}
#'  \item{NAME_1}{String, level 1 name, state for the USA}
#'  \item{HASC_1}{String, two letter country code, and two letter state code, eg. US.CA}
#'  \item{CCN_1}{ NA in the case of USA}
#'  \item{CCA_1}{'' in the case of USA}
#'  \item{TYPE_1}{what type is the level 1 - State or Federal District in the case of USA}
#'  \item{ENGTYPE_1}{what type is the level 1 - State or Federal District in the case of USA}
#'  \item{NL_NAME_1}{'' in the case of USA}
#'  \item{VARNAME_1}{String, another state name descriptor in the case of USA}

#'}

"port1"

#' Parameters coming from Gaussian fits to posterior distributions of fits to the CDC's ILI data
#'

#' Commute Rate by Country Based on Flight Volume Data
#'
#'
#' A data frame with 166 rows and 167 columns (first column should be ignored) denoting the normalized flight volume from origin (row) to all other countries,
#' Diagonal elements are \emph{NOT} included.  The three-letter ISO3 name convention is used to describe the countries.
#' @section Warning:
#' \bold{This Data Set is Not Currently Supported}

"flights"

#' Pearson Correlation Table for \% ILI 
#'
#' A numeric correlation table with 132 rows and 6 columns showing for each season and each HHS region the season that has the highest correlation with it.
#' This season is referred to as the 'prior season'
#' Correlation is also shown for the national data- which is denoted as region 11.
#' \describe{
#'   \item{Region}{Starting year for the prior flu season for regions 1-10 and national - 11}
#'   \item{year.start}{Starting year for the flu season}
#'   \item{year.end}{Ending year for the flu season}
#'   \item{year.start.prior}{Starting year for the prior flu season}
#'   \item{year.end.prior}{Ending year for the prior flu season}
#'   \item{Cor}{The Pearson correlation between the two seasons}
#'   ...
#' }

"corTable"

#' Pearson Correlation Table for \% ILI - New version
#'
#' A numeric correlation table with 132 rows and 6 columns showing for each season and each HHS region the season that has the highest correlation with it.
#' This season is referred to as the 'prior season'
#' Correlation is also shown for the national data- which is denoted as region 11.
#' \describe{
#'   \item{Region}{Region name,  RegionX and national - United.States}
#'   \item{start}{Starting year for the flu season}
#'   \item{end}{Ending year for the flu season}
#'   \item{R0-mean}{Mean value for R0}
#'   \item{R0-SD}{Standard deviation for R0}
#'   \item{pC-mean}{Mean value for percent clinical, pC}
#'   \item{pC-SD}{Standard deviation for percent clinical, pC}
#'   \item{SH-mean}{Mean value for SH term}
#'   \item{SH-SD}{Standard deviation for SH term}
#'   \item{SV-mean}{Mean value for SV term}
#'   \item{SV-SD}{Standard deviation for SV term}
#'   \item{Cor}{The Pearson correlation between the two seasons}
#'   ...
#' }

"priorTable"



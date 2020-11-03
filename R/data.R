#' Medical data from 35 patients
#'
#' A dataset containing three variables (creatinine clearance C; digoxin clearance D; urine flow U) from 35 patients.
#'
#' @format A data frame with 35 rows and 3 variables:
#' \describe{
#'   \item{C}{creatinine clearance, in ml/min/1.73m^2}
#'   \item{D}{digoxin clearance, in ml/min/1.73m^2}
#'   \item{U}{urine flow, in ml/min}
#' }
#' @source Edwards, D. (2012). Introduction to graphical modelling, Section 3.1.4, Springer Science & Business Media.
"med"

#' 2017 Korea presidential election data
#'
#' A dataset containing 14 variables, consists of the voting results earned by the top five candidates from 250 electoral districts in Korea.
#'
#' @format A data frame with 1250 rows and 14 variables:
#' \describe{
#'   \item{PrecinctCode}{250 precinct codes designated by the election committee (4 digits)}
#'   \item{CityCode}{250 city codes of administrative standard code management system (5 digits)}
#'   \item{Municipalities}{250 municipalities}
#'   \item{City}{17 metropolitan cities}
#'   \item{Name}{Symbols 1-5, corresponding to Moon Jae-in, Hong Jun-pyo, Ahn Cheol-soo, Yoo Seung-min, Shim Sang-jung}
#'   \item{PoliticalParty}{The political parties of the first three candidates correspond to rogressivism, conservatism and centrism respectively.}
#'   \item{CandidateName}{The names of the candidates}
#'   \item{FirstPlace}{No.1 candidate by city, county and district}
#'   \item{AveAge}{Average age of voters in 17 years: statistics on resident registration population of the Ministry of Government Administration and Home Affairs}
#'   \item{AveYearEdu}{Average number of years of education for voters}
#'   \item{AveHousePrice}{Average price per square meter in 17 years}
#'   \item{AveInsurance}{The average insurance premium for each city, county, district}
#'   \item{VoteRate}{Vote rate by candidate}
#'   \item{NumVote}{Number of votes by candidate}
#' }
#' @source \url{https://github.com/OhmyNews/2017-Election}
"ElecData"

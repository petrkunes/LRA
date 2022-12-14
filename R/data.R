#' Pollen counts from Černé jezero in Šumava Mts., Czech Republic
#'
#' @source Carter, V.A., Chiverrell, R.C., Clear, J.L., Kuosmanen, N., Moravcová, A., Svoboda, M., Svobodová-Svitavská, H., Leeuwen, V., Van Leeuwen, J., van der Knaap, W.O., & Kuneš, P. 2018. Quantitative palynology informing conservation ecology in the Bohemian/Bavarian Forests of Central Europe. Frontiers in Plant Science 8: 1–14.
#'
#' @format A data frames with columns as time windows and rows:
#' \describe{
#'  \item{Radius}{Site radius in m.}
#'  \item{Model}{model type: 1 = Prentice bog model, 2 = Sugita lake model.}
#'  \item{Abies...}{Pollen counts for each taxon.}
#' }
#' @examples
#' \dontrun{
#'  PC_Cerne
#' }
"PC_Cerne"

#' Pollen counts from Prášilské jezero in Šumava Mts., Czech Republic
#'
#' @source Carter, V.A., Chiverrell, R.C., Clear, J.L., Kuosmanen, N., Moravcová, A., Svoboda, M., Svobodová-Svitavská, H., Leeuwen, V., Van Leeuwen, J., van der Knaap, W.O., & Kuneš, P. 2018. Quantitative palynology informing conservation ecology in the Bohemian/Bavarian Forests of Central Europe. Frontiers in Plant Science 8: 1–14.
#' @source Carter, V.A., Moravcová, A., Chiverrell, R.C., Clear, J.L., Finsinger, W., Dreslerová, D., Halsall, K., & Kuneš, P. 2018. Holocene-scale fire dynamics of central European temperate spruce-beech forests. Quaternary Science Reviews 191: 15–30.
#' @format A data frames with columns as time windows and rows:
#' \describe{
#'  \item{Radius}{Site radius in m.}
#'  \item{Model}{model type: 1 = Prentice bog model, 2 = Sugita lake model.}
#'  \item{Abies...}{Pollen counts for each taxon.}
#' }
#' @examples
#' \dontrun{
#'  PC_Prasilske
#' }
"PC_Prasilske"

#' Pollen counts from Rybárenská slať in Šumava Mts., Czech Republic
#'
#' @source Carter, V.A., Chiverrell, R.C., Clear, J.L., Kuosmanen, N., Moravcová, A., Svoboda, M., Svobodová-Svitavská, H., Leeuwen, V., Van Leeuwen, J., van der Knaap, W.O., & Kuneš, P. 2018. Quantitative palynology informing conservation ecology in the Bohemian/Bavarian Forests of Central Europe. Frontiers in Plant Science 8: 1–14.
#' @source Svobodová, H., Soukupová, L., & Reille, M. 2002. Diversified development of mountain mires, Bohemian Forest, Central Europe, in the last 13,000 years. Quaternary International 91: 123–135.
#' @format A data frames with columns as time windows and rows:
#' \describe{
#'  \item{Radius}{Site radius in m.}
#'  \item{Model}{model type: 1 = Prentice bog model, 2 = Sugita lake model.}
#'  \item{Abies...}{Pollen counts for each taxon.}
#' }
#' @examples
#' \dontrun{
#'  PC_Rybarenska
#' }
"PC_Rybarenska"


#' Set of relative pollen productivity estimates for the Gaussian plume model
#'
#' @source Abraham, V., Oušková, V., & Kuneš, P. 2014. Present-day vegetation helps quantifying past land cover in selected regions of the Czech Republic. PLoS ONE 9: e100117.
#'
#' @format A data frame with columns:
#' \describe{
#'  \item{alpha}{Relative pollen productivity to Poaceae = 1.}
#'  \item{vg}{Fall speed of pollen in m/s.}
#' }
#' @examples
#' \dontrun{
#'  avg
#' }
"avg"

#' Set of error estimates (variance-covariance matrix) for the relative pollen productivity estimates for the Gaussian plume model.
#'
#' @source Abraham, V., Oušková, V., & Kuneš, P. 2014. Present-day vegetation helps quantifying past land cover in selected regions of the Czech Republic. PLoS ONE 9: e100117.
#'
#' @format A matrix with identical columns and rows. The diagonal represents error estimates. Other values are covariances, here empty (data not available).
#'
#' @examples
#' \dontrun{
#'  alvc
#' }
"alvc"

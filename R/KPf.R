#' Pollen dispersal-deposition coefficient
#'
#' Calculates the taxon-specific pollen dispersal-deposition coefficient required for \code{\link{REVEALS}} function.
#'
#' @param vg fall speed of pollen in m/s
#' @param u wind speed in m/s
#' @param Zmax maximum extent of the region in metres
#' @param model type of the deposition model for the type of sedimentation basin: "peatland" (Prentice, 1985), "lake" (Sugita, 1993)
#' @param dwm type of dispersal model used: "gpm neutral" (Prentice, 1985), "lsm unstable" - from \code{\link{DispersalFactorK}} function in \code{DISQOVER} package (Theuerkauf, et al. 2016)
#'
#'
#' @references Prentice, I.C. 1985. Pollen representation, source area, and basin size: Toward a unified theory of pollen analysis. Quaternary Research 23: 76–86.
#' @references Sugita, S. 1993. A Model of Pollen Source Area for an Entire Lake Surface. Quaternary Research 39: 239–244.
#' @references Sugita, S. 2007. Theory of quantitative reconstruction of vegetation I: pollen from large sites REVEALS regional vegetation composition. Holocene 17: 229–241.
#' @references Theuerkauf, M., Couwenberg, J., Kuparinen, A., & Liebscher, V. 2016. A matter of dispersal: REVEALSinR introduces state-of-the-art dispersal models to quantitative vegetation reconstruction. Vegetation History and Archaeobotany 25: 541–553.

#'
#'
#'
#' @export

KPf <- function(vg, u, Zmax, radius, model, dwm) {
  b <- 75.2 * vg / u

  #Prentice bog model

  if (dwm == "gpm neutral") {
    if (model == "peatland") {
      KP <-
        (exp(-1 * b * radius ^ 0.125) - exp(-1 * b * (Zmax) ^ 0.125))
    }
    if (model == "lake")
    {
      #Sugita lake model

      xa <- b * (Zmax - radius) ^ (1 / 8)
      xb <- b * (Zmax + radius) ^ (1 / 8)
      xc <- b * (radius + radius) ^ (1 / 8)

      KP <-
        (4 * pi * radius / (b ^ 8)) * (Igamma(8, xa) - Igamma(8, xb) + Igamma(8, xc)) ## NEED TO CHECK!

    }
  }
  if (dwm == "lsm unstable") {
    KP <- DispersalFactorK(vg, model, radius, dwm, Zmax)
  }

  return(KP)
}

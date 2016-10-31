########################################################
### This is a function for Landcover Reconstruction Algorithm (Sugita 2007a,b)
### Release version: 1 (30 October 2016)
### Version includes: REVEALS - calculation for single and multiple sites, bog and lake model
### Gaussian plume model, Lagrangian stochastic model from DISQOVER package
########################################################

### DEFINE Sugita's KP
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
      
      #KP<-(4*pi*radius/(b^8))*(pgamma(xa,8)-pgamma(xb,8)+pgamma(xc,8))
      KP <-
        (4 * pi * radius / (b ^ 8)) * (Igamma(8, xa) - Igamma(8, xb) + Igamma(8, xc)) ## NOT SURE IF THIS IS CORRECT FUNCTION - NEED CHECK!
      
    }
  }
  if (dwm == "lsm unstable") {
    KP <- DispersalFactorK(vg, model, radius, dwm, Zmax)
  }
  
  return(KP)
}


### MAIN REVEALS FUNCTION WITH ERROR ESTIMATES FOR MULTIPLE SITES

REVEALS <-
  function(pollen_counts, avg, alvc, u, Zmax, dwm) {
    # lst <- read.table(pollen_counts, sep = ",", colClasses = "character")
    if (is.data.frame(pollen_counts)) {maxk <- 1} else {maxk <- length(pollen_counts)}
    # maxk <- length(rownames(lst)) ## NO OF SITES
    # if (length(colnames(lst))>1) {maxk <- 1}
    ### READ TIME INTERVALS
    timeint <- c()
    if (maxk > 1) {
    for (k in 1:maxk) {
      
      polcount <- pollen_counts[[k]]
      timeint <- c(timeint, colnames(polcount))
      }
      
      timeint <- unique(timeint)
      } else {timeint <- colnames(pollen_counts)}
    
    
#     avg <-
#       read.table(
#         data_avg, row.names = 1, header = T, sep = ","
#       )
    # alvc <- read.table(file_alvc, sep = ",", row.names = 1, header = T)
    maxt <- length(timeint) # NO OF TIMESLICES
    Vmean <- matrix(0, nrow = length(rownames(avg)), ncol = maxt)
    rownames(Vmean) <- rownames(avg)
    colnames(Vmean) <- timeint
    varVmean <- matrix(0, nrow = length(rownames(avg)), ncol = maxt)
    rownames(varVmean) <- rownames(avg)
    colnames(varVmean) <- timeint
    covVmean <-
      matrix(0, nrow = length(rownames(avg)), ncol = length(rownames(avg)) * maxt)
    Vsites <- matrix(0, nrow = 1, ncol = maxt)
    colnames(Vsites) <- timeint
    REVEALS_V <- REVEALS_E <- list()
    
    for (k in 1:maxk) {
      
      if (maxk > 1) {polcount <- pollen_counts[[k]]} else {polcount <- pollen_counts}
      
      polcount[is.na(polcount)] <- 0
      matrix_V <- polcount[-1:-2,,drop=FALSE] * 0
      REVEALS_V[[k]] <- matrix_V
      REVEALS_E[[k]] <- matrix_V
    }
    
    
    #rownames(REVEALS_V) <- rownames(avg)
    #colnames(REVEALS_V) <- timeint
    
    
    for (kt in 1:maxt) {
      ### LOOP FOR TIMESLICES
      ### u_calc
      cat("Estimating time window ",timeint[kt],": ",kt, " out of ",maxt,"\n", sep =
            "")
      varqa <- matrix(0, nrow = length(rownames(avg)), ncol = maxk)
      covqa <-
        matrix(0, nrow = length(alvc), ncol = length(alvc) * maxk)
      na <- matrix(0, nrow = length(rownames(avg)), ncol = maxk)
      qa <- matrix(0, nrow = length(rownames(avg)), ncol = maxk)
      
      for (k in 1:maxk) {
        if (maxk > 1) {polcount <- pollen_counts[[k]]} else {polcount <- pollen_counts}
        polcount[is.na(polcount)] <- 0
        sumpol <- sum(polcount[-1:-2, timeint[kt]])
        if (sumpol > 0) {
          for (i in 1:(length(rownames(polcount)) - 2)) {
            # LOOP FOR SPECIES i
            na[i, k] <- polcount[i + 2, timeint[kt]] / avg[i, 1]
            qa[i, k] <- na[i, k] / sumpol #pollen percentage
          }
        }
      }
      for (k in 1:maxk) {
        if (maxk > 1) {polcount <- pollen_counts[[k]]} else {polcount <- pollen_counts}
        polcount[is.na(polcount)] <- 0
        sumpol <- sum(polcount[-1:-2, timeint[kt]])
        if (sumpol > 0) {
          q <- polcount[-1:-2, timeint[kt]] / sumpol
          for (i in 1:(length(rownames(polcount)) - 2)) {
            if (q[i] > 0) {
              varqa[i, k] <-
                (qa[i, k] * qa[i, k]) * (((1 - q[i]) / (sumpol * q[i])) + (alvc[i, i] /
                                                                             avg[i, 1] * avg[i, 1])) # NEED TO CHECK FOR INPUT DATA (VARCOV or SE?) SEEMS THAT INPUT DATA IS VARIANCE (NO SE)
            }
          }
        }
      }
      for (k in 1:maxk) {
        if (maxk > 1) {polcount <- pollen_counts[[k]]} else {polcount <- pollen_counts}
        polcount[is.na(polcount)] <- 0
        sumpol <- sum(polcount[-1:-2, timeint[kt]])
        if (sumpol > 0) {
          q <- polcount[-1:-2, timeint[kt]] / sumpol
          for (i in 1:(length(rownames(polcount)) - 2)) {
            for (j in 1:(length(rownames(polcount)) - 2)) {
              # LOOP FOR SPECIES j
              if (i != j) {
                covqa[i, (length(rownames(polcount)) - 2) * (k - 1) + j] <-
                  (qa[i, k] * qa[j, k]) * (-1 / sumpol + alvc[i, j] / (avg[i, 1] * avg[j, 1]))
              }
            }
          }
        }
        
      }
      
      
      
      
      ### V_calc
      varqak <-
        matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
      varsumqak <- matrix(0, nrow = 1, ncol = maxk)
      covqak <-
        matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
      qak <-
        matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
      sumqak <- matrix(0, nrow = 1, ncol = maxk)
      nak <-
        matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
      sumnak <- matrix(0, nrow = 1, ncol = maxk)
      
      
      for (k in 1:maxk) {
        if (maxk > 1) {polcount <- pollen_counts[[k]]} else {polcount <- pollen_counts}
        polcount[is.na(polcount)] <- 0
        sumpol <- sum(polcount[-1:-2, timeint[kt]])
        radius <- polcount[1, timeint[kt]]
        if (sumpol > 0) {
          if (polcount[2, timeint[kt]] == 1) {
            model <- "peatland"
          } else {
            model <- "lake"
          }
          
          
          for (i in 1:(length(rownames(polcount)) - 2)) {
            qak[i, k] <- qa[i, k] / KPf(avg[i, 2], u, Zmax, radius, model, dwm)
            nak[i, k] <-
              na[i, k] / KPf(avg[i, 2], u, Zmax, radius, model, dwm)
          }
          sumqak[k] <- sum(qak[, k])
          sumnak[k] <- sum(nak[, k])
        }
      }
      ### V
      
      V <-
        matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
      rownames(V) <- rownames(avg)
      if (maxk>1) {colnames(V) <- names(pollen_counts)} else {colnames(V) <- "V"}
      
      for (k in 1:maxk) {
        if (maxk > 1) {polcount <- pollen_counts[[k]]} else {polcount <- pollen_counts}
        polcount[is.na(polcount)] <- 0
        sumpol <- sum(polcount[-1:-2, timeint[kt]])
        if (sumpol > 0) {
          for (i in 1:(length(rownames(polcount)) - 2)) {
            V[i, k] <- qak[i, k] / sumqak[k]
          }
        }
      }
      
      ### var V
      varV <-
        matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
      rownames(varV) <- rownames(avg)
      if (maxk>1) {colnames(varV) <- names(pollen_counts)} else {colnames(varV) <- "varV"}
      
      for (k in 1:maxk) {
        if (maxk > 1) {polcount <- pollen_counts[[k]]} else {polcount <- pollen_counts}
        polcount[is.na(polcount)] <- 0
        sumpol <- sum(polcount[-1:-2, timeint[kt]])
        radius <- polcount[1, timeint[kt]]
        if (sumpol > 0) {
          if (polcount[2, timeint[kt]] == 1) {
            model <- "peatland"
          } else {
            model <- "lake"
          }
          
          
          for (i in 1:(length(rownames(polcount)) - 2)) {
            varqak[i, k] <-
              varqa[i, k] / (KPf(avg[i, 2], u, Zmax, radius, model, dwm) ^ 2)
          }
          ccc1 <- ccc2 <- 0
          for (i in 1:(length(rownames(polcount)) - 2)) {
            ccc1 <- ccc1 + varqak[i, k]
            for (j in 1:(length(rownames(polcount)) - 2)) {
              if (i != j) {
                ccc2 <-
                  ccc2 + covqa[i, (length(rownames(polcount)) - 2) * (k - 1) + j] / (KPf(avg[i, 2], u, Zmax, radius, model, dwm) * KPf(avg[j, 2], u, Zmax, radius, model, dwm))
              }
            }
          }
          varsumqak[1, k] <- ccc1 + ccc2
          ccc3 <- 0
          for (i in 1:(length(rownames(polcount)) - 2)) {
            for (j in 1:(length(rownames(polcount)) - 2)) {
              ccc3 <-
                ccc3 + covqa[i, (length(rownames(polcount)) - 2) * (k - 1) + j] / (KPf(avg[i, 2], u, Zmax, radius, model, dwm) * KPf(avg[j, 2], u, Zmax, radius, model, dwm))
            }
            covqak[i, k] <- ccc3
          }
          for (i in 1:(length(rownames(polcount)) - 2)) {
            if (qak[i, k] > 0) {
              ccc1 <- varqak[i, k] / (qak[i, k] * qak[i, k])
              ccc2 <- varsumqak[k] / (sumqak[k] * sumqak[k])
              ccc3 <- 2 * covqak[i, k] / (qak[i, k] * sumqak[k])
              ccc <- V[i, k] ^ 2
              varV[i, k] <- ccc * (ccc1 + ccc2 - ccc3)
            }
          }
        }
      }
      
      ### OUTPUT FILES FOR SINGLE SITES
      for (m in 1:maxk) {
        REVEALS_V[[m]][, kt] <- V[, m]
        REVEALS_E[[m]][, kt] <- sqrt(varV)[, m]
      }
      
      st <-
        which(colSums(V) > 0) # sites available with pollen counts for time period
      nst <- length(st) # number of sites available
      Vsites[1,kt] <- nst
      if (maxk > 1 && nst > 0) {
        ### meanV_calc
        cat("   Calculating mean V\n", sep = "")
        
        nka <-
          vector(mode = "numeric", length = length(rownames(avg)))
        sumnorg <- vector(mode = "numeric", length = 1000)
        sumnk <- vector(mode = "numeric", length = 1000)
        sumnka <- vector(mode = "numeric", length = 1000)
        Q <-
          vector(mode = "numeric", length = length(rownames(avg)))
        Qa <- matrix(0, nrow = length(rownames(avg)), ncol = 1000)
        Vx <- matrix(0, nrow = length(rownames(avg)), ncol = 1000)
        
        
        if (nst > 1) {
          btns <- matrix(0, nrow = 1000, ncol = nst)
          for (k in 1:1000) {
            btns[k,] <- sample(st, replace = TRUE)
          }
          
          for (bs in 1:1000) {
            nk <- vector(mode = "numeric", length = length(rownames(avg)))
            norg <-
              vector(mode = "numeric", length = length(rownames(avg)))
            for (k in 1:nst) {
              polcnt <- pollen_counts[[btns[bs, k]]]
              if (polcnt[2, timeint[kt]] == 1) {
                model <- "peatland"
              } else {
                model <- "lake"
              }
              radius <- polcnt[1, timeint[kt]]
              for (i in 1:length(rownames(avg))) {
                nk[i] <- nk[i] + polcnt[i + 2, kt] / KPf(avg[i, 2], u, Zmax, radius, model, dwm)
                norg[i] <- norg[i] + polcnt[i + 2, kt]
              }
            }
            for (i in 1:length(rownames(avg))) {
              nka[i] <- nk[i] / avg[i, 1]
            }
            for (i in 1:length(rownames(avg))) {
              sumnorg[bs] <- sumnorg[bs] + norg[i]
            }
            for (i in 1:length(rownames(avg))) {
              sumnk[bs] <- sumnk[bs] + nk[i]
              sumnka[bs] <- sumnka[bs] + nka[i]
            }
            if (sumnk[bs] > 0) {
              for (i in 1:length(rownames(avg))) {
                Q[i] <- nk[i] / sumnk[bs]
              }
            }
            for (i in 1:length(rownames(avg))) {
              Qa[i, bs] <- Q[i] / avg[i, 1]
            }
            sumQa <- 0
            for (i in 1:length(rownames(avg))) {
              sumQa <- sumQa + Qa[i, bs]
            }
            if (sumQa > 0) {
              for (i in 1:length(rownames(avg))) {
                Vx[i, bs] <- Qa[i, bs] / sumQa
              }
            }
          }
          
          ### mean_var_calc
          cat("   Estimating errors\n", sep = "")
          sumQa <-
            vector(mode = "numeric", length = length(rownames(avg)))
          
          total <- 0
          for (i in 1:length(rownames(avg))) {
            for (bs in 1:1000) {
              sumQa[i] <- sumQa[i] + Qa[i, bs]
              total <- total + Qa[i, bs]
            }
          }
          if (total > 0) {
            for (i in 1:length(rownames(avg))) {
              Vmean[i, kt] <- sumQa[i] / total
            }
          }
          
          for (i in 1:length(rownames(avg))) {
            sx <- sxx <- 0
            for (bs in 1:1000) {
              sx <- sx + Vx[i, bs]
              sxx <- sxx + Vx[i, bs] ^ 2
            }
            sss <- sxx - sx * sx / 1000
            if (sss >= 0) {
              varVmean[i, kt] <- sss / 999
            }
          }
          mx <- mean(sumnka)
          for (i in 1:length(rownames(avg))) {
            for (j in 1:length(rownames(avg))) {
              if (i != j) {
                if (mx > 0) {
                  covVmean[i, length(rownames(avg)) * (kt - 1) + j] <-
                    -Vmean[i, kt] * Vmean[j, kt] / mx
                  
                }
              }
            }
          }
          
        } else {
          Vmean[,kt] <- V[,st]
          varVmean[,kt] <- varV[,st]
        }
      }
    } ### END OF kt LOOP (TIMESLICES)
    if (maxk > 1) {
      return(
        list(
          V_sites = REVEALS_V, SE_sites = REVEALS_E, Mean_V = Vmean, Mean_SE = sqrt(varVmean), No_sites = Vsites
        )
      )
    } else {
      return(list(V_sites = REVEALS_V, SE_sites = REVEALS_E))
    }
  } ### END OF FUNCTION REVEALS



### libraries
if (!require("reshape")) {install.packages("reshape")}
if (!require("compositions")) {install.packages("compositions")}
if (!require("DISQOVER")) {download.file("http://disqover.botanik.uni-greifswald.de/wp-content/uploads/2016/04/DISQOVER_0.7-5.tar.gz", "boot")
install.packages("boot", repos = NULL, type = "source")}

if (!require("zipfR")) install.packages("zipfR", dependencies = "Depends")

require(DISQOVER)
require(zipfR)

cat("Welcome to REVEALS program - based on v 4.9 with calculation of error estimates for multiple sites\n\n")





### DEFINE Sugita's KP
KPf <- function(vg, u, Zmax, radius, model) {

	b <- 75.2 * vg/u

	#Prentice bog model
	if (model == "Prentice") {
		KP <- (exp(-1 * b * radius^0.125) - exp(-1 * b * (Zmax * 1000)^0.125))
	} else {

		#Sugita lake model
		
		xa <- b * (Zmax * 1000 - radius)^(1/8)
		xb <- b * (Zmax * 1000 + radius)^(1/8)
		xc <- b * (radius + radius)^(1/8)

		#KP<-(4*pi*radius/(b^8))*(pgamma(xa,8)-pgamma(xb,8)+pgamma(xc,8))
		KP <- (4 * pi * radius/(b^8)) * (Igamma(8, xa) - Igamma(8, xb) + Igamma(8, xc)) ## NOT SURE IF THIS IS CORRECT FUNCTION - NEED CHECK!

	}

	return(KP)
}

### MAIN REVEALS FUNCTION OLD

reveals_old <- function(file_name_list, file_name_avg, u, Zmax) {
	avg <- read.table(file_name_avg, row.names = 1, header = T, sep = ",")
	lst <- read.table(file_name_list)

	paldatasample <- read.table(as.character(lst[1, ]), sep = ",", row.names = 1, check.names = F, header = T)
	veg <- matrix(nrow = length(row.names(paldatasample)) - 2, ncol = length(paldatasample))
	rownames(veg) <- rownames(paldatasample)[3:length(rownames(paldatasample))]
	colnames(veg) <- colnames(paldatasample)

	## LOOP FOR TIMELAYERS
	for (w in 1:(length(paldatasample))) {

		## LOOP FOR ALL SITES
		allsites <- matrix(nrow = nrow(avg), ncol = length(rownames(lst)))
		allsitesprop <- matrix(nrow = nrow(avg), ncol = length(rownames(lst)))
		for (m in 1:(length(rownames(lst)))) {
			polcount <- read.table(as.character(lst[m, ]), sep = ",", row.names = 1, check.names = F, header = T)
			polcount[is.na(polcount)] <- 0

			## pollen sum
			sumv <- 0
			for (j in 1:(length(rownames(polcount)) - 2)) {
				vg <- avg[j, 2]
				radius <- polcount[1, w]
				if (polcount[2, w] == 1) {
					model <- "Prentice"
				} else {
					model <- "Sugita"
				}

				v <- polcount[j + 1, w]/(avg[j, 1] * KPf(vg, u, Zmax, radius, model))
				sumv <- sumv + v
			}

			## vegetation proportion for 1 species
			for (i in 1:(length(rownames(polcount)) - 2)) {
				vg <- avg[i, 2]
				radius <- polcount[1, w]
				if (polcount[2, w] == 1) {
					model <- "Prentice"
				} else {
					model <- "Sugita"
				}
				v1 <- polcount[i + 2, w]/(avg[i, 1] * KPf(vg, u, Zmax, radius, model))
				allsitesprop[i, m] <- v1/sumv
				allsites[i, m] <- v1
			}
		} ## END OF ALL SITES LOOP
		veg[, w] <- rowSums(allsites)/sum(rowSums(allsites))


	} ## END OF TIMELAYERS LOOP

	return(veg)

}


### REVEALS WITH ERROR ESTIMATES - WORKING ONLY FOR ONE LOCALITY

REVEALS <- function(file_name_list, file_name_avg, alvc, u, Zmax) {


	maxk <- length(file_name_list) ## NO OF SITES
	avg <- read.table(file_name_avg, row.names = 1, header = T, sep = ",")
	alvc <- read.table("alpha.varcov.csv", sep = ",", row.names = 1, header = T)
	lst <- read.table(file_name_list)
	maxt <- length(colnames(read.table(as.character(lst[1, ]), sep = ",", row.names = 1, check.names = F, header = T))) # NO OF TIMESLICES
	REVEALS_V <- read.table(as.character(lst[1, ]), sep = ",", row.names = 1, check.names = F, header = T)
	REVEALS_V[, ] <- 0
	REVEALS_V <- REVEALS_V[-1:-2, ]
	REVEALS_E <- REVEALS_V

	for (t in 1:maxt) {
		### u_calc
		
		varqa <- matrix(0, nrow = length(rownames(avg)), ncol = maxk)
		covqa <- matrix(0, nrow = length(alvc), ncol = length(alvc) * maxk)
		na <- matrix(0, nrow = length(rownames(avg)), ncol = maxk)
		qa <- matrix(0, nrow = length(rownames(avg)), ncol = maxk)

		for (k in 1:maxk) {
			polcount <- read.table(as.character(lst[k, ]), sep = ",", row.names = 1, check.names = F, header = T)
			polcount[is.na(polcount)] <- 0
			sumpol <- sum(polcount[-1:-2, t])
			if (sumpol > 0) {
				for (i in 1:(length(rownames(polcount)) - 2)) { # LOOP FOR SPECIES i
					na[i, k] <- polcount[i + 2, t]/avg[i, 1]
					qa[i, k] <- na[i, k]/sumpol #pollen percentage
				}
			}
		}
		for (k in 1:maxk) {
			polcount <- read.table(as.character(lst[k, ]), sep = ",", row.names = 1, check.names = F, header = T)
			polcount[is.na(polcount)] <- 0
			sumpol <- sum(polcount[-1:-2, t])
			if (sumpol > 0) {
				q <- polcount[-1:-2, t]/sumpol
				for (i in 1:(length(rownames(polcount)) - 2)) {

					if (q[i] > 0) {
						varqa[i, k] <- (qa[i, k] * qa[i, k]) * (((1 - q[i])/(sumpol * q[i])) + (alvc[i, i]/avg[i, 1] * avg[i, 1])) # NEED TO CHECK FOR INPUT DATA (VARCOV or SE?) SEEMS THAT INPUT DATA IS VARIANCE (NO SE)
					}
				}
			}
		}
		for (k in 1:maxk) {
			polcount <- read.table(as.character(lst[k, ]), sep = ",", row.names = 1, check.names = F, header = T)
			polcount[is.na(polcount)] <- 0
			sumpol <- sum(polcount[-1:-2, t])
			if (sumpol > 0) {
				q <- polcount[-1:-2, t]/sumpol
				for (i in 1:(length(rownames(polcount)) - 2)) {
					for (j in 1:(length(rownames(polcount)) - 2)) { # LOOP FOR SPECIES j
						if (i != j) {
							covqa[i, j + (k - 1)] <- (qa[i, k] * qa[j, k]) * (-1/sumpol + alvc[i, j]/(avg[i, 1] * avg[j, 1]))
						}
					}
				}
			}

		}




		### V_calc
		varqak <- matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
		varsumqak <- matrix(0, nrow = 1, ncol = maxk)
		covqak <- matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
		qak <- matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
		sumqak <- matrix(0, nrow = 1, ncol = maxk)
		nak <- matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
		sumnak <- matrix(0, nrow = 1, ncol = maxk)


		for (k in 1:maxk) {
			polcount <- read.table(as.character(lst[k, ]), sep = ",", row.names = 1, check.names = F, header = T)
			polcount[is.na(polcount)] <- 0
			sumpol <- sum(polcount[-1:-2, t])
			radius <- polcount[1, t]
			if (polcount[2, t] == 1) {
				model <- "Prentice"
			} else {
				model <- "Sugita"
			}

			if (sumpol > 0) {
				for (i in 1:(length(rownames(polcount)) - 2)) {
					qak[i, k] <- qa[i, k]/KPf(avg[i, 2], u, Zmax, radius, model)
					nak[i, k] <- na[i, k]/KPf(avg[i, 2], u, Zmax, radius, model)
				}
				sumqak[k] <- sum(qak[, k])
				sumnak[k] <- sum(nak[, k])
			}
		}
		### V
		
		V <- matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)

		for (k in 1:maxk) {
			polcount <- read.table(as.character(lst[k, ]), sep = ",", row.names = 1, check.names = F, header = T)
			polcount[is.na(polcount)] <- 0
			sumpol <- sum(polcount[-1:-2, t])
			if (sumpol > 0) {
				for (i in 1:(length(rownames(polcount)) - 2)) {
					V[i, k] <- qak[i, k]/sumqak[k]
				}
			}
		}

		### var V
		varV <- matrix(0, nrow = length(rownames(polcount)) - 2, ncol = maxk)
		for (k in 1:maxk) {
			polcount <- read.table(as.character(lst[k, ]), sep = ",", row.names = 1, check.names = F, header = T)
			polcount[is.na(polcount)] <- 0
			sumpol <- sum(polcount[-1:-2, t])
			if (polcount[2, t] == 1) {
				model <- "Prentice"
			} else {
				model <- "Sugita"
			}
			if (sumpol > 0) {

				for (i in 1:(length(rownames(polcount)) - 2)) {
					varqak[i, k] <- varqa[i, k]/(KPf(avg[i, 2], u, Zmax, radius, model)^2)
				}
				ccc1 <- ccc2 <- 0
				for (i in 1:(length(rownames(polcount)) - 2)) {
					ccc1 <- ccc1 + varqak[i, k]
					for (j in 1:(length(rownames(polcount)) - 2)) {
						if (i != j) {
							ccc2 <- ccc2 + covqa[i, j]/(KPf(avg[i, 2], u, Zmax, radius, model) * KPf(avg[j, 2], u, Zmax, radius, model))
						}
					}
				}
				varsumqak[1, k] <- ccc1 + ccc2
				ccc3 <- 0
				for (i in 1:(length(rownames(polcount)) - 2)) {
					for (j in 1:(length(rownames(polcount)) - 2)) {
						ccc3 <- ccc3 + covqa[i, j]/(KPf(avg[i, 2], u, Zmax, radius, model) * KPf(avg[j, 2], u, Zmax, radius, model))
					}
					covqak[i, k] <- ccc3
				}
				for (i in 1:(length(rownames(polcount)) - 2)) {
					if (qak[i, k] > 0) {
						ccc1 <- varqak[i, k]/(qak[i, k] * qak[i, k])
						ccc2 <- varsumqak[k]/(sumqak[k] * sumqak[k])
						ccc3 <- 2 * covqak[i, k]/(qak[i, k] * sumqak[k])
						ccc <- V[i, k]^2
						varV[i, k] <- ccc * (ccc1 + ccc2 - ccc3)
					}
				}
			}
		}

		### TEMPORARY OUTPUT FILES
		
		REVEALS_V[, t] <- V
		REVEALS_E[, t] <- sqrt(varV)

	} ### END OF t LOOP (TIMESLICES)
	return(list(REVEALS_V, REVEALS_E))
} ### END OF FUNCTION REVEALS








### DATA INPUT
## files should be in the same format as in LRA.REVEALS.v4.2.2
## first row includes ages/or depths /or sample names
## second row includes radius in m
## third row includes model type - 1 = Prentice, 2 = Sugita
## fourth row and onwards include taxa


### commands
library(zipfR)
file_name_avg <- "alpha.vg.csv"
file_name_list <- "filename_list.csv"
alvc <- read.table("alpha.varcov.csv", sep = ",", row.names = 1, header = T)
u <- 3
Zmax <- 100
REVEALS("filename_list.csv", "alpha.vg.csv", "alpha.varcov.csv", 3, 100)


rm(list = ls())
################################################################################ part 1 functions and transitions#####################################
# Gillespie fuction with Tau-leaping correction
gillesp <- function(start, ratefun, trans, pars, times = 0:50, deltaT) {
  tfn <- function(Y) {
    if (tevents > 1) {
      Z <- apply(as.matrix(Y, byrow = T), 2, sum)
    }
    if (tevents == 1) {
      Z <- Y
    }
    return(Z)
  }

  t0 <- times[1] ## set time to starting time
  ntimes <- length(times)
  X <- start ## set state to starting state
  res <- matrix(
    nrow = length(times), ncol = length(start),
    dimnames = list(times, names(start))
  )
  ## matrix for results
  for (ctr in 1:(ntimes - 1)) {
    res[ctr, ] <- X
    while (t0 < times[ctr + 1]) {
      rates <- ratefun(X, pars, t0)
      if (all(rates == 0)) break
      totrate <- sum(rates)
      elapsed <- deltaT
      tevents <- rpois(1, totrate * deltaT)
      t0 <- t0 + deltaT ## update time

      if (tevents > 0) {
        which.trans <-
          sample(1:nrow(trans), size = tevents, prob = rates, replace = TRUE) ## pick transition

        while (
          sum((X + tfn(trans[which.trans, ])) < 0) > 0) {
          tevents <- max(c(1, rpois(1, totrate * deltaT)))
          which.trans <-
            sample(1:nrow(trans), size = tevents, prob = rates, replace = TRUE)
        }

        X <- X + tfn(trans[which.trans, ])
      }
    }
  }
  cbind(times, res)
}


###
statenames <- c("Sc", " Shw", " Srh", " Sv", " Ec", " Ehw", " Erh", " Ev", " Ecd", " Ic1", " Ih1", " Ihw1", " Irh1", " Iv1", " Iccc1", " Inccc2", " Inh2", " Ic2", " Ih2", " Iccc2", " Fc", " Fh", " Rh", " Rccc", " R", " CC", " CD", "CP")
transnames <- c("ScToEc", "ScToInccc2", "ScToInh2", "ScToSv", "ShwToEhw", "ShwToSrh", "SrhToErh", "SrhToShw", "SvToEv", "SvToSc", "EcToIc1", "EcToEcd", "EcToEv", "EhwToIhw1", "EhwToErh", "ErhToIrh1", "ErhToEhw", "EvToIv1", "EvToEc", "EcdToIc1", "EcdToIh1", "Ic1ToIccc1", "Ic1ToIh1", "Ic1ToIc2", "Ic1ToR", "Ic1ToIh2", "Ic1ToFc", "Ic1ToIccc2", "Ic1ToIv1", "Ih1ToIc1", "Ih1ToIh2", "Ih1ToRh", "Ihw1ToIh2", "Ihw1ToFh", "Ihw1ToRh", "Irh1ToIhw1", "Iv1ToIc1", "Iccc1ToIc1", "Iccc1ToFh", "Iccc1ToIccc2", "Iccc1ToRccc", "Inccc2ToEc", "Inccc2ToSc", "Inh2ToEc", "Inh2ToSc", "Ic2ToIccc2", "Ic2ToIh2", "Ic2ToIFc", "Ic2ToR", "Ih2ToIc2", "Ih2ToFh", "Ih2ToR", "Ih2ToRh", "Iccc2ToIc2", "Iccc2ToFh", "Iccc2ToRccc", "FcToR", "FhToR", "RhToR", "RcccToR")
trans <- matrix(c
(
  -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
  1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 1, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, -1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0
),
byrow = TRUE, nrow = length(transnames), ncol = length(statenames), dimnames = list(transnames, statenames)
)
ratefun <- function(X, pars, time) {
  vals <- c(as.list(pars), as.list(X)) # X is the initial status
  rates <- with(vals, c(

    # those leave Sc
    ScToEc = (betaIc1 * Ic1 + betaIc2 * Ic2 + betaIh1 * Ih1 + betaIh2 * Ih2 + betaIccc1 * Iccc1 + betaIccc2 * Iccc2 + betaIrh1 * Irh1 + betaFc * Fc + betaFh * Fh) * Sc / N,
    ScToInccc2 = pnccc * (betaIc1 * Ic1 + betaIc2 * Ic2 + betaIh1 * Ih1 + betaIh2 * Ih2 + betaIccc1 * Iccc1 + betaIccc2 * Iccc2 + betaIrh1 * Irh1 + betaFc * Fc + betaFh * Fh) * Sc / N * alpha * leg * gamma1 * ccc1 * fcccmax2 * Ic1,
    ScToInh2 = pnh * (betaIc1 * Ic1 + betaIc2 * Ic2 + betaIh1 * Ih1 + betaIh2 * Ih2 + betaIccc1 * Iccc1 + betaIccc2 * Iccc2 + betaIrh1 * Irh1 + betaFc * Fc + betaFh * Fh) * Sc / N * alpha * leg * gamma1 * (1 - ccc1) * h1 * fhmax2 * Ic1,
    ScToSv = rhoV * (Ih1 + Ih2 + Ihw1) * (Sc / (Sc + Ec + Ic1)),
    # tgise keave Shw
    ShwToEhw = (psi * betaIc1 * Ic1 + psi * betaIc2 * Ic2 + betaIhwo1 * Ihw1 + betaIhw1 * Ih1 + betaIhw2 * Ih2 + betaIv1 * Iv1 + betaFhw * Fh + psi * betaFc * Fc) * Shw / N,
    ShwToSrh = rhoRh * Shw,
    # those leave Srh
    SrhToErh = (betaIc1 * Ic1 + betaIc2 * Ic2 + betaIrh1 * Irh1 + betaFc * Fc) * Srh / N,
    SrhToShw = rhoH * Srh,
    # those leave Sv
    SvToEv = (betaIhwo1 * Ihw1 + betaIhw1 * Ih1 + betaIhw2 * Ih2 + betaIv1 * Iv1 + betaFhw * Fh) * Sv / N,
    SvToSc = rhoRv * Sv,
    # those leave Ec
    EcToIc1 = alpha * Ec,
    EcToEcd = alphaCho * Ec,
    EcToEv = rhoV * (Ih1 + Ih2 + Ihw1) * (Ec / (Sc + Ec + Ic1)),
    # those leave Ehw
    EhwToIhw1 = alphaHw * Ehw,
    EhwToErh = rhoRh * Ehw,
    # those leave Erh
    ErhToIrh1 = alphaRh * Erh,
    ErhToEhw = rhoH * Erh,
    # those leave Ev
    EvToIv1 = alphaV * Ev,
    EvToEc = rhoRv * Ev,
    # those leave Ecd
    EcdToIc1 = gammaE1 * (1 - he1) * Ecd,
    EcdToIh1 = gammaE1 * he1 * Ecd,
    # those leave Ic1
    Ic1ToIccc1 = gamma1 * ccc1 * fcccmax1 * Ic1,
    Ic1ToIh1 = agu * els * gamma1 * (1 - ccc1) * h1 * fhmax1 * Ic1,
    Ic1ToIc2 = els * gamma1 * (1 - ccc1) * (1 - h1) * delta1 * Ic1,
    Ic1ToR = gamma1 * (1 - ccc1) * (1 - h1) * (1 - delta1) * Ic1,
    Ic1ToIh2 = leg * gamma1 * (1 - ccc1) * h1 * fhmax2 * Ic1,
    Ic1ToFc = agu * leg * gamma1 * (1 - ccc1) * (1 - h1) * delta1 * Ic1,
    Ic1ToIccc2 = leg * gamma1 * ccc1 * fcccmax2 * Ic1,
    Ic1ToIv1 = rhoV * (Ih1 + Ih2 + Ihw1) * (Ic1 / (Sc + Ec + Ic1)),
    # those leave Ih1
    Ih1ToIc1 = els * gammaH1 * phiH1 * Ih1,
    Ih1ToIh2 = agu * els * gammaH1 * (1 - phiH1) * deltaH1 * Ih1,
    Ih1ToRh = agu * els * gammaH1 * (1 - phiH1) * (1 - deltaH1) * Ih1,
    # those leave Ihw1
    Ihw1ToIh2 = gammaHw1 * hw1 * Ihw1,
    Ihw1ToFh = gammaHw1 * (1 - hw1) * deltaHw1 * Ihw1,
    Ihw1ToRh = gammaHw1 * (1 - hw1) * (1 - deltaHw1) * Ihw1,
    # those leave Irh1
    Irh1ToIhw1 = omegaRh * Irh1,
    # those leave Iv1
    Iv1ToIc1 = rhoRv * Iv1,
    # those leave Iccc1
    Iccc1ToIc1 = gammaCcc1 * phiCcc1 * Iccc1,
    Iccc1ToFh = mel * gammaCcc1 * (1 - phiCcc1) * deltaCcc1 * Iccc1,
    Iccc1ToIccc2 = gammaCcc1 * (1 - phiCcc1) * deltaCcc1 * Iccc1,
    Iccc1ToRccc = gammaCcc1 * (1 - phiCcc1) * (1 - deltaCcc1) * Iccc1,
    # those leave Inccc2
    Inccc2ToEc = gammaNccc * eppsilonNccc * Inccc2,
    Inccc2ToSc = gammaNccc * (1 - eppsilonNccc) * Inccc2,
    # those leave Inh2
    Inh2ToEc = gammaNh * eppsilonNh * Inh2,
    Inh2ToSc = gammaNh * (1 - eppsilonNh) * Inh2,
    # those leave Ic2
    Ic2ToIccc2 = els * gamma2 * ccc2 * fcccmax2 * Ic2,
    Ic2ToIh2 = els * gamma2 * (1 - ccc2) * h2 * fhmax2 * Ic2,
    Ic2ToIFc = els * gamma2 * (1 - ccc2) * (1 - h2) * delta2 * Ic2,
    Ic2ToR = els * gamma2 * (1 - ccc2) * (1 - h2) * (1 - delta2) * Ic2,
    # those leave Ih2
    Ih2ToIc2 = els * gammaH2 * phiH2 * Ih2,
    Ih2ToFh = gammaH2 * (1 - phiH2) * deltaH2 * deltaB2 * Ih2,
    Ih2ToR = gammaH2 * (1 - phiH2) * deltaH2 * (1 - deltaB2) * Ih2,
    Ih2ToRh = gammaH2 * (1 - phiH2) * (1 - deltaH2) * Ih2,
    # those leave Iccc2
    Iccc2ToIc2 = gammaCcc2 * phiCcc2 * Iccc2,
    Iccc2ToFh = gammaCcc2 * (1 - phiCcc2) * deltaCcc2 * Iccc2,
    Iccc2ToRccc = gammaCcc2 * (1 - phiCcc2) * (1 - deltaCcc2) * Iccc2,
    # those leave Fc;
    FcToR = gammaFc * Fc,
    # those leave Fh
    FhToR = gammaFh * Fh,
    # those leave Rh
    RhToR = gammaRh * Rh,
    # those leave Rccc
    RcccToR = gammaRccc * Rccc
  ))
}
################################################################################################################################################################################################################
################################################################################ part 2 pars of 37 models ######################################################################################################
################################################################################################################################################################################################################

################################ Agusto 2015
### use the parameters under Low-effectiveness level of the basic public health control strategy {###assume the worst scenario}
#### Pratmeters
pars.agu15 <- c(
  psi = 0, gammaHw1 = 0.5366,
  hw1 = 0,
  deltaHw1 = 1 - 0.42,
  agu = 0,
  leg = 1,
  els = 1,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  rhoRv = (1 / 10) * 24,
  rhoV = (1 / 10) * 24,
  rhoRh = (1 / 8) * 24,
  rhoH = (1 / 8) * 24,
  omegaRh = 0.5,
  alpha = 0.5239,
  alphaHw = 0.5239,
  alphaRh = 0.5239,
  alphaV = 0.5239,
  betaIc1 = 1.5 * 0.3045,
  betaIc2 = 1.5 * 0.3045,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 1 * 0.3045 * 1,
  betaIv1 = 1 * 0.3045,
  betaIhwo1 = 1 * 0.3045,
  betaIhw1 = 1 * 0.3045,
  betaIhw2 = 1 * 0.3045,
  betaFc = 1.5 * 0.3045 * 1,
  betaFh = 0,
  betaFhw = 0.1 * 1 * 0.3045 * 1,
  delta1 = 1,
  delta2 = 0.5366 * (1 - 0.48) / (0.5366 + 0.21),
  deltaH1 = 0,
  deltaH2 = 0.5366 * (1 - 0.42) / ((1 - 0.21) * 1 + 0.5366),
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  gamma1 = 0.5472,
  gammaH1 = 0,
  gamma2 = 0.5366 + 0.21,
  gammaH2 = (1 - 0.21) * 1 + 0.5366,
  phiH1 = 0,
  phiH2 = (1 - 0.21) * 1 / ((1 - 0.21) * 1 + 0.5366),
  phiCcc1 = 0,
  phiCcc2 = 0,
  h1 = 0.5,
  h2 = 0.21 / (0.5366 + 0.21),
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 0.5,
  gammaFh = 0.5,
  N = 10000
)
############################################################## parameters in Camacho 2014
# hospital transmission was based on the begining when there is no measrues taken, very high.


gammaH <- 1 / 3
gammaD <- 1 / 7.49
gammaI <- 1 / 10
gammaDH <- 1 / (1 / gammaD - 1 / gammaH)
gammaIH <- 1 / (1 / gammaI - 1 / gammaH)
delta1 <- 0.88
delta1o <- delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD)
delta2o <- delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH)
delta2 <- 0
deltaH1 <- 0
deltaH2 <- 0.88
deltaB2 <- 1
deltaCcc1 <- 0
deltaCcc2 <- 0
h1 <- 0.21
h2 <- 0
h1o <- h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH)
gamma1 <- gammaH * h1o / h1
gamma2 <- 0
gammaH1 <- 0
gammaH2 <- gammaDH * delta2o / delta1
pars.cam14 <- c(
  psi = 0, gammaHw1 = 0,
  hw1 = 0,
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 1 / 5.99,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.1 * (1 - 0.3 * 0.00005),
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 3.24 * (1 - 0),
  betaIh1 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 0.78 * (1 - 0.3 * 0.00005),
  betaFh = 0.78 * (1 - 0.3 * 0.00005),
  betaFhw = 0,
  gammaH = 1 / 3,
  gammaD = 1 / 7.49,
  gammaI = 1 / 10,
  gammaDH = 1 / (1 / gammaD - 1 / gammaH),
  gammaIH = 1 / (1 / gammaI - 1 / gammaH),
  delta1 = 0.88,
  delta1o = delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD),
  delta2o = delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH),
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0.88,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0.21,
  h2 = 0,
  h1o = h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH),
  gamma1 = gammaH * h1o / h1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = gammaDH * delta2o / delta1,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 1 / 0.99,
  gammaFh = 1 / 0.99,
  N = 10000
)
############################################################# parameters in Elisenberg 2014 full model SEIHFR
################################# data are based on All 3 countries together
#### parameters
pars.eli14_full_all <- c(
  psi = 0, gammaHw1 = 0,
  hw1 = 0,
  deltaHw1 = 0,
  agu = 1,
  leg = 0,
  els = 1,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  he1 = 0,
  gammaE1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 1 / 9,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.24,
  betaIc2 = 3 * 0.24,
  betaIh = 0.1,
  betaIh1 = 0.24 * 0.1,
  betaIh2 = 3 * 0.24 * 0.1,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 3 * 0.24,
  betaFh = 3 * 0.24,
  betaFhw = 0,
  delta1 = 0.9,
  delta2 = 0.95,
  deltaH1 = 0.6,
  deltaH2 = 0.8,
  deltaB2 = 0.1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1 / 10,
  gammaRccc = 0,
  gammaFc = 0.5,
  gammaFh = 0.5,
  gamma1o = 1 / 6,
  gamma2o = 1 / 2,
  h1 = 0.1667 / (1 / 6 + 0.1667),
  gamma1 = 0.1667 + 1 / 6,
  h2 = 0.75 / (1 / 2 + 0.75),
  gamma2 = 1 / 2 + 0.75,
  gammaH1 = 1 / 6,
  gammaH2 = 1 / 2,
  N = 10000
)
############################################################## parameters in Gomes 2014, SEIHFR
## this paper combines epidemic and mobility, here only consider compartment model, and beta are estimated based on the functions and parameters provided in the paper
psi <- 0
gammaHw1 <- 0
hw1 <- 0
deltaHw1 <- 0
agu <- 1
leg <- 1
els <- 0
mel <- 0
fhmax1 <- 1
fhmax2 <- 1
fcccmax1 <- 0
fcccmax2 <- 0
pnccc <- 0
pnh <- 0
alphaCho <- 0
gammaE1 <- 0
he1 <- 0
eppsilonNccc <- 0
eppsilonNh <- 0
ccc1 <- 0
ccc2 <- 0
phiH1 <- 0
phiH2 <- 0
phiCcc1 <- 0
phiCcc2 <- 0
rhoRv <- 0
rhoV <- 0
rhoRh <- 0
rhoH <- 0
omegaRh <- 0
alpha <- 1 / 7
alphaHw <- 0
alphaRh <- 0
alphaV <- 0
betaIc1 <- 0.146844
betaIc2 <- 0
betaIh1 <- 0
betaIh2 <- 0.1146355
betaIccc1 <- 0
betaIccc2 <- 0
betaIrh1 <- 0
betaIv1 <- 0
betaIhwo1 <- 0
betaIhw1 <- 0
betaIhw2 <- 0
betaFc <- 0.5454545
betaFh <- 0.5454545
betaFhw <- 0
gammaH <- 1 / 5
gammaD <- 1 / 9.6
gammaI <- 1 / 10
gammaDH <- 1 / (1 / gammaD - 1 / gammaH)
gammaIH <- 1 / (1 / gammaI - 1 / gammaH)
delta1 <- 0.55
delta1o <- delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD)
delta2o <- delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH)
delta2 <- 0
deltaH1 <- 0
deltaH2 <- 0.55
deltaB2 <- 1
deltaCcc1 <- 0
deltaCcc2 <- 0
h1 <- 0.8
h2 <- 0
h1o <- h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH)
gamma1 <- gammaH * h1o / h1
gamma2 <- 0
gammaH1 <- 0
gammaH2 <- gammaDH * delta2o / delta1
gammaCcc0 <- 0
gammaH0 <- 0
gammaCcc1 <- 0
gammaCcc2 <- 0
gammaNccc <- 0
gammaNh <- 0
gammaRh <- 1
gammaRccc <- 0
gammaFc <- 1 / 2
gammaFh <- 1 / 2
pars.gom14_SEIHFR <- c(
  psi = 0, gammaHw1 = 0,
  hw1 = 0,
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 1 / 7,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.146844,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0.1146355,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 0.5454545,
  betaFh = 0.5454545,
  betaFhw = 0,
  gammaH = 1 / 5,
  gammaD = 1 / 9.6,
  gammaI = 1 / 10,
  gammaDH = 1 / (1 / gammaD - 1 / gammaH),
  gammaIH = 1 / (1 / gammaI - 1 / gammaH),
  delta1 = 0.55,
  delta1o = delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD),
  delta2o = delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH),
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0.55,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0.8,
  h2 = 0,
  h1o = h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH),
  gamma1 = gammaH * h1o / h1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = gammaDH * delta2o / delta1,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 1 / 2,
  gammaFh = 1 / 2,
  N = 10000
)
############################################################## parameters in Legrand 2007
################################# legrand07.congo95
psi <- 0
gammaHw1 <- 0
hw1 <- 0
deltaHw1 <- 0
agu <- 1
leg <- 1
els <- 0
mel <- 0
fhmax1 <- 1
fhmax2 <- 1
fcccmax1 <- 0
fcccmax2 <- 0
pnccc <- 0
pnh <- 0
alphaCho <- 0
gammaE1 <- 0
he1 <- 0
eppsilonNccc <- 0
eppsilonNh <- 0
ccc1 <- 0
ccc2 <- 0
phiH1 <- 0
phiH2 <- 0
phiCcc1 <- 0
phiCcc2 <- 0
rhoRv <- 0
rhoV <- 0
rhoRh <- 0
rhoH <- 0
omegaRh <- 0
alpha <- 1 / 7
alphaHw <- 0
alphaRh <- 0
alphaV <- 0
betaIc1 <- 0.588 / 7
betaIc2 <- 0
betaIh1 <- 0
betaIh2 <- 0.794 / 7
betaIccc1 <- 0
betaIccc2 <- 0
betaIrh1 <- 0
betaIv1 <- 0
betaIhwo1 <- 0
betaIhw1 <- 0
betaIhw2 <- 0
betaFc <- 7.653 / 7
betaFh <- 7.653 / 7
betaFhw <- 0
gammaH <- 1 / 5
gammaD <- 1 / 9.6
gammaI <- 1 / 10
gammaDH <- 1 / (1 / gammaD - 1 / gammaH)
gammaIH <- 1 / (1 / gammaI - 1 / gammaH)
delta1 <- 0.81
delta1o <- delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD)
delta2o <- delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH)
delta2 <- 0
deltaH1 <- 0
deltaH2 <- 0.81
deltaB2 <- 1
deltaCcc1 <- 0
deltaCcc2 <- 0
h1 <- 0.8
h2 <- 0
h1o <- h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH)
gamma1 <- gammaH * h1o / h1
gamma2 <- 0
gammaH1 <- 0
gammaH2 <- gammaDH * delta2o / delta1
gammaCcc0 <- 0
gammaH0 <- 0
gammaCcc1 <- 0
gammaCcc2 <- 0
gammaNccc <- 0
gammaNh <- 0
gammaRh <- 1
gammaRccc <- 0
gammaFc <- 1 / 2
gammaFh <- 1 / 2
pars.leg07_con95 <- c(
  psi = 0, gammaHw1 = 0,
  hw1 = 0,
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 1 / 7,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.588 / 7,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0.794 / 7,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 7.653 / 7,
  betaFh = 7.653 / 7,
  betaFhw = 0,
  gammaH = 1 / 5,
  gammaD = 1 / 9.6,
  gammaI = 1 / 10,
  gammaDH = 1 / (1 / gammaD - 1 / gammaH),
  gammaIH = 1 / (1 / gammaI - 1 / gammaH),
  delta1 = 0.81,
  delta1o = delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD),
  delta2o = delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH),
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0.81,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0.8,
  h2 = 0,
  h1o = h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH),
  gamma1 = gammaH * h1o / h1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = gammaDH * delta2o / delta1,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 1 / 2,
  gammaFh = 1 / 2,
  N = 10000
)
################################# legrand07.uganda00
psi <- 0
gammaHw1 <- 0
hw1 <- 0
deltaHw1 <- 0
agu <- 1
leg <- 1
els <- 0
mel <- 0
fhmax1 <- 0
fhmax2 <- 1
fcccmax1 <- 0
fcccmax2 <- 0
pnccc <- 0
pnh <- 0
alphaCho <- 0
gammaE1 <- 0
he1 <- 0
eppsilonNccc <- 0
eppsilonNh <- 0
ccc1 <- 0
ccc2 <- 0
phiH1 <- 0
phiH2 <- 0
phiCcc1 <- 0
phiCcc2 <- 0
rhoRv <- 0
rhoV <- 0
rhoRh <- 0
rhoH <- 0
omegaRh <- 0
alpha <- 1 / 12
alphaHw <- 0
alphaRh <- 0
alphaV <- 0
betaIc1 <- 3.532 / 7
betaIc2 <- 0
betaIh1 <- 0
betaIh2 <- 0.012 / 7
betaIccc1 <- 0
betaIccc2 <- 0
betaIrh1 <- 0
betaIv1 <- 0
betaIhwo1 <- 0
betaIhw1 <- 0
betaIhw2 <- 0
betaFc <- 0.462 / 7
betaFh <- 0.462 / 7
betaFhw <- 0
gammaH <- 1 / 4.2
gammaD <- 1 / 8
gammaI <- 1 / 10
gammaDH <- 1 / (1 / gammaD - 1 / gammaH)
gammaIH <- 1 / (1 / gammaI - 1 / gammaH)
delta1 <- 0.53
delta1o <- delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD)
delta2o <- delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH)
delta2 <- 0
deltaH1 <- 0
deltaH2 <- 0.53
deltaB2 <- 1
deltaCcc1 <- 0
deltaCcc2 <- 0
h1 <- 0.8
h2 <- 0
h1o <- h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH)
gamma1 <- gammaH * h1o / h1
gamma2 <- 0
gammaH1 <- 0
gammaH2 <- gammaDH * delta2o / delta1
gammaCcc0 <- 0
gammaH0 <- 0
gammaCcc1 <- 0
gammaCcc2 <- 0
gammaNccc <- 0
gammaNh <- 0
gammaRh <- 1
gammaRccc <- 0
gammaFc <- 1 / 2
gammaFh <- 1 / 2
pars.leg07_uga00 <- c(
  psi = 0, gammaHw1 = 0,
  hw1 = 0,
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 1 / 12,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 3.532 / 7,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0.012 / 7,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 0.462 / 7,
  betaFh = 0.462 / 7,
  betaFhw = 0,
  gammaH = 1 / 4.2,
  gammaD = 1 / 8,
  gammaI = 1 / 10,
  gammaDH = 1 / (1 / gammaD - 1 / gammaH),
  gammaIH = 1 / (1 / gammaI - 1 / gammaH),
  delta1 = 0.53,
  delta1o = delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD),
  delta2o = delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH),
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0.53,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0.8,
  h2 = 0,
  h1o = h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH),
  gamma1 = gammaH * h1o / h1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = gammaDH * delta2o / delta1,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 1 / 2,
  gammaFh = 1 / 2,
  N = 10000
)
############################################################## parameters in River 2014
################################# River_Liberia 2014
## list single parameters
agu <- 1
leg <- 1
els <- 0
mel <- 0
fhmax1 <- 0
fhmax2 <- 1
fcccmax1 <- 0
fcccmax2 <- 0
pnccc <- 0
pnh <- 0
alphaCho <- 0
gammaE1 <- 0
he1 <- 0
eppsilonNccc <- 0
eppsilonNh <- 0
ccc1 <- 0
ccc2 <- 0
phiH1 <- 0
phiH2 <- 0
phiCcc1 <- 0
phiCcc2 <- 0
rhoRv <- 0
rhoV <- 0
rhoRh <- 0
rhoH <- 0
omegaRh <- 0
alphaCho <- 0
alpha <- 1 / 12
alphaHw <- 0
alphaRh <- 0
alphaV <- 0
betaIc1 <- 0.16
betaIc2 <- 0
betaIh1 <- 0
betaIh2 <- 0.062
betaIccc1 <- 0
betaIccc2 <- 0
betaIrh1 <- 0
betaIv1 <- 0
betaIhwo1 <- 0
betaIhw1 <- 0
betaIhw2 <- 0
betaFc <- 0.489
betaFh <- 0.489
betaFhw <- 0
gammaH <- 1 / 3.24
gammaD <- 1 / 13.31
gammaI <- 1 / 15
gammaDH <- 1 / (1 / gammaD - 1 / gammaH)
gammaIH <- 1 / (1 / gammaI - 1 / gammaH)
delta1 <- 0.5
delta1o <- delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD)
delta2o <- delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH)
delta2 <- 0
deltaH1 <- 0
deltaH2 <- 0.5
deltaB2 <- 1
deltaCcc1 <- 0
deltaCcc2 <- 0
h1 <- 0.197
h2 <- 0
h1o <- h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH)
gamma1 <- gammaH * h1o / h1
gamma2 <- 0
gammaH1 <- 0
gammaH2 <- gammaDH * delta2o / delta1
gammaCcc0 <- 0
gammaH0 <- 0
gammaCcc1 <- 0
gammaCcc2 <- 0
gammaNccc <- 0
gammaNh <- 0
gammaRh <- 1
gammaRccc <- 0
gammaFc <- 1 / 2.01
gammaFh <- 1 / 2.01
##### Parameter vector
pars.riv14_lib <- c(
  psi = 0, gammaHw1 = 0,
  hw1 = 0,
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alphaCho = 0,
  alpha = 1 / 12,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.16,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0.062,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0.489,
  betaFh = 0.48,
  betaFhw = 0,
  gammaH = 1 / 3.24,
  gammaD = 1 / 13.31,
  gammaI = 1 / 15,
  gammaDH = 1 / (1 / gammaD - 1 / gammaH),
  gammaIH = 1 / (1 / gammaI - 1 / gammaH),
  delta1 = 0.5,
  delta1o = delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD),
  delta2o = delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH),
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0.5,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0.197,
  h2 = 0,
  h1o = h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH),
  gamma1 = gammaH * h1o / h1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = gammaDH * delta2o / delta1,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 1 / 2.01,
  gammaFh = 1 / 2.01,
  N = 10000
)
################################# River_Sierra Leone 2014
## list
psi <- 0
gammaHw1 <- 0
hw1 <- 0
deltaHw1 <- 0
agu <- 1
leg <- 1
els <- 0
mel <- 0
fhmax1 <- 0
fhmax2 <- 1
fcccmax1 <- 0
fcccmax2 <- 0
pnccc <- 0
pnh <- 0
alphaCho <- 0
gammaE1 <- 0
he1 <- 0
eppsilonNccc <- 0
eppsilonNh <- 0
ccc1 <- 0
ccc2 <- 0
phiH1 <- 0
phiH2 <- 0
phiCcc1 <- 0
phiCcc2 <- 0
rhoRv <- 0
rhoV <- 0
rhoRh <- 0
rhoH <- 0
omegaRh <- 0
alphaCho <- 0
alpha <- 1 / 10
alphaHw <- 0
alphaRh <- 0
alphaV <- 0
betaIc1 <- 0.128
betaIc2 <- 0
betaIh1 <- 0
betaIh2 <- 0.080
betaIccc1 <- 0
betaIccc2 <- 0
betaIrh1 <- 0
betaIv1 <- 0
betaIhwo1 <- 0
betaIhw1 <- 0
betaIhw2 <- 0
betaFc <- 0.111
betaFh <- 0.111
betaFhw <- 0
gammaH <- 1 / 4.12
gammaD <- 1 / 10.38
gammaI <- 1 / 20
gammaDH <- 1 / (1 / gammaD - 1 / gammaH)
gammaIH <- 1 / (1 / gammaI - 1 / gammaH)
delta1 <- 0.75
# propotion of motality; originally writen as delta
delta1o <- delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD)
delta2o <- delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH)
delta2 <- 0
deltaH1 <- 0
deltaH2 <- 0.75
deltaB2 <- 1
deltaCcc1 <- 0
deltaCcc2 <- 0
h1 <- 0.197
h2 <- 0
h1o <- h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH)
gamma1 <- gammaH * h1o / h1
gamma2 <- 0
gammaH1 <- 0
gammaH2 <- gammaDH * delta2o / delta1
gammaCcc0 <- 0
gammaH0 <- 0
gammaCcc1 <- 0
gammaCcc2 <- 0
gammaNccc <- 0
gammaNh <- 0
gammaRh <- 1
gammaRccc <- 0
gammaFc <- 1 / 4.5
gammaFh <- 1 / 4.5
pars.riv14_sie <- c(
  psi = 0, gammaHw1 = 0,
  hw1 = 0,
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alphaCho = 0,
  alpha = 1 / 10,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.128,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0.080,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 0.111,
  betaFh = 0.111,
  betaFhw = 0,
  gammaH = 1 / 4.12,
  gammaD = 1 / 10.38,
  gammaI = 1 / 20,
  gammaDH = 1 / (1 / gammaD - 1 / gammaH),
  gammaIH = 1 / (1 / gammaI - 1 / gammaH),
  delta1 = 0.75,
  delta1o = delta1 * gammaI / (delta1 * gammaI + (1 - delta1) * gammaD),
  delta2o = delta1 * gammaIH / (delta1 * gammaIH + (1 - delta1) * gammaDH),
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0.75,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0.197,
  h2 = 0,
  h1o = h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) / (h1 * (gammaI * (1 - delta1o) + gammaD * delta1o) + (1 - h1) * gammaH),
  gamma1 = gammaH * h1o / h1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = gammaDH * delta2o / delta1,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 1 / 4.5,
  gammaFh = 1 / 4.5,
  N = 10000
)
###############################################################################################################################################################
################################################################################# SEIHR MODEL, WITHOUT F  ##############################################################
###############################################################################################################################################################


########################################################### parameters in Chowell 2015
################################# Chowell15.westafrica
#### Pratmeters
pars.cho15 <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 0,
  els = 1,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 1 / 4,
  he1 = 0.5,
  gammaE1 = 1 / 3 + 1 / 3,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 0,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.3056,
  betaIc2 = 0,
  betaIh1 = 0.6 * 0.3056,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,
  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 6,
  gamma2 = 0,
  gammaH1 = 1 / 7,
  gammaH2 = 0,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Fasina 2014 (SEIHR, no F), deduced from R=betaI/gammaI, BetaIc = betaIh,  () see Weitz2015 method.
#### Pratmeters
pars.fas14 <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 0,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 1 / 9,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 1.75 * (1 / 6),
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 1.75 * (1 / 6),
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,
  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,
  gammaH = 1 / 5,
  gammaD = 0,
  gammaI = 1 / 6,
  gammaDH = 0,
  gammaIH = 1 / (6 - 5),
  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = (1 / 5) / (1 / 5 + 1 / 6),
  h2 = 0,
  gamma1 = 1 / 6 + 1 / 5,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 1 / 1,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,

  N = 10000
)
############################################################## parameters in Khan 2015 ; this model do not account feuneral transmission explicily.

pars.kha14_lib <- c(
  psi = 1.6, gammaHw1 = 0.16 + 0.1 + 0.1,
  hw1 = 0.16 / (0.16 + 0.1 + 0.1),
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,
  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,
  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,
  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,
  alpha = 0.1,
  alphaHw = 0.1,
  alphaRh = 0,
  alphaV = 0,
  betaIc1 = 0.3906,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0.3906 * 0.7,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 1.6 * 0.3906 * 0.7,
  betaIhw1 = 1.6 * 0.3906 * 0.7,
  betaIhw2 = 1.6 * 0.3906 * 0.7,
  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,
  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  h1 = 0.16 / (0.16 + 0.1 + 0.1),
  h2 = 0,
  gamma1 = 0.16 + 0.1 + 0.1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0.2 + 0.5,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
## Kha sierra leone
pars.kha14_sie <- c(
  psi = 1.6, gammaHw1 = 0.16 + 0.1 + 0.1,
  hw1 = 0.16 / (0.16 + 0.1 + 0.1),
  deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 0.1,
  alphaHw = 0.1,
  alphaRh = 0,
  alphaV = 0,


  betaIc1 = 0.344,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0.7 * 0.344,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 1.6 * 0.344 * 0.7,
  betaIhw1 = 1.6 * 0.344 * 0.7,
  betaIhw2 = 1.6 * 0.344 * 0.7,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0.16 / (0.16 + 0.1 + 0.1),
  h2 = 0,
  gamma1 = 0.16 + 0.1 + 0.1,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0.2 + 0.5,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Kucharski  2015
############### Sierra Leone 2014 December;
#### Pratmeters
pars.kuc15 <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,

  agu = 0,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 1,
  pnccc = 1,
  pnh = 1,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0.09,
  eppsilonNh = 0,
  ccc1 = (1 / 3) / (1 / 3.0 + 1 / 4.6 + 1 / 9.5),
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alphaCho = 0,
  alpha = 1 / 9.4,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,


  gamma1 = 1 / 3.0 + 1 / 4.6 + 1 / 9.5,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 1 / (9.5 - 4.6),

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 1 / (9.5 - 3),
  gammaNccc = 1 / 7,
  gammaNh = 1 / 7,

  gammaRh = 1,
  gammaRccc = 1,
  gammaFc = 0,
  gammaFh = 0,

  h1 = (1 / 4.6) / (1 / 3 + 1 / 4.6 + 1 / 9.5),
  h2 = 0,

  betaIc1 = 1.94 * (1 / 9.5 + 1 / 4.6),
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = (1 - 0.63) * 1.94 * (1 / 9.5 + 1 / 4.6),
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  N = 10000
)
###############################################################################################################################################################################################
################################################################################### SEIHR MODEL WITHOUT H   ##################################################################################
############################################################################################################################################################################################


########################################### Elisenberg part model (SEIIFR), not hospitalization , PART MODEL ALL COUNTRY

#### parameters
pars.eli14_part_all <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 0,
  els = 1,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  he1 = 0,
  gammaE1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 10,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.11,
  betaIc2 = 3.38 * 0.11,
  betaIh = 0,
  betaIh1 = 0,
  betaIh2 = 0,

  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 3.44 * 0.11,
  betaFh = 3.44 * 0.11,
  betaFhw = 0,

  delta1 = 0.56 / 0.99,
  delta2 = 0.99,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0.36,
  gammaFh = 0.36,

  gamma1 = 0.19,
  gamma2 = 0.53,

  h1 = 0,
  h2 = 0,

  gammaH1 = 0,
  gammaH2 = 0,
  N = 10000
)
################## Elisenberg part model (SEIIFR), Guinea

#### parameters
pars.eli14_part_gui <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 0,
  els = 1,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  he1 = 0,
  gammaE1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 0.11,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.16,
  betaIc2 = 2.04 * 0.16,
  betaIh = 0,
  betaIh1 = 0.16 * 0,
  betaIh2 = 0.16 * 2.04 * 0,

  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 2.63 * 0.16,
  betaFh = 2.63 * 0.16,
  betaFhw = 0,

  delta1 = 0.65 / 0.99,
  delta2 = 0.99,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0.4,
  gammaFh = 0.4,

  gamma1 = 0.19,
  gamma2 = 0.92,

  h1 = 0,
  h2 = 0,

  gammaH1 = 0,
  gammaH2 = 0,
  N = 10000
)
################## Elisenberg part model (SEIIFR), liberia
#### parameters
pars.eli14_part_lib <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 0,
  els = 1,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  he1 = 0,
  gammaE1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 0.11,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.12,
  betaIc2 = 4.9 * 0.12,
  betaIh = 0,
  betaIh1 = 0.12 * 0,
  betaIh2 = 0.12 * 4.9 * 0,

  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 3.84 * 0.12,
  betaFh = 3.84 * 0.12,
  betaFhw = 0,

  delta1 = 0.63 / 0.93,
  delta2 = 0.93,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0.35,
  gammaFh = 0.35,

  gamma1 = 0.2,
  gamma2 = 0.96,

  h1 = 0,
  h2 = 0,

  gammaH1 = 0,
  gammaH2 = 0,

  N = 10000
)
################## Elisenberg part model (SEIIFR), sierra

#### parameters
pars.eli14_part_sie <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 0,
  els = 1,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  he1 = 0,
  gammaE1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 0.1,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.12,
  betaIc2 = 1.73 * 0.12,
  betaIh = 0,
  betaIh1 = 0.12 * 0,
  betaIh2 = 0.12 * 1.73 * 0,

  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 2.97 * 0.12,
  betaFh = 2.97 * 0.12,
  betaFhw = 0,

  delta1 = 0.38 / 0.96,
  delta2 = 0.96,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,
  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,
  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0.38,
  gammaFh = 0.38,

  gamma1 = 0.15,
  gamma2 = 0.73,

  h1 = 0,
  h2 = 0,

  gammaH1 = 0,
  gammaH2 = 0,

  N = 10000
)
################################# Weitz 2015 (SEIFR, no   H)###### the transmission rate needs recalculation!
###### Guinea
#### Pratmeters
pars.wei15_gui <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 11,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = exp(8.8 * 0.011) / 9.08,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 2.2 * exp(8.8 * 0.011) / 9.08,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0.7,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 6,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 1 / 2,
  gammaFh = 0,

  N = 10000
)
###### Liberia
#### Pratmeters
pars.wei15_lib <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 11,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = exp(8.8 * 0.048) / 9.08,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 2.2 * exp(8.8 * 0.048) / 9.08,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0.7,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 6,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 1 / 2,
  gammaFh = 0,

  N = 10000
)
###### Sierra Leone
#### Pratmeters
pars.wei15_sie <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 11,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = exp(8.8 * 0.032) / 9.08,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 2.2 * exp(8.8 * 0.032) / 9.08,
  betaFh = 0,
  betaFhw = 0,


  delta1 = 0.7,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 6,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 1 / 2,
  gammaFh = 0,

  N = 10000
)
############################################################################################# F#########################################################################################################
######################################################################## SEIR MODEL, WITHOUT H, F #########################################################################################################
######################################################################################################################################################################################################



################################# Althaus 2014 ( R0 and gamma beta were provided in the table in the paper
################################# Guinea

#### Pratmeters
pars.alt14_gui <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 5.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.27,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 5.61,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Librea



#### Pratmeters
pars.alt14_lib <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 5.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.28,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 5.61,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Sierra leone



#### Pratmeters
pars.alt14_sie <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 5.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.45,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 5.61,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Bashar 2014 (beta and gamma was not provided, but R0 is provided, so beta0 (before intervention could be calculated)
################################# Guinea
#### Pratmeters
pars.bas14_gui <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 6.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 1.145 * (1 / 11),
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 11,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Sierra Leone  ########delet this one, in the original paper, this one is not explored either, as the R0 do not simulated the real one.

#### Pratmeters
pars.bas14_sie <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 6.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 1.24 * (1 / 11),
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 11,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
########################################################### parameters in Chowell 2004
################################# Chowell04.congo95

#### Pratmeters
pars.cho04_con95 <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 5.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.33,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 5.61,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Chowell04.uganda00


#### Pratmeters
pars.cho04_uga00 <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 3.35,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.38,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,


  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 3.5,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 1,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,


  N = 10000
)
############################################################## parameters in Gomes 2014,SEIR model
pars.gom14_SEIR <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 1,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 7,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 2.1 * 0.1,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,

  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,

  gamma1 = 1 / 10,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,

  N = 10000
)
################################# Lekone 2006 cong95
################################# Vague prior

#### Pratmeters
pars.lek06_vag <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 9.431,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.243,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 5.712,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# informative prior
#### Pratmeters
pars.lek06_inf <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 10.11,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.209,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 6.523,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
#################################################### Meltzer 2014

########################### Libria
#### Pratmeters
pars.mel14_lib <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 6.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.1 * 0.02 + 0.075 * 0.03 + 0.825 * 0.3,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 6,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
#################################################### Meltzer 2014

########################### sieria
#### Pratmeters
pars.mel14_sie <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 6.3,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.15 * 0.02 + 0.25 * 0.03 + 0.6 * 0.3,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 6,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Ferrary 2004, SIR model, only R0 and gamma were provided, beta was extimated from R0*gamma.E is missing. I assume there is E, but E to I is very fast rate.
################################# cong95

#### Pratmeters
pars.fer04_con95 <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 3.65 * (1 / 14),
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 14,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# uganda 2000

#### Pratmeters
pars.fer04_uga00 <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 0,
  leg = 0,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 1,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 1.79 * (1 / 14),
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 1,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 14,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Shamen 2014 (SEIFR, no H, but F is not transmisinable ) ## time variant, based on the reported number on 27/sep 2014
################################# Guinea

#### Pratmeters
pars.sha14_gui <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,

  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 10.75,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.17,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 7.08,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0,

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,

  N = 10000
)
################################# Liberia
#### Pratmeters
pars.sha14_lib <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0,
  agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,


  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 5.87,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.17,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,


  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 9.99,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0, #

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################# Siera Leone

#### Pratmeters
pars.sha14_sie <- c(
  psi = 0, gammaHw1 = 0, hw1 = 0, deltaHw1 = 0, agu = 1,
  leg = 1,
  els = 0,
  mel = 0,
  fhmax1 = 0,
  fhmax2 = 0,
  fcccmax1 = 0,
  fcccmax2 = 0,
  pnccc = 0,
  pnh = 0,

  alphaCho = 0,
  gammaE1 = 0,
  he1 = 0,

  eppsilonNccc = 0,
  eppsilonNh = 0,
  ccc1 = 0,
  ccc2 = 0,
  phiH1 = 0,
  phiH2 = 0,
  phiCcc1 = 0,
  phiCcc2 = 0,

  rhoRv = 0,
  rhoV = 0,
  rhoRh = 0,
  rhoH = 0,
  omegaRh = 0,

  alpha = 1 / 12.61,
  alphaHw = 0,
  alphaRh = 0,
  alphaV = 0,

  betaIc1 = 0.34,
  betaIc2 = 0,
  betaIh1 = 0,
  betaIh2 = 0,
  betaIccc1 = 0,
  betaIccc2 = 0,
  betaIrh1 = 0,
  betaIv1 = 0,
  betaIhwo1 = 0,
  betaIhw1 = 0,
  betaIhw2 = 0,

  betaFc = 0,
  betaFh = 0,
  betaFhw = 0,

  delta1 = 0,
  delta2 = 0,
  deltaH1 = 0,
  deltaH2 = 0,
  deltaB2 = 0,
  deltaCcc1 = 0,
  deltaCcc2 = 0,

  h1 = 0,
  h2 = 0,
  gamma1 = 1 / 2.91,
  gamma2 = 0,
  gammaH1 = 0,
  gammaH2 = 0, #

  gammaCcc0 = 0,
  gammaH0 = 0,
  gammaCcc1 = 0,
  gammaCcc2 = 0,
  gammaNccc = 0,
  gammaNh = 0,

  gammaRh = 0,
  gammaRccc = 0,
  gammaFc = 0,
  gammaFh = 0,
  N = 10000
)
################################################################################################################################################################################################################
################################################################################ part 3 run models ######################################################################################################
################################################################################################################################################################################################################
#### the common part for all managemnts

nrun <- 1000
times <- seq(0, 2050, by = 0.25)
### the intervention efficacy with the cost function of "Cheap and effective" under low, intermediate and high budget levels is c(0.50,0.99,1.00);
### the intervention efficacy with the cost function of "Expensive and effective" under low, intermediate and high budget levels is c(0.08,0.50,0.99);
### the intervention efficacy with the cost function of "Cheap and partly effective" under low, intermediate and high budget levels is c(0.38,0.74,0.74);

intensity <- c(0.5, 0.99, 1, 0.08, 0.5, 0.99, 0.38, 0.74, 0.75)
betaIc.m <- 0.24340284
# community transmission based on the mean model
betaIh.m <- 0.103692478
# hospital transmission based on the mean model
betaIf.m <- 0.53268485
# funeral transmission based on the mean model
betaIchf.m <- 0.293571905
betaIch.m <- 0.193721753951408
betaIcf.m <- 0.357152112361342
deltaC.m <- 0.667967113
deltaH.m <- 0.631825725
theta.m <- 0.551832569
dC <- (1 - theta.m) * deltaC.m / ((1 - theta.m) * deltaC.m + theta.m * deltaH.m)
dH <- theta.m * deltaH.m / ((1 - theta.m) * deltaC.m + theta.m * deltaH.m)
wC <- 1 / (1 + theta.m + (1 - theta.m) * deltaC.m + theta.m * deltaH.m)
wH <- theta.m / (1 + theta.m + (1 - theta.m) * deltaC.m + theta.m * deltaH.m)
wF <- ((1 - theta.m) * deltaC.m + theta.m * deltaH.m) / (1 + theta.m + (1 - theta.m) * deltaC.m + theta.m * deltaH.m)
wCH_C <- 1 / (1 + theta.m)
wCH_H <- theta.m / (1 + theta.m)
wCF_C <- 1 / (1 + (1 - theta.m) * deltaC.m + theta.m * deltaH.m)
wCF_F <- ((1 - theta.m) * deltaC.m + theta.m * deltaH.m) / (1 + (1 - theta.m) * deltaC.m + theta.m * deltaH.m)

pC <- wC * betaIc.m / betaIchf.m
pH <- wH * betaIh.m / betaIchf.m
pF <- wF * betaIf.m / betaIchf.m
pCH_C <- wCH_C * betaIc.m / betaIch.m
pCH_H <- wCH_H * betaIh.m / betaIch.m
pCF_C <- wCF_C * betaIc.m / betaIcf.m
pCF_F <- wCF_F * betaIf.m / betaIcf.m


##### this part is intervention and model type specific######################################################################################
##### Simulation under the intervenion of reducing community transmission ###################################################################
# full model
runModel.full <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaIc1"] <- pars.t["betaIc1"] * (1 - intensity [i])
    pars.t["betaIc2"] <- pars.t["betaIc2"] * (1 - intensity [i])
    pars.t["betaIrh1"] <- pars.t["betaIrh1"] * (1 - intensity [i])

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "comtra.rda", sep = ""))
}


# H model
runModel.h <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaIc1"] <- (1 - pCF_C) * pars.t["betaIc1"] + (1 - intensity [i]) * pCF_C * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (1 - pCF_C) * pars.t["betaIc2"] + (1 - intensity [i]) * pCF_C * pars.t["betaIc2"]

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "comtra.rda", sep = ""))
}


# F model
runModel.f <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaIc1"] <- (1 - pCH_C) * pars.t["betaIc1"] + (1 - intensity [i]) * pCH_C * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (1 - pCH_C) * pars.t["betaIc2"] + (1 - intensity [i]) * pCH_C * pars.t["betaIc2"]

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)


      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")


      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "comtra.rda", sep = ""))
}


# SEIR model
runModel <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaIc1"] <- (1 - pC) * pars.t["betaIc1"] + (1 - intensity [i]) * pC * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (1 - pC) * pars.t["betaIc2"] + (1 - intensity [i]) * pC * pars.t["betaIc2"]

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "comtra.rda", sep = ""))
}


##### Simulation under the intervenion of improving hospitalization ###################################################################
# full model
runModel.full <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model

    pars.t["he1"] <- ifelse(pars.t["he1"] * (1 + intensity[i]) < 1, pars.t["he1"] * (1 + intensity[i]), 1)
    pars.t["h1"] <- ifelse(pars.t["h1"] * (1 + intensity[i]) < 1, pars.t["h1"] * (1 + intensity[i]), 1)
    pars.t["h2"] <- ifelse(pars.t["h2"] * (1 + intensity[i]) < 1, pars.t["h2"] * (1 + intensity[i]), 1)
    pars.t["hw1"] <- ifelse(pars.t["hw1"] * (1 + intensity[i]) < 1, pars.t["hw1"] * (1 + intensity[i]), 1)

    pars.t["betaIh1"] <- pars.t["betaIh1"] * (1 - intensity[i])
    pars.t["betaIh2"] <- pars.t["betaIh2"] * (1 - intensity[i])
    pars.t["betaIhwo1"] <- pars.t["betaIhwo1"] * (1 - intensity[i])
    pars.t["betaIhw1"] <- pars.t["betaIhw1"] * (1 - intensity[i])
    pars.t["betaIhw2"] <- pars.t["betaIhw2"] * (1 - intensity[i])
    pars.t["betaIv1"] <- pars.t["betaIv1"] * (1 - intensity[i])

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "hospital.rda", sep = ""))
}


# H model
runModel.h <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model ## for H model , kha15, it also has betaIhawo1

    pars.t["he1"] <- ifelse(pars.t["he1"] * (1 + intensity[i]) < 1, pars.t["he1"] * (1 + intensity[i]), 1)
    pars.t["h1"] <- ifelse(pars.t["h1"] * (1 + intensity[i]) < 1, pars.t["h1"] * (1 + intensity[i]), 1)
    pars.t["h2"] <- ifelse(pars.t["h2"] * (1 + intensity[i]) < 1, pars.t["h2"] * (1 + intensity[i]), 1)
    pars.t["hw1"] <- ifelse(pars.t["hw1"] * (1 + intensity[i]) < 1, pars.t["hw1"] * (1 + intensity[i]), 1)

    pars.t["betaIh1"] <- pars.t["betaIh1"] * (1 - intensity[i])
    pars.t["betaIh2"] <- pars.t["betaIh2"] * (1 - intensity[i])
    pars.t["betaIhwo1"] <- pars.t["betaIhwo1"] * (1 - intensity[i])
    pars.t["betaIhw1"] <- pars.t["betaIhw1"] * (1 - intensity[i])
    pars.t["betaIhw2"] <- pars.t["betaIhw2"] * (1 - intensity[i])
    pars.t["betaIv1"] <- pars.t["betaIv1"] * (1 - intensity[i])


    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "hospital.rda", sep = ""))
}


# F model
runModel.f <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model

    theta.new <- ifelse(theta.m * (1 + intensity[i]) < 1, theta.m * (1 + intensity[i]), 1)
    wCH_C.new <- 1 / (1 + theta.new)
    wCH_H.new <- theta.new / (1 + theta.new)
    pars.t["betaIc1"] <- (wCH_C.new / wCH_C) * pCH_C * pars.t["betaIc1"] + (wCH_H.new / wCH_H) * pCH_H * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (wCH_C.new / wCH_C) * pCH_C * pars.t["betaIc2"] + (wCH_H.new / wCH_H) * pCH_H * pars.t["betaIc2"]
    pars.t["delta1"] <- (1 - theta.new) / (1 - theta.m) * dC * pars.t["delta1"] + theta.new / theta.m * dH * pars.t["delta1"]
    pars.t["delta2"] <- (1 - theta.new) / (1 - theta.m) * dC * pars.t["delta2"] + theta.new / theta.m * dH * pars.t["delta2"]

    pars.t["betaIc1"] <- (1 - pCH_H) * pars.t["betaIc1"] + (1 - intensity[i]) * pCH_H * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (1 - pCH_H) * pars.t["betaIc2"] + (1 - intensity[i]) * pCH_H * pars.t["betaIc2"]


    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "hospital.rda", sep = ""))
}


# SEIR model
runModel <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model

    theta.new <- ifelse(theta.m * (1 + intensity[i]) < 1, theta.m * (1 + intensity[i]), 1)
    wC.new <- 1 / (1 + theta.new + (1 - theta.new) * deltaC.m + theta.new * deltaH.m)
    wH.new <- theta.new / (1 + theta.new + (1 - theta.new) * deltaC.m + theta.new * deltaH.m)
    wF.new <- ((1 - theta.new) * deltaC.m + theta.new * deltaH.m) / (1 + theta.new + (1 - theta.new) * deltaC.m + theta.new * deltaH.m)
    pars.t["betaIc1"] <- (wC.new / wC) * pC * pars.t["betaIc1"] + (wH.new / wH) * pH * pars.t["betaIc1"] + (wF.new / wF) * pF * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (wC.new / wC) * pC * pars.t["betaIc2"] + (wH.new / wH) * pH * pars.t["betaIc2"] + (wF.new / wF) * pF * pars.t["betaIc2"]

    pars.t["betaIc1"] <- (1 - pH) * pars.t["betaIc1"] + (1 - intensity[i]) * pH * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (1 - pH) * pars.t["betaIc2"] + (1 - intensity[i]) * pH * pars.t["betaIc2"]


    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)


      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "hospital.rda", sep = ""))
}

##### Simulation under the intervenion of reducing funeral transmission ###################################################################

# full model
runModel.full <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaFc"] <- pars.t["betaFc"] * (1 - intensity [i])
    pars.t["betaFh"] <- pars.t["betaFh"] * (1 - intensity [i])
    pars.t["betaFhw"] <- pars.t["betaFhw"] * (1 - intensity [i])

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "funtra.rda", sep = ""))
}


# H model
runModel.h <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))


  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaIc1"] <- (1 - pCF_F) * pars.t["betaIc1"] + (1 - intensity[i]) * pCF_F * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (1 - pCF_F) * pars.t["betaIc2"] + (1 - intensity[i]) * pCF_F * pars.t["betaIc2"]

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "funtra.rda", sep = ""))
}


# F model
runModel.f <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))


  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaFc"] <- pars.t["betaFc"] * (1 - intensity[i])
    pars.t["betaFh"] <- pars.t["betaFh"] * (1 - intensity[i])
    pars.t["betaFhw"] <- pars.t["betaFhw"] * (1 - intensity[i])

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")

      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "funtra.rda", sep = ""))
}


# SEIR model
runModel <- function(pars.model, modelName, nrun, intensity) {
  output <- array(NA, dim = c(length(seq(1, 2002 * 4 + 1, by = 28)), length(statenames) + 1, nrun, length(intensity)))
  dur <- matrix(NA, nrow = nrun, ncol = length(intensity))

  for (i in 1:length(intensity)) {
    pars.t <- pars.model
    pars.t["betaIc1"] <- (1 - pF) * pars.t["betaIc1"] + (1 - intensity[i]) * pF * pars.t["betaIc1"]
    pars.t["betaIc2"] <- (1 - pF) * pars.t["betaIc2"] + (1 - intensity[i]) * pF * pars.t["betaIc2"]

    for (r in 1:nrun) {
      pars.sto <- pars.t
      pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")] <-
        rlnorm(14, log(pars.sto[c("betaIc1", "betaIc2", "betaIh1", "betaIh2", "betaIccc1", "betaIccc2", "betaIrh1", "betaIv1", "betaIhwo1", "betaIhw1", "betaIhw2", "betaFc", "betaFh", "betaFhw")]), 0.1)

      Gmodel <- gillesp(
        start = c(
          Sc = 9989, Shw = 10, Srh = 0, Sv = 0, Ec = 0, Ehw = 0, Erh = 0, Ev = 0, Ecd = 0, Ic1 = 1, Ih1 = 0, Ihw1 = 0, Irh1 = 0, Iv1 = 0, Iccc1 = 0, Inccc2 = 0, Inh2 = 0, Ic2 = 0, Ih2 = 0, Iccc2 = 0, Fc = 0, Fh = 0, Rh = 0, Rccc = 0, R = 0, CC = 1, CD = 0, CP = 1
        ),
        times = times, ratefun = ratefun, trans = trans, pars = pars.sto, deltaT = 0.1
      )
      cat(date(), " - ", i, "\n")
      dur[r, i] <- sort(Gmodel[Gmodel[, length(statenames) + 1] == 0, 1])[1]
      output[, , r, i] <- Gmodel[seq(1, 2002 * 4 + 1, by = 28), ]
    }
  }
  save(output, dur, file = paste(modelName, "funtra.rda", sep = ""))
}

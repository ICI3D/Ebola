##################################################################################
## Ebola vaccine trial equipoise model
##
## Accompanying:
##
## Bellan SE, JRC Pulliam, J Dushoff, and LA Meyers. (2014) Ebola Vaccine
## Trials: The Ethical Mandate for a Therapeutic Safety Net. _The British
## Medical Journal_ 349:g7518.
##
## Steve Bellan, October 2014

###############
### LICENSE
###
### This code is made available under a Creative Commons Attribution 4.0
### International License. You are free to reuse this code provided that you
### give appropriate credit, provide a link to the license, and indicate if
### changes were made.
### You may do so in any reasonable manner, but not in any way that suggests
### the licensor endorses you or your use. Giving appropriate credit includes
### citation of the above publication *and* providing a link to this repository:
###
### https://github.com/ICI3D/Ebola
###############


## CONTROL PARAMETERS
SAVEPLOTS <- FALSE # Change this to TRUE to save figures 
if(SAVEPLOTS){
	try(system('mkdir Figures'))
	if(file.exists("Figures")){
		plotPath <- 'Figures'
	}else{
		warning("Figures directory does not exist; plotPath will not be defined, and figures will be plotted but not saved.")
		SAVEPLOTS <- FALSE
	}
}

## THERAPEUTIC EFFICACY

## Is there observational evidence that experimental therapies are effective?
## Ongoing outbreak, 3/14 CFR amongst cases known to have been given experimental therapies and with definitive outcomes (Table S1).
## Ongoing outbreak, .642 CFR amongst 1303 hopsitalized cases with definitive outcomes.

## Restrict analysis to 20-40ish age groups to control for age effect
expther <- c(8,0) ## 0/8 died w/ experimental therapy (19-40)
schif <- c(8,27) ## 27/35 died in SL intensive care (Schieffelin et al.) 21-40 yrs
bah <- c(15,4) ## 4/19 died in Conakry intensive care (Schieffelin et al.) 19-40 yrs
schifbah <- c(23,31) ## 31/54 died in SL & Conakry intensive care (pooled) 19-40 yrs
wa <- c(261,577)## 577/838 in all West Africa (WHO NEJM) 15-44 yrs
fisher.test(cbind(expther,schif)) ## sig
fisher.test(cbind(expther,bah)) ## nsig
fisher.test(cbind(expther,schifbah)) ## sig
fisher.test(cbind(expther,wa)) ## sig
oratio <- function(x) x[1,1]*x[2,2] / (x[1,2]*x[2,1])
lapply(list(schif,bah,schifbah,wa), function(x) oratio(cbind(expther,x))) ## all Inf because 8/8 survival

pefxn <- function(x) 1-  x[2,1]/sum(x[,1]) / (x[2,2]/sum(x[,2]))
lapply(list(schif,bah,schifbah,wa), function(x) pefxn(cbind(expther,x))) ## all 1 because 8/8 survival

protEffic <- 1-(3/16)/.57 ## 67% between all experimental therapy (not just low age groups, and so conservative) and pooled Bah and Schieffelin results
protEffic <- .70 ## assume 70% protective efficacy for figure. This is conservative based on

## available data for above reason. Yet, any protective efficacy above 0 would
## reduce risk diferrences between arms. As argued in the main text, the fact that
## these therapies are consisently used in the developed world suggests that the
## anticipated protective efficacy (including risk of adverse effects) is expected
## to be >0.

## Use binom.test to get confidence intervals
## NOTE: This ignores variance around .57, which we feel is justifibale since
## the sample size is so big)

protEffic95CI <- 1 - binom.test(3,16)$conf[2:1]/.57
protEffic95CI ## 95% CI: 21-93%
## Note: for small sample sizes, binom.test() gives more accurate confidence
## intervals than approximate methods that acounts for the sample size in
## comparison group (1303):
library(epicalc)
csi(3,13,837,466)
## 95% CI here would be 33-89%; we use the more conservative, wider values
## above in the manuscript text


## EXPLORING EQUIPOISE

## Expected number of deaths among participants by trial arm before stopping
## criteria is met in a hypothetical vaccine trial, with and without a therapy
## that reduces the case fatality rate by 70%. Clinical equipoise is fully
## achieved when the expected number of deaths is identical for the two arms
## (or close enough to be counterbalanced by benefits from an active control
## or adverse vaccine side effects). The therapeutic safety net reduces the
## overall mortality and substantially bridges the gap between the two arms.
## However, it has no impact on the number of infections, and thus does not
## reduce power to detect vaccine efficacy. This scenario assumes that the
## experimental vaccine has a 50% chance of having 80% efficacy and a 50% 
## chance of being completely ineffective. Following vaccine trial proposals,
## GSK we further assume that study participants have a 10% annual risk of
## Ebola infection, and that 30 participants will become infected before the
## trial reaches its stopping criteria. In the absence of a therapeutic safety
## net, 57% of these 30 cases are expected to be fatal (Table 1).

## Calculate mortality risk in both arms of a vaccine trial
mortRisk <- function(exposureHazard = .1, cfr = .57, probVaccWorks = .5, VaccEfficacy = .8, therapEfficacy = .75) {
	riskVacc <- (probVaccWorks * (1-VaccEfficacy) + (1-probVaccWorks)) * exposureHazard * cfr * (1-therapEfficacy)
	riskCont <- exposureHazard * cfr * (1-therapEfficacy)
	return(data.frame(riskVacc=riskVacc, riskCont=riskCont))
}

cfr <- .57 ## amongst definitive hospitalized cases in west Africa
numInfectAtStop <- cfr*30 ## GSK estimated # infections at stopping criteria * cfr = stopping # of deaths
therapEfficacy <- c(0,round(protEffic,2)) ## Example therapeutic efficacy
nArms <- 1000 ## number patients in each arm

morts <- mortRisk(therapEfficacy = therapEfficacy)
if(SAVEPLOTS) png(file.path(plotPath, 'Mortality Balance.png'), w = 4, h= 3.5, units='in', res = 300)
par(lwd=2, 'ps'=12, mar = c(4,6,4,.5))
plot(0, 0, xlab = '', ylab = 'Expected # Ebola Deaths Before \nStopping Criteria Met',
		 bty = 'n', xlim = c(.7,2.3), ylim = c(0,.07), lty = 1, las = 1, type = 'l', axes=F)
points(1:2, morts$riskVacc, pch = 19, cex = 1.5)
points(1:2, morts$riskCont, pch = 21, cex = 1.5)
axis(1, at = 1:2, labels = c('no therapy',paste0(therapEfficacy[2]*100,'% effective\n therapy')), padj = 1)
yticksRisk <- seq(0,.07, by = .01)
## Stopping criteria assumed to be met at 30 infections by GSK.
## 30 inf * 64% cfr = 19.26 deaths. Let's figure out how long the trial has gone until this happens.
## 19.4 deaths = AnnualMortalityRisk * people * time period in yrs
timeToStop <- numInfectAtStop/sum(morts[1,])/(nArms*2) ## in years
yticksDeaths <- pretty(yticksRisk*numInfectAtStop/sum(morts[1,]), 5) ## Convert annual risk to # of deaths by the end of trial
par(xpd=T)
axis(2, at = yticksDeaths*sum(morts[1,])/numInfectAtStop, labels = yticksDeaths, las = 1)
par(xpd=T, mar = rep(0,4))
legend(.5, max(yticksRisk)*1.4, leg = c('control arm','vaccine arm'), pch = c(21,19), ncol = 2)
if(SAVEPLOTS) graphics.off()

## Look at it over full range of therapeutic efficacies
therapEfficacy <- seq(0, 1, by = .01)
morts <- mortRisk(therapEfficacy = therapEfficacy)
if(SAVEPLOTS) png(file.path(plotPath, 'Mortality Balance vs therapeutic efficacy.png'), w = 3.5, 4, units='in', res = 300)
par(lwd=2, 'ps'=12, mar = c(5,4,1,.5))
plot(therapEfficacy, morts$riskCont, xlab = 'Therapeutic Efficacy', ylab = 'Annual Ebola Mortality Risk',
		 bty = 'n', xlim = c(0,1), ylim = c(0,max(morts)), lty = 1, las = 1, type = 'l')
lines(therapEfficacy, morts$riskVacc, lty = 2)
legend('topright', leg = c('control arm','vaccine arm'), lty = 1:2, bty = 'n')
if(SAVEPLOTS) graphics.off()

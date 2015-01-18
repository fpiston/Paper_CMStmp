# Read the data
cmstmp  <- read.csv('./fertility_environment.csv')

# Data transformation to days from 1970 
cmstmp$Heading.date <- as.numeric(as.Date(cmstmp$Heading.date, format="%d/%m/%y"))

## Codification of new variables
cmstmp$No.cuaj <- 20 - cmstmp$N.gr20 # number of sterile flowers
cmstmp$Trial20 <- 20 # number of observed flowers
cmstmp$Fertility1 <- cmstmp$Fertility/100 # transformation to fertility proportion
cmstmp$planta  <- as.factor(paste(cmstmp$genotipo, cmstmp$planta, sep="-"))
cmstmp$Clone <- as.factor(cmstmp$Clone)
cmstmp$Tallo.emer <- as.factor(cmstmp$Tallo.emer)
cmstmp$gr.spik <- (cmstmp$N.gr.total/cmstmp$N.espigui.espi)
cmstmp$flores.total <- (cmstmp$N.espigui.espi * 5) # We assume 5 viable flowers per spikelet
cmstmp$gr.flores <- (cmstmp$N.gr.total/cmstmp$flores.total)
cmstmp$flores.fallos <- (cmstmp$flores.total-cmstmp$N.gr.total)

################################################################################
## Data explanation
################################################################################

## genotipo      : Factor; plant genotype
## planta        : Factor; individual plants of each genotype
## Clone         : Factor, clone id. individual plants were cloned to have
##### replications in both locations
## Ubicacion     : Factor; locations where the plants were grown
## Heading.date  : num; date at heading stage
## Lenght.day    : num; day lenght at anthesis
## Tmean         : num; average day temperature at anthesis date
## Tmin          : num; minimum day temperature at anthesis date
## Tmax          : num; maximum day temperature at anthesis date
## N.espigui.espi: int; number of spikelets per spike
## N.gr20        : int; number of grains of the 20 central flowers
## N.gr.total    : int; number of total grains of the spike
## Fertility     : int; fertility measures expressed as the percentage of grains of the 20 central
#### flowers
## Tallo.emer    : Factor; order of emergence
## No.cuaj       : num; number of grains
## Trial20       : num; number of total flowers
## Fertility1    : num; fertility estimator as proportion of number of grains of the 20
#### central flowers
## gr.spik       : num; total grains per spike
## flores.total  : num; total number of flowers per spike
## gr.flores     : num; fertility estimator expressed as total grains per total number of flowers
## flores.fallos : num; number of sterile flowers
## Te1H          : Factor; plant with addition of 1HchS telo
## Te6H          : Factor; plant with addition of 6HchS telo
## Tr1HB         : Factor; plant with translocation of 1HchS·1B
## Tr6HD         : Factor; plant with translocation of 1HchS·1B

################################################################################
## Correlations between the two fertility estimators                           #
################################################################################
summary(fercoT <- lm(cmstmp$Fertility1 ~ cmstmp$gr.flores))
summary(fercoTj <- lm(Fertility1 ~ gr.flores, data = cmstmp[cmstmp$Ubicacion=='jaulon',]))
summary(fercoTi <- lm(Fertility1 ~ gr.flores, data = cmstmp[cmstmp$Ubicacion=='invernadero',]))

## Corresponds to paper Figure 1
pdf('../paper/img/correlations-fertility-est.pdf', width=9, height=9)

par(mfrow=c(2,2))
plot(cmstmp$Fertility1 ~ cmstmp$gr.flores, ylab='Relative fertility', xlab='Total fertility', type='n')
abline(fercoT, col='red', lwd=2, lty=2)
points(cmstmp$Fertility1 ~ cmstmp$gr.flores)
text(0.6,0.1,paste('R² =', sprintf('%.2f',summary(fercoT)$r.squared)),cex=2)

plot(Fertility1 ~ gr.flores, data = cmstmp[cmstmp$Ubicacion=='jaulon',], main='Field', ylab='Relative fertility', xlab='Total fertility', type='n')
abline(fercoT, col='red', lwd=2, lty=2)
points(Fertility1 ~ gr.flores, data = cmstmp[cmstmp$Ubicacion=='jaulon',])
text(0.6,0.1,paste('R² =', sprintf('%.2f',summary(fercoTj)$r.squared)),cex=2)

plot(Fertility1 ~ gr.flores, data = cmstmp[cmstmp$Ubicacion=='invernadero',], main='Greenhouse', ylab='Relative fertility', xlab='Total fertility', xlim=c(0,0.9), type='n')
abline(fercoT, col='red', lwd=2, lty=2)
points(Fertility1 ~ gr.flores, data = cmstmp[cmstmp$Ubicacion=='invernadero',])
text(0.6,0.1,paste('R² =', sprintf('%.2f',summary(fercoTi)$r.squared)),cex=2)
dev.off()

############################################################
## Fit linear model                     ####################
############################################################

## As both fertility estimators are highly correlated we used
## the 20 flower fertility proportion values ('Fertility1')
## We decided to use this estimator because it represents well the fertility and it is
## less time consuming for future works

## Saturated model
lm1 <- lm(Fertility1 ~ planta:genotipo + (genotipo + Ubicacion + Heading.date + Lenght.day + Tmean + Tmin + Tmax)^2, data = cmstmp)

## Check linear model assumptions using plots
plot(lm1)

## the heteroscedasticity of residuals has some pattern in the 0 and close values
## Too many 0 is difficult to resolve with data transformation. We tried the log
## and arc sine transformation but the problem was not resolved . Moreover the resulting R² model 
## is better without any data transformation.
## Other posibility was to use the generalized linear model with a binomial or quasibinomial
## residual struture. We also did it, but  model complexity was increased with very similar
## results. So, we decided to use the simplest method, in our case, linear model with
## gaussian residual struture (the normal one). 

##### Model reduction based on significant test
## We removed the non-significant elements (see the anova table of model 'lm1')

lm1r1 <- lm(Fertility1 ~ (genotipo + Ubicacion + Tmean)^2, data = cmstmp)

## In addition:
## The interaction 'planta:genotipo' was removed because we did not have
## enough data to support this additional model element (the 'planta' inside 'genotipo').

## 'Lenght.day' and 'Heading.date' highly correlated with the factor 'Ubicacion'
## so they were removed them from the model. Moreover, they did not explain too much proportion of the
## variance.


##### Linear model assumption check by plot
par(mfrow=c(2,2))
plot(lm1r1, main='lm1r1, Fertilidad20') # It fits relatively well , except for the abundance of zeros
##dev.off()

#####  ANOVA tables and model summary
## Corresponds to paper Table 3
anova(lm1r1)

summary(lm1r1)



################################################################################
## Paper plots                                                                 #
################################################################################

## Corresponds to paper Figure 2
pdf('../paper/img/fertility1-genotype-place.pdf', width=11, height=7)


m <- rbind(c(1,1,1,1,1,2))
print(m)
layout(m)

par(mar = c(0, 0, 0, 0), oma = c(4, 6, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))



boxplot(
  formula = Fertility1 ~ genotipo,
  data    = cmstmp,
  boxwex  = 0.25,
  at      = 1:8 - 0.15,
  subset  = Ubicacion == "invernadero",
  col     = "bisque",
##  main    = "blah",
  xlab    = "",
  xaxt    = 'n',
  yaxt    = 'n',
  ylab    = "",
  ylim    = c( -0.01, 1.01),
  yaxs    = "i"
)
mtext('Fertility', side = 2, line = 3, outer = FALSE, at = NA,
           adj = NA, padj = NA, cex = 1.3, col = NA, font = NA)
axis(2, cex.axis=1.4)

boxplot(
  formula = Fertility1 ~ genotipo,
  data    = cmstmp,
  boxwex  = 0.25,
  at      = 1:8 + 0.15,
  subset  = Ubicacion == "jaulon",
  col     = "gray80",
  add     = TRUE,
  xlab    = "",
  xaxt    = 'n',
  yaxt    = 'n'
)

labels <- levels(cmstmp$genotipo)
axis(1, labels = labels, at = c(1:8), cex.axis=1.2)

legend("topleft", c("Greenhouse", "Field"),fill = c('bisque', 'gray80'),bty = "n", cex=1.5)

boxplot(
  formula = Fertility1 ~ Ubicacion,
  data    = cmstmp,
  yaxt    = 'n',
  xaxt    = 'n',
  col     = c('bisque','gray80'),
  ylim    = c( -0.01, 1.01),
  yaxs    = "i",
  boxwex  = 0.5
)


axis(1, labels = 'Total', at = 1.5, cex.axis=1.5)
dev.off()



################################################################################
## Evaluation of the chromosomal configuration using lsmeans                    #
################################################################################

library(car)
## Genotypes with the telo 1HchS addition
cmstmp$Te1H <- recode(cmstmp$genotipo, "c('T218T21AH1S','T593T21AH1S') = 1; else = 0")
## Genotypes with the telo 6HchS addition
cmstmp$Te6H <- recode(cmstmp$genotipo, "c('T593T21','T593T21AH1S','T593T552') = 1; else = 0")
## Genotypes with 1HchS.1BL translocation
cmstmp$Tr1HB <- recode(cmstmp$genotipo, "c('T218T552','T236T552','T593T552','T650T552') = 1; else = 0")
## Genotypes with 6HchS.6DL translocation
cmstmp$Tr6HD <- recode(cmstmp$genotipo, "c('T218T650','T650T552') = 1; else = 0")

## Genotypes with 1HchS translocation or addition
cmstmp$A1H <- recode(cmstmp$genotipo, "c('T218T21AH1S','T593T21AH1S','T218T552','T236T552','T593T552','T650T552') = 1; else = 0")
## Genotypes with 6HchS translocation or addition
cmstmp$A6H <- recode(cmstmp$genotipo, "c('T593T21','T593T21AH1S','T593T552','T218T650','T650T552') = 1; else = 0")


################################################################################
## Fit model with the different chromosomal configurations like factors        #
################################################################################

library(lsmeans)

lm1r1C <- lm(Fertility1 ~ Ubicacion + Tmean + Te1H + Te6H + Tr1HB + Tr6HD + Te1H:Te6H + Te6H:Tr1HB + Tr1HB:Tr6HD + Ubicacion:Tmean + Ubicacion:Te1H + Ubicacion:Te6H + Ubicacion:Tr1HB + Ubicacion:Tr6HD + Tmean:Te1H + Tmean:Te6H + Tmean:Tr1HB + Tmean:Tr6HD, data = cmstmp)

anova(lm1r1C)

## Tr1HB1:Tr6HD1, we do not have enough data to calcule if the interaction is significant 
## but we can calculate de additive effect of Tr1HB1 and Tr6HD1 together. 

means <- na.omit(as.data.frame(summary(lsmeans(lm1r1C, ~ Te1H + Te6H + Tr1HB + Tr6HD))))
colnames(means) <- c('Te1H','Te6H','Tr1HB','Tr6HD','lsmean','SE','df','lower.CL','upper.CL')

means <- means[-c(1,6,8,10,11,12,14:16),]
means$Chr <- c('Te1H','Te6H','Te1H+\nTe6H','Tr1HB','Te6H+\nTr1HB','Tr6HD','Tr1HB+\nTr6HD')

## Corresponds to paper Figure 3
pdf('../paper/img/lsmean-cromosomes.pdf', width=6, height=7)

x <- barplot(means$lsmean, ylim = c(0,1.01), ylab='Fertility')
arrows(x,means$lower.CL,x,means$upper.CL, length = 0.1, angle = 90, code = 3)
text(x-0.4, par("usr")[3] - 0.04, srt = 0, adj = c(0,0.5), labels = means$Chr, xpd = TRUE, cex = 1)

box()

dev.off()

lm1r1CA <- lm(Fertility1 ~ (Ubicacion + Tmean + A1H + A6H)^2, data = cmstmp)
means <- lsmeans(lm1r1CA, ~ Te1H + Te6H + Tr1HB + Tr6HD)


################################################################################
## Evaluation of the factor Genotype by lsmeans                                #
################################################################################
lm1r1 <- lm(Fertility1 ~ (genotipo + Ubicacion + Tmean)^2, data = cmstmp)
means <- lsmeans(lm1r1, pairwise ~  genotipo)
means <- na.omit(as.data.frame(summary(means[[1]])))

x11()
x <- barplot(means$lsmean, ylim = c(0,1.01), ylab='Fertility')
arrows(x,means$lower.CL,x,means$upper.CL, length = 0.1, angle = 90, code = 3)
text(x, par("usr")[3] - 0.03, srt = 45, adj = c(1,0), labels = means$genotipo, xpd = TRUE, cex = 1)

box()

################################################################################
## Evaluation of the mean temperature and Genotype effect by lsmeans           #
################################################################################
means <- lsmeans(lm1r1, list( ~ genotipo*Tmean*Ubicacion), at=list(Tmean=c(10,12,14,16,18,20,22,24)))

meansT <- as.data.frame(means[[1]])

## Corresponds to paper Figure 4
pdf('../paper/img/fitted-fertility-genotype.pdf', width=12, height=9)
m <- rbind(c(1,2,3,4),c(5,6,7,8))
print(m)
layout(m)

par(mar = c(2, 2, 2, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
par(cex = 1)
geno <- c("T218T21AH1S", "T218T552", "T218T650", "T236T552", "T593T21", "T593T21AH1S", "T593T552", "T650T552")

for (i in geno){


prueba <- meansT[which(meansT$genotipo==i),]
pruebaJ <- meansT[which(meansT$genotipo==i & meansT$Ubicacion=='jaulon'),]
pruebaI <- meansT[which(meansT$genotipo==i & meansT$Ubicacion=='invernadero'),]
plot(lsmean ~ Tmean, data = prueba, col=as.numeric(prueba$Ubicacion), type='n')
  if (i == 'T218T21AH1S') {
    legend("topright", c("Greenhouse", "Field"),fill = c(1, 2),bty = "n", cex=1.2)
  }
  if (i == 'T218T650' || i == 'T593T21AH1S' || i == 'T593T552' || i == 'T650T552') {
    u <- par("usr")
    rect(u[1], u[3], u[2], u[4], col = "grey90")
  }
arrows(pruebaI$Tmean,pruebaI$lower.CL,pruebaI$Tmean,pruebaI$upper.CL, length = 0.05, angle = 90, code = 3, col=1)
##polygon(c(pruebaI$Tmean,rev(pruebaI$Tmean)),c(pruebaI$upper.CL,rev(pruebaI$lower.CL)),col="grey",border="red") 
lines(lsmean ~ Tmean, pruebaI,col=1)
points(lsmean ~ Tmean, pruebaI, pch=20,col=1, cex=2)
arrows(pruebaJ$Tmean,pruebaJ$lower.CL,pruebaJ$Tmean,pruebaJ$upper.CL, length = 0.05, angle = 90, code = 3, col=2)
##polygon(c(pruebaJ$Tmean,rev(pruebaJ$Tmean)),c(pruebaJ$upper.CL,rev(pruebaJ$lower.CL)),col="grey",border="red") 
lines(lsmean ~ Tmean, pruebaJ,col=2)
points(lsmean ~ Tmean, pruebaJ, pch=20,col=2, cex=2)

#rect(par('usr')[1],par('usr')[4] , par('usr')[2], par('usr')[4]+((par('usr')[4]-par('usr')[3])/12), col = "grey90", xpd = TRUE)

title(main=i)

}


mtext('Average temperature (°C)', side = 1, line = 2, outer = TRUE, at = NA, 
          adj = NA, padj = NA, cex = 1.5, col = NA, font = NA)
mtext('lsmeans of Fertility values', side = 2, line = 2, outer = TRUE, at = NA,
           adj = NA, padj = NA, cex = 1.5, col = NA, font = NA)
dev.off()

################################################################################
## Relation of each chromosomal configuration and Tmean with fertility          #
## by lsmeans                                                                  #
################################################################################

meansCT <- lsmeans(lm1r1C, ~ (Te1H + Te6H + Tr1HB + Tr6HD)*Tmean*Ubicacion, at=list(Tmean=c(10,12,14,16,18,20,22,24)))


## Te1H Te6H Tr1HB Tr6HD Tmean   Ubicacion        lsmean         SE  df
##    0    0     0     0    10 invernadero -0.2744001244 0.10885329 958
##    1    0     0     0    10 invernadero  0.3744387872 0.07799682 958
##    0    1     0     0    10 invernadero  0.2946693489 0.08330872 958
##    1    1     0     0    10 invernadero  0.9435082605 0.08156856 958
##    0    0     1     0    10 invernadero  0.1993088959 0.08481251 958
##    1    0     1     0    10 invernadero  0.8481478075 0.10570712 958
##    0    1     1     0    10 invernadero  0.7683783692 0.07696964 958
##    1    1     1     0    10 invernadero  1.4172172808 0.12403137 958
##    0    0     0     1    10 invernadero  0.4598011493 0.08915566 958
##    1    0     0     1    10 invernadero  1.1086400609 0.10763823 958
##    0    1     0     1    10 invernadero  1.0288706226 0.10656613 958

##    1    1     0     1    10 invernadero  1.6777095342 0.14308116 958
##    0    0     1     1    10 invernadero  0.9335101695 0.08777818 958
##    1    0     1     1    10 invernadero  1.5823490811 0.14521926 958
##    0    1     1     1    10 invernadero  1.5025796428 0.12145978 958
##    1    1     1     1    10 invernadero  2.1514185544 0.18334319 958

meansCT <- as.data.frame(meansCT[[1]])
meansCT <- meansCT[!(meansCT$Te1H==0 & meansCT$Te6H==0 & meansCT$Tr1HB==0 & meansCT$Tr6HD==0),]

meansCT <- meansCT[!(meansCT$Te1H==1 & meansCT$Te6H==0 & meansCT$Tr1HB==1 & meansCT$Tr6HD==0),]

meansCT <- meansCT[!(meansCT$Te1H==1 & meansCT$Te6H==1 & meansCT$Tr1HB==1 & meansCT$Tr6HD==0),]

meansCT <- meansCT[!(meansCT$Te1H==1 & meansCT$Te6H==0 & meansCT$Tr1HB==0 & meansCT$Tr6HD==1),]

meansCT <- meansCT[!(meansCT$Te1H==0 & meansCT$Te6H==1 & meansCT$Tr1HB==0 & meansCT$Tr6HD==1),]

meansCT <- meansCT[!(meansCT$Te1H==1 & meansCT$Te6H==1 & meansCT$Tr1HB==0 & meansCT$Tr6HD==1),]

meansCT <- meansCT[!(meansCT$Te1H==1 & meansCT$Te6H==0 & meansCT$Tr1HB==1 & meansCT$Tr6HD==1),]

meansCT <- meansCT[!(meansCT$Te1H==0 & meansCT$Te6H==1 & meansCT$Tr1HB==1 & meansCT$Tr6HD==1),]

meansCT <- meansCT[!(meansCT$Te1H==1 & meansCT$Te6H==1 & meansCT$Tr1HB==1 & meansCT$Tr6HD==1),]


meansCT$Chro <- rep(c("Te1H","Te6H","Te1H+Te6H","Tr1HB","Te6H+Tr1HB","Tr6HD","Tr1HB+Tr6HD"), times=8)


## Corresponds to paper Figure 5

## c(bottom, left, top, right) defaultc(5, 4, 4, 2) +  0.1
pdf('../paper/img/fitted-fertility-chro-Tmean.pdf', width=9, height=6)
par(mar = c(4, 4, 3, -0.1)+0.1)
#par(tcl = -0.25)
#par(mgp = c(2, 0.6, 0))
par(mfrow=c(1,2))
## convert factor to numeric for convenience
meansCT$Chro2 <- as.numeric(as.factor(meansCT$Chro))
nChro <- max(meansCT$Chro2)


# get the range for the x and y axis
yrange <- range(meansCT$lsmean)
xrange <- range(meansCT$Tmean)

# set up the plot
plot(xrange, yrange, type="n", xlab="Temperature (°C)",
   ylab="LSmeans Fertility" )

colors <- rainbow(nChro)
linetype <- c(1:nChro)
plotchar <- seq(18,18+nChro,1)

## Field ##

## add lines
for (i in 1:nChro) {
  tree <- meansCT[(meansCT$Ubicacion=='jaulon' & meansCT$Chro2==i),]
  lines(tree$Tmean, tree$lsmean, type="b", lwd=1.5,lty=linetype[i], col=colors[i], pch=plotchar[i])
}

# add a title and subtitle
title("Field")

# add a legend
legend(xrange[1], yrange[2],levels(as.factor(meansCT$Chro)), cex=0.8, col=colors,
   pch=plotchar, lty=linetype)

## Greenhouse ##
par(mar = c(4, -0.1, 3, 4)+0.1)
plot(xrange, yrange, type="n", xlab="Temperature (°C)",   yaxt    = 'n',
   ylab="" )
## add lines
for (i in 1:nChro) {
  tree <- meansCT[(meansCT$Ubicacion=='invernadero' & meansCT$Chro2==i),]
  lines(tree$Tmean, tree$lsmean, type="b", lwd=1.5,lty=linetype[i], col=colors[i], pch=plotchar[i])
}

# add a title and subtitle
title("Greenhouse")

# add a legend
##legend(xrange[1], yrange[2],levels(as.factor(meansCT$Chro)), cex=0.8, col=colors, pch=plotchar, lty=linetype)

dev.off()

################################################################################

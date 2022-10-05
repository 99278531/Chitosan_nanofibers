# Chitosan_nanofibers
Repositorio para el estudio de Efectos terapéuticos de nanofibras electrohiladas de quitosano en heridas cutáneas de animales: una revisión sistemática y metanálisis.


Este es el repositorio del script utilizado para los calculos del proyecto: Efectos terapéuticos de nanofibras electrohiladas de quitosano en heridas cutáneas de animales: una revisión sistemática y metanálisis.

library(readxl)
file.choose()
ruta <-"C:\\Users\\alfin\\Desktop\\R77.xlsx"
S16403 <-read_excel(ruta)

library(metafor)
library(boot)
library(ggplot2)
library(boot)
library(readxl)

meta=S16403
rm(S16403)
meta[,12]<-NULL
str(meta)

ADGmd <-escalc( measure="MD", 
                m1i=Treated.mean, 
                sd1i=Treated.SD, 
                n1i=Treated.N, 
                m2i=Control.Mean, 
                sd2i=Control.SD, 
                n2i=Control.N, 
                data=meta)
ADGmd

ADGmetamd<-rma(yi, vi, data=ADGmd) 


summary(ADGmetamd)

ADGmetamd$I2

confint(ADGmetamd)

forest(ADGmetamd, order="obs") 



forest(ADGmetamd, xlim=c(-1.6,1.6),order="obs" , digits=c(2,1), cex=.8)



b_res <- rma(yi, vi, data=ADGmd, slab=ADGmd$Source)

qqnorm(b_res, main = "Random-Effects Model")

radial(b_res, main = "Random-Effects Model")

baujat(b_res)
inf <- influence(b_res)
print(inf)
plot(inf)

forest(b_res, xlim=c(-100.0,60), showweights=TRUE, order="obs",digits=c(2,1), cex=.8, 
psize=1, col="blue",alim=c(-60.0,25.0), xlab= "Efecto del grupo con nanofibras electrohiladas de quitosano")
text(-94,39,"Autor, año", cex=.8,font=2)
text(38,39,"Pesos", cex=.8,font=2)
text(50,39," Alfa [ IC - 95% ]", cex=.8,font=2)
par(mar=c(4,4,1,2))


inf_b_res<- influence(b_res)
dev.off()
par(mar=rep(2, 4)) 
plot(inf_b_res)


regtest(b_res)


### fit mixed-effects model with absolute latitude as predictor
res <- rma(yi, vi, mods = ~ ablat, data=dat)

### draw plot
regplot(res, xlim=c(10,60), predlim=c(10,60), xlab="Absolute Latitude", refline=0,
        atransf=exp, at=log(seq(0.2,1.6,by=0.2)), digits=1, las=1, bty="l",
        label=c(4,7,12,13), offset=c(1.6,0.8), labsize=0.9)

######################################
#########################################

funnel.meta(b_res,
            xlim = c(-0.5, 2),
            studlab = TRUE)

funnel(b_res, xlab = "Correlation coefficient")
regtest(b_res)
ranktest(b_res)



######sensibilidad 

leave1out(b_res, digits = 3)
plot(leave)



###################################
####metaregresion





dat=ADGmd
res.modq <- rma(yi, vi, mods = ~ concentration, data=dat) 
res.modq 
 
res.modpol <- rma(yi, vi, mods = ~ factor(polymer_sintetic), data=dat) 
res.modpol

res.modnat <- rma(yi, vi, mods = ~factor(natural), data=dat) 
res.modnat

res.modant <- rma(yi, vi, mods = ~factor(antibacterial), data=dat) 
res.modant

res.modmetal <- rma(yi, vi, mods = ~factor(metal), data=dat) 
res.modmetal

res.modantiox<- rma(yi, vi, mods = ~factor(antioxidant), data=dat) 
res.modantiox

res.modother <- rma(yi, vi, mods = ~factor(other), data=dat) 
res.modother


res.modnum<- rma(yi, vi, mods = ~ num_comb, data=dat) 
res.modnum


### draw plot


regplot(res.modpol,  xlab="NFQS + Polímero sintetico", ylab="Efecto en herida cutanea (DM)", labsize=11)
par(mar=c(4,4,1,2))



qqnorm(res.modpol, main = "Random-Effects Model")
################################

res.modantiox<- rma(yi, vi, mods = ~factor(antioxidant), data=dat)

res.modantiox

### fit mixed-effects model with absolute latitude as predictor
res <- rma(yi, vi, mods = ~ ablat, data=dat)

### draw plot
regplot(res.modantiox, xlim=c(10,60), predlim=c(10,60), xlab="NFQS + Antioxidante", refline=0,
        atransf=exp, at=log(seq(0.2,1.6,by=0.2)), digits=1, las=1, bty="l",
        label=c(4,7,12,13), offset=c(1.6,0.8), labsize=0.9)


regplot(res.modantiox, xlim=c(10,60), predlim=c(10,60), xlab="NFQS + Antioxidante")
######################################

### random-effects model using the log relative risks and then transformation
### of the estimated average log relative risk to the average relative risk

res <- rma(yi, vi, data = dat)
predict(res, transf = exp, digits = 2)

### externally standardized residuals

res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)
rstudent(b_res)

### influential case diagnostics
res <- rma(yi, vi, mods = cbind(ablat, year), data = dat)
inf <- influence(b_res)
plot(inf, plotdfb = TRUE)



###################################################
regplot(res.modnum, xlim=c(10,60), predlim=c(10,60), xlab="Number of combinations", refline=0,
        atransf=exp, at=log(seq(0.2,1.6,by=0.2)), digits=1, las=1, bty="l",
        label=c(4,7,12,13), offset=c(1.6,0.8), labsize=0.9)
res1 <- rma(yi, vi, data = dat, method = "FE")

radial(res1, main = "Fixed-Effects Model")

res2 <- rma(yi, vi, data = dat, method = "REML")

radial(res2, main = "Random-Effects Model")


res <- rma(yi, vi, data = dat)
qqnorm(res, main = "Random-Effects Model")

res <- rma(yi, vi, mods = cbind(ablat), data = dat)
qqnorm(res, main = "Mixed-Effects Model")



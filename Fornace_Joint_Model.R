##################################################################################################
################### Joint Model of Infection and Detection #######################################
##################################################################################################

library(fields)
library(INLA)


## Load data
# Dataset including GPS coordinates, covariates, detected and infected status by location

## Z = probability detected
z <- df$screened

## Y = probability positive
y <- ifelse(df$pcr > 0, 1, df$pcr)

## Set intercepts
z.b0 <- rep(1, length(z))
y.b0 <- rep(1, length(y))

## Mean-centre and scale selected environmental covariates
# z = pop + fric + fr + asp
# y = rds + ups + twi + bc7 + cfd 
covar <- df[c("pop", "fric", "fr", "asp", "rds", "ups", "twi", "bc7", "cfd")]
mean_covariates = apply(covar,2,mean)
sd_covariates = apply(covar,2,sd)
covar =
  scale(covar,
        mean_covariates, sd_covariates)
covar <- data.frame(covar)

## Extract covariates

# Pop
z.pop <- covar$pop
y.pop <- covar$pop
# Fric
z.fric <- covar$fric
y.fric <- covar$fric
# Forest
z.fr <- covar$fr
y.fr <- covar$fr
# Aspect
z.asp <- covar$asp
y.asp <- covar$asp
# Roads
z.rds <- covar$rds
y.rds <- covar$rds
# Upslope area
z.ups <- covar$ups
y.ups <- covar$ups
# TWI
z.twi <- covar$twi
y.twi <- covar$twi
# Climate
z.bc7 <- covar$bc7
y.bc7 <- covar$bc7

# Closed canopy forest
z.cfd <- covar$cfd
y.cfd <- covar$cfd

## Spatial component
## Scale coordinates
xmin0 <- min(df$x)
ymin0 <- min(df$y)
df$x <- df$x - xmin0
df$y <- df$y - ymin0
df$x <- df$x/1000
df$y <- df$y/1000

## Set domain
xmin <- min(df$x)
ymin <- min(df$y)
xmax <- max(df$x)
ymax <- max(df$y)
x <- c(xmin, xmin, xmax, xmax)
y <- c(ymin, ymax, ymin, ymax)
boundary <- cbind(x,y)

## Create mesh
pts <- cbind(df$x, df$y)
mesh1 <- inla.mesh.create.helper(points=pts, max.edge=c(1,2), cut=0.5)
plot(mesh1)
points(pts[,1], pts[,2], pch=19, cex=.5, col="red")

## Create observation matrix
A <- inla.spde.make.A(mesh1, loc=pts)

## SPDE without penalised priors
spde <- inla.spde2.matern(mesh=mesh1, alpha=2)

## Y = probability positive
y <- ifelse(df$pcr > 0, 1, df$pcr)

############################## Model detection probability (z) ###############################

# Data
dat.z <- list(z.b0=z.b0, z.fric=z.fric, z.fr=z.fr, z.asp=z.asp, z.pop=z.pop, z.rds=z.rds, z.bc7=z.bc7,
              z.ups=z.ups, z.twi=z.twi, z.cfd=z.cfd)

# Priors
fixed.priors.z <- list(mean = list(z.b0=0, z.pop=0, z.fric=0, z.fr=0, z.asp=0, z.rds= 0, z.bc7=0, z.ups=0,
                                   z.twi=0, z.cfd=0), 
                     prec=list(z.b0 = 1/100, z.pop=1/100, z.fric=1/100, z.fr=1/100, z.asp=1/100, z.rds= 1/100, 
                               z.bc7=1/100, z.ups=1/100, z.twi=0, z.cfd=1/100))

## Create data stack 
stk.z <- inla.stack(data=list(z=z, y=cbind(z, NA)), A=list(A,1), tag="est.z",
                   effects=list(list(i.z=1:spde$n.spde, i.zc=1:spde$n.spde), 
                   list(data.frame(dat.z))))

## Fit model without spatial component
form0.z <- z ~ 0 + z.b0 + z.pop + z.fr + z.asp + z.fric 
res0.z  <- inla(form0.z, family = "binomial", data=inla.stack.data(stk.z),
                control.predictor = list(A=inla.stack.A(stk.z), compute=TRUE),
                control.fixed=fixed.priors.z, control.family=list(link="logit"),
                control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE))
summary(res0.z)

## Model with spatial component
form.z <- z ~ 0 + z.b0 + z.pop + z.fr + z.asp + z.fric + f(i.z, model=spde)
res.z <- inla(form.z, family = "binomial", data=inla.stack.data(stk.z),
              control.predictor = list(A=inla.stack.A(stk.z), compute=TRUE), 
              control.fixed=fixed.priors.z, control.family = list(link="logit"),
              control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE))
summary(res.z)

############################## Model infection probability (y) ###############################

# Data
dat.y <- list(y.b0=y.b0, y.pop=y.pop, y.fric=y.fric, y.fr=y.fr, y.asp=y.asp, y.rds=y.rds, y.bc7=y.bc7,
              y.ups=y.ups, y.twi=y.twi, y.cfd=y.cfd)

# Priors
fixed.priors.y <- list(mean = list(y.b0=0, y.pop=0, y.fric=0, y.fr=0, y.asp=0, y.rds= 0, y.bc7=0, y.ups=0,
                                   y.twi=0, y.cfd=0),
                       prec=list(y.b0 = 1/100, y.pop=1/100, y.fric=1/100, y.fr=1/100, y.asp=1/100, y.rds= 1/100,
                                 y.bc7=1/100, y.ups=1/100, y.twi=1/100, y.cfd=1/100))

## Create data stack
stk.y <- inla.stack(data=list(r=y, y=cbind(NA, y)), A=list(A,1), tag="est.y", 
                    effects=list(i.y=1:spde$n.spde, list(data.frame(dat.y))))

## Fit model without spatial component lo + bc7 + rds + fric
form0.y <- r ~ 0 + y.b0 + y.rds + y.bc7 + y.ups + y.twi + y.cfd
res0.y  <- inla(form0.y, family = "binomial", data=inla.stack.data(stk.y),
                control.predictor = list(A=inla.stack.A(stk.y), compute=TRUE),
                control.fixed=fixed.priors.y, control.family=list(link="logit"),
                control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE))
summary(res0.y)

## Model with spatial component
form.y <- r ~ 0 + y.b0 + y.rds + y.bc7 + y.ups + y.twi + y.cfd + f(i.y, model=spde)
res.y <- inla(form.y, family = "binomial", data=inla.stack.data(stk.y),
              control.predictor = list(A=inla.stack.A(stk.y), compute=TRUE),
              control.fixed=fixed.priors.y, control.family = list(link="logit"),
              control.compute=list(dic=TRUE, cpo=TRUE, config=TRUE))
summary(res.y)

## Calculate spatial range
spde.est.y <- inla.spde2.result(inla = res.y, name = "i.y", spde = spde, do.transf = TRUE)
inla.zmarginal(spde.est.y$marginals.range.nominal[[1]])

######################################### Joint modelling ################################################

# z = pop + fric + fr + asp
# y = rds + ups + twi + bc7 + cfd 

## Joint data stack
stk.zy <- inla.stack(stk.z, stk.y)

## Joint priors
fixed.priors <- list(mean = list(y.b0=0, y.pop=0, y.fric=0, y.fr=0, y.asp=0, y.rds=0, y.ups=0, y.twi=0, y.bc7=0, y.cfd=0,
                                 z.b0=0, z.pop=0, z.fric=0, z.fr=0, z.asp=0, z.rds=0, z.ups=0, z.twi=0, z.bc7=0, z.cfd=0),
                     prec=list(z.b0 = 1/100, z.pop=0, z.fric=1/100, z.fr=1/100, z.asp=1/100, z.rds= 1/100, z.ups=1/100,
                               z.twi=1/100, z.bc7=1/100, z.cfd=1/100,
                               y.b0 = 1/100, y.pop=1/100, y.fric=1/100, y.fr=1/100, y.asp=1/100, y.rds= 1/100, y.ups=1/100,
                               y.twi=1/100, y.bc7=1/100, y.cfd=1/100))

## Joint non-spatial model
form.zy0 <- y ~ 0 + z.b0 + y.b0 + z.pop + z.fric + z.fr + z.asp + y.rds + y.ups + y.twi + y.bc7 + y.cfd 
res.zy0 <- inla(form.zy0, family=c('binomial', 'binomial'), data=inla.stack.data(stk.zy),
               control.compute=list(dic=TRUE, config=TRUE), control.fixed = fixed.priors,
               control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE))
summary(res.zy0)

## Joint model
form.zy <- y ~ 0 + z.b0 + y.b0 + z.pop + z.fric + z.fr + z.asp + y.rds + y.ups + y.twi + y.bc7 + y.cfd + f(i.z, model=spde) + f(i.y, model=spde)
res.zy <- inla(form.zy, family=c('binomial', 'binomial'), data=inla.stack.data(stk.zy),
               control.compute=list(dic=TRUE, config=TRUE), control.fixed = fixed.priors,
               control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE))
summary(res.zy)

## Create data stack
stk.z <- inla.stack(data=list(z=z, y=cbind(z, NA)), A=list(A,1), tag="est.z",
                    effects=list(list(i.z=1:spde$n.spde, i.zc=1:spde$n.spde),
                    list(data.frame(dat.z))))

## Joint data stack
stk.zy <- inla.stack(stk.z, stk.y)

## Joint model with shared spatial component
form.zyc <- y ~ 0 +z.b0 + y.b0 + z.pop + z.fric + z.fr + z.asp + y.rds + y.ups + y.twi + y.bc7 + y.cfd +
  f(i.z, model=spde) + f(i.y, model=spde) + f(i.zc, copy="i.y", fixed=FALSE)
res.zyc <- inla(form.zyc, family=c('binomial', 'binomial'), data=inla.stack.data(stk.zy),
                control.compute=list(dic=TRUE, config=TRUE), control.fixed = fixed.priors,
                control.family=list(list(link="logit"),list(link="logit")),
                control.predictor=list(A=inla.stack.A(stk.zy), compute=TRUE))

summary(res.zyc)

## Compare with combined DICs
rbind(separate=c(z=res.z$dic$dic, y=res.y$dic$dic),
      joint=tapply(res.zyc$dic$local.dic, res.zyc$dic$family, sum))

## Compare with combined DICs

# Combined with joint
rbind(separate=c(z=res.z$dic$dic, y=res.y$dic$dic),
      joint=tapply(res.zyc$dic$local.dic, res.zyc$dic$family, sum))

# Combined separate spatial
rbind(separate=c(z=res.z$dic$dic, y=res.y$dic$dic),
      joint=tapply(res.zy$dic$local.dic, res.zy$dic$family, sum))



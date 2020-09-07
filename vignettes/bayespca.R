## ------------------------------------------------------------------------
set.seed(141)
I <- 100
J <- 20
V1 <- rnorm(I, 0, 50)
V2 <- rnorm(I, 0, 30)
V3 <- rnorm(I, 0, 10)
X <- matrix(c(rep(V1, 7), rep(V2, 7), rep(V3, 6)), I, J)
X <- X + matrix(rnorm(I * J, 0, 1), I, J)

## ------------------------------------------------------------------------
# Install and load package
# devtools::install_github("davidevdt/bayespca")
library(bayespca)

# De-activate data center and scaling;
ctrl <- vbpca_control(center = FALSE, scalecorrection = -1,
						plot.lowerbound = FALSE)
# Estimate vbpca with fixed prior variances (equal to 1)
# for the elements of W
mod1 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'fixed',
				control = ctrl, verbose = FALSE )

# Test the class of mod1:
is.vbpca(mod1)



## ------------------------------------------------------------------------
mod1$muW

## ------------------------------------------------------------------------
mod1$P

## ------------------------------------------------------------------------
mod1$elbo

mod1$time

## ------------------------------------------------------------------------
mod2 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'fixed',
				updatetau = TRUE, control = ctrl, verbose = FALSE )

mod2$muW

## ------------------------------------------------------------------------
mod2$Tau

## ------------------------------------------------------------------------
# Estimate the model
mod3 <- vbpca(X, D = 3, maxIter = 1e+03,
				priorvar = 'jeffrey', control = ctrl, verbose = FALSE )
mod3$muW


mod3$Tau 

## ------------------------------------------------------------------------
# Set hyperparameter values 
ctrl2 <- vbpca_control(center = FALSE, scalecorrection = -1,
                       plot.lowerbound = FALSE, 
					   alphatau = 2, betatau = .5)
					   
					   
					   
# Estimate the model 
mod4 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'invgamma', 
              control = ctrl2, verbose = FALSE )
			  
			  
mod4$muW 


mod4$Tau

 

## ------------------------------------------------------------------------
# Set hyperparameter values 
ctrl3 <- vbpca_control(center = FALSE, scalecorrection = -1,
                       plot.lowerbound = FALSE, 
					   alphatau = c(.5, 50, 3), betatau = c(.5, .01, 10), 
					   hypertype = 'component')
					   
					   
# Estimate the model 
mod5 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'invgamma', 
              control = ctrl3, verbose = FALSE )

			  
mod5$muW 



mod5$Tau 

## ------------------------------------------------------------------------
# Specify component-specific Gamma(.01, 10) hyperpriors on betatau 
ctrl4 <- vbpca_control(center = FALSE, scalecorrection = -1, 
                       plot.lowerbound = FALSE, 
					   alphatau = 1, betatau = 1,
					   gammatau = .01, deltatau = 10, 
					   hypertype = 'component')

					   
# Estimate the model 
mod6 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'invgamma', 
              control = ctrl4, verbose = FALSE )
			  
			  
mod6$muW 


mod6$Tau 

## ------------------------------------------------------------------------
mod6$priorBeta 

## ---- fig.cap = "Prior variances for the first 3 components."------------
# Fixed prior global variances, updated via Type-II maximum likelihood: 
mod7 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'fixed',
              updatetau = TRUE, control = ctrl, verbose = FALSE, 
			  global.var = TRUE)
			  
			  
mod7$muW 

mod7$Tau 

## ---- fig.cap = "Scree-plot for 10 components. "-------------------------
mod8 <- vbpca(X, D = 10, maxIter = 1e+03, priorvar = 'fixed',
				updatetau = TRUE, global.var = TRUE,
				control = ctrl, verbose = FALSE )

## ------------------------------------------------------------------------
# SVS, fixed priorInclusion and InverseGamma(5, 1) for tau, v0 = .005
ctrl5 <- vbpca_control(center = FALSE, scalecorrection = -1,
						plot.lowerbound = FALSE,
						alphatau = 5, betatau = 1,
						beta1pi = -1, v0 = 5e-03)
# Estimate the model with priorInclusion = 0.5
mod9 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'invgamma',
				SVS = TRUE, priorInclusion = 0.5, control = ctrl5,
				verbose = FALSE )
				
mod9$muW

## ------------------------------------------------------------------------
# SVS, priorInclusion with Beta(1,1) priors and InverseGamma(5, 1) for tau, v0 = .005
ctrl6 <- vbpca_control(center = FALSE, scalecorrection = -1,
						plot.lowerbound = FALSE, alphatau = 5,
						betatau = 1, beta1pi = 1, beta2pi = 1,
						v0 = 5e-03)
						
# Estimate the model
mod10 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'invgamma',
				SVS = TRUE, priorInclusion = 0.5, control = ctrl6,
				verbose = FALSE )
				
mod10$muW

## ------------------------------------------------------------------------
mod9$inclusionProbabilities

mod10$inclusionProbabilities

## ---- fig.show='hold', fig.cap = "True and Estimated inclusion probabilities.", fig.width=6, fig.height=4----
trueInclusions <- matrix(0, J, 3)
trueInclusions[1:7, 1] <- 1
trueInclusions[8:14, 2] <- 1
trueInclusions[15:20, 3] <- 1

par(mfrow=c(1,2))
image(1:ncol(trueInclusions), 1:nrow(trueInclusions),
		t(trueInclusions[J:1, ]), ylab = "", axes = FALSE,
		main = "True Inclusions", xlab = "",
		col = RColorBrewer::brewer.pal(9, "Blues"))
axis(side = 1, at = 1:3, labels = paste("Component ", 1:3 ))
axis(side = 2, at = 1:20, labels = paste("Var ", J:1 ))

fields::image.plot(1:ncol(trueInclusions), 1:nrow(trueInclusions),
					t(mod9$inclusionProbabilities[J:1, ]), ylab = "", axes = FALSE,
					main = "Estimated Inclusions", xlab = "",
					col = RColorBrewer::brewer.pal(9, "Blues"))
axis(side = 1, at = 1:3, labels = paste("Component ", 1:3 ))

## ------------------------------------------------------------------------
mod10$priorInclusion

## ------------------------------------------------------------------------
# Beta priors with different degrees of sparsity for each component
ctrl7 <- vbpca_control(center = FALSE, scalecorrection = -1,
						plot.lowerbound = FALSE,
						alphatau = 5, betatau = 1,
						beta1pi = c(0.01, 1, 10), beta2pi = 1,
						v0 = 5e-03)
# Estimate the model
mod11 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'invgamma', SVS = TRUE,
				priorInclusion = rep(0.5, 3), control = ctrl7, verbose = FALSE )
				
mod11$muW

mod11$priorInclusion

mod11$inclusionProbabilities

## ---- fig.show='hold', fig.cap = "High posterior density intervals. "----
# Set hyperparameter values and require 50% probability density intervals 
ctrl8 <- vbpca_control(center = FALSE, scalecorrection = -1, 
                        plot.lowerbound = FALSE, 
					    alphatau = 2, betatau = .5, 
					    hpdi = TRUE, probHPDI = 0.5)

# Estimate the model 
mod12 <- vbpca(X, D = 3, maxIter = 1e+03, priorvar = 'invgamma',
              control = ctrl8, verbose = TRUE )

# Plot HPD intervals for variables 1:10, component 1 
plothpdi(mod12, d = 1, vars = 1:10)


## ------------------------------------------------------------------------
PCs <- X %*% mod1$muW 
head(PCs, 15)


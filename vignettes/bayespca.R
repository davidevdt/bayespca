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

# Estimate vbpca with fixed prior precisions (equal to 1)
# for the elements of W
mod1 <- vbpca(X, D = 3, maxIter = 1e+03, alphatau=0,
              center = FALSE, scalecorrection = -1,
              plot.lowerbound = FALSE,
			  verbose = FALSE )

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
mod2 <- vbpca(X, D = 3, maxIter = 1e+03, alphatau=0,
             updatetau = TRUE, center = FALSE,
             scalecorrection = -1,
             plot.lowerbound = FALSE,
             verbose = FALSE )

mod2$muW

## ------------------------------------------------------------------------
mod2$Tau

## ------------------------------------------------------------------------
# Estimate the model
mod3 <- vbpca(X, D = 3, maxIter = 1e+03,
              alphatau = 2, betatau = .5,
              center = FALSE, scalecorrection = -1,
              plot.lowerbound = FALSE,
              verbose = FALSE )
mod3$muW


mod3$Tau

## ------------------------------------------------------------------------
# Estimate the model
mod4 <- vbpca(X, D = 3, maxIter = 1e+03,
              center = FALSE, scalecorrection = -1,
              alphatau = c(.5, 50, 3), betatau = c(.5, .01, 10),
              plot.lowerbound = FALSE,
              verbose = FALSE )

mod4$muW

mod4$Tau

## ---- fig.cap = "Prior precisions for the first 3 components."-----------
# Fixed prior global variances, updated via Type-II maximum likelihood:
mod5 <- vbpca(X, D = 3, maxIter = 1e+03, alphatau=0,
              updatetau = TRUE,  
              center = FALSE, scalecorrection = -1,
              plot.lowerbound = FALSE,
              verbose = FALSE, global.var = TRUE)

mod5$muW


mod5$Tau

## ---- fig.cap = "Scree-plot for 10 components. "-------------------------
mod6 <- vbpca(X, D = 10, maxIter = 1e+03, alphatau=0,
              updatetau = TRUE,  
              center = FALSE, scalecorrection = -1,
              plot.lowerbound = FALSE,
              verbose = FALSE, global.var = TRUE)

mod6$Tau

## ------------------------------------------------------------------------
mod7 <- vbpca(X, D = 3, maxIter = 1e+03, alphatau=0,
              updatetau = TRUE, center = FALSE,
              scalecorrection = -1,
              plot.lowerbound = FALSE,
              verbose = FALSE)

mod7$muW


mod7$Tau

## ---- fig.cap = "Heatmap of Tau. ", fig.width=5, fig.height=5.5----------
mat_mod_7 <- plotheatmap(mod7, matrix_type="Tau", bound_tau=50)
mat_mod_7$W
mat_mod_7$Tau 

## ---- fig.show='hold', fig.cap = "High posterior density intervals. ", fig.width=4, fig.height=6----
# Set hyperparameter values and require 90% probability density intervals
# Estimate the model
mod8 <- vbpca(X, D = 3, maxIter = 1e+03,
                alphatau = .001, betatau = .001,
                center = FALSE, scalecorrection = -1,
                hpdi = TRUE, probHPDI = 0.9,
                plot.lowerbound = FALSE,
                verbose = TRUE )

# Plot HPD intervals for variables 1:10, component 1
plothpdi(mod8, d = 1, vars = 1:20)

## ------------------------------------------------------------------------
PCs <- X %*% mod1$muW 
head(PCs, 15)


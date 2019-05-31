context("bayespca")

library(bayespca)

# Set seed
set.seed(4321)

# Tolerance
tol <- 1e-08


# Generate data 
I <- 1e+02
V1 <- rnorm(I, 0, 50)
V2 <- rnorm(I, 0, 30)
V3 <- rnorm(I, 0, 10)

X <- cbind(V1, V1, V1, V2, V2, V2, V3, V3)
X <- X + matrix(rnorm(I * 8, 0, 1), I, 8)

##################### bayespca vs. Ordinary PCA - With and without scaling 
pca_mod <- prcomp(X, center = TRUE, scale. = FALSE)
vbpca_ctrl <- vbpca_control(center = TRUE, scalecorrection = -1, plot.lowerbound = FALSE)
vbpca_mod <- vbpca(X, D = 3, tau = 1e+05, priorvar = "fixed", updatetau = FALSE, control = vbpca_ctrl)

# Test 1: VBPCA recovers PCA for low regularization (large tau)
testthat::test_that("Test 1: vbpca vs. pca, no scaling", {
    testthat::expect_equivalent(pca_mod[[2]][,1:3], vbpca_mod[[1]])
    testthat::expect_equivalent(t(vbpca_mod[[1]]) %*% vbpca_mod[[1]], diag(1, 3) )
})




pca_mod <- prcomp(X, center = TRUE, scale. = TRUE)
vbpca_ctrl <- vbpca_control(center = TRUE, scalecorrection = 1, plot.lowerbound = FALSE)
vbpca_mod <- vbpca(X, D = 3, tau = 1e+05, priorvar = "fixed", updatetau = FALSE, control = vbpca_ctrl)

# Test 2: VBPCA recovers PCA for low regularization (large tau)
testthat::test_that("Test 2: vbpca vs. pca, scaling", {
    testthat::expect_equivalent(pca_mod[[2]][,1:3], vbpca_mod[[1]])
    testthat::expect_equivalent(t(vbpca_mod[[1]]) %*% vbpca_mod[[1]], diag(1, 3) )
})






##################### bayespca vs. Ordinary PCA - Reconstructed data 
pred_pca <- pca_mod$x[,1:3] %*% t(pca_mod[[2]][,1:3])
pred_vbpca <- scale(X) %*% vbpca_mod$muW[,1:3] %*% t(vbpca_mod$P[,1:3])  

testthat::test_that("Test 3: vbpca vs. pca, reconstructed matrix", {
    testthat::expect_equivalent(pred_pca, pred_vbpca)
})

 
 
 
##################### bayespca - Regularization on tau 
vbpca_ctrl <- vbpca_control(center = TRUE, scalecorrection = -1, plot.lowerbound = FALSE)
vbpca_mod <- vbpca(X, D = 3, tau = 1e-05, priorvar = "fixed", updatetau = FALSE, control = vbpca_ctrl)
vbpca_mod2 <- vbpca(X, D = 3, tau = 1e+05, priorvar = "fixed", updatetau = FALSE, control = vbpca_ctrl)
 
testthat::test_that("Test 4: vbpca, effect of tau on absolute values of weights", {
    testthat::expect_equivalent( sum(abs(vbpca_mod[[1]]) < abs(vbpca_mod[[2]])), 3 * 8    )
})



##################### bayespca - Regularization with Inverse Gamma
vbpca_ctrl <- vbpca_control(center = TRUE, scalecorrection = -1, plot.lowerbound = FALSE,
                            alphatau = 100, betatau = 1e-02)
vbpca_ctrl2 <- vbpca_control(center = TRUE, scalecorrection = -1, plot.lowerbound = FALSE, 
                             alphatau = 1e-02, betatau = 1e+02)
vbpca_mod <- vbpca(X, D = 3, tau = 1, priorvar = "invgamma",, control = vbpca_ctrl)
vbpca_mod2 <- vbpca(X, D = 3, tau = 1, priorvar = "invgamma", control = vbpca_ctrl2)

testthat::test_that("Test 5: effect of InvGamma parameters", {
	testthat::expect_equal( sum(abs(vbpca_mod2[[1]])) -  sum(abs(vbpca_mod[[1]])) >= tol, TRUE )
}) 




##################### bayespca - Stochastic Variable Selection 
vbpca_ctrl <- vbpca_control(center = TRUE, scalecorrection = -1, plot.lowerbound = FALSE,
                            alphatau = 1, betatau = 1, beta1pi = -1)
vbpca_mod <- vbpca(X, D = 3, tau = 1, priorvar = "invgamma", SVS = TRUE, 
                   priorInclusion = 0.1, control = vbpca_ctrl)
vbpca_mod2 <- vbpca(X, D = 3, tau = 1, priorvar = "invgamma", SVS = TRUE,
                    priorInclusion = 0.9, control = vbpca_ctrl)

testthat::test_that("Test 6: effect of prior inclusion probability (SVS)", {
	testthat::expect_equal( sum( vbpca_mod[[9]] < vbpca_mod2[[9]] ), 3 * 8 )
}) 







vbpca_ctrl <- vbpca_control(center = TRUE, scalecorrection = -1, plot.lowerbound = FALSE,
                            alphatau = 1, betatau = 1, beta1pi = 1, beta2pi = 50, v0 = 1e-04)
vbpca_ctrl2 <- vbpca_control(center = TRUE, scalecorrection = -1, plot.lowerbound = FALSE,
                            alphatau = 1, betatau = 1, beta1pi = 320, beta2pi = 50, v0 = 1e-04)
							
vbpca_mod <- vbpca(X, D = 3, tau = 1, priorvar = "invgamma", SVS = TRUE, 
                   priorInclusion = 0.5, control = vbpca_ctrl)
vbpca_mod2 <- vbpca(X, D = 3, tau = 1, priorvar = "invgamma", SVS = TRUE,
                    priorInclusion = 0.5, control = vbpca_ctrl2)

					
testthat::test_that("Test 7: effect of Beta priors on inclusion probabilities (SVS)", {
	testthat::expect_equal( sum( vbpca_mod[[9]] < vbpca_mod2[[9]] ), 3 * 8 )
	testthat::expect_equal( colSums( vbpca_mod2[[9]] > .5 ), c(3, 3, 2) )
	testthat::expect_equal( colSums( vbpca_mod[[9]] > .5 ), c(0, 0, 0) )
}) 
		




		




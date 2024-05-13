
### --- Data simulation --- ###
n <- 10000
theta <- c(-4, 1, 2)

x <- cbind(rep(1, n), rnorm(n), rnorm(n))
score <-  x %*% theta
y <- ifelse(1/(1 + exp(-score) ) > runif(n), 1, 0)

w <- rchisq(n, df = 4) * rnorm(n, mean = mu_s, sd = sd_s) * 100

### --- Cost specification --- ###
a_cost <- .1
b_cost <- 10
c_cost <- 1

### --- Optimization --- ###
# --- Optimization parameters --- #
d <- 3
R <- 50
A <- diag(d)
B <- rep(R/d, d)
A <- rbind(A, - diag(d))
B <- c(B, rep(R/d, d))

# --- Obective function --- #
loss_theta <- function(theta, lambda = 0){
  n <- length(y)
  p <-  plogis(as.matrix(x)%*%theta )
  aec <- (c_cost*y*w + p*(-c_cost*y*w + (1-y)*w*a_cost + b_cost) ) / n 
  loss <- sum(aec + lambda * sum(abs(theta)))
  return( -loss )
}

# --- Starting value --- #
st.val <- glm(y ~ x - 1, family=binomial(link = "logit"))$coefficients

# --- Optimization process --- #
logit_lf_mod <- maxLik(loss_theta, 
                       method = "BFGS", 
                       start = st.val, 
                       control=list(tol=-1, reltol=1e-15, gradtol=1e-15),
                       constraints = list(ineqA = A, ineqB = B))

# --- Cost-sensitive estimator --- #
theta_hat <- logit_lf_mod$estimate
theta_hat

# --- Prediction --- #
pred <- plogis(x%*%c(theta_hat)) 



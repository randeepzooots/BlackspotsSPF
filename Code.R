library(readxl)
library(summarytools)
library(magrittr)
library(dplyr)
library(zipfR)
library(ggplot2)
library(kableExtra)
I_80_Compiled <- read_excel("Blackspot Data/I-80 Compiled.xlsx", 
                            sheet = "I-80 Compiled")
View(I_80_Compiled)
I_80_Compiled
colnames(I_80_Compiled)


cont_var <- read_excel("Blackspot Data/I-80 Compiled.xlsx", 
                            sheet = "Continuous ")
cat_var <- read_excel("Blackspot Data/I-80 Compiled.xlsx", 
                       sheet = "Categorical")


library(MASS)
exp <- crash_fi ~ lanes + lane_width + OSW + ISW + m_width + m_barrier +
  cz_width + r_barrier + en_ramp_i + ex_ramp_i + w_type_i + en_ramp_d + 
  ex_ramp_d + w_type_d + aadt

summary(m1 <- glm.nb(exp,data = I_80_Compiled))
m1$fitted.values
plot(I_80_Compiled$crash_fi,type = 'l',col = 'blue')
lines(m1$fitted.values,col = 'red')
m1$coefficients

hist(I_80_Compiled$crash_fi)



nll <- function(beta)
{
  -sum(Y*log(alpha)+Y*(dm%*%beta)-(Y+(1/alpha))*log(1+alpha*exp(dm%*%beta))+
         log(gamma(Y+(1/alpha)))-log(gamma(Y+1))-log(gamma(1/alpha)))
}

dm <- model.matrix(~lanes + lane_width + OSW + ISW + m_width + m_barrier +
                     cz_width + r_barrier + en_ramp_i + ex_ramp_i + w_type_i + 
                     en_ramp_d + ex_ramp_d + w_type_d + aadt,data = X)

I_80_Compiled$lanes <- as.factor(I_80_Compiled$lanes)

b_est <- optim(beta,nll)
b_est

X <- subset(I_80_Compiled,select = -c(segment,yr,RSOW,RSIW,crash_pdo))



colnames(dm)
alpha <- m1$theta
Y <- X['crash_fi']

cbind(r$par,m1$coefficients)



op <- exp(dm%*%r$par)

plot(op,type = 'l',col = 'red')
lines(m1$fitted.values, col = 'green')

X_p <- mutate(X,Indicator = ifelse(crash_fi > 3,1,0))
X_p <- mutate(X_p,Ibeta1 = c(Ibeta((1/alpha)/((1/alpha)+exp(dm%*%beta)),1/alpha,crash_fi+1)))
X_p <- mutate(X_p,Ibeta2 = beta(1/alpha,crash_fi+1))

I <- X_p$Indicator
nll_c <- function(beta)
{
  -sum(I*(Y*log(alpha)+Y*(dm%*%beta)-(Y+(1/alpha))*log(1+alpha*exp(dm%*%beta))+
            log(gamma(Y+(1/alpha)))-log(gamma(Y+1))-log(gamma(1/alpha)))+(1-I)*
         (log(Ibeta1)-log(Ibeta2)))
}
b_est_c <- optim(beta,nll_c)
Ibeta(1/alpha/(1/alpha+0.6),1/alpha,Y+1)
Ibeta1 <- c(Ibeta((1/alpha)/((1/alpha)+exp(dm%*%beta)),1/alpha,X$crash_fi+1))
Ibeta2 <- beta(1/alpha,X$crash_fi+1)
cbind(r_c$par,r$par)

params <- cbind(m1$coefficients,b_est$par,b_est_c$par)
colnames(params) <- c('glm.nb','MLE','MLE Censored')

Y_hat <- data.frame(cbind(m1$fitted.values,y_hat,y_hat_c))
ggplot(Y_hat)+geom_line(aes(1:NROW(Y_hat),X1))+
  geom_line(aes(1:NROW(Y_hat),X2),color = 'blue')+
  geom_line(aes(1:NROW(Y_hat),X3),color = 'red')+
  xlab('')+ylab('Crashes')+theme_bw()



op_c <- exp(dm%*%r_c$par)
plot(op,type = 'l',col = 'red')
lines(op_c, col = 'green')
lines(Y,col = 'light blue')

jpeg('Comparison1.jpg',width = 1080,height = 720,res = 600)



function (mu.link = "log",
          sigma.link = "log")
{
  mstats <- checklink(
    "mu.link",
    "Negative Binomial type I",
    substitute(mu.link),
    c("inverse", "log", "identity",
      "sqrt")
  )
  dstats <- checklink(
    "sigma.link",
    "Negative Binomial type I",
    substitute(sigma.link),
    c("inverse", "log", "identity",
      "sqrt")
  )
  structure(
    list(
      family = c("NBIlc", "left censored Negative Binomial type I"),
      parameters = list(mu = TRUE, sigma = TRUE),
      nopar = 2,
      type = "Discrete",
      mu.link = as.character(substitute(mu.link)),
      sigma.link = as.character(substitute(sigma.link)),
      mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
      mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
      mu.dr = mstats$mu.eta,
      sigma.dr = dstats$mu.eta,
      dldm = function (y, mu, sigma)
        attr(
          gamlss::numeric.deriv(dNBIlc(y, mu, sigma, log = TRUE),
                                "mu", delta = NULL),
          "gradient"
        ),
      d2ldm2 = function (mu,
                         sigma)
      {
        -1 / (mu * (1 + mu * sigma))
      },
      dldd = function (y, mu, sigma)
        attr(
          gamlss::numeric.deriv(dNBIlc(y, mu, sigma, log = TRUE),
                                "sigma", delta = NULL),
          "gradient"
        ),
      d2ldd2 = function (y,
                         mu, sigma)
      {
        dldd <- -((1 / sigma) ^ 2) * (digamma(y[, 1] + (1 / sigma)) -
                                        digamma(1 /
                                                  sigma) - log(1 + mu * sigma) - (y[,
                                                                                    1] - mu) * sigma /
                                        (1 + mu * sigma))
        d2ldd2 <-
          -dldd ^ 2
        d2ldd2 <-
          ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)
        d2ldd2
      },
      d2ldmdd = function (y)
        rep(0, length(y[, 1])),
      G.dev.incr = function (y, mu,
                             sigma, ...)
        - 2 * dNBIlc(
          y,
          mu = mu,
          sigma = sigma,
          log = TRUE
        ),
      rqres = expression(
        rqres(
          pfun = "pNBIlc ",
          censored = "left",
          type = "Discrete",
          ymin = 0,
          y = y,
          mu = mu,
          sigma = sigma
        )
      ),
      mu.initial = expression(mu <- (y[, 1] + mean(y[, 1])) / 2),
      sigma.initial = expression(sigma <- rep(max(((var(y[,
                                                          1]) - mean(y[, 1])) /
                                                     (mean(y[, 1]) ^ 2)
      ), 0.1), length(y[,
                        1]))),
      mu.valid = function(mu)
        all(mu > 0),
      sigma.valid = function(sigma)
        all(sigma >
              0),
      y.valid = function (y)
        all(y[, 1] >= 0),
      mean = function(mu, sigma)
        mu,
      variance = function(mu,
                          sigma)
        mu + sigma * mu ^ 2
    ),
    class = c("gamlss.family",
              "family")
  )
}
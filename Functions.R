
#' Fit degpd and zidegpd

p.G <- function(u, type = 1, prob, kappa, delta) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (type == 0) {
    return(u)
  } else if (type == 1) {
    return(u^kappa)
  } else if (type == 2) {
    F.min <- 1-pnorm(1, mean = 0, sd = 1/sqrt(kappa)) #Phi(sqrt(kappa)*(a-u)--->Phi(sqrt(kappa)*(0-1)= 1-Phi(sqrt(kappa)) cdf of standard normal dist, denonirator Eq.7 A flexible extended generalized Pareto distribution for tail estimation
    F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))#Phi(sqrt(kappa)*(b-u)--->Phi(sqrt(kappa)*(1-1)=1/2
    p <- (pnorm(u, mean = 1, sd = 1/sqrt(kappa)) - F.min)/(F.max -F.min)
    #p <- (pnorm(q, mean = 1, sd = 1/sqrt(kappa)) - F.min)/(F.max -F.min)
    return(p)
  }else if (type == 3) {
    lower=1/32
    upper=1/2
    stopifnot(lower <= upper,kappa > 0)
    aa <- rep(lower, length(u))
    bb <- rep(upper, length(u))
    normalize.factor <- pbeta(bb, kappa, kappa) - pbeta(aa, kappa, kappa)
    #tt <- pbeta(apply(cbind(apply(cbind((upper-lower)*u +lower , bb), 1, min), aa), 1, max), kappa, kappa)
    tt <- pbeta((upper-lower)*u +lower, kappa, kappa)
    tt <- tt / normalize.factor
    return(tt)
  }else if (type == 4) {
    return(1 - pbeta((1 - u)^delta, 1/delta, 2))
  } else if (type == 5) {
    return((1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2))
  } else if (type == 6) {
    return(prob * u^kappa + (1 - prob) * u^delta)
  }
}

#p.G(runif(10), type = 5, delta = 2, kappa = 1, prob = 0.1)

d.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA, log = FALSE) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (log == FALSE) {
    if (type == 0) {
      return(1)
    } else if (type == 1) {
      return(kappa * u^(kappa - 1))
    } else if (type == 2) {
    F.min <-1- pnorm(1, mean = 0, sd = 1/sqrt(kappa))
    F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
    d<- dnorm(u,mean = 1, sd = 1/sqrt(kappa))
    den <- (sqrt(kappa)*d)/(F.max -F.min)
    return(den)
    }else if (type == 3) {
      lower=1/32
      upper=1/2
      stopifnot(lower <= upper,kappa > 0)
      tt <- rep(0, length(u))
      normalize.factor <- pbeta(upper, kappa, kappa) - pbeta(lower, kappa, kappa)
      tt[u >= lower & u <= upper] <- dbeta(u[u >= lower & u <= upper],
                                           kappa, kappa) / normalize.factor
      return(tt)
    }
    else if (type == 4) {
      return(dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 - u)^(delta - 1))
    } else if (type == 5) {
      return((kappa/2) * (1 - pbeta((1 - u)^delta, 1/delta, 2))^(kappa/2 - 1) * dbeta((1 - u)^delta, 1/delta, 2) * delta * (1 -
                                                                                                                              u)^(delta - 1))
    } else if (type == 6) {
      return(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1))
    }
  } else {
    if (type == 0) {
      return(0)
    } else if (type == 1) {
      return(log(kappa) + (kappa - 1) * log(u))
    }else if (type == 2) {
      F.min <-1- pnorm(1, mean = 0, sd = 1/sqrt(kappa), log.p = TRUE)
      F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa), log.p = TRUE)
      d<- dnorm(u,mean = 1, sd = 1/sqrt(kappa), log.p = TRUE)
      den <- (sqrt(kappa)*d)/(F.max -F.min)
      return(den)
    }else if (type == 3) {
      lower=1/32
      upper=1/2
      stopifnot(lower <= upper,kappa > 0)
      tt <- rep(0, length(u))
      normalize.factor <- pbeta(upper, kappa, kappa, log = TRUE) - pbeta(lower, kappa, kappa, log = TRUE)
      tt[u >= lower & u <= upper] <- dbeta(u[u >= lower & u <= upper],
                                           kappa, kappa,log = TRUE) / normalize.factor
      return(tt)
    }else if (type == 4) {
      return(dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) + log(delta) + (delta - 1) * log(1 - u))
    } else if (type == 5) {
      return(log(kappa/2) + (kappa/2 - 1) * log(1 - pbeta((1 - u)^delta, 1/delta, 2)) + dbeta((1 - u)^delta, 1/delta, 2, log = TRUE) +
               log(delta) + (delta - 1) * log(1 - u))
    } else if (type == 6) {
      return(log(prob * kappa * u^(kappa - 1) + (1 - prob) * delta * u^(delta - 1)))
    }
  }
}
# x<-runif(10)
# d.G(x,  type=2,kappa=2)
##
q.G <- function(u, type = 1, prob = NA, kappa = NA, delta = NA) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (type == 0) {
    return(u)
  } else if (type == 1) {
    return(u^(1/kappa))
  } else if (type == 2) {
    F.min <- 1-pnorm(1, mean = 0, sd = 1/sqrt(kappa))
    F.max <- pnorm(1, mean = 1, sd = 1/sqrt(kappa))
    q <- (qnorm(u* (F.max - F.min) + (F.min), mean = 0, 
                sd = 1))/sqrt(kappa) +1
    return(q)
  }else if (type == 3) {
    lower=1/32
    upper=1/2
    tt <- u
    pin <- pbeta(lower, kappa, kappa) + u * (pbeta(upper, kappa, kappa) - pbeta(lower, kappa, kappa))
    tt <- (qbeta(pin, kappa, kappa)-lower)/(upper-lower)
    return(tt)
  }else if (type == 4) {
    return(1 - qbeta(1 - u, 1/delta, 2)^(1/delta))
  } else if (type == 5) {
    return(1 - qbeta(1 - u^(2/kappa), 1/delta, 2)^(1/delta))
  } else if (type == 6) {
    dummy.func <- function(u, p, prob = NA, kappa = NA, delta = NA) {
      return(p.G(u = u, prob = prob, kappa = kappa, delta = delta, type = 4) - p)
    }
    find.root <- function(u, prob = NA, kappa = NA, delta = NA) {
      return(uniroot(dummy.func, interval = c(0, 1), p = u, prob = prob, kappa = kappa, delta = delta)$root)
    }
    return(sapply(u, FUN = find.root, prob = prob, kappa = kappa, delta = delta))
  }
}
#q.G(runif(10), kappa = 2, type = 4)
##
r.G <- function(n, prob = NA, kappa = NA, delta = NA, type = 1, unifsamp = NULL, direct = FALSE) {
  if (!type %in% 1:6) {
    stop("Invalid `type' argument")
  }
  type <- type[1]
  if (type %in% c(1,2, 3, 5,6) && missing(kappa)) {
    stop("Argument `kappa' missing.")
  }
  if (type %in% c(4, 5, 6) && missing(delta)) {
    stop("Argument `delta' missing.")
  }
  if (type == 6 && missing(prob)) {
    stop("Argument `prob' missing.")
  }
  if (is.null(unifsamp)) {
    unifsamp <- runif(n)
  }
  if (type != 6 | (type == 6 & direct)) {
      return(q.G(unifsamp, prob = prob, kappa = kappa, delta = delta, type = type))
  } else if (type == 6 & !direct) {
      components <- sample(x = c(1, 2), size = n, replace = TRUE, prob = c(prob, 1 - prob))
      res <- c()
      res[components == 1] <- q.G(unifsamp[components == 1], prob = NA, kappa = kappa, delta = NA, type = 1)
      res[components == 2] <- q.G(unifsamp[components == 2], prob = NA, kappa = delta, delta = NA, type = 1)
      return(res)
  }
}

#r.G(10, type = 6, delta = 1, kappa = 1, prob = 0.1)
############################################
##    EXTENDED GPD TYPE 1 to 4           ##
############################################

pegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(p.G(extraDistr::pgpd(q, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type))
}

pegpd(runif(10), type = 3, kappa = 1, sigma = 1, xi=0.1)
##
degpd <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, log = FALSE) {
  if (log == FALSE) {
    return(d.G(extraDistr::pgpd(x, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type) * extraDistr::pgpd(x, sigma = sigma,
                                                                                                                               xi = xi))
  } else {
    return(d.G(extraDistr::pgpd(x, sigma = sigma, xi = xi), prob = prob, kappa = kappa, delta = delta, type = type, log = TRUE) + extraDistr::dgpd(x,
                                                                                                                                           sigma= sigma, xi = xi, log = TRUE))
  }
}

degpd(runif(10), type = 2, kappa = 1, sigma = 1, xi=0.3)
#
qegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(extraDistr::qgpd(q.G(p, prob = prob, kappa = kappa, delta = delta, type = type), sigma = sigma, xi = xi))
}

qegpd(runif(10), type = 3, kappa = 1, sigma = 1, xi=0.3)


regpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL, censoring = c(0, Inf)) {
    return(extraDistr::qgpd(r.G(n, prob = prob, kappa = kappa, delta = delta, type = type, unifsamp), sigma = sigma,
                            xi = xi))
}

regpd(10, type = 2, kappa = 1, sigma = 1, xi=0.3)

############################################
##   Discrete EXTENDED GPD TYPE 1 to 4    ##
############################################
##################Density function##########

ddiscegpd <- function(x, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(pegpd(x+1,prob = prob , kappa =kappa, delta = delta, sigma = sigma, xi = xi, type = type)- pegpd(x,prob = prob , kappa =kappa, delta = delta, sigma = sigma, xi = xi, type = type))
}
ddiscegpd(runif(10),type = 3, kappa = 1, sigma = 1, xi=0.3)
##################CDF################################
pdiscegpd <- function(q, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  return(pegpd(q+1,prob = prob , kappa =kappa, delta = delta, sigma = sigma, xi = xi, type = type))
}
################Quantile function###############################
qdiscegpd <- function(p, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1) {
  q<-ceiling(qegpd(p, prob = prob, kappa = kappa, delta = delta, type = type, sigma = sigma, xi=xi))-1
  q[q < 0] <- 0
  if(is.matrix(p)) {
    q <- matrix(q, ncol = ncol(p), nrow = nrow(p))
    colnames(q) <- colnames(p)
    rownames(q) <- rownames(p)
  }
  return(q)
}


rdiscegpd <- function(n, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL) {
  return(floor(regpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi,type = type, unifsamp)))
}

# x<-rdiscegpd(10000, type = 3, kappa = 10, sigma = 1, xi=0.3)
#################################################################################################
################Zero Inflated Discrete EGPD  ###################################################
################################################################################################
#--------------------------------------------------------------------------------------------------
dzidiscegpd<-function(x,  pi = NA,prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1)
{
  if (type==1) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==2) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==3) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==4) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
  }
  if (type==5) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==6) {
    if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (any(pi <= 0) | any(pi >= 1) )  stop(paste("pi must be between 0 and 1", "\n", ""))
  #if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
  #if (any(kappa <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  #if (any(delta <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", ""))
  if (any(xi <= 0) )  stop(paste("xi must be greater than 0", "\n", ""))
  if (any(x < 0) )  stop(paste("x must be 0 or greater than 0", "\n", ""))
  ly <- max(length(x),length(pi),length(prob),length(kappa),length(delta),length(sigma),length(xi))
  x <- rep(x, length = ly)
  pi <- rep(pi, length = ly)
  prob <- rep(prob, length = ly)
  kappa <- rep(kappa, length = ly)
  delta <- rep(delta, length = ly)
  sigma <- rep(sigma, length = ly)
  xi <- rep(xi, length = ly)
  fy <- rep(0, ly)
  fy <- ifelse((x==0), pi+(1-pi)*ddiscegpd(0, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type), (1-pi)*ddiscegpd(x, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type))
  fy
}

dzidiscegpd(x, type = 3, pi=0.1, kappa = 1, sigma = 1, xi=0.3)
##################################
# CDF ---------------------------------------------------------------------------------------------
pzidiscegpd<-function(q, pi = NA,prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1)
{
  if (type==1) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==2) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==3) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==4) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
  }
  if (type==5) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==6) {
    if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (any(pi <= 0) | any(pi >= 1) )  stop(paste("pi must be between 0 and 1", "\n", ""))
  #if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
  #if (any(kappa <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  #if (any(delta <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", ""))
  if (any(xi <= 0) )  stop(paste("xi must be greater than 0", "\n", ""))
  
  if (any(q < -1) )  stop(paste("y must be 0 or greater than 0", "\n", ""))
  ly <- max(length(q),length(pi),length(prob),length(kappa),length(delta),length(sigma),length(xi))
  q <- rep(q, length = ly)
  pi <- rep(pi, length = ly)
  prob <- rep(prob, length = ly)
  kappa <- rep(kappa, length = ly)
  delta <- rep(delta, length = ly)
  sigma <- rep(sigma, length = ly)
  xi <- rep(xi, length = ly)
  cdf <- rep(0,ly)
  cdf <- pdiscegpd(q, prob = prob, kappa = kappa, delta =delta, sigma = sigma, xi = xi, type = type)
  cdf <- ifelse(q < 0, 0, pi + (1 - pi) * cdf)
  cdf
}

################################################################################################
#Quantile function-----------------------------------------------------------------------------------------

qzidiscegpd <- function(p,  pi = NA,prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1)
{
  
  if (type==1) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==2) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==3) {
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==4) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
  }
  if (type==5) {
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (type==6) {
    if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
    if (any(delta <= 0) )  stop(paste("delta must be greater than 0", "\n", ""))
    if (any(kappa <= 0) )  stop(paste("kappa must be greater than 0", "\n", ""))
  }
  if (any(pi <= 0) | any(pi >= 1) )  stop(paste("pi must be between 0 and 1", "\n", ""))
  #if (any(prob <= 0) | any(prob >= 1) )  stop(paste("prob must be between 0 and 1", "\n", ""))
  #if (any(kappa <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  #if (any(delta <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(sigma <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(xi <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))
  if (any(p <= 0) | any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))
  # if (log.p == TRUE) p <- exp(p)   else p <- p
  #if (lower.tail == TRUE)  p <- p  else p <- 1 - p
  ly <- max(length(p),length(pi),length(prob),length(kappa),length(delta),length(sigma),length(xi))
  p <- rep(p, length = ly)
  pi <- rep(pi, length = ly)
  prob <- rep(prob, length = ly)
  kappa <- rep(kappa, length = ly)
  delta <- rep(delta, length = ly)
  sigma <- rep(sigma, length = ly)
  xi <- rep(xi, length = ly)
  pnew <- ((p-pi)/(1-pi))-(1e-7)
  pnew <- ifelse(pnew > 0, pnew,0 )
  q<- qdiscegpd(pnew, prob = prob, kappa = kappa, delta =delta, sigma = sigma, xi = xi, type = type)
  q
}

###Zero Inflated discrete degpd ###
rzidiscegpd <- function(n,pi=NA, prob = NA, kappa = NA, delta = NA, sigma = NA, xi = NA, type = 1, unifsamp = NULL ) {
  z <- rbinom(n,size=1,prob=pi)
  y <- (1-z)*rdiscegpd(n, prob = prob, kappa = kappa, delta = delta, sigma = sigma, xi = xi, type = type, unifsamp = NULL)
  return(y)
}


# x<-rzidiscegpd(1000, pi=0.2, kappa = 15, sigma = 1, xi=0.5, type = 2)
# plot(table(x))
###################################################################################################
#Generator-----------------------------------------------------------------------------------------
##################################
discegpd.nll <- function(theta, x, type){
  if (type==1){
    if (theta[1] <= 0 | theta[2] <= 0 | theta[3] <= 10^(-6) | theta[3] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(ddiscegpd(x, kappa=theta[1], sigma=theta[2], xi=theta[3], type = 1)))
      return(v)
    }
  }else if (type==2){
    if (theta[1] <= 0 | theta[2] <= 0 | theta[3] <= 10^(-6) | theta[3] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(ddiscegpd(x, kappa=theta[1], sigma=theta[2], xi=theta[3], type = 2)))
      return(v)
    }
  }else if (type==3){
    if (theta[1] <= 0 | theta[2] <= 0 | theta[3] <= 10^(-6) | theta[3] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(ddiscegpd(x, kappa=theta[1], sigma=theta[2], xi=theta[3], type = 3)))
      return(v)
    }
  }else if (type == 4) {
    if (theta[1] <= 0 | theta[1] > 100 | theta[2] <= 0 | theta[3] <= 10^(-6) | theta[3] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(ddiscegpd(x, delta=theta[1], sigma = theta[2], xi=theta[3], type = 4)))
      return(v)
    }
  }else if (type == 5) {
    if (theta[1] <= 0 | theta[2] <= 0 | theta[2] > 100 | theta[3] <= 0 | theta[4] <= 10^(-6) | theta[4] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(ddiscegpd(x, kappa=theta[1],delta=theta[2], sigma = theta[3], xi=theta[4], type = 5)))
      return(v)
    }
  } else if (type == 6) {
    if (theta[1] < 0 | theta[1] > 1 | theta[2] <= 0 | theta[3] <= 0 | theta[4] <= 0 | theta[5] <= 10^(-6) | theta[5] > 0.99 |
        theta[3] < theta[2]) {
      return(Inf)
    } else {
      v<- -sum(log(ddiscegpd(x, prob=theta[1],kappa=theta[2],delta=theta[3], sigma = theta[4], xi=theta[5], type = 6)))
      return(v)
    }
  }
}
# x<-rdiscegpd(100, delta = 1, sigma = 1, xi=0.2, type = 4)
# discegpd.nll(c(0.1,1,0.3), x, type = 3)

discegpd.fit.ml <- function(x, type = 1, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, print = TRUE) {
  if (type == 1) {
    theta0 <- c(kappa0, sigma0, xi0)
    opt <- optim(par = theta0, fn = discegpd.nll, x = x, type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("kappa", "sigma", "xi")
    return(thetahat)
  } else if (type == 2) {
    theta0 <- c(kappa0, sigma0, xi0)
    opt <- optim(par = theta0, fn = discegpd.nll, x = x, type = type, method = "Nelder-Mead", control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("kappa", "sigma", "xi")
    return(thetahat)
  }else if (type == 3) {
    theta0 <- c(kappa0, sigma0, xi0)
    opt <- optim(par = theta0, fn = discegpd.nll, x = x, type = type, method = "Nelder-Mead", control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("kappa", "sigma", "xi")
    return(thetahat)
  }else if (type == 4) {
    theta0 <- c(delta0, sigma0, xi0)
    opt <- optim(par = theta0, fn = discegpd.nll, x = x, type = type, method = "Nelder-Mead", control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("delta", "sigma", "xi")
    return(thetahat)
  } else if (type == 5) {
    theta0 <- c(kappa0, delta0, sigma0, xi0)
    opt <- optim(par = theta0, fn = discegpd.nll, x = x,  type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c( "kappa" ,"delta", "sigma", "xi")
    return(thetahat)
  } else if (type == 6) {
    theta0 <- c(prob0, kappa0, delta0, sigma0, xi0)
    opt <- optim(par = theta0, fn = discegpd.nll, x = x, type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("prob", "kappa", "delta", "sigma", "xi")
    return(thetahat)
  }
}


##NLL
zi.discegpd.nll <- function(theta, x, type){
  if (type==1){
    if (theta[1] <=0 | theta[1] >= 1 | theta[2] <= 0 |theta[2] >20 | theta[3] <= 0 | theta[4] <= 10^(-6) | theta[4] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(dzidiscegpd(x,pi=theta[1], kappa=theta[2], sigma=theta[3], xi=theta[4], type = 1)))
      return(v)
    }
  }else if (type == 2){
    if (theta[1] <=0 | theta[1] >= 1 | theta[2] <= 0 |theta[2] >20 | theta[3] <= 0 | theta[4] <= 10^(-6) | theta[4] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(dzidiscegpd(x,pi=theta[1], kappa=theta[2], sigma=theta[3], xi=theta[4], type = 2)))
      return(v)
    }
  } else if (type == 3){
    if (theta[1] <=0 | theta[1] >= 1 | theta[2] <= 0 |theta[2] >20 | theta[3] <= 0 | theta[4] <= 10^(-6) | theta[4] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(dzidiscegpd(x,pi=theta[1], kappa=theta[2], sigma=theta[3], xi=theta[4], type = 3)))
      return(v)
    }
  } else if (type == 4) {
    if (theta[1] <=0 | theta[1] >= 1 | theta[2] <= 0 | theta[2] > 100 | theta[3] <= 0 | theta[4] <= 10^(-6) | theta[4] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(dzidiscegpd(x,pi=theta[1], delta=theta[2], sigma=theta[3], xi=theta[4], type = 4)))
      return(v)
    }
  }else if (type == 5) {
    if (theta[1] <=0 | theta[1] >= 1 | theta[2] <= 0 | theta[2] > 50 | theta[3] <= 0 | theta[3] > 60 | theta[4] <= 0 | theta[5] <= 10^(-6) | theta[5] > 0.99) {
      return(Inf)
    } else {
      v<- -sum(log(dzidiscegpd(x, pi=theta[1], kappa=theta[2],delta=theta[3], sigma = theta[4], xi=theta[5], type = 5)))
      return(v)
    }
  } else if (type == 6) {
    if (theta[1] <=0 | theta[1] >= 1 | theta[2] < 0 | theta[2] > 1 | theta[3] <= 0 | theta[4] <= 0 | theta[5] <= 0 | theta[6] <= 10^(-6) | theta[6] > 0.99 |
        theta[4] < theta[3]) {
      return(Inf)
    } else {
      v<- -sum(log(dzidiscegpd(x, pi=theta[1],prob=theta[2],kappa=theta[3],delta=theta[4], sigma = theta[5], xi=theta[6], type = 6)))
      return(v)
    }
  }
}
###
zidiscegpd.fit.ml <- function(x, type = 1,pi0 = NA, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, print = TRUE) {
  if (type == 1) {
    theta0 <- c(pi0,kappa0, sigma0, xi0)
    opt <- optim(par = theta0, fn = zi.discegpd.nll, x = x, type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("pi","kappa", "sigma", "xi")
    return(thetahat)
  } else if (type == 2) {
    theta0 <- c(pi0,kappa0, sigma0, xi0)
    opt <- optim(par = theta0, fn = zi.discegpd.nll, x = x, type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("pi","kappa", "sigma", "xi")
    return(thetahat)
  } else if (type == 3) {
    theta0 <- c(pi0,kappa0, sigma0, xi0)
    opt <- optim(par = theta0, fn = zi.discegpd.nll, x = x, type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("pi","kappa", "sigma", "xi")
    return(thetahat)
  } else if (type == 4) {
    theta0 <- c(pi0, delta0, sigma0, xi0)
    opt <- optim(par = theta0, fn = zi.discegpd.nll, x = x, type = type, method = "Nelder-Mead", control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("pi","delta", "sigma", "xi")
    return(thetahat)
  } else if (type == 5) {
    theta0 <- c(pi0, kappa0, delta0, sigma0, xi0)
    opt <- optim(par = theta0, fn = zi.discegpd.nll, x = x,  type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("pi", "kappa", "delta", "sigma", "xi")
    return(thetahat)
  } else if (type == 6) {
    theta0 <- c(pi0, prob0, kappa0, delta0, sigma0, xi0)
    opt <- optim(par = theta0, fn = zi.discegpd.nll, x = x, type = type, method = "Nelder-Mead",control = list(maxit = 1000), hessian = FALSE)
    if (print) {
      print(opt)
    }
    thetahat <- opt$par
    names(thetahat) <- c("pi","prob", "kappa", "delta", "sigma", "xi")
    return(thetahat)
  }
}


#######################BOOSTRAP####
discegpd.fit.ml.boot <- function(data, i, type = 1, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, print = FALSE) {
  return(discegpd.fit.ml(data[i], type = type, prob0 = prob0, kappa0 = kappa0, delta0 = delta0, sigma0 = sigma0, xi0 = xi0, print = print))
}

##################

###############Bootstrap estimates for confirdence Intervals #####################
zidiscegpd.fit.ml.boot <- function(data, i, type = 1,pi0=NA, prob0 = NA, kappa0 = NA, delta0 = NA, sigma0 = NA, xi0 = NA, print = FALSE) {
  return(zidiscegpd.fit.ml(data[i], type = type,pi0=pi0, prob0 = prob0, kappa0 = kappa0, delta0 = delta0, sigma0 = sigma0, xi0 = xi0, print = print))
}


##########################


############

fit.model <- function(data, model = 1, family= c("discegpd", "zidiscegpd"), init, confint = FALSE, R = 1000, ncpus = 1,plots = TRUE) {
  if(family=="discegpd"){
    # Sanity checks
    initsize <- switch(model, 3, 3,3,3, 4, 5)
    if (length(init) != initsize) {
      stop("Invalid starting values in `init'; incorrect length.")
    }
    data = data
    x <- seq(0, max(data, na.rm = TRUE), by = 1)
    ddiscegpds <- c()
    qdiscegpds <- c()
    if (model == 1) {
      fit.mle <- discegpd.fit.ml(x = data, type = 1,kappa0 = init[1], sigma0 = init[2], xi0 = init[3])
      if (confint) {
        fit.mle.boot <- boot::boot(data = data, statistic = discegpd.fit.ml.boot, R = R, type = 1, kappa0 = init[1], sigma0 = init[2],
                                   xi0 = init[3], parallel = "multicore", ncpus = ncpus)
        CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
        CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
        CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
        CIs.mle <- cbind(CI.mle.kappa, CI.mle.sigma, CI.mle.xi)
      }
      ddiscegpd.mle <- ddiscegpd(x = x, type = 1, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
      ddiscegpds <- c(ddiscegpds, ddiscegpd.mle)
      if (plots) {
        qdiscegpd.mle <- qdiscegpd(p = c(1:length(data))/(length(data) + 1), type = 1, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
        qdiscegpds <- c(qdiscegpds, qdiscegpd.mle)
        if (confint) {
          q.mle.boot <- mapply(FUN = qdiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(1),kappa= as.list(fit.mle.boot$t[,1]), sigma = as.list(fit.mle.boot$t[, 2]), xi = as.list(fit.mle.boot$t[, 3]))
          q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
          q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
        }
      }
    } else if (model == 2) {
      fit.mle <- discegpd.fit.ml(x = data, type = 2,kappa0 = init[1], sigma0 = init[2], xi0 = init[3])
      if (confint) {
        fit.mle.boot <- boot::boot(data = data, statistic = discegpd.fit.ml.boot, R = R, type = 2, kappa0 = init[1], sigma0 = init[2],
                                   xi0 = init[3], parallel = "multicore", ncpus = ncpus)
        CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
        CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
        CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
        CIs.mle <- cbind(CI.mle.kappa, CI.mle.sigma, CI.mle.xi)
      }
      ddiscegpd.mle <- ddiscegpd(x = x, type = 2, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
      ddiscegpds <- c(ddiscegpds, ddiscegpd.mle)
      if (plots) {
        qdiscegpd.mle <- qdiscegpd(p = c(1:length(data))/(length(data) + 1), type = 2, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
        qdiscegpds <- c(qdiscegpds, qdiscegpd.mle)
        if (confint) {
          q.mle.boot <- mapply(FUN = qdiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(2),kappa= as.list(fit.mle.boot$t[,1]), sigma = as.list(fit.mle.boot$t[, 2]), xi = as.list(fit.mle.boot$t[, 3]))
          q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
          q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
        }
      }
    } else if (model == 3) {
      fit.mle <- discegpd.fit.ml(x = data, type = 3,kappa0 = init[1], sigma0 = init[2], xi0 = init[3])
      if (confint) {
        fit.mle.boot <- boot::boot(data = data, statistic = discegpd.fit.ml.boot, R = R, type = 3, kappa0 = init[1], sigma0 = init[2],
                                   xi0 = init[3], parallel = "multicore", ncpus = ncpus)
        CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
        CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
        CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
        CIs.mle <- cbind(CI.mle.kappa, CI.mle.sigma, CI.mle.xi)
      }
      ddiscegpd.mle <- ddiscegpd(x = x, type = 3, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
      ddiscegpds <- c(ddiscegpds, ddiscegpd.mle)
      if (plots) {
        qdiscegpd.mle <- qdiscegpd(p = c(1:length(data))/(length(data) + 1), type = 3, kappa = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
        qdiscegpds <- c(qdiscegpds, qdiscegpd.mle)
        if (confint) {
          q.mle.boot <- mapply(FUN = qdiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(3),kappa= as.list(fit.mle.boot$t[,1]), sigma = as.list(fit.mle.boot$t[, 2]), xi = as.list(fit.mle.boot$t[, 3]))
          q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
          q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
        }
      }
    } else if (model == 4) {
      fit.mle <- discegpd.fit.ml(x = data, type = 4,  delta0 = init[1], sigma0 = init[2], xi0 = init[3])
      if (confint) {
        fit.mle.boot <- boot::boot(data = data, statistic = discegpd.fit.ml.boot, R = R, type = 4, delta0 = init[1], sigma0 = init[2],
                                   xi0 = init[3], parallel = "multicore", ncpus = ncpus)
        CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
        CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
        CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
        CIs.mle <- cbind(CI.mle.delta, CI.mle.sigma, CI.mle.xi)
      }
      ddiscegpd.mle <- ddiscegpd(x = x, type = 4, delta = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
      ddiscegpds <- c(ddiscegpds, ddiscegpd.mle)
      if (plots) {
        qdiscegpd.mle <- qdiscegpd(p = c(1:length(data))/(length(data) + 1), type = 4, delta = fit.mle[1], sigma = fit.mle[2], xi = fit.mle[3])
        qdiscegpds <- c(qdiscegpds, qdiscegpd.mle)
        if (confint) {
          q.mle.boot <- mapply(FUN = qdiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(4), delta = as.list(fit.mle.boot$t[, 1]), sigma = as.list(fit.mle.boot$t[, 2]), xi = as.list(fit.mle.boot$t[, 3]))
          q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
          q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
        }
      }
    } else if (model == 5) {
      fit.mle <- discegpd.fit.ml(x = data, type = 5, kappa0 = init[1], delta0 = init[2], sigma0 = init[3], xi0 = init[4])
      if (confint) {
        fit.mle.boot <- boot::boot(data = data, statistic = discegpd.fit.ml.boot, R = R, type = 5, kappa0 = init[1], delta0 = init[2],
                                   sigma0 = init[3], xi0 = init[4], parallel = "multicore", ncpus = ncpus)
        CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
        CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
        CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
        CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
        CIs.mle <- cbind(CI.mle.kappa, CI.mle.delta, CI.mle.sigma, CI.mle.xi)
      }
      ddiscegpd.mle <- ddiscegpd(x = x, type = 5, kappa = fit.mle[1], delta = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
      ddiscegpds <- c(ddiscegpds, ddiscegpd.mle)
      if (plots) {
        qdiscegpd.mle <- qdiscegpd(p = c(1:length(data))/(length(data) + 1), type = 5,kappa = fit.mle[1], delta = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
        qdiscegpds <- c(qdiscegpds, qdiscegpd.mle)
        if (confint) {
          q.mle.boot <- mapply(FUN = qdiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(5),  kappa = as.list(fit.mle.boot$t[, 1]),delta = as.list(fit.mle.boot$t[, 2]), sigma = as.list(fit.mle.boot$t[, 3]), xi = as.list(fit.mle.boot$t[,4]))
          q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
          q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
        }
      }
    } else if (model == 6) {
      fit.mle <- discegpd.fit.ml(x = data, type = 6,prob0 = init[1], kappa0 = init[2], delta0 = init[3], sigma0 = init[4], xi0 = init[5])
      if (confint) {
        fit.mle.boot <- boot::boot(data = data, statistic = discegpd.fit.ml.boot, R = R, type = 6, prob0 = init[1], kappa0 = init[2], delta0 = init[3], sigma0 = init[4], xi0 = init[5], parallel = "multicore",
                                   ncpus = ncpus)
        CI.mle.prob <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
        CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index =2, type = "perc")$perc[4:5]
        CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
        CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
        CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 5, type = "perc")$perc[4:5]
        CIs.mle <- cbind(CI.mle.prob, CI.mle.kappa, CI.mle.delta, CI.mle.sigma, CI.mle.xi)
      }
      ddiscegpd.mle <-  ddiscegpd(x = x, type = 6,prob = fit.mle[1], kappa = fit.mle[2], delta = fit.mle[3], sigma = fit.mle[4], xi = fit.mle[5])
      ddiscegpds <- c(ddiscegpds, ddiscegpd.mle)
      if (plots) {
        qdiscegpd.mle <- qdiscegpd(p = c(1:length(data))/(length(data) + 1), type = 6, prob = fit.mle[1], kappa = fit.mle[2], delta = fit.mle[3], sigma = fit.mle[4], xi = fit.mle[5])
        qdiscegpds <- c(qdiscegpds, qdiscegpd.mle)
        if (confint) {
          q.mle.boot <- mapply(FUN =qdiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(6),  prob = as.list(fit.mle.boot$t[, 1]),kappa = as.list(fit.mle.boot$t[, 2]), delta = as.list(fit.mle.boot$t[, 3]), sigma = as.list(fit.mle.boot$t[,4]), xi = as.list(fit.mle.boot$t[, 5]))
          q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
          q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
        }
      }
    }
    fits<- list(mle=fit.mle)
    if (confint){
      CIs<- list(mle=CIs.mle)
    }else{
      CIs<- NULL
    }
    if (plots) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      mat <- matrix(c(1:(1 )), nrow = 1, ncol = 1 , byrow = TRUE)
      layout(mat)
      par(mar = c(4, 4, 1, 1))
      
      # plot(table(data)/length(data), ylab="Density",xlab="Data", main="", col = "red", lwd= 1.5)
      # lines(x+0.5, ddiscegpd.mle, col = "blue",type="h",lwd= 1.5 )
      # legend("topright", legend=c("Empirical", "Fitted"), col=c("red", "blue"), lty=1, cex=1,lwd= 1.5)
      # grid()
      plot(data, qdiscegpd.mle, asp = 1, xlab = "Empirical quantiles", ylab = "Fitted quantiles", ylim = range(qdiscegpds, na.rm = TRUE),
           type = "n")
      if (confint) {
        polygon(x = c(sort(data), sort(data, decreasing = TRUE)), y = c(q.mle.L, q.mle.U[length(q.mle.U):1]), col = rgb(0,  0, 1, alpha = 0.1), lty = 1, border = rgb(0, 0, 1, alpha = 0.5))
      }
      #grid()
      lines(sort(data), qdiscegpd.mle, lty = 1, type = "b", pch = 20, col = "gray")
      abline(0, 1, col = "red")
    }
    return(list(fit=fits, confint=CIs))
  }else if (family=="zidiscegpd"){
      # Sanity checks
      initsize <- switch(model, 4, 4,4 ,4,5, 6)
      if (length(init) != initsize) {
        stop("Invalid starting values in `init'; incorrect length.")
      }
      x <- seq(0, max(data, na.rm = TRUE), by = 1)
      dzidiscegpds <- c()
      qzidiscegpds <- c()
      if (model == 1) {
        fit.mle <- zidiscegpd.fit.ml(x = data, type = 1, pi0 = init[1],kappa0 = init[2], sigma0 = init[3], xi0 = init[4])
        if (confint) {
          fit.mle.boot <- boot::boot(data = data, statistic = zidiscegpd.fit.ml.boot, R = R, type = 1,pi0 = init[1], kappa0 = init[2], sigma0 = init[3],
                                     xi0 = init[4], parallel = "multicore", ncpus = ncpus)
          CI.mle.pi <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
          CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
          CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
          CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
          CIs.mle <- cbind(CI.mle.pi, CI.mle.kappa, CI.mle.sigma, CI.mle.xi)
        }
        zidiscegpd.mle <- dzidiscegpd(x = x, type = 1,pi=fit.mle[1], kappa = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
        dzidiscegpds <- c(dzidiscegpds, zidiscegpd.mle)
        if (plots) {
          qzidiscegpd.mle <- qzidiscegpd(p = c(1:length(data))/(length(data) + 1), type = 1,pi=fit.mle[1], kappa = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
          qzidiscegpds <- c(qzidiscegpds, qzidiscegpd.mle)
          if (confint) {
            q.mle.boot <- mapply(FUN = qzidiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(1), pi = as.list(fit.mle.boot$t[,1]),kappa= as.list(fit.mle.boot$t[,2]), sigma = as.list(fit.mle.boot$t[, 3]), xi = as.list(fit.mle.boot$t[, 4]))
            q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
            q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
          }
        }
      } else if (model == 2) {
        fit.mle <- zidiscegpd.fit.ml(x = data, type = 2, pi0 = init[1],kappa0 = init[2], sigma0 = init[3], xi0 = init[4])
        if (confint) {
          fit.mle.boot <- boot::boot(data = data, statistic = zidiscegpd.fit.ml.boot, R = R, type = 2,pi0 = init[1], kappa0 = init[2], sigma0 = init[3],
                                     xi0 = init[4], parallel = "multicore", ncpus = ncpus)
          CI.mle.pi <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
          CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
          CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
          CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
          CIs.mle <- cbind(CI.mle.pi, CI.mle.kappa, CI.mle.sigma, CI.mle.xi)
        }
        zidiscegpd.mle <- dzidiscegpd(x = x, type = 2,pi=fit.mle[1], kappa = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
        dzidiscegpds <- c(dzidiscegpds, zidiscegpd.mle)
        if (plots) {
          qzidiscegpd.mle <- qzidiscegpd(p = c(1:length(data))/(length(data) + 1), type = 2,pi=fit.mle[1], kappa = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
          qzidiscegpds <- c(qzidiscegpds, qzidiscegpd.mle)
          if (confint) {
            q.mle.boot <- mapply(FUN = qzidiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(2), pi = as.list(fit.mle.boot$t[,1]),kappa= as.list(fit.mle.boot$t[,2]), sigma = as.list(fit.mle.boot$t[, 3]), xi = as.list(fit.mle.boot$t[, 4]))
            q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
            q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
          }
        }
      } else if (model == 3) {
        fit.mle <- zidiscegpd.fit.ml(x = data, type = 3, pi0 = init[1],kappa0 = init[2], sigma0 = init[3], xi0 = init[4])
        if (confint) {
          fit.mle.boot <- boot::boot(data = data, statistic = zidiscegpd.fit.ml.boot, R = R, type = 3,pi0 = init[1], kappa0 = init[2], sigma0 = init[3],
                                     xi0 = init[4], parallel = "multicore", ncpus = ncpus)
          CI.mle.pi <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
          CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
          CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
          CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
          CIs.mle <- cbind(CI.mle.pi, CI.mle.kappa, CI.mle.sigma, CI.mle.xi)
        }
        zidiscegpd.mle <- dzidiscegpd(x = x, type = 3,pi=fit.mle[1], kappa = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
        dzidiscegpds <- c(dzidiscegpds, zidiscegpd.mle)
        if (plots) {
          qzidiscegpd.mle <- qzidiscegpd(p = c(1:length(data))/(length(data) + 1), type = 3,pi=fit.mle[1], kappa = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
          qzidiscegpds <- c(qzidiscegpds, qzidiscegpd.mle)
          if (confint) {
            q.mle.boot <- mapply(FUN = qzidiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(3), pi = as.list(fit.mle.boot$t[,1]),kappa= as.list(fit.mle.boot$t[,2]), sigma = as.list(fit.mle.boot$t[, 3]), xi = as.list(fit.mle.boot$t[, 4]))
            q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
            q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
          }
        }
      } else if (model == 4) {
        fit.mle <- zidiscegpd.fit.ml(x = data, type = 4,  pi0 = init[1],delta0 = init[2], sigma0 = init[3], xi0 = init[4])
        if (confint) {
          fit.mle.boot <- boot::boot(data = data, statistic = zidiscegpd.fit.ml.boot, R = R, type = 4,pi0 = init[1], delta0 = init[2], sigma0 = init[3],
                                     xi0 = init[4], parallel = "multicore", ncpus = ncpus)
          CI.mle.pi <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
          CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
          CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
          CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
          CIs.mle <- cbind(CI.mle.pi, CI.mle.delta, CI.mle.sigma, CI.mle.xi)
        }
        zidiscegpd.mle <- dzidiscegpd(x = x, type = 4, pi = fit.mle[1], delta = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
        dzidiscegpds <- c(dzidiscegpds, zidiscegpd.mle)
        if (plots) {
          qzidiscegpd.mle <- qzidiscegpd(p = c(1:length(data))/(length(data) + 1), type = 4, pi = fit.mle[1], delta = fit.mle[2], sigma = fit.mle[3], xi = fit.mle[4])
          qzidiscegpds <- c(qzidiscegpds, qzidiscegpd.mle)
          if (confint) {
            q.mle.boot <- mapply(FUN = qzidiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(4), pi = as.list(fit.mle.boot$t[,
                                                                                                                                             1]),delta = as.list(fit.mle.boot$t[, 2]), sigma = as.list(fit.mle.boot$t[, 3]), xi = as.list(fit.mle.boot$t[, 4]))
            q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
            q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
          }
        }
      } else if (model == 5) {
        fit.mle <- zidiscegpd.fit.ml(x = data, type = 5,pi0 = init[1], kappa0 = init[2], delta0 = init[3], sigma0 = init[4], xi0 = init[5])
        if (confint) {
          fit.mle.boot <- boot::boot(data = data, statistic = zidiscegpd.fit.ml.boot, R = R, type = 5, pi0 = init[1],kappa0 = init[2], delta0 = init[3],
                                     sigma0 = init[4], xi0 = init[5], parallel = "multicore", ncpus = ncpus)
          CI.mle.pi <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
          CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
          CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
          CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
          CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 5, type = "perc")$perc[4:5]
          CIs.mle <- cbind(CI.mle.pi,CI.mle.kappa, CI.mle.delta, CI.mle.sigma, CI.mle.xi)
        }
        zidiscegpd.mle <- dzidiscegpd(x = x, type = 5, pi = fit.mle[1],kappa = fit.mle[2], delta = fit.mle[3], sigma = fit.mle[4], xi = fit.mle[5])
        dzidiscegpds <- c(dzidiscegpds, zidiscegpd.mle)
        if (plots) {
          qzidiscegpd.mle <- qzidiscegpd(p = c(1:length(data))/(length(data) + 1), type = 5,pi = fit.mle[1],kappa = fit.mle[2], delta = fit.mle[3], sigma = fit.mle[4], xi = fit.mle[5])
          qzidiscegpds <- c(qzidiscegpds, qzidiscegpd.mle)
          if (confint) {
            q.mle.boot <- mapply(FUN = qzidiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(5), pi = as.list(fit.mle.boot$t[,
                                                                                                                                             1]), kappa = as.list(fit.mle.boot$t[, 2]),delta = as.list(fit.mle.boot$t[, 3]), sigma = as.list(fit.mle.boot$t[, 4]), xi = as.list(fit.mle.boot$t[,
                                                                                                                                                                                                                                                                                               5]))
            q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
            q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
          }
        }
      } else if (model == 6) {
        fit.mle <- zidiscegpd.fit.ml(x = data, type = 6, pi0 = init[1],prob0 = init[2], kappa0 = init[3], delta0 = init[4], sigma0 = init[5], xi0 = init[6])
        if (confint) {
          fit.mle.boot <- boot::boot(data = data, statistic = zidiscegpd.fit.ml.boot, R = R, type = 6, pi0 = init[1],prob0 = init[2], kappa0 = init[3], delta0 = init[4], sigma0 = init[5], xi0 = init[6], parallel = "multicore",
                                     ncpus = ncpus)
          CI.mle.pi <- boot::boot.ci(boot.out = fit.mle.boot, index = 1, type = "perc")$perc[4:5]
          CI.mle.prob <- boot::boot.ci(boot.out = fit.mle.boot, index = 2, type = "perc")$perc[4:5]
          CI.mle.kappa <- boot::boot.ci(boot.out = fit.mle.boot, index = 3, type = "perc")$perc[4:5]
          CI.mle.delta <- boot::boot.ci(boot.out = fit.mle.boot, index = 4, type = "perc")$perc[4:5]
          CI.mle.sigma <- boot::boot.ci(boot.out = fit.mle.boot, index = 5, type = "perc")$perc[4:5]
          CI.mle.xi <- boot::boot.ci(boot.out = fit.mle.boot, index = 6, type = "perc")$perc[4:5]
          CIs.mle <- cbind(CI.mle.pi,CI.mle.prob, CI.mle.kappa, CI.mle.delta, CI.mle.sigma, CI.mle.xi)
        }
        zidiscegpd.mle <-  dzidiscegpd(x = x, type = 6, pi = fit.mle[1],prob = fit.mle[2], kappa = fit.mle[3], delta = fit.mle[4], sigma = fit.mle[5], xi = fit.mle[6])
        dzidiscegpds <- c(dzidiscegpds, zidiscegpd.mle)
        if (plots) {
          qzidiscegpd.mle <- qzidiscegpd(p = c(1:length(data))/(length(data) + 1), type = 6, pi = fit.mle[1],prob = fit.mle[2], kappa = fit.mle[3], delta = fit.mle[4], sigma = fit.mle[5], xi = fit.mle[6])
          qzidiscegpds <- c(qzidiscegpds, qzidiscegpd.mle)
          if (confint) {
            q.mle.boot <- mapply(FUN =qzidiscegpd, p = list(c(1:length(data))/(length(data) + 1)), type = list(6), pi = as.list(fit.mle.boot$t[,
                                                                                                                                            1]), prob = as.list(fit.mle.boot$t[, 2]),kappa = as.list(fit.mle.boot$t[, 3]), delta = as.list(fit.mle.boot$t[, 4]), sigma = as.list(fit.mle.boot$t[,
                                                                                                                                                                                                                                                                                                5]), xi = as.list(fit.mle.boot$t[, 6]))
            q.mle.L <- apply(q.mle.boot, 1, quantile, 0.025, na.rm = TRUE)
            q.mle.U <- apply(q.mle.boot, 1, quantile, 0.975, na.rm = TRUE)
          }
        }
      }
      fits<- list(mle=fit.mle)
      if (confint){
        CIs<- list(mle=CIs.mle)
      }else{
        CIs<- NULL
      }
      if (plots) {
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        mat <- matrix(c(1:(1 )), nrow = 1, ncol = 1 , byrow = TRUE)
        layout(mat)
        par(mar = c(4, 4, 1, 1))
        
        # plot(table(data)/length(data), ylab="Density",xlab="Data", main="", col = "red", lwd= 1.5)
        # lines(x+0.5, zidiscegpd.mle, col = "blue", type="h", lwd= 1.5)
        # legend("topright", legend=c("Empirical", "Fitted"), col=c("red", "blue"), lty=1, cex=1,lwd= 1.5)
        
        #grid()
        plot(data, qzidiscegpd.mle, asp = 1, xlab = "Empirical quantiles", ylab = "Fitted quantiles", ylim = range(qzidiscegpds, na.rm = TRUE),
             type = "n")
        if (confint) {
          polygon(x = c(sort(data), sort(data, decreasing = TRUE)), y = c(q.mle.L, q.mle.U[length(q.mle.U):1]), col = rgb(0,  0, 1, alpha = 0.1), lty = 1, border = rgb(0, 0, 1, alpha = 0.5))
        }
        #grid()
        lines(sort(data), qzidiscegpd.mle, lty = 1, type = "b", pch = 20, col = "grey")
        abline(0, 1, col = "red")
        #abline(1, 1, col = "lightgrey")
      }
      return(list(fit=fits, confint=CIs))
      
    }
  
}

##Return level plots
# returnlevelplot <- function(y, pd) {
#   library(ggplot2)
#   
#   n <- length(y)
#   q <- sort(y)
#   
#   p <- (1:n) / (n + 1)
#   
#   t <- 1 / (1 - p)
#   
#   data <- data.frame(Return_Period = t, Return_Level = q)
#   
#   ggplot(data, aes(x = Return_Period, y = Return_Level)) +
#     geom_point(col = "grey", pch = 16) +
#     geom_line(aes(y = quantile(pd, p)), linetype = "dashed", color = "black") +
#     scale_x_log10(breaks = 10^(-2:4), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#     labs(x = "m-observation return period", y = "Return Level") +
#     theme_minimal()
# }


returnlevelplot <- function(y, pd) {
  n <- length(y)
  q <- sort(y)
  
  # Calculate probabilities
  p <- (1:n) / (n + 1)
  
  # Return periods
  t <- 1 / (1 - p)
  
  # Create the plot using base R plotting, with log scale on x-axis
  plot(t, q, log = "x", pch = 1, col = "grey",
       xlab = "m-observation return period", ylab = "Return Level", 
       main = "", xaxt = "n")
  
  # Add a dashed line for the quantiles of pd
  lines(t, quantile(pd, p), lty = 2, col = "red")
  
  # Set custom x-axis ticks with powers of 10 only
  axis(1, at = 10^(-2:4), labels = parse(text = paste0("10^", -2:4)))
}



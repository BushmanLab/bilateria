
fitpowbend2<-function (x, trunc, start.value, ...) {
	dots <- list(...)
	if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
	if (!missing(trunc)) {
		if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
	}
	if (missing(start.value)) {
		gamahat <- function(ga, xvec) eq <- -sum(log(xvec)) - VGAM::zeta(ga, deriv = 1) * length(xvec)/VGAM::zeta(ga)
		shat <- uniroot(gamahat, interval = c(1.01, 20), xvec = x)$root
		oMhat <- 3
	}
	else {
		shat <- start.value[1]
		oMhat <- start.value[2]
	}
	if (!"method" %in% names(dots)) {
		dots$method <- "L-BFGS-B"
		if (!"lower" %in% names(dots)) 
			dots$lower = c(s = 0.1, oM = 1)
		if (!"upper" %in% names(dots)) 
			dots$upper = c(s = 2.999, oM = 10000)
	}
	if (missing(trunc)) {
		LL <- function(s, oM) -sum(dpowbend(x, s = s, oM = oM, log = TRUE))
	}
	else {
		LL <- function(s, oM) -sum(dtrunc2("powbend", x = x, coef = list(s = s, oM = oM), trunc = trunc, log = TRUE))
	}
	result <- c(
    lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(do.call("mle2", c(list(LL, start = list(s = shat, oM = oMhat*mult), data = list(x = x)), dots)),silent=TRUE)),
	  lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(do.call("mle2", c(list(LL, start = list(s = shat*mult, oM = oMhat), data = list(x = x)), dots)),silent=TRUE)),
	  lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(do.call("mle2", c(list(LL, start = list(s = shat*mult, oM = oMhat*mult), data = list(x = x)), dots)),silent=TRUE)),
	  lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(do.call("mle2", c(list(LL, start = list(s = shat/mult, oM = oMhat*mult), data = list(x = x)), dots)),silent=TRUE)),
	  lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(do.call("mle2", c(list(LL, start = list(s = mult, oM = 1), data = list(x = x)), dots)),silent=TRUE)),
	  lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(do.call("mle2", c(list(LL, start = list(s = 1, oM = mult), data = list(x = x)), dots)),silent=TRUE)),
	  lapply(c(1000),function(mult)try(do.call("mle2", c(list(LL, start = list(s = coef(fitpower(xx,trunc=trunc)), oM = mult), data = list(x = x)), dots)),silent=TRUE))
  )
  result<-result[[which.min(sapply(result,function(xx){if(inherits(xx,'try-error'))return(Inf);return(xx@min)}))]]
  if(inherits(result,'try-error'))stop(result)
	new("fitsad", result, sad = "powbend", trunc = ifelse(missing(trunc), NaN, trunc))
}


fitnbinom2<- function (x, trunc = 0, start.value, ...) {
  dots <- list(...)
  if (any(x <= 0) | any(!sads:::is.wholenumber(x))) 
    stop("All x must be positive integers")
  if (!is.null(trunc)) {
    if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)) {
    muhat <- length(x)/(length(x) + mean(x))
    sizehat <- muhat * mean(x)
  }
  else {
    sizehat <- start.value[[1]]
    muhat <- start.value[[2]]
  }
  if (is.null(trunc)) {
    LL <- function(size, mu) -sum(dnbinom(x, size = size, mu = mu, log = TRUE))
  }
  else {
    LL <- function(size, mu) -sum(dtrunc2("nbinom", x = x, coef = list(size = size, mu = mu), trunc = trunc, log = TRUE))
  }
  result <- lapply(c(.01,.1,1,10,100),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(size = mult, mu = muhat), data = list(x = x)),method='Nelder-Mead',skip.hessian=TRUE, dots))),silent=TRUE))
  result<-result[[which.min(sapply(result,function(xx){if(inherits(xx,'try-error'))return(Inf);return(xx@min)}))]]
  if(inherits(result,'try-error'))stop(result)
  new("fitsad", result, sad = "nbinom", trunc = ifelse(is.null(trunc), NaN, trunc))
}

fitls2<-function (x, trunc, start.value, upper = length(x), ...) {
    dots <- list(...)
    if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
    S <- length(x)
    N <- sum(x)
    if (missing(start.value)) {
        f1 <- function(a) S + a * log((a/(a + N)))
        sol <- uniroot(f1, interval = c(1/N, N))
        alfa <- sol$root
        X <- N/(N + alfa)
    } else {
        alfa <- start.value
    }
    if (!missing(trunc)) {
        if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
        else LL <- function(N, alpha) -sum(dtrunc2("ls", x = x, coef = list(N = N, alpha = alpha), trunc = trunc, log = TRUE))
    }
    if (missing(trunc)) LL <- function(N, alpha) -sum(dls2(x, N, alpha, log = TRUE))
    result <- do.call("mle2", c(list(LL, start = list(alpha = alfa), 
        data = list(x = x), fixed = list(N = N), method = "Brent", 
        lower = 0, upper = upper), dots))
    if (abs(as.numeric(result@coef) - upper) < 1e-07) warning("mle equal to upper bound provided. \\n Try new value for the 'upper' argument")
    new("fitsad", result, sad = "ls", trunc = ifelse(missing(trunc), NaN, trunc))
}



dls2<-function (x, N, alpha, log = FALSE) {
    N[!is.finite(N) | N <= 0] <- NaN
    alpha[!is.finite(alpha) | alpha <= 0] <- NaN
    X <- N/(N + alpha)
    gama <- function(y) 1/log(1/(1 - y))
    y <- log(gama(X)) + x*log(X)-log(x)
    if (any(is.nan(y))) warning("NaNs produced")
    if (any(!sads:::is.wholenumber(x))) warning("non integer values in x")
    y[!sads:::is.wholenumber(x) | x < 1] <- 0
    if (log) return(y)
    else return(exp(y))
}


fitbs2<-function (x, trunc, ...) {
    dots <- list(...)
    if (any(x <= 0)) stop("All x must be positive")
    s <- length(x)
    n <- sum(x)
    if (!missing(x)) {
        if (!missing(trunc)) {
            if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
            else LL <- function(N, S) -sum(dtrunc2("bs", x = x, coef = list(N = N, S = S), trunc = trunc, log = TRUE))
        }
        if (missing(trunc)) LL <- function(N, S) -sum(dbs2(x = x, N = N, S = S, log = TRUE))
        result <- do.call("mle2", c(list(minuslogl = LL, data = list(x = x), fixed = list(N = n, S = s), eval.only = TRUE), dots))
        result@details$convergence = 0
        new("fitsad", result, sad = "bs", trunc = ifelse(missing(trunc), NaN, trunc))
    }
}


dbs2<-function(x,N,S,log=FALSE){
    N[!is.finite(N) | N <= 0] <- NaN
    S[!is.finite(S) | S <= 0] <- NaN
    y <- log(S - 1) + (S - 2)*log(1 - x/N)-log(N)
    if (any(is.nan(y))) warning("NaNs produced")
    y[x < 0 | x > N] <- 0
    if (log) return(y)
    else return(exp(y))
}


fitmzsm2 <-function (x, trunc, start.value, upper = length(x), ...) {
    dots <- list(...)
    if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
    if (sum(x) < 100) warning("\\n small sample size (J<100); \\n mzsm may not be a good approximation")
    if (!missing(trunc)) {
        if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
    }
    if (missing(start.value)) thetahat <- length(x)
    else thetahat <- start.value
    if (missing(trunc)) {
        LL <- function(J, theta) -sum(dmzsm2(x, J = J, theta = theta, log = TRUE))
    }
    else {
        LL <- function(J, theta) -sum(dtrunc2("mzsm", x = x, coef = list(J = J, theta = theta), trunc = trunc, log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = list(theta = thetahat), 
        fixed = list(J = sum(x)), data = list(x = x), method = "Brent", 
        lower = 0.001, upper = upper), dots))
    if (abs(as.numeric(result@coef) - upper) < 1e-07) warning("mle equal to upper bound provided. \\n Try value for the 'upper' arguent")
    new("fitsad", result, sad = "mzsm", trunc = ifelse(missing(trunc), NaN, trunc))
}

dmzsm2<-function (x, J, theta, log = FALSE) {
    J[!is.finite(J) | J <= 0] <- NaN
    theta[!is.finite(theta) | theta <= 0] <- NaN
    mzsm <- function(y, J, theta) log(theta)-log(y) + (theta - 1)*log(1 - y/J)
    sn <- mzsm(y = x, J = J, theta = theta)
    mu <- mzsm(y = 1:J, J = J, theta = theta)
    lpn <- suppressWarnings(sn - log(sum(exp(mu))))
    if (any(!sads:::is.wholenumber(x))) warning("non integer values in x")
    lpn[x <= 0 | x > J | !sads:::is.wholenumber(x)] <- -Inf
    if (any(is.nan(lpn))) warning("NaNs produced")
    if (log) return(lpn)
    else return(exp(lpn))
}


fitgamma2<- function (x, trunc, start.value, ...) {
  dots <- list(...)
  if (any(x <= 0)) stop("All x must be positive")
  if (!missing(trunc)) {
    if (min(x) <= trunc) stop("truncation point should, be lower than the lowest data value")
  }
  if (missing(start.value)) {
    if (missing(trunc)) {
      ka <- (mean(x)/sd(x))^2
      theta <- var(x)/mean(x)
      kahat <- function(k, dados) {
        eq <- length(dados) * (log(k) - log(mean(dados)) - digamma(k)) + sum(log(dados))
      }
      ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x)$root
      theta <- mean(x)/ka
    }
    else {
      xh <- hist(x, plot = FALSE)
      xbr <- xh$breaks
      Eh <- matrix(ncol = 2, nrow = length(xbr) - 1)
      for (i in 1:(length(xbr) - 1)) {
        m1 <- matrix(c(1, 1, -1, 1), ncol = 2)
        Eh[i, ] <- solve(m1, xbr[i:(i + 1)])
      }
      P <- xh$counts/sum(xh$counts)
      P[P == 0] <- min(P[P > 0])
      Y <- log(P[-length(P)]) - log(P[-1]) - (log(Eh[-nrow(Eh), 2]) - log(Eh[-1, 2]))
      X1 <- Eh[-1, 1] - Eh[-nrow(Eh), 1]
      X2 <- log(Eh[-nrow(Eh), 1]) - log(Eh[-1, 1])
      st1 <- unname(coef(lm(Y ~ X1 + X2 - 1)))
      ka <- st1[2] + 1
      theta <- 1/st1[1]
    }
  }
  else {
    ka <- start.value[1]
    theta <- start.value[2]
  }
  if (missing(trunc)) {
    LL <- function(shape, rate) -sum(dgamma(x, shape, rate, log = TRUE))
  }
  else {
    LL <- function(shape, rate) -sum(dtrunc2("gamma", x = x, coef = list(shape = shape, rate = rate), trunc = trunc, log = TRUE))
  }
  result <- c(
    lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(shape = ka, rate = 1/theta*mult), data = list(x = x)), dots))),silent=TRUE)),
    lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(shape = ka*mult, rate = 1/theta), data = list(x = x)), dots))),silent=TRUE)),
    lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(shape = ka*mult, rate = 1/theta*mult), data = list(x = x)), dots))),silent=TRUE)),
    lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(shape = ka*mult, rate = 1/theta/mult), data = list(x = x)), dots))),silent=TRUE)),
    lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(shape = 1, rate = mult), data = list(x = x)), dots))),silent=TRUE)),
    lapply(c(.001,.01,.1,1,10,100,1000),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(shape = mult, rate = 1), data = list(x = x)), dots))),silent=TRUE))
  )
  result<-result[[which.min(sapply(result,function(xx){if(inherits(xx,'try-error'))return(Inf);return(xx@min)}))]]
  if(inherits(result,'try-error'))stop(result)
  new("fitsad", result, sad = "gamma", trunc = ifelse(missing(trunc), NaN, trunc))
}


fitgeom2<-function (x, trunc = 0, start.value, ...) {
    dots <- list(...)
    if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
    if (!is.null(trunc)) {
        if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
    }
    if (missing(start.value)) phat <- 1/(mean(x))
    else phat <- start.value
    if (!"method" %in% names(dots)) {
        dots$method <- "Brent"
        if (!"lower" %in% names(dots)) dots$lower = max(c(phat/10, 1e-08))
        if (!"upper" %in% names(dots)) dots$upper = min(c(phat * 10, 0.99))
    }
    if (is.null(trunc)) {
        LL <- function(prob) -sum(dgeom(x, prob, log = TRUE))
    }
    else {
        LL <- function(prob) -sum(dtrunc2("geom", x = x, coef = prob, trunc = trunc, log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = list(prob = phat), data = list(x = x)), dots))
    new("fitsad", result, sad = "geom", trunc = ifelse(is.null(trunc), NaN, trunc))
}


dtrunc2<-function (f, x, trunc, coef, log = FALSE) {
    pf <- get(paste("p", f, sep = ""), mode = "function")
    df <- get(paste("d", f, sep = ""), mode = "function")
    tt <- rep(0, length(x))
    if (!missing(trunc)) tt[x > trunc] <- do.call(df, c(list(x = x[x > trunc],log=TRUE), coef))-do.call(pf, c(list(q = trunc,lower.tail=FALSE,log=TRUE), coef))
    else tt <- do.call(df, c(list(x = x), coef,log=TRUE))
    if (log) return(tt)
    return(exp(tt))
}

dpoilog2 <- function (x, mu, sig, log = FALSE) {
    if (length(mu) > 1 | length(sig) > 1) stop("Vectorization of parameters not implemented")
    if (any(!sads:::is.wholenumber(x))) warning("non integer values in x")
    to.zero <- !is.wholenumber(x) | x < 0
    to.NaN <- NULL
    if (!is.finite(mu)) to.NaN <- 1:length(x)
    if (!is.finite(sig) | sig <= 0) to.NaN <- 1:length(x)
    x[!is.wholenumber(x) | x < 0] <- 0
    mu[!is.finite(mu)] <- 1
    sig[!is.finite(sig) | sig <= 0] <- 1
    y <- poilog::dpoilog(x, mu, sig,log=TRUE)
    y[to.NaN] <- NaN
    y[to.zero] <- 0
    if (any(is.nan(y))) warning("NaNs produced")
    if (log) return(y)
    else return(exp(y))
}

fitpoilog2<- function (x, trunc = 0, ...) {
    dots <- list(...)
    if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
    if (!is.null(trunc)) {
        if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
        else {
            if (trunc == 0) {
                pl.par <- poilog::poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))), zTrunc = TRUE)$par
            }
            else pl.par <- poilog::poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))))$par
            LL <- function(mu, sig) -sum(dtrunc("poilog", x = x, coef = list(mu = mu, sig = sig), trunc = trunc, log = TRUE))
        }
    }
    if (is.null(trunc)) {
        pl.par <- poilogMLE(x, startVals = c(mu = mean(log(x)) + log(0.5), sig = sd(log(x))), zTrunc = FALSE)$par
        LL <- function(mu, sig) -sum(dpoilog2(x, mu, sig, log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = as.list(pl.par), data = list(x = x)), dots))
    new("fitsad", result, sad = "poilog", trunc = ifelse(is.null(trunc), NaN, trunc))
}

qfinder2<-function (dist, want, coef, init = 0) {
    if (want < 0) return(0)
    if (want >= 1) return(Inf)
    q0 <- sum(do.call(dist, c(list(x = 0:init), coef)))
    if (is.nan(q0)) return(NaN)
    if (q0 >= want) return(init)
    step <- 1
    guess <- init + 2
    last <- init + 1
    cum <- q0
    repeat {
        my.q <- do.call(dist, c(list(x = last:guess), coef))
        my.sq <- sum(my.q)
        if (cum + my.sq > want) break
        if (my.sq < 1e-14) stop("quantile function did not converge!")
        last <- guess + 1
        step <- step * 2
        guess <- guess + step
        cum <- cum + my.sq
    }
    my.sq <- cum + cumsum(my.q)
    add <- which(my.sq < want - .Machine$double.eps)
    if (length(add)) last <- last + max(add)
    return(last)
}

lazyRadPred<-function(object){
  sad<-object@sad
  coef<-coef(object)
  trunc<-object@trunc
  S<-length(object@data$x)
  N<-sum(object@data$x)
	distribution <- distr(sad)
	if (distribution == "discrete"){
		## Approximates the [q] function instead of calling it directly to save some
		## computational time (as [q] is inneficiently vectorized)
		y <- 1:(2*N)
		Y <- ppoints(S)
		if(!is.nan(trunc))
			X <- do.call(ptrunc, list(sad, q = y, coef = coef, lower.tail=F, trunc = trunc))
		else {
			psad <- get(paste("p", sad, sep=""), mode = "function")
			qsad <- get(paste("q", sad, sep=""), mode = "function")
			X <- do.call(psad, c(list(q = y, lower.tail = F), coef))
		}
		ab <- approx(x=c(1, X), y=c(0, y), xout=Y, method="constant")$y
		if (any(is.na(ab)))warning('Extreme values generated by radpred. Skipping since they will not appear on plot')
		## Extreme values of abundance are out of bounds for approx. Explicit form:
		#for (i in 1:length(ab)) {
			#if (!is.na(ab[i])) break;
			#cat("Note: extreme values generated by radpred. Calculations will take a while...\n")
			#if(! is.nan(trunc))
				##ab[i] <- do.call(qtrunc, list(sad, p = Y[i], coef = coef, lower.tail=FALSE, trunc = trunc))
			#else
				#ab[i] <- do.call(qsad, c(list(p = Y[i], lower.tail=FALSE), coef))
		#}
    ab[is.na(ab)]<-Inf
	}
	else if(distribution == "continuous"){
		Y <- ppoints(S)
		if(!is.nan(trunc))
			ab <- do.call(qtrunc, list(sad, p = Y, coef = coef, lower.tail=F, trunc = trunc))
		else{
			qsad <- get(paste("q", sad, sep=""), mode = "function")
			ab <- do.call(qsad, c(list(p = Y, lower.tail = F), coef))
		}
	} else{
		stop("Please provide a valid distribution")
  }
  new("rad", data.frame(rank=1:S, abund=ab))
}

#http://andrewgelman.com/2016/06/11/log-sum-of-exponentials/
log_sum_exp<-function(u, v) max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))

pZM<-function(N,q,s){
  #fs<-1/(1:N+q)^s
  logFs<- -s*log(1:N+q)
  H<-Reduce(log_sum_exp,logFs)
  exp(logFs-H)
}
pRadZM<-function(ns,q,s,log=TRUE){
  if(q<=-1)return(ifelse(log,-Inf,0))
  N<-sum(ns)
  rad<-table(factor(ns,levels=1:N))
  ps<-pZM(N,q,s)
  out<-dmultinom(rad,sum(rad),ps,log=log)
  return(out)
}
fitZM<-function(ns){
  fit<-withCallingHandlers(
    nlm(function(qs,ns)-pRadZM(ns,qs[1],qs[2]),c(0,1),ns),
    warning=function(w)if(grepl('Inf replaced by maximum positive value',w))invokeRestart('muffleWarning')
  )
  return(c('q'=fit$estimate[1],'s'=fit$estimate[2],'log-likelihood'=-fit$minimum))
}



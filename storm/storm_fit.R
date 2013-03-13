library(MASS)

fit.skellam.points <- function(data, bins) {
    mcd <- tryCatch(mcd <- cov.mcd(data), error=function(e){return(NULL)})
    if(is.null(mcd) || is.na(mcd)) {
        return(NULL)
    }

    data <- data.frame(data, weight=1)

    dmin <- min(data[,1])
    dmax <- max(data[,1])
    interval <- (dmax - dmin) / bins

    for (i in seq(from=dmin, to=dmax - interval, by=interval)) {
        cmax <- i + interval
        currbin <- data[data$V1 >= i & data$V1 <= cmax,2]
        weight <- 1 / var(currbin)
        if (is.na(weight) || weight == Inf) {
            weight <- 0
        }
        data[data$V1 >= i & data$V1 <= cmax,3] <- weight
    }
    fit <- lm(data[mcd$best,2] ~ data[mcd$best,1], weights=data[mcd$best,3])

    return(coefficients(fit))
}

fit.filter <- function(img) {
    dims <- dim(img)
    mux <- ceiling(dims[2] / 2)
    muy <- ceiling(dims[1] / 2)
    sigmaxstart = dims[2] / 10
    sigmaystart = dims[1] / 10
    img[muy, mux] <- 0

    offset <- min(img)
    scaling <- max(img)

    data <- data.frame(y=1:dims[1], x=rep(1:dims[2], each=dims[1]), val=img[1:length(img)], weights=1)
    data$weights[(mux-1) * dims[1] + muy] <- 0
    fit <- nls(val ~ a+(b-a)*exp(-((x - mux)^2)/(2*sigmax^2)-((y - muy)^2)/(2*sigmay^2)),start=list(sigmax=sigmaxstart,sigmay=sigmaystart, a=offset,b=scaling), data=data, weights=weights, trace=F, control=list(warnOnly=T))

    return(c(coefficients(fit)[1], coefficients(fit)[2]))
}

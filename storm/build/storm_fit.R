library(MASS)

fit.storm.points <- function(data, bins) {
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

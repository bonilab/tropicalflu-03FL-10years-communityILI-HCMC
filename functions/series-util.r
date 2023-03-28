library(compiler)

## —————————————————————————————————————————————————————————————————————————————
#' Convert daily data to weekly data.
#' method can be 'mean', 'median' or 'sum'
## —————————————————————————————————————————————————————————————————————————————
series.toNWeekly <- function( dailySeries, dates = NULL, method,
                              nWeeksPerPeriod = 1, breakMonday = "2009-07-27" ) {
    nWeeklySeries <- NULL

    if (is.null(dates)) {
        daily.dates <- as.Date(names(dailySeries))
    } else {
        daily.dates <- as.Date(dates)
    }

    monday <- as.Date(breakMonday)
    while ( monday > min(daily.dates) ) {
        monday <- monday - (7 * nWeeksPerPeriod)
    }
    nextMonday <- monday + (7 * nWeeksPerPeriod)

    daily.length <- length(dailySeries)

    lastPost <- 0
    for (post in 1 : daily.length) {
        if (post == daily.length || daily.dates[post + 1] >= nextMonday) {

            if ( lastPost < post ) {
                if ( (lastPost+1) < post ) {
                    if (method == "mean") {
                        curNWeekData <- mean(dailySeries[ (lastPost+1) : post ], na.rm = T)
                    } else if (method == "median") {
                        curNWeekData <- median(dailySeries[ (lastPost+1) : post ], na.rm = T)
                    } else if (method == "sum") {
                        curNWeekData <- sum(dailySeries[ (lastPost+1) : post ], na.rm = T)
                    } else {
                        stop("Aggregating method must be 'mean' or 'median' or 'sum'.")
                    }

                } else {
                    curNWeekData <- dailySeries[post]
                }

                lastPost <- post

            } else {
                curNWeekData <- NA
            }

            nWeeklySeries <- append(nWeeklySeries, curNWeekData)
            names(nWeeklySeries)[ length(nWeeklySeries) ] <- as.character(monday)

            monday <- nextMonday
            nextMonday <- monday + (7 * nWeeksPerPeriod)
        }
    }


    return( nWeeklySeries )
}

## Compile the function
series.toNWeekly <- cmpfun(series.toNWeekly)



## —————————————————————————————————————————————————————————————————————————————
## Generate a white-noise series that has the same mean and the same standard deviation of the
## given series.
## —————————————————————————————————————————————————————————————————————————————
series.simulateWhiteNoise <- function( dataSeries ) {
    dataMean <- mean(dataSeries, na.rm = T)
    dataSd <- sd(dataSeries, na.rm = T)

    availables <- which(!is.na(dataSeries))
    nAvailables <- length(availables)

    dataSeries[availables] <- rnorm(nAvailables, dataMean, dataSd)

    return(dataSeries)
}

## Compile the function
series.simulateWhiteNoise <- cmpfun(series.simulateWhiteNoise)



## —————————————————————————————————————————————————————————————————————————————
## Convert daily data to monthly data.
## —————————————————————————————————————————————————————————————————————————————
series.toMonthly <- function( dailySeries, dates = NULL, representativeDayOfMonth = 15 ) {


    firstOfNextMonth <- function(curDate) {
        curDate <- as.character(as.Date(curDate))
        substr(curDate, 9, 10) <- "01"

        nextMonth <- as.Date(curDate) + 31
        nextMonth <- as.character(nextMonth)
        substr(nextMonth, 9, 10) <- "01"

        return( as.Date(nextMonth) )
    }


    monthlySeries <- NULL

    thisMonth <- as.Date("2009-08-01")
    nextMonth <- firstOfNextMonth(thisMonth)

    daily.length <- length(dailySeries)

    if (is.null(dates)) {
        daily.dates <- as.Date(names(dailySeries))
    } else {
        daily.dates <- as.Date(dates)
    }

    lastPost <- 1
    for (post in 1 : daily.length) {
        while (daily.dates[post] >= nextMonth) {
            if ( lastPost < post ) {
                if ( (post - 1 - lastPost) > 1 ) {
                    curMonthData <- mean(dailySeries[ lastPost : (post-1) ], na.rm = T)
                } else {
                    curMonthData <- dailySeries[lastPost]
                }
                lastPost <- post

            } else {
                curMonthData <- NA
            }

            if (post > 1) {
                monthlySeries <- append(monthlySeries, curMonthData)
                names(monthlySeries)[ length(monthlySeries) ] <- as.character(thisMonth + representativeDayOfMonth - 1)
            }

            thisMonth <- nextMonth
            nextMonth <- firstOfNextMonth(thisMonth)
        }
    }

    if ( (daily.length - lastPost) > 1 ) {
        curMonthData <- mean(dailySeries[ lastPost : daily.length ], na.rm = T)
    } else {
        curMonthData <- dailySeries[ lastPost ]
    }
    monthlySeries <- append(monthlySeries, curMonthData)
    names(monthlySeries)[ length(monthlySeries) ] <- as.character(thisMonth + representativeDayOfMonth - 1)

    return( monthlySeries )
}

## Compile the function
series.toMonthly <- cmpfun(series.toMonthly)



## —————————————————————————————————————————————————————————————————————————————
#'
## —————————————————————————————————————————————————————————————————————————————
series.movingAverage <- function(dataSeries, dates = NULL, windowSize, preserveMissing = F) {
    result <- dataSeries
    halfWindowSize <- round (windowSize - 1) / 2

    len <- length(dataSeries)

    if ( is.null(dates) ) {
        dates <- as.Date(names(dataSeries))
    } else if ( length(dates) <= 1 && is.na(dates) ) {
        dates <- 1 : length(dataSeries)
    } else {
        dates <- as.Date(dates, origin = "1970-01-01")
    }

    for (i in 1 : len) {
        curDate <- dates[i]
        startDate <- curDate - halfWindowSize
        endDate <- curDate + halfWindowSize

        startIndex <- i - halfWindowSize
        endIndex <- i + halfWindowSize
        if (startIndex < 1) startIndex <- 1
        if (endIndex > len) endIndex <- len

        while (dates[startIndex] < startDate) startIndex <- startIndex + 1
        while (dates[endIndex] > endDate) endIndex <- endIndex - 1

        result[i] <- mean(dataSeries[startIndex : endIndex], na.rm = T)
    }

    if (preserveMissing) {
        result[ which(is.na(dataSeries)) ] <- NA    # preserve missing data points
    }

    return(result)
}

## Compile the function
series.movingAverage <- cmpfun(series.movingAverage)



## —————————————————————————————————————————————————————————————————————————————
#'
## —————————————————————————————————————————————————————————————————————————————
series.movingPastAverage <- function(dataSeries, dates = NULL, windowSize = 7) {
    result <- dataSeries

    len <- length(dataSeries)

    if (is.null(dates)) {
        dates <- as.Date(names(dataSeries))
    } else {
        dates <- as.Date(dates)
    }

    for (i in 1 : len) {
        curDate <- dates[i]
        startDate <- curDate - windowSize + 1
        startIndex <- i - windowSize + 1

        if (startIndex < 1) startIndex <- 1
        while (dates[startIndex] < startDate) startIndex <- startIndex + 1

        result[i] <- mean(dataSeries[startIndex : i], na.rm = T)
    }

    result[ which(is.na(dataSeries)) ] <- NA    # preserve missing data points

    return(result)
}

## Compile the function
series.movingPastAverage <- cmpfun(series.movingPastAverage)



## —————————————————————————————————————————————————————————————————————————————
#'
## —————————————————————————————————————————————————————————————————————————————
series.discretize <- function(input, nOutSteps = 7, inputIsDiscrete = F) {
    sNames <- names(input)

    input <- as.numeric(input)

    maxInput <- max(input)
    minInput <- min(input)

    inputRange <- maxInput - minInput
    if (inputIsDiscrete) {
        inputRange <- inputRange + 1
    }
    stepSize <- inputRange / nOutSteps

    inputLen <- length(input)
    for (i in 1 : inputLen) {
        input[i] <- floor( (input[i] - minInput) / stepSize ) + 1

        if (input[i] > nOutSteps) {
            input[i] <- nOutSteps
        }
    }

    input <- as.factor(input)
    names(input) <- sNames
    return(input)
}



## —————————————————————————————————————————————————————————————————————————————
#'
## —————————————————————————————————————————————————————————————————————————————
series.movingZscore <- function( dataSeries, dates = NULL,
                                 windowSize = 365, minAvailDays = NULL ) {
    if (is.null(dates)) {
        dates <- names(dataSeries)
    }
    dates <- as.Date( dates, origin = "1970-01-01" )

    if (is.null(minAvailDays)) {
        minAvailDays <- windowSize / 2
    }
    if (minAvailDays > windowSize) {
        minAvailDays <- windowSize
    }


    wHalfSize <- floor(windowSize / 2)

    seriesLen <- length(dataSeries)
    zscores <- rep( NA, seriesLen )
    names(zscores) <- names(dataSeries)

    for ( pos in 1 : seriesLen ) {
        if ( is.na(dataSeries[pos]) ) {
            next
        }

        curDate <- dates[pos]
        sDate <- curDate - wHalfSize
        eDate <- curDate + wHalfSize

        sPos <- max( 1, pos - wHalfSize )
        ePos <- min( seriesLen, pos + wHalfSize )

        while ( dates[sPos] < sDate ) { sPos <- sPos + 1 }
        while ( dates[ePos] > eDate ) { ePos <- ePos - 1 }

        nAvailDays <- sum( !is.na(dataSeries[sPos:ePos]) )
        if (nAvailDays >= minAvailDays && nAvailDays > 0) {
            rangeMean <- mean( dataSeries[sPos:ePos], na.rm = T )
            rangeSd <- sd( dataSeries[sPos:ePos], na.rm = T )
            zscores[pos] <- ( dataSeries[pos] - rangeMean ) / rangeSd
        }
    }

    return(zscores)
}

## Compile the function
series.movingZscore <- cmpfun(series.movingZscore)


## —————————————————————————————————————————————————————————————————————————————
series.generateLagSeries <- function(series, lag) {
    seriesLength <- length(series)
    lagSeries <- append( rep(NA, lag), series[1 : (seriesLength - lag)] )
    return(lagSeries)
}
## —————————————————————————————————————————————————————————————————————————————


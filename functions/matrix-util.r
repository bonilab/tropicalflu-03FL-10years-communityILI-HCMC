library(compiler)
#source("./scripts/utils/series-util.r")

####################################################################################################
vectorSum <- function(vecA, vecB) {
    len <- length(vecA)
    newVec <- c()

    for (i in 1 : len) {
        if ( is.na(vecA[i]) ) {
            newVec[i] <- vecB[i]
        } else if ( is.na(vecB[i]) ) {
            newVec[i] <- vecA[i]
        } else {
            newVec[i] <- vecA[i] + vecB[i]
        }
    }
    names(newVec) <- names(vecA)
    return(newVec)
}


####################################################################################################
matrix.removeNaColumns <- function( dataMatrix ) {
    naCols <- which(colSums( !is.na(dataMatrix) ) == 0)
    if (length(naCols) > 0) {
        newMatrix <- dataMatrix[, -naCols]
    }
    return(newMatrix)
}



####################################################################################################
#' Calculate the Pearson's correlation between each individual row and the column mean of the
#' matrix.
#' @param rowsReplacedByWhiteNoise The list of indices of rows that must be replaced by a
#'                                 white-noise series before testing the correlations.
#'                                 By default, this parameter receives NULL value, meaning that
#'                                 all the rows will be preserved as they are given.
#'                                 If this parameter receives "all" for its value, all the rows in
#'                                 the given matrix will be replaced by white-noise series.
#' @param rowsCalculated    The list of indices of rows of which the Pearson's correlation will be
#'                          calculated and returned. Specifying this parameter may help improving
#'                          the computational speed.
#'                          By default, the correlations of all rows will be calculated.
####################################################################################################
matrix.correlRowVsColMean <- function( dataMatrix,
                                       aggregatedSeries = NULL,
                                       smthWindowSize = 0,
                                       smthAggreSeries = T,
                                       rowsReplacedByWhiteNoise = NULL,
                                       rowsCalculated = "all",
                                       returnFullStats = F ) {

    if (!is.null(rowsReplacedByWhiteNoise)) {
        if (rowsReplacedByWhiteNoise == "all") {
            rowsReplacedByWhiteNoise <- c(1 : nrow(dataMatrix))
        }
        for (row in rowsReplacedByWhiteNoise) {
            dataMatrix[row, ] <- series.simulateWhiteNoise( dataMatrix[row, ] )
        }
    }

    if (is.null(aggregatedSeries)) {
        aggregatedSeries <- colMeans(dataMatrix, na.rm = T)
    }

    if (smthAggreSeries && smthWindowSize > 1) {
        aggregatedSeries <- series.movingAverage( aggregatedSeries,
                                                  windowSize = smthWindowSize, preserveMissing = F )
    }

    if (is.null(rowsCalculated) || rowsCalculated == "all") {
        rowsCalculated <- c( 1 : nrow(dataMatrix) )
    }

    corFrame <- NULL
    for (row in rowsCalculated) {
        rowData <- dataMatrix[row, ]
        if (smthWindowSize > 1) {
            rowData <- series.movingAverage( rowData,
                                             windowSize = smthWindowSize, preserveMissing = F )
        }

        if (returnFullStats) {
            corr <- cor.test( aggregatedSeries, rowData, conf.level = 0.95,
                              method = "pearson", na.action = "na.omit" )
            corr <- data.frame( rowName = rownames(dataMatrix)[row],
                                correlation = corr$estimate,
                                ciLo_95 = corr$conf.int[1],
                                ciHi_95 = corr$conf.int[2],
                                pVal = corr$p.value)
            corFrame <- rbind(corFrame, corr)

        } else {
            corr <- cor( aggregatedSeries, rowData,
                         use = "complete.obs", method = "pearson" )
            corFrame <- append(corFrame, corr)
        }
    }

    if (returnFullStats) {
        rownames(corFrame) <- NULL
    } else {
        names(corFrame) <- rownames(dataMatrix)[rowsCalculated]
    }

    # acf(aggregatedSeries, lag.max = 0.7 * length(aggregatedSeries), na.action = na.pass, col = "grey")

    return(corFrame)
}

## Compile the function
matrix.correlRowVsColMean <- cmpfun(matrix.correlRowVsColMean)



####################################################################################################
## Convert the values on each row of the given matrix to z-scores.
## minAvailCol :  If the number of available columns on a row is less than this number,
##                the returned z-scores for all the columns of that row will be NAs.
## getOneCol   :  If this parameter is 0 or FALSE, the z-scores of each single cell of the
##                given matrix will be calculated and returned.
##                If this parameter receives a positive number, only the z-scores of the
##                column that has the index of this given value will be calculated and returned.
####################################################################################################
# rowZScore <- function( dataMatrix, minAvailCol = 30, getOneCol = 0 ) {
#     sd <- apply(dataMatrix, 1, sd, na.rm = T)
#     mean <- rowMeans(dataMatrix, na.rm = T)
#
#     nRow <- nrow(dataMatrix)
#     nCol <- ncol(dataMatrix)
#
#     if (getOneCol == 0) {
#         zscore <- matrix( data = NA, nrow = nRow, ncol = nCol,
#                           dimnames = list(rownames(dataMatrix), colnames(dataMatrix)) )
#     } else {
#         zscore <- rep( NA, nRow )
#         names(zscore) <- rownames(dataMatrix)
#     }
#
#     for (row in 1 : nRow) {
#         if ( sum(!is.na(dataMatrix[row, ])) >= minAvailCol ) {
#             if (getOneCol == 0) {
#                 for (col in 1 : nCol) {
#                     if ( !is.na(dataMatrix[row, col]) ) {
#                         zscore[row, col] <- ( dataMatrix[row, col] - mean[row] ) / sd[row]
#                     }
#                 }
#
#             } else {
#                 if ( !is.na(dataMatrix[row, getOneCol]) ) {
#                     zscore[row] <- ( dataMatrix[row, getOneCol] - mean[row] ) / sd[row]
#                 }
#             }
#         }
#     }
#
#     return(zscore)
# }



# Calculate yearly detrended z-score for each row (clinic)
# rowZScore.yearDetrend <- function( dataMatrix ) {
#
#     findYearBreakColumns <- function( dataMatrix ) {
#         tets <- c('2010-02-14', '2011-02-03', '2012-01-23', '2013-02-10', '2014-01-31',
#                   '2015-02-19', '2016-02-08')
#         tets <- as.Date(tets)
#
#         dates <- as.Date(colnames(dataMatrix))
#         nDates <- length(dates)
#
#         breaks <- 1
#         i <- 1
#         for (k in 2 : nDates) {
#             if (dates[k - 1] < tets[i] && tets[i] <= dates[k]) {
#                 breaks <- c(breaks, k)
#
#                 i <- i + 1
#                 if (i > length(tets)) {
#                     break
#                 }
#             }
#         }
#         breaks <- c(breaks, nDates + 1)
#
#         return(breaks)
#     }
#
#
#     breaks <- findYearBreakColumns(dataMatrix)
#     nBreaks <- length(breaks)
#
#     nRow <- nrow(dataMatrix)
#     nCol <- ncol(dataMatrix)
#     zscore <- matrix( data = NA, nrow = nRow, ncol = nCol,
#                       dimnames = list(rownames(dataMatrix), colnames(dataMatrix)) )
#
#     for ( i in 1 : (nBreaks - 1) ) {
#         sCol <- breaks[i]
#         eCol <- breaks[i + 1] - 1
#
#         zscore[ , sCol:eCol] <- rowZScore( dataMatrix[ , sCol:eCol] )
#     }
#
#     return(zscore)
# }




# rowZScore.movingLocal <- function( dataMatrix, dates = NULL, windowSize = 365 ) {
#     ## FIXME: This fuction implicitly assumes that the consecutive columns in
#     ##        the dataMatrix represent consecutive days.
#     windowHalfSize <- floor(windowSize / 2)
#
#     nRow <- nrow(dataMatrix)
#     nCol <- ncol(dataMatrix)
#     zscore <- matrix( data = NA, nrow = nRow, ncol = nCol,
#                       dimnames = list(rownames(dataMatrix), colnames(dataMatrix)) )
#
#     for ( col in 1 : nCol ) {
#         sCol <- max( 1, col - windowHalfSize )
#         eCol <- min( nCol, col + windowHalfSize )
#         localColIndex <- col - sCol + 1
#
#         zscore[ , col] <- rowZScore( dataMatrix[ , sCol:eCol], getOneCol = localColIndex )
#     }
#
#     return(zscore)
# }



# Convert the values in each row to yearly percentiles.
# rowPScore.yearDetrend <- function(dataMatrix) {
#
#     toPScore <- function(series) {
#         pscore <- series
#
#         names(series) <- NULL
#         series <- series[which(!is.na(series))]
#         series <- sort(series)
#         seriesLen <- length(series)
#
#         for (i in 1 : length(pscore)) {
#             if (is.na(pscore[i])) {
#                 next
#             }
#
#             minIndex <- 1
#             while (series[minIndex] < pscore[i]) {
#                 minIndex <- minIndex + 1
#             }
#             maxIndex <- minIndex
#             while (maxIndex < seriesLen && series[maxIndex + 1] == pscore[i]) {
#                 maxIndex <- maxIndex + 1
#             }
#             pscore[i] <- ((minIndex - 1) + (maxIndex - 1)) / (2 * (seriesLen - 1))
#         }
#
#         return(pscore)
#     }
#
#
#     breaks <- findYearBreakColumns(dataMatrix)
#     nBreaks <- length(breaks)
#
#     nRow <- nrow(dataMatrix)
#     pscore <- dataMatrix
#
#     for ( i in 1 : (nBreaks - 1) ) {
#         sCol <- breaks[i]
#         eCol <- breaks[i + 1] - 1
#
#         for (row in 1 : nRow) {
#             pscore[row, sCol:eCol] <- toPScore(dataMatrix[row, sCol:eCol])
#         }
#     }
#
#     return(pscore)
# }


#' Convert daily data to weekly data.
#' @return A matrix with clinics in rows and weeks in columns (column names are the first days - Monday - of the weeks).
#'         Weekly data are weekly means of the daily data.
# toWeeklyData <- function( dailyData ) {
#     weeklyData <- NULL
#
#     monday <- as.Date("2009-07-27")
#     nextMonday <- monday + 7
#
#     daily.nCols <- ncol(dailyData)
#     daily.dates <- as.Date(colnames(dailyData))
#
#     lastCol <- 1
#     for (col in 1 : daily.nCols) {
#         if (daily.dates[col] >= nextMonday) {
#             if ( lastCol <= (col - 1) ) {
#                 if ( (col - 1 - lastCol) > 1 ) {
#                     curWeekData <- rowMeans(dailyData[ , lastCol : (col-1) ], na.rm = T)
#                 } else {
#                     curWeekData <- dailyData[ , lastCol]
#                 }
#                 weeklyData <- cbind(weeklyData, curWeekData)
#                 colnames(weeklyData)[ ncol(weeklyData) ] <- as.character(monday)
#
#                 lastCol <- col
#             }
#
#             monday <- nextMonday
#             nextMonday <- monday + 7
#         }
#     }
#
#     if ( (daily.nCols - lastCol) > 1 ) {
#         curWeekData <- rowMeans(dailyData[ , lastCol : daily.nCols ], na.rm = T)
#     } else {
#         curWeekData <- dailyData[ , lastCol]
#     }
#     weeklyData <- cbind(weeklyData, curWeekData)
#     colnames(weeklyData)[ ncol(weeklyData) ] <- as.character(monday)
#
#     return( as.matrix(weeklyData) )
# }



## Each column is a day, each row is an independent time series
# matrix.movingAverage <- function(dailyMatrix, windowSize = 7, preserveMissing = T) {
#     result <- dailyMatrix
#     halfWindowSize <- round (windowSize - 1) / 2
#
#     nCols <- ncol(dailyMatrix)
#     dates <- as.Date( colnames(dailyMatrix) )
#
#     for (col in 1 : nCols) {
#         curDate <- dates[col]
#         startDate <- curDate - halfWindowSize
#         endDate <- curDate + halfWindowSize
#
#         startCol <- col - halfWindowSize
#         endCol <- col + halfWindowSize
#         if (startCol < 1) startCol <- 1
#         if (endCol > nCols) endCol <- nCols
#
#         while (dates[startCol] < startDate) startCol <- startCol + 1
#         while (dates[endCol] > endDate) endCol <- endCol - 1
#
#         result[ , col] <- rowMeans(dailyMatrix[ , startCol : endCol], na.rm = T)
#     }
#
#     if (preserveMissing) {
#         ## preserve missing data points
#         result[ which(is.na(dailyMatrix)) ] <- NA
#     }
#
#     return(result)
# }




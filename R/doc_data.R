#' Radomized data for testing models
#'
#' Contains 5 Variables, one dependent, 4 independent. The fourth independent is correlated with the dependent
#'
#' @docType data
#'
#' @usage data(data)
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references K.T.Krahl (2023)
#'
#' @source \href{https://www.hanseatic-statistics.de/}
#'
#' @examples
#' data(data)
#' \donttest{LinReg('dp',c('iv_1','iv_2','iv_3',iv_4),data, bootstrap = FALSE, Number_Bootstrapps=1000, outlier_controll = FALSE, plot=TRUE)}
"data"

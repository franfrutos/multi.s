#' Computes split-half internal consistency measures via permutation for human behavior tasks.
#'
#' @description Computes split-half reliability indexes ussing a permutation approach for the difference in outcome.
#'
#' @param data A data frame containing at least an outcome variable, the variable you want to analyze and a variable for subject's identificators.
#'
#' @param outcome The name of your outcome variable in strings. The default is "RT".
#'
#' @param average The type of averaging that you want to apply for computing split-half. The default is "mean".
#'
#' @param type The type of score that will be applied. Default is "difference".
#'
#' @param q Quantile to apply if mean or median is not used as averaging method. Default is .80.
#'
#' @param permutations Number of times that you want to compute the split-half. Default is 10, but more is recommended.
#'
#' @param variable Name of the variable to be analyzed in strings (e.g. "Congruency"). Default is set to NULL, and is required to be filled.
#'
#' @param condition A variable's name containing different conditions for splitting analysis. This argument is optional.
#'
#' @param subject A variable's name for subject's indentification. The default is "subject".
#'
#' @param include_block A logical that determine if each half computed should apply constrains depending on blocking. If TRUE, split-half will be computed splitting by block. Default is FALSE.
#'
#' @param block The name of the variable that you want to apply constrains in strings. If include_block is TRUE, block is required to be filled. Default is null.
#'
#' @param return_iteration A logical controlling if te split-half estimation for each iteration is returned in the output. Default is FALSE.
#'
#' @return A data frame containing split-half and spearman-brown estimators for n permutations.
#'
#' @import tidyr
#' @import Rcpp
#' @import grid
#' @import ggplot2
#' @importFrom stats complete.cases cor median na.omit quantile sd
#' @importFrom robustbase colMedians
#' @importFrom matrixStats colQuantiles
#' @importFrom dplyr select summarise group_by mutate n_distinct
#' @importFrom plyr arrange
#' @useDynLib multi.s, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
#'
splith <- function(data, outcome = "RT", average = "mean", type = "difference", q = .80, permutations = 10,
                   variable = NULL, condition = NULL, subject = "subject", include_block = FALSE,
                   block = NULL, return_iterations = FALSE) {


  if (is.data.frame(data) == FALSE) {
    stop(paste0("you have to provide a data frame, not a", typeof(data)))
  }

  # set data to data frame for avoidding tibble issues:
  data <- as.data.frame(data)

  #save function settings:
  call = list(Permutations = permutations,
              outcome = outcome,
              average = average,
              permutations = permutations,
              variable = variable,
              condition = condition,
              include_block = include_block,
              block = block,
              return_iterations = return_iterations)

  #creating a list for output
  output <- list(call = call,
                 data = data)

  iteration <- 1:permutations

  #how many participants?
  n_par <- n_distinct(data[, subject])

  #create a vector for participants
  plist <- sort(unique(data[, subject]))

  #create a vector for condition (if given)
  if (!is.null(condition)){
    clist <- unique(data[, condition])
  } else {
    clist <- 1
    condition <- "condition"
    data$condition <- rep(1, each = nrow(data))
  }

  #create a vector for Blocks
  if (include_block){
    if(is.null(block)) {
      warning("If include_block is set to TRUE, and no blocking variable is provided, there will be no stratifications in splits.")
    } else {
      blist <- sort(unique(data[, block]))
    }
  } else {
    blist <- 1
    block <- "block"
    data$block <- rep(1, each = nrow(data))
  }

  #create a vector with variables
  if (!is.null(variable)){
    if (!length(unique(data[, variable])) %in% 1:2){
      stop("your variable needs one or two values values")
    } else if (length(unique(data[, variable])) == 2) {
      if (type != "difference") stop('If your variable has two values, "type" must be set to difference.')
      vlist <- sort(unique(data[, variable]))
    } else {
      if (type == "difference") stop('A difference cannot be computed if there is only one level in "variable".')
      vlist <- 1
      variable <- "variable"
      data$variable <- rep(1, each = nrow(data))    }
  } else {
    if (type == "average") {
      vlist <- 1
      variable <- "variable"
      data$variable <- rep(1, each = nrow(data))
    } else {
      stop('If you set "type" to "average", you need to select a variable with two values.')
    }
  }


  #  function to detect if a vector is binary (used for automaticly check if outcome is RT or ACC type)
  # is.binary = function(v, naVal="NA") {
  #   if (!is.numeric(v)) stop("Only numeric vectors are accepted.")
  #   vSet = unique(v)
  #   if (!missing(naVal)) vSet[vSet == naVal] = NA
  #   vSet = vSet[!is.na(vSet)]
  #
  #   !(any(as.integer(vSet) != vSet) || length(vSet) > 2)
  # }


  # checks whether user difference score is based on means or medians
  #if (!is.binary(data[, outcome])){
    if (average == "mean") {
      ave_fun <- function(val) {
        tryCatch(
          error = function(err) NaN,
          colMeans(val)
        )
      }
      ave_fun_basic <- function(val) {
        mean(val)
      }
    } else if (average == "median") {
      ave_fun <- function(val) {
        tryCatch(
          error = function(err) NaN,
          robustbase::colMedians(val)
        )
      }
      ave_fun_basic <- function(val) {
        median(val)
      }
    } else {
      ave_fun <- function(val, q) {
        tryCatch(
          error = function(err) NaN,
          matrixStats::colQuantiles(val, probs = q)
        )
      }
      ave_fun_basic <- function(val, q) {
        quantil(val, q)
      }
    }
  #} else {
    #stop("You need to provide RT data")
  #}

  # create the data.frame to populate
  findata <-
    data.frame(
      Condition = rep(clist, each = (length(plist) * length(iteration))),
      Participant = rep(plist, each = length(iteration)),
      Iteration = rep(iteration, times = (
        length(clist) * length(plist)
      )),
      bias1 = NA,
      bias2 = NA
    )

  # loop counter
  l <- 1

  # participant loop counter for progress bar
  ppt <- 1

  # create vectors to contain both halfs to be compared
  bias1v <-
    vector(length = (length(clist) * length(plist) *
                       length(iteration)))
  bias2v <-
    vector(length = (length(clist) * length(plist) *
                       length(iteration)))

  for(j in clist){

    pb <- txtProgressBar(min = 0,
                         max = n_par,
                         style = 3)
    setTxtProgressBar(pb, 0)

    for(i in plist){

      tempcm.1 <- list()
      tempcm.2 <- list()

      if (type == "difference") {
        tempim.1 <- list()
        tempim.2 <- list()
      }

      for(k in blist){

        tempcvb <-
          subset(data[, outcome],
                 data[, subject] == i &
                   data[, block] == k &
                   data[, condition] == j &
                   data[, variable] == vlist[1])

        if (length(tempcvb) != 0){
          midtrial.con <- sum(!is.na(tempcvb)) / 2

          tempcm <- samploop(a = matrix(nrow = length(tempcvb), ncol = permutations, 0),
                             b = tempcvb)


          tempcm.1[[k]] <- tempcm[1:floor(midtrial.con), ]
          tempcm.2[[k]] <- tempcm[(floor(midtrial.con) + 1):length(tempcvb), ]
        }

        if (type == "difference") {
          tempivb   <-
            subset(data[, outcome],
                   data[, subject] == i &
                     data[, block] == k &
                     data[, condition] == j &
                     data[, variable] == vlist[2])

          if (length(tempivb) != 0){
            midtrial.incon <- sum(!is.na(tempivb)) / 2

            tempim <- samploop(a = matrix(nrow = length(tempivb), ncol = permutations, 0),
                               b = tempivb)

            tempim.1[[k]] <- tempim[1:floor(midtrial.incon), ]
            tempim.2[[k]] <- tempim[(floor(midtrial.incon) + 1):length(tempivb), ]
            }
          }
        }

      tempcm.1 <- do.call(rbind, tempcm.1)
      tempcm.2 <- do.call(rbind, tempcm.2)
      if (type == "difference") {
        tempim.1 <- do.call(rbind, tempim.1)
        tempim.2 <- do.call(rbind, tempim.2)
      }

      if (type == "difference") {
        if (average %in% c("mean", "median")) {
          bias1v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.1)  - ave_fun(tempim.1)
          bias2v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.2)  - ave_fun(tempim.2)
        } else {
          bias1v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.1, q)  - ave_fun(tempim.1, q)
          bias2v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.2, q)  - ave_fun(tempim.2, q)
        }
      } else {
        if (average %in% c("mean", "median")) {
          bias1v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.1)
          bias2v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.2)
        } else {
          bias1v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.1, q)
          bias2v[l:(l + permutations - 1)] <-
            ave_fun(tempcm.2, q)
        }
      }


      l <- l + permutations
      ppt <- ppt + 1
      setTxtProgressBar(pb, ppt)

    }
    ppt <- 1 # reset the progress bar

    print(paste("condition", j, "complete"))
  }

  print("Calculating split half estimates")

  findata$bias1 <- bias1v
  findata$bias2 <- bias2v

  if (sum(is.na(findata$bias1)) +
      sum(is.na(findata$bias2)) > 0){
    print("the following are participants/conditions with missing data")
    omitted <- findata[!complete.cases(findata),]
    output$omitted <- omitted
    print(unique(omitted[c("Condition", "Participant")]))
    print(
      "note: these iterations will be removed from the split half
            reliability calculations, in that condition"
    )
    warning(
      "Bias indices missing:
              at least one participant has missing data from at one condition
              These cases are removed from calculating reliability estimates
              $omitted contains the missing cases"
    )
  }

  # remove NA rows
  findata2 <-  na.omit(findata)
  findata2$Iteration <- as.factor(findata2$Iteration)

  out <- findata2 %>%
    dplyr::group_by(Condition, Iteration) %>%
    dplyr::summarise(
      n = round(sum(!is.na(bias1)), 2),
      splithalf = cor(bias1, bias2, use = "pairwise.complete"),
      spearmanbrown = (2 * cor(bias1, bias2,
                               use = "pairwise.complete")) /
        (1 + (2 - 1) * abs(
          cor(bias1, bias2,
              use = "pairwise.complete")
        ))
    )
  round.to = 2
  out2 <- out %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(
      n = mean(n),
      splithalf_estimate = round(mean(splithalf), round.to),
      splithalf95CI_lower = round(quantile(splithalf, c(.025), names = F), round.to),
      splithalf95CI_upper = round(quantile(splithalf, c(.975), names = F), round.to),
      spearmanbrown_estimate = round(mean(spearmanbrown), round.to),
      spearmanbrown95CI_lower = round(quantile(spearmanbrown, c(.025), names = F), round.to),
      spearmanbrown95CI_upper = round(quantile(spearmanbrown, c(.975), names = F), round.to)
    ) %>%
    as.data.frame()

  colnames(out2) <- c(
    "Condition",
    "n",
    "splithalf",
    "95_low",
    "95_high",
    "spearmanbrown",
    "SB_low",
    "SB_high"
  )


  print(paste(
    "split half estimates for",
    permutations,
    "random splits",
    sep = " "
  ))

  if (return_iterations == TRUE) {
    output$estimates <- out
    output$raw_bias <- findata
  }

  output$final_estimates <- out2

  return(output)
}

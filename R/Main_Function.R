
utils::globalVariables(c(
  "Alpha", "Gamma", "Gvariable", "Hier", "Dataset", "Xvariable",
  "coef", "intercept", "lm", "pred", "slope",
  "type", "value", "x_max", "x_min"
))

#' Calculate stability for a single time series.
#'
#' \code{iSTAY_Single} computes the stability of order q for a single time series.
#'
#' @param data A \code{vector} of time series data, or a \code{data.frame} with sampling units as rows and time points as columns.
#' @param order.q a numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime Logical (\code{TRUE} or \code{FALSE}), indicating whether to use all time points in the data.
#' @param start_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the starting column (time point) for the analysis interval.
#' @param end_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the ending column (time point) for the analysis interval.
#'
#'
#'
#' @return a dataframe with columns: \cr
#' Dataset--the input dataset \cr
#' Order_q--order of stability \cr
#' Stability--stability measures of order q
#'
#'
#' @examples
#' # Compute the stability of individual plots
#' data("Data_Jena_20_metacommunities")
#' individual_plots <- do.call(rbind, Data_Jena_20_metacommunities)
#' output_individual_plots <- iSTAY_Single(data = individual_plots, order.q=c(1,2), Alltime = TRUE)
#' output_individual_plots
#'
#' # Compute the stability of individual populations
#' data("Data_Jena_76_metapopulations")
#' individual_populations <- do.call(rbind, Data_Jena_76_metapopulations)
#' output_individual_populations <- iSTAY_Single(data = individual_populations, order.q = c(1,2), Alltime = TRUE)
#' output_individual_populations
#'
#' @export

iSTAY_Single <- function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL){
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  NA_num <- sum(is.na(data))
  if (NA_num != 0) {
    stop("There are some NA in the data.")
  }
  if (sum(which(order.q <= 0))) {
    stop("Order q need larger than 0.")
  }
  if (Alltime != TRUE) {
    if (is.null(start_T) | is.null(end_T)) {
      stop("Need to set the length of time series for calculating.")
    }
    if ((start_T >= end_T) | length(start_T) != 1 | length(end_T) != 1) {
      stop("Starting and ending time need to be a number, and ending time needs larger than starting time.")
    }
  }
  stability <- function(vector, q) {
    if (sum(vector != 0) == 0) {
      out <- NA
    }
    else {
      K <- length(vector)
      if (q == 1) {
        H <- sum((vector[vector > 0] / sum(vector)) * log(vector[vector > 0] / sum(vector))) * (-1)
        out <- (exp(H) - 1) / (K - 1)
      }
      else {
        up <- 1 - (sum((vector[vector > 0] / sum(vector))^q))^(1 / (1-q))
        out <- up / (1 - K)
      }
    }
    return(out)
  }
  if (Alltime == TRUE) {
    stab <- as.matrix(apply(data, 1, function(vec) sapply(order.q, function(qq) stability(vec, q = qq))))
  }
  else {
    subdata <- data[, c(start_T : end_T)]
    stab <- as.matrix(apply(subdata, 1, function(vec) sapply(order.q, function(qq) stability(vec, q = qq))))
  }
  result <- data.frame(Assemblage = rep(rownames(as.data.frame(data)), length(order.q)), 
                       Order_q = rep(order.q, each = nrow(data)), 
                       Stability = as.vector(t(stab)))
  colnames(result)[1] <- c("Dataset")
  return(result)
}





#' Calculate stability and synchrony for multiple time series.
#'
#' \code{iSTAY_Multiple} computes gamma, alpha, and beta stability, as well as synchrony, for multiple time-series data.
#'
#' @param data A \code{data.frame} containing multiple time series data, with sampling units as rows and time points as columns, or a \code{list} of \code{data.frames} with each data frame representing multiple time series.
#' @param order.q A numerical vector specifying the orders of stability and synchrony. Default is c(1,2).
#' @param Alltime Logical (\code{TRUE} or \code{FALSE}), indicating whether to use all time points in the data.
#' @param start_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the starting column (time point) for the analysis interval.
#' @param end_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the ending column (time point) for the analysis interval.
#'
#'
#' @return a data frame with the following columns: \cr
#'  Dataset--the input dataset \cr
#'  Order_q--order of stability or synchrony \cr
#'  Gamma, Alpha, Beta--stability measures of order q \cr 
#'  Synchrony--synchrony measure of order q
#'
#' @examples
#' # Stability of multiple time series
#' data("Data_Jena_20_metacommunities")
#' metacommunities <- Data_Jena_20_metacommunities
#' output_metacommunities <- iSTAY_Multiple(data = metacommunities, order.q = c(1,2), Alltime = TRUE)
#' output_metacommunities
#'
#' # Stability of metapopulations
#' data("Data_Jena_76_metapopulations")
#' metapopulations <- Data_Jena_76_metapopulations
#' output_metapopulations <- iSTAY_Multiple(data = multiple_species, order.q = c(1,2), Alltime = TRUE)
#' output_metapopulations
#'
#'
#' @export

iSTAY_Multiple <- function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL){
  NA_num <- sum(is.na(data))
  if (NA_num != 0) {
    stop("There are some NA in the data.")
  }
  if (sum(which(order.q <= 0))) {
    stop("Order q need larger than 0.")
  }
  if (Alltime != TRUE) {
    if (is.null(start_T) | is.null(end_T)) {
      stop("Need to set the length of time series for calculating.")
    }
    if ((start_T >= end_T) | length(start_T) != 1 | length(end_T) != 
        1) {
      stop("Starting and ending time need to be a number, and ending time needs larger than starting time.")
    }
  }
  Stabillity_Multiple <- function(ZZ, q) {
    K <- ncol(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    if (length(which(z_iplus == 0)) != 0) {
      ZZ <- ZZ[-which(z_iplus == 0), ]
    }
    ZZ <- as.matrix(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    z_plusk <- apply(ZZ, 2, sum)
    z_plusplus <- sum(ZZ)
    p_i <- as.data.frame(apply(ZZ, 2, function(w) w / z_iplus))
    p_pool <- z_plusk / z_plusplus
    ww <- z_iplus / z_plusplus
    if (q == 1) {
      p_i_new <- p_i
      p_pool_new <- p_pool
      alpha <- (exp(-sum(ZZ[ZZ > 0] / z_plusplus * log(p_i_new[p_i_new > 0]))) - 1) / (K - 1)
      gamma <- (exp(-sum(p_pool[p_pool > 0] * log(p_pool_new[p_pool_new > 0]))) - 1) / (K - 1)
    } else {
      alpha <- (sum(z_iplus / z_plusplus * (ZZ / z_iplus)^q)^(1 / (1-q)) - 1) / (K - 1)
      gamma <- (sum(p_pool^q)^(1 / (1 - q)) - 1) / (K - 1)
    }
    return(c(gamma, alpha, (gamma / alpha), (gamma - alpha)))
  }
  Synchrony <- function(ZZ, q) {
    M <- nrow(ZZ)
    K <- ncol(ZZ)
    if (M == 1) {
      value <- 1
    }
    else {
      z_iplus <- apply(ZZ, 1, sum)
      if (length(which(z_iplus == 0)) != 0) {
        ZZ <- ZZ[-which(z_iplus == 0), ]
      }
      ZZ <- as.matrix(ZZ)
      z_iplus <- apply(ZZ, 1, sum)
      if (length(z_iplus) <= 1) {
        value <- NA
      }
      else {
        z_plusk <- apply(ZZ, 2, sum)
        z_plusplus <- sum(ZZ)
        p_i <- apply(ZZ, 2, function(w) w / z_iplus)
        ww <- z_iplus / z_plusplus
        if (q == 1) {
          J <- exp(-sum(ZZ[ZZ > 0] / z_plusplus * log(ZZ[ZZ > 0] / z_plusplus)))
          J <- ifelse(J < K, J, J / M + K * (M - 1) / M)
          pool <- z_plusk/z_plusplus
          pool[which(pool == 0)] <- 10^(-15)
          G <- exp(-sum(pool * log(pool)))
          A <- exp(-sum(ZZ[ZZ > 0] / z_plusplus * log(p_i[p_i > 0])))
          value <- (J - G)/(J - A)
        }
        else {
          J <- sum((ZZ/z_plusplus)^q)^(1/(1-q))
          J <- ifelse(J < K, J, J / M + K * (M - 1) / M)
          G <- sum((z_plusk / z_plusplus)^q) ^ (1 / (1-q))
          A <- sum(apply(ZZ, 2, function(w) (w / z_iplus)^q * ww))^(1 / (1-q))
          value <- (J - G) / (J - A)
        }
      }
    }
    return(value)
  }
  if (is.data.frame(data) | is.matrix(data)) {
    if (Alltime == TRUE) {
      subdata <- data
    }
    else {
      subdata <- data[, c(start_T : end_T)]
    }
    out <- as.matrix(sapply(order.q, function(qq) c(Stabillity_Multiple(subdata, q = qq), Synchrony(subdata, q = qq))))
    result <- data.frame(Site = rep(1, length(order.q)),
                         Order_q = order.q, t(out))
    colnames(result)[3:7] <- c("Stab_Gamma", "Stab_Alpha", 
                               "Stab_Beta_multiplicative", 
                               "Stab_Beta_additive", 
                               "Synchrony")
    result <- result[, -5]
  }
  else {
    out <- lapply(order.q, function(qq) {
      cal <- lapply(data, function(ZZ) {
        if (Alltime == T) {
          subZZ <- ZZ
        }
        else {
          subZZ <- ZZ[, c(start_T:end_T)]
        }
        outout <- c(Stabillity_Multiple(subZZ, q = qq), 
                    Synchrony(subZZ, q = qq))
        result <- data.frame(Order_q = qq, t(outout))
        colnames(result)[2:6] <- c("Stab_Gamma", "Stab_Alpha", 
                                   "Stab_Beta_multiplicative", 
                                   "Stab_Beta_additive", 
                                   "Synchrony")
        result <- result[, -4]
        return(result)
      })
      cal2 <- do.call(rbind, cal)
      calcal <- data.frame(Site = names(data), cal2)
      return(calcal)
    })
    result <- do.call(rbind, out)
  }
  colnames(result) <- c("Dataset", "Order_q", "Gamma", "Alpha", "Beta", "Synchrony")
  return(result)
}

#' Calculate stability and synchrony at each hierarchical level
#'
#' \code{iSTAY_Hier} computes gamma, alpha, and beta stability, as well as synchrony, at each hierarchical level for time series of biomass or other variables.
#'
#' @param data A \code{data.frame} containing the hierarchical data, with sampling units as rows and time points as columns.
#' @param structure The hierarchical structure of the input data.
#' @param order.q A numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime Logical (\code{TRUE} or \code{FALSE}), indicating whether to use all time points in the data.
#' @param start_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the starting column (time point) for the analysis interval.
#' @param end_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the ending column (time point) for the analysis interval.
#'
#' @import dplyr
#'
#' @return a data frame with the following columns: \cr
#'  Hier_level--hierarchical level(e.g. Level 1: population, Level 2: community, Level 3: block, and level 4: overall data)) \cr
#'  Order_q--order of stability or synchrony \cr
#'  Gamma, Alpha, Beta--stability measures of order q \cr
#'  Synchrony--synchrony measure of order q
#'
#' @examples
#'
#' data("Data_Jena_462_populations")
#' data("Data_Jena_hierarchical_structure")
#' output_hier <- iSTAY_Hier(data = Data_Jena_462_populations, structure = Data_Jena_hierarchical_structure,
#'                          order.q=c(1,2), Alltime=TRUE)
#' output_hier
#'
#'
#' @export

iSTAY_Hier <- function (data, structure, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL){
  if (Alltime == FALSE) {
    data <- data[, start_T:end_T]
  }
  TT <- ncol(data)
  del <- which(apply(data, 1, sum) == 0)
  if (length(del) != 0) {
    data <- data[-del, ]
    structure <- structure[-del, ]
  }
  data <- as.matrix(data)
  
  N <- sum(data)
  G <- apply(data, 2, sum) / N
  
  S_gamma <- function(q){
    if(q == 1){
      (exp(-sum(G[G > 0] * log(G[G > 0]))) - 1) / (TT - 1)
    }
    else{
      ((sum(G^q))^(1 / (1 - q)) - 1) / (TT - 1)
    }
  }
  out_gamma <- sapply(order.q, S_gamma)
  out <- data.frame(order_q = order.q, 
                    gamma_value = out_gamma)
  
  S_alpha <- function(q){
    if(q == 1){
      (exp(sum((apply(p, 2, function(x) -sum(x[x>0] * log(x[x>0])))) * w)) - 1) / (TT - 1)
    }
    else{
      (sum((apply(p, 2, function(x) sum((x[x>0])^q))) * w)^(1 / (1 - q)) - 1) / (TT - 1)
    }
  }
  
  nStruc <- ncol(structure) + 1
  
  for (h in 1:(ncol(structure) - 1)) {
    J <- unique(structure[, h])
    NJ <- length(J)
    MJ <- matrix(0, NJ, TT)
    for (j in 1 : NJ) {
      JJ <- which(structure[, h] == J[j])
      if(length(JJ) == 1) {MJ[j, ] <- data[JJ, ]}
      else{MJ[j, ] <- colSums(data[JJ, ])} 
    }
    p <- (apply(MJ, 1, function(x) {
      
      x / max(sum(x), 10^(-15))
      
    }))
    w <- rowSums(MJ) / N
    
    out <- cbind(out, sapply(order.q, S_alpha))
  }
  
  p <- (apply(data, 1, function(x) {
    
    x / max(sum(x), 10^(-15))
    
  }))
  w <- rowSums(data) / N
  
  out <- cbind(out, sapply(order.q, S_alpha))
  
  # Highiest Level Synchrony
  beta <- out[,2] - out[,3]
  
  J <- unique(structure[, 1])
  NJ <- length(J)
  MJ <- matrix(0, NJ, TT)
  for (j in 1 : NJ) {
    JJ <- which(structure[, 1] == J[j])
    if(length(JJ) == 1) {MJ[j, ] <- data[JJ, ]}
    else{MJ[j, ] <- colSums(data[JJ, ])} 
  }
  p <- (apply(MJ, 1, function(x) {
    
    x / max(sum(x), 10^(-15))
    
  }))
  w <- rowSums(MJ) / N
  
  Max_h <- function(q){
    if(q == 1){
      A <- exp(-sum(w * apply(p, 2, function(x) sum(x[x > 0] * log(x[x > 0])))))
      WP <- t(p) * w
      J <- exp(-sum(WP[WP > 0] * log(WP[WP > 0])))
      J <- ifelse(J > TT, 1 / NJ * J + (1 - 1 / NJ) * TT, J)
      (J - A) / (TT - 1)
    }
    else{
      A <- sum(w * apply(p, 2, function(x) sum(x[x>0]^q)))^(1 / (1 - q))
      J <- sum((w^q) * apply(p, 2, function(x) sum(x[x>0]^q)))^(1 / (1 - q))
      J <- ifelse(J > TT, 1 / NJ * J + (1 - 1 / NJ) * TT, J)
      (J - A) / (TT - 1)
    }
  }
  
  beta_max <- sapply(order.q, Max_h)
  
  out <- cbind(out, 1 - beta / beta_max)
  
  for (h in 2 : (ncol(structure))) {
    beta <- out[, (h + 1)] - out[, (h + 2)]
    
    K <- unique(structure[, h - 1])
    NK <- length(K)
    MK <- matrix(0, NK, TT)
    for (k in 1 : NK) {
      KK <- which(structure[, h - 1] == K[k])
      if(length(KK) == 1) {MK[k,] <- data[KK,]}
      else{MK[k,] <- colSums(data[KK,])} 
    }
    
    wk.temp <- rep(0, nrow(data))
    nk.temp <- rep(0, nrow(data))
    for (k in 1:NK) {
      KK <- which(structure[, h - 1] == K[k])
      wk.temp[KK] <- rowSums(MK)[k] / N
      nk.temp[KK] <- rowSums(MK)[k]
    }
    
    J <- unique(structure[, h])
    NJ <- length(J)
    MJ <- matrix(0, NJ, TT)
    wk <- rep(0, NJ)
    nk <- rep(0, NJ)
    for (j in 1 : NJ) {
      JJ <- which(structure[, h] == J[j])
      
      if(length(JJ) == 1) {
        MJ[j,] <- data[JJ,]
        wk[j] <- wk.temp[JJ]
        nk[j] <- nk.temp[JJ]
      }
      else{
        MJ[j,] <- colSums(data[JJ,])
        wk[j] <- (wk.temp[JJ[1]])
        nk[j] <- nk.temp[JJ[1]]
      } 
    }
    p <- (apply(MJ, 1, function(x) {
      
      x / max(sum(x), 10^(-15))
      
    }))
    w <- rowSums(MJ) / N
    nn <- apply(MJ, 1, sum)
    
    Max_d <- function(q){
      if(q == 1){
        A <- exp(-sum((wk * (nn / nk)) * apply(p, 2, function(x) sum(x[x>0] * log(x[x>0])))))
        WP <- t(p) * (nn / nk)
        J <- exp(-sum(wk * apply(WP, 1, function(x) sum(x[x>0] * log(x[x>0])))))
        J <- ifelse(J > TT, 1 / NJ * J + (1 - 1 / NJ) * TT, J)
        (J - A) / (TT - 1)
      }
      else{
        A <- sum((wk * (nn / nk)) * apply(p, 2, function(x) sum(x[x>0]^q)))^(1 / (1 - q))
        J <- sum((wk * (nn / nk)^q) * apply(p, 2, function(x) sum(x[x>0]^q)))^(1 / (1 - q))
        J <- ifelse(J > TT, 1 / NJ * J + (1 - 1 / NJ) * TT, J)
        (J - A) / (TT - 1)
      }
    }
    
    beta_max <- sapply(order.q, Max_d)
    
    out <- cbind(out, (1 - beta / beta_max))
    
    
  }
  
  colnames(out) <- c("Order.q" ,"S_gamma", paste0("S_alpha",(nStruc - 1) : 1), paste0("Synchrony", (nStruc-1):1))
  
  return(data.frame(Hier_level = rep(nStruc : 1, each = length(order.q)),
                    Order_q = rep(order.q, nStruc),
                    Gamma = c(out[, 2], out[, 2 : nStruc] |> unlist()),
                    Alpha = c(rep(NA, length(order.q)), out[, 3 : (nStruc + 1)] |> unlist()),
                    Synchrony = c(rep(NA, length(order.q)), out[, (nStruc + 2) : (nStruc * 2)] |> unlist())) |>
           dplyr::mutate(Beta = Gamma - Alpha, .before = Synchrony))
  
}



#' ggplot2 extension for plotting stability and synchrony profiles.
#'
#' \code{ggiSTAY_qprofile} is a graphical function based on the output from the function \code{iSTAY_Single}, \code{iSTAY_Multiple} or \code{iSTAY_Hier}. It generates stability (and synchrony, if multiple time series are included) profiles that depict how stability and synchrony vary with the order q > 0.
#'
#' @param output the output obtained from \code{iSTAY_Single}, \code{iSTAY_Multiple} or \code{iSTAY_Hier}.
#'
#' @import ggpubr
#'
#'
#' @return For an \code{iSTAY_Single} object, this function return a figure showing the stability profile.\cr
#' For an \code{iSTAY_Multiple} object, it returns a figure displaying the profiles for gamma, alpha, and beta stability, as well as synchrony.\cr
#' For an \code{iSTAY_Hier} object, it returns a figure displaying the profiles for gamma, alpha, and beta stability, as well as synchrony.
#'
#'
#' @examples
#' data("Data_Jena_20_metacommunities")
#' data("Data_Jena_76_metapopulations")
#' data("Data_Jena_462_populations")
#' data("Data_Jena_hierarchical_structure")
#'
#' ## Single time series analysis
#' # Plot the stability profiles of two selected plots
#' individual_plots <- do.call(rbind, Data_Jena_20_metacommunities)
#' output_individual_plots_q <- iSTAY_Single(data = individual_plots[which(rownames(individual_plots) %in% c("B1_4.B1A04", "B4_2.B4A14")),],
#'                                     order.q = seq(0.1,2,0.1), Alltime = TRUE)
#' ggiSTAY_qprofile(output = output_individual_plots_q)
#'
#' # Plot the stability profiles of two selected populations
#' individual_populations <- do.call(rbind, Data_Jena_76_metapopulations)
#' output_individual_populations_q <- iSTAY_Single(data = individual_populations[which(rownames(individual_populations) %in% c("B1A06_B1_16.BM_Ant.odo", "B1A06_B1_16.BM_Cam.pat")),],
#'                                        order.q = seq(0.1,2,0.1), Alltime = TRUE)
#' ggiSTAY_qprofile(output = output_individual_populations_q)
#'
#'
#' ## Multiple time series analysis
#' # Plot the gamma, alpha and beta stability profiles, as well as synchrony profiles of two selected metacommunities
#' metacommunities <- Data_Jena_20_metacommunities
#' output_metacommunities_q <- iSTAY_Multiple(data = metacommunities[which(names(metacommunities) %in% c("B1_1",  "B3_2"))],
#'                                          order.q = seq(0.1,2,0.1), Alltime = TRUE)
#' ggiSTAY_qprofile(output = output_metacommunities_q)
#'
#' # Plot the gamma, alpha and beta stability profiles, as well as synchrony profiles of two selected metapopulations
#' metapopulations <- Data_Jena_76_metapopulations
#' output_metapopulations_q <- iSTAY_Multiple(data = metapopulations[which(names(metapopulations) %in% c("B1A04_B1_4", "B4A14_B4_2"))],
#'                                             order.q = seq(0.1,2,0.1), Alltime = TRUE)
#' ggiSTAY_qprofile(output = output_metapopulations_q)
#'
#'
#' ## Hierarchical time series analysis
#' output_hier_q <- iSTAY_Hier(data = Data_Jena_462_populations,
#'                            structure = Data_Jena_hierarchical_structure,
#'                            order.q = seq(0.1,2,0.1), Alltime = TRUE)
#' ggiSTAY_qprofile(output = output_hier_q)
#'
#' @export

ggiSTAY_qprofile <- function(output){
  
  # local helper: safe color generator (no extra deps)
  ggplotColors <- function(n) {
    grDevices::hcl(h = seq(15, 375, length.out = n + 1), l = 65, c = 100)[1:n]
  }
  
  if (length(which(colnames(output) == "Stability")) != 0) {
    if (length(which(colnames(output) == "Dataset")) == 0 ||
        length(which(colnames(output) == "Order_q")) == 0) {
      stop('Please put the complete output of "iSTAY_Single", "iSTAY_Multiple" or "iSTAY_Hier" function.')
    } else {
      outtype <- "single"
    }
  } else if (length(which(colnames(output) == "Hier_level")) != 0) {
    if (length(which(colnames(output) == "Order_q")) == 0 ||
        length(which(colnames(output) == "Gamma")) == 0 ||
        length(which(colnames(output) == "Alpha")) == 0 ||
        length(which(colnames(output) == "Beta")) == 0 ||
        length(which(colnames(output) == "Synchrony")) == 0) {
      stop('Please put the complete output of "iSTAY_Single", "iSTAY_Multiple" or "iSTAY_Hier" function.')
    } else {
      outtype <- "hier"
    }
  } else {
    if (length(which(colnames(output) == "Dataset")) == 0 ||
        length(which(colnames(output) == "Order_q")) == 0 ||
        length(which(colnames(output) == "Gamma")) == 0 ||
        length(which(colnames(output) == "Alpha")) == 0 ||
        length(which(colnames(output) == "Beta")) == 0 ||
        length(which(colnames(output) == "Synchrony")) == 0) {
      stop('Please put the complete output of "iSTAY_Single", "iSTAY_Multiple" or "iSTAY_Hier" function.')
    } else {
      outtype <- "multiple"
    }
  }
  
  if (outtype == "single") {
    
    if (length(unique(output$`Dataset`)) <= 4) {
      cbPalette <- c("#EA0000","#0066CC","#64A600","#fb8500")
    } else {
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#fb8500"),
                     ggplotColors(length(unique(output$`Dataset`)) - 4))
    }
    
    output$`Dataset` <- factor(output$`Dataset`,
                                      levels = unique(output$`Dataset`))
    
    plotout <- ggplot2::ggplot(
      data = output,
      mapping = ggplot2::aes(x = Order_q, y = Stability, color = `Dataset`)
    ) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::ylab("Stability") +
      ggplot2::xlab("Order of q") +
      ggplot2::labs(color = "Dataset") +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = cbPalette) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 13),
        legend.text  = ggplot2::element_text(size = 12),
        legend.key.size = grid::unit(0.8, "cm"),
        axis.title   = ggplot2::element_text(size = 16)
      )
    
  } else if (outtype == "hier") {
    
    maxhier  <- max(output$Hier_level)
    hier_num <- unique(output$Hier_level)
    qq       <- unique(output$Order_q)
    type_name   <- c(paste0("Gamma(", maxhier, ")"),
                     paste0("Alpha(", hier_num[-1], ")"))
    type_name2_2 <- c(paste0("Beta(", hier_num[-1], ")"), "Alpha(1)")
    type_diff    <- paste0("Hier", hier_num[-1])
    
    plotdat1 <- data.frame(
      Order_q  = rep(qq, length(hier_num)),
      Stability = c(dplyr::filter(output, Hier_level == maxhier)$Gamma,
                    dplyr::filter(output, Hier_level != maxhier)$Alpha),
      type = rep(type_name, each = length(qq))
    )
    plotdat1$type <- factor(plotdat1$type, levels = type_name)
    
    plotdat2 <- data.frame(
      Order_q  = rep(qq, (length(hier_num) - 1)),
      Synchrony = dplyr::filter(output, Hier_level != maxhier)$Synchrony,
      type = rep(type_diff, each = length(qq))
    )
    plotdat2$type <- factor(plotdat2$type, levels = type_diff)
    
    plotdat2_2 <- data.frame(
      Order_q  = rep(qq, length(hier_num)),
      Stability = c(dplyr::filter(output, Hier_level != maxhier)$Beta,
                    dplyr::filter(output, Hier_level == 1)$Alpha),
      type = rep(type_name2_2, each = length(qq))
    )
    plotdat2_2$type <- factor(plotdat2_2$type, levels = type_name2_2)
    
    if (length(unique(plotdat1$type)) <= 4) {
      cbPalette <- c("#EA0000","#0066CC","#64A600","#fb8500")
    } else {
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#fb8500"),
                     ggplotColors(length(unique(plotdat1$type)) - 4))
    }
    
    if (length(unique(plotdat2$type)) <= 4) {
      cbPalette_2 <- c("#EA0000","#0066CC","#64A600","#fb8500")
    } else {
      cbPalette_2 <- c(c("#EA0000","#0066CC","#64A600","#fb8500"),
                       ggplotColors(length(unique(plotdat2$type)) - 4))
    }
    
    plotout1 <- ggplot2::ggplot(plotdat1, ggplot2::aes(x = Order_q, y = Stability, color = type)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::ylab("Hierarchical Stability") +
      ggplot2::labs(color = "") +
      ggplot2::scale_colour_manual(values = cbPalette) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 16),
        plot.margin = grid::unit(c(1,1,1,1), "cm"),
        legend.key.size = grid::unit(0.8, "cm"),
        legend.text = ggplot2::element_text(size = 12)
      )
    
    plotout2 <- ggplot2::ggplot(plotdat2, ggplot2::aes(x = Order_q, y = Synchrony, color = type)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::ylab(expression(paste("Hierarchical Synchrony"))) +
      ggplot2::labs(color = "Hierarchical level") +
      ggplot2::scale_colour_manual(values = cbPalette_2) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 16),
        plot.margin = grid::unit(c(1,1,1,1), "cm"),
        legend.key.size = grid::unit(0.8, "cm"),
        legend.text = ggplot2::element_text(size = 12)
      )
    
    plotout2_2 <- ggplot2::ggplot(plotdat2_2, ggplot2::aes(x = Order_q, y = Stability, color = type)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::ylab(expression(paste("Hierarchical Stability"))) +
      ggplot2::labs(color = "") +
      ggplot2::scale_colour_manual(values = cbPalette) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text  = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 16),
        plot.margin = grid::unit(c(1,1,1,1), "cm"),
        legend.key.size = grid::unit(0.8, "cm"),
        legend.text = ggplot2::element_text(size = 12)
      )
    
    plotout <- list()
    plotout[[1]] <- plotout1
    plotout[[2]] <- ggpubr::ggarrange(plotout2_2, plotout2, ncol = 2)
    
  } else {  # multiple
    
    output$Dataset <- factor(output$Dataset, levels = unique(output$Dataset))
    
    stab_plotdat <- data.frame(
      Dataset    = rep(output$Dataset, 4),
      Order_q = rep(output$Order_q, 4),
      value   = c(output$Gamma, output$Alpha, output$Beta, output$Synchrony),
      type    = rep(c("Gamma", "Alpha", "Beta", "Synchrony"), each = nrow(output))
    )
    stab_plotdat$type <- factor(stab_plotdat$type, levels = c("Gamma","Alpha","Beta","Synchrony"))
    
    if (length(unique(stab_plotdat$Dataset)) <= 4) {
      cbPalette <- c("#EA0000","#0066CC","#64A600","#fb8500")
    } else {
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#fb8500"),
                     ggplotColors(length(unique(stab_plotdat$Dataset)) - 4))
    }
    
    plotout <- ggplot2::ggplot(stab_plotdat, ggplot2::aes(x = Order_q, y = value, color = Dataset)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::facet_wrap(~ type, nrow = 2, scales = "free") +
      ggplot2::ylab("Stability and Synchrony") +
      ggplot2::xlab("Order of q") +
      ggplot2::labs(color = "Dataset") +
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = cbPalette) +
      ggplot2::theme(
        strip.text   = ggplot2::element_text(size = 13),
        legend.title = ggplot2::element_text(size = 13),
        legend.text  = ggplot2::element_text(size = 12),
        legend.key.size = grid::unit(0.8, "cm"),
        axis.title   = ggplot2::element_text(size = 16)
      )
  }
  
  return(plotout)
}



#' ggplot2 extension for plotting diversity–stability and diversity–synchrony relationships.
#'
#' \code{ggiSTAY_analysis} is a graphical function based on the output from the functions \code{iSTAY_Single} or \code{iSTAY_Multiple}. It generates plots showing the relationships between stability (and synchrony if multiple time series are included) and an additional variable, such as diversity or another relevant factor.
#'
#' @param output The output obtained from \code{iSTAY_Single} or \code{iSTAY_Multiple}. It must include (or be combined with) a column corresponding to the variable specified in \code{x_variable}. If \code{by_group} is not \code{NULL},  must also include a column corresponding the variable specified in \code{by_group}.
#' @param x_variable The name of the column representing the diversity (or other) variable to be used as the x-axis in the plot.
#' @param by_group The name of the column representing a categorical variable used to color points by group. The argument is required if \code{model = "LMM"},  as the model uses it as random effect for both intercept and slope. Default is \code{NULL}.
#' @param model Specifies the fitting model. Use \code{model = "lm"} for a linear model; or \code{model = "LMM"} for a linear mixed model with random effects for both intercept and slope. Default is \code{model = "LMM"}.
#'
#'
#' @import stringr
#' @import lme4
#' @import dplyr
#'
#'
#' @return For an \code{iSTAY_Single} object, this function return a figure showing the relationship between the diversity (or other) variable and stability. \cr
#' For an \code{iSTAY_Multiple} object, this function returns a figure showing the diversity (or other) variable and gamma, alpha, and beta stability, as well as synchrony.
#'
#' @examples
#' data("Data_Jena_20_metacommunities")
#' data("Data_Jena_76_metapopulations")
#' data("Data_Jena_462_populations")
#' data("Data_Jena_hierarchical_structure")
#'
#' ## Single time series analysis
#' # Analyze the stability of individual plots and diversity-stability relationship
#' individual_plots <- do.call(rbind, Data_Jena_20_metacommunities)
#' output_individual_plots_div <- iSTAY_Single(data=individual_plots, order.q=c(1,2), Alltime=TRUE)
#' output_individual_plots_div <- data.frame(output_individual_plots_div,
#'                                      log2_sowndiv = log2(as.numeric(do.call(rbind,
#'                                        strsplit(output_individual_plots_div[,1],"[._]+"))[,2])),
#'                                      block = do.call(rbind,
#'                                        strsplit(output_individual_plots_div[,1],"[._]+"))[,1])
#'
#' ggiSTAY_analysis(output = output_individual_plots_div, x_variable = "log2_sowndiv",
#'                     by_group="block", model="LMM")
#'
#' # Analyze the stability of individual populations and diversity-stability relationships
#' individual_populations <- do.call(rbind, Data_Jena_76_metapopulations)
#' output_individual_populations_div <- iSTAY_Single(data = individual_populations,
#'                                          order.q=c(1,2), Alltime=TRUE)
#' output_individual_populations_div <- data.frame(output_individual_populations_div,
#'                               log2_sowndiv = log2(as.numeric(do.call(rbind,
#'                                       strsplit(output_individual_populations_div[,1],"[._]+"))[,3])),
#'                               block=do.call(rbind,
#'                                     strsplit(output_individual_populations_div[,1],"[._]+"))[,2])
#'
#' ggiSTAY_analysis(output = output_individual_populations_div, x_variable = "log2_sowndiv",
#'                     by_group = "block", model = "LMM")
#'
#'
#' ## Multiple time series analysis
#' # Analyze the stability and synchrony within each metacommunity and their relationships with diversity
#' metacommunities <- Data_Jena_20_metacommunities
#' output_metacommunities_div <- iSTAY_Multiple(data = metacommunities, order.q = c(1,2), Alltime = TRUE)
#' output_metacommunities_div <- data.frame(output_metacommunities_div, 
#'                                     log2_sowndiv = log2(as.numeric(do.call(rbind, strsplit(output_metacommunities_div[, 1], "_"))[, 2])),
#'                                     block = do.call(rbind, strsplit(output_metacommunities_div[, 1], "_"))[, 1])
#'
#' ggiSTAY_analysis(output = output_metacommunities_div, x_variable = "log2_sowndiv",
#'                     by_group = "block", model = "LMM")
#'
#' # Analyze the stability and synchrony within each metapopulation and their relationships with diversity
#' metapopulations <- Data_Jena_76_metapopulations
#' output_metapopulations_div <- iSTAY_Multiple(data = metapopulations,
#'                                               order.q=c(1,2), Alltime=TRUE)
#' output_metapopulations_div <- data.frame(output_metapopulations_div,
#'                              log2_sowndiv = log2(as.numeric(do.call(rbind,
#'                              strsplit(output_metapopulations_div[,1],"[._]+"))[,3])),
#'                              block = do.call(rbind, strsplit(output_metapopulations_div[,1],"_"))[,2])
#'
#' ggiSTAY_analysis(output = output_metapopulations_div, x_variable = "log2_sowndiv",
#'                     by_group = "block", model = "LMM")
#'
#' @export

ggiSTAY_analysis <- function(output, x_variable, by_group=NULL, model="LMM"){ 
  
  # local helper: safe color generator
  ggplotColors <- function(n) {
    grDevices::hcl(h = seq(15, 375, length.out = n + 1), l = 65, c = 100)[1:n]
  }
  
  ## --- Basic checks ---------------------------------------------------------
  
  if(length(which(colnames(output)=="Stability"))!=0){ 
    if((length(which(colnames(output)=="Dataset"))==0 ) | 
       length(which(colnames(output)=="Order_q"))==0){ 
      stop('Please put the complete output of "iSTAY_Single" or "iSTAY_Multiple" function.') }
  }
  else{ 
    if(length(which(colnames(output)=="Dataset"))==0 | 
       length(which(colnames(output)=="Order_q"))==0 | 
       length(which(colnames(output)=="Gamma"))==0 | 
       length(which(colnames(output)=="Alpha"))==0 | 
       length(which(colnames(output)=="Beta"))==0 | 
       length(which(colnames(output)=="Synchrony"))==0){ 
      stop('Please put the complete output of "iSTAY_Single" or "iSTAY_Multiple" function.') } 
  } 
  if(length(which(colnames(output)==x_variable))==0){ 
    stop('The output data need to combine a column setting as x_variable.') }
  else{ colnames(output)[which(colnames(output)==x_variable)] <- c("Xvariable") } 
  if(is.null(by_group)==FALSE){ if(length(which(colnames(output)==by_group))==0){ 
    stop('The output data need to combine a column setting as by_group.') }
  else{ colnames(output)[which(colnames(output)==by_group)] <- c("Gvariable") } } 
  if(model=="LMM"){ 
    if(is.null(by_group)==TRUE){ 
      stop('For linear mixed model, you need to set by_group variable as random effect.') } } 
  if(is.null(by_group)==FALSE){ output$Gvariable <- as.factor(output$Gvariable) } 
  

  ## Single ------------------------------------------------------------------


  if(length(which(colnames(output)=="Stability"))!=0){ 
    
    # fit model 
    lm_sign <-c() 
    lm_slope <- c() 
    plotdata <- c() 
    part_fit <- c() 
    plotdata_text_part <- c() 
    
    for (qq in unique(output$Order_q)) {
      subdata <- dplyr::filter(output, Order_q == qq)
      
      if (identical(model, "lm")) {
        MODEL  <- stats::lm(Stability ~ Xvariable, subdata)
        smry   <- summary(MODEL)
        sign   <- smry$coefficients[2, 4]
        pred_value <- stats::predict(MODEL, newdata = subdata)
        
      } else {
        MODEL  <- lmerTest::lmer(Stability ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata)
        smry   <- summary(MODEL)
        sign   <- smry$coefficients[2, 5]
        pred_value <- stats::predict(MODEL, newdata = subdata, re.form = NA)
        
        ran <- c()
        for (bb in rownames(stats::coef(MODEL)$Gvariable)) {
          ran <- rbind(ran, range(dplyr::filter(subdata, Gvariable == bb)$Xvariable))
        }
        slope_intercept <- data.frame(
          Order_q  = rep(paste0("q = ", qq), nrow(stats::coef(MODEL)$Gvariable)),
          Gvariable = rownames(stats::coef(MODEL)$Gvariable),
          intercept = stats::coef(MODEL)$Gvariable[, 1],
          slope     = stats::coef(MODEL)$Gvariable[, 2],
          ran
        )
        colnames(slope_intercept)[5:6] <- c("x_min", "x_max")
        part_fit <- rbind(part_fit, slope_intercept)
        
        plotdata_text_part <- rbind(
          plotdata_text_part,
          data.frame(
            slope   = paste0("slope = ", format(round(stats::coef(MODEL)$Gvariable[, 2], 4), nsmall = 4)),
            Order_q = rep(paste0("q = ", qq), nrow(stats::coef(MODEL)$Gvariable)),
            Gvariable = rownames(stats::coef(MODEL)$Gvariable)
          )
        )
      }
      
      lm_sign  <- rbind(lm_sign,  c(sign = ifelse(sign < 0.05, "significant", "non-significant")))
      lm_slope <- rbind(lm_slope, c(slope = smry$coefficients[2, 1]))
      plotdata <- rbind(plotdata, data.frame(subdata, pred = pred_value))
    }
    
    lm_sign <- as.data.frame(lm_sign) 
    lm_slope <- as.data.frame(lm_slope) 
    rownames(lm_sign) <- paste("q = ", unique(output$Order_q), sep="") 
    rownames(lm_slope) <- paste("q = ", unique(output$Order_q), sep="") 
    
    plotdata$Order_q <- paste("q = ", plotdata$Order_q, sep="") 
    plotdata$sign <- sapply(plotdata$Order_q, function(yy){ 
      sign <- lm_sign[which(rownames(lm_sign)==yy),1] 
      return(sign) 
    }) 
    plotdata$sign <- factor(plotdata$sign, levels=c("significant", "non-significant")) 
    slope_text <- data.frame(slope = paste("slope = ",round(lm_slope[,1],4),sep=""), 
                             Order_q = rownames(lm_slope)) 
    
    if(!any(plotdata$sign == "non-significant")) { 
      dummy <- plotdata[1, ] 
      dummy$Xvariable <- NA 
      dummy$pred <- NA 
      dummy$Stability <- NA 
      dummy$sign <- "non-significant" 
      plotdata <- rbind(plotdata, dummy) 
    } 
    if(!any(plotdata$sign == "significant")) { 
      dummy <- plotdata[1, ] 
      dummy$Xvariable <- NA 
      dummy$pred <- NA 
      dummy$Stability <- NA 
      dummy$sign <- "significant" 
      plotdata <- rbind(plotdata, dummy) 
    } 
    
    if(model=="LMM"){ 
      plotdata_text_part$Order_q <- as.factor(plotdata_text_part$Order_q) 
      plotdata_text_part$Gvariable <- as.factor(plotdata_text_part$Gvariable) 
      part_fit$Order_q <- as.factor(part_fit$Order_q) 
      part_fit$Gvariable <- as.factor(part_fit$Gvariable) 
    } 
    
    
    
    if (!is.null(by_group)) {
      # palette
      if (length(unique(plotdata$Gvariable)) <= 4) {
        cbPalette <- c("#EA0000", "#64A600", "#0066CC", "#fb8500")
      } else {
        cbPalette <- c(c("#EA0000", "#64A600", "#0066CC", "#fb8500"),
                       ggplotColors(length(unique(plotdata$Gvariable)) - 4))
      }
      
      ### 0717 revise 
      if(length(unique(plotdata$Gvariable)) == 4){ 
        hhjust1 <- 1 
        vvjust1 <- -3.5 
        hhjust2 <- rep(c(2.1,1),4) 
        vvjust2 <- rep(c(-1.7,-1.7,0.1,0.1),2) 
      }
      else{ 
        gnum <- length(unique(plotdata$Gvariable)) 
        if(gnum%%2==0){ create <- gnum/2 
        }
        else{ 
          create <- (gnum + 1)/2 
        }
        hhjust1 <- 1 
        vcreate <- 0.1+(-1.8)*c(create:0) 
        vvjust1 <- vcreate[1] 
        if(gnum%%2==1){ 
          hhjust2 <- rep(c(2.1,1),create)[-length(rep(c(2.1,1),create))] 
          hhjust2 <- rep(hhjust2,2) 
          vvjust2 <- rep(vcreate[-1], each=2)[-length(rep(vcreate[-1], each=2))] 
          vvjust2 <- rep(vvjust2, 2) 
        }
        else{ 
          hhjust2 <- rep(rep(c(2.1,1),create),2) 
          vvjust2 <- rep(rep(vcreate[-1], each=2),2) 
        } 
      } 
      

      
      if (identical(model, "LMM")) {
        plotout <- ggplot2::ggplot() +
          ggplot2::geom_point(data = plotdata,
                              ggplot2::aes(x = Xvariable, y = Stability, color = Gvariable), size = 2.7) +
          ggplot2::geom_segment(data = part_fit,
                                ggplot2::aes(x = x_min, xend = x_max,
                                             y = intercept + slope * x_min,
                                             yend = intercept + slope * x_max,
                                             color = Gvariable)) +
          ggplot2::geom_line(data = plotdata,
                             ggplot2::aes(x = Xvariable, y = pred, linetype = sign),
                             linewidth = 1.2, color = "black") +
          ggplot2::scale_linetype_manual(values = c("solid", "dashed"), drop = FALSE) +
          ggplot2::scale_color_manual(values = cbPalette) +
          ggplot2::facet_wrap(~ Order_q, scales = "fixed",
                              ncol = min(5, length(unique(plotdata$Order_q)))) +
          ggplot2::labs(linetype = "", color = by_group) +
          ggplot2::xlab(x_variable) + ggplot2::ylab("Stability") +
          ggplot2::theme_bw() +
          ggplot2::geom_text(data = slope_text, ggplot2::aes(x = -Inf, y = -Inf, label = slope),
                             x = max(plotdata$Xvariable, na.rm = TRUE),
                             y = min(plotdata$Stability, na.rm = TRUE),
                             size = 3.5, hjust=hhjust1, vjust=vvjust1) +
          ggplot2::geom_text(data = plotdata_text_part,
                             ggplot2::aes(x = -Inf, y = -Inf, label = slope, color = Gvariable),
                             x = max(plotdata$Xvariable, na.rm = TRUE),
                             y = min(plotdata$Stability, na.rm = TRUE),
                             size = 3.5, hjust=hhjust2, vjust=vvjust2,
                             key_glyph = ggplot2::draw_key_path) +
          ggplot2::theme(
            strip.text   = ggplot2::element_text(size = 13),
            legend.title = ggplot2::element_text(size = 13),
            legend.text  = ggplot2::element_text(size = 12),
            legend.key.size = grid::unit(0.8, "cm"),
            axis.title   = ggplot2::element_text(size = 16)
          ) +
          ggplot2::guides(linetype = ggplot2::guide_legend(keywidth = 2.5))
      } else {
        
        plotout <- ggplot2::ggplot() +
          ggplot2::geom_point(data = plotdata,
                              ggplot2::aes(x = Xvariable, y = Stability, color = Gvariable), size = 2.7) +
          ggplot2::geom_line(data = plotdata,
                             ggplot2::aes(x = Xvariable, y = pred, linetype = sign),
                             linewidth = 1.2, color = "black") +
          ggplot2::scale_linetype_manual(values = c("solid", "dashed"), drop = FALSE) +
          ggplot2::scale_color_manual(values = cbPalette) +
          ggplot2::facet_wrap(~ Order_q, scales = "fixed",
                              ncol = min(5, length(unique(plotdata$Order_q)))) +
          ggplot2::labs(linetype = "", color = by_group) +
          ggplot2::xlab(x_variable) + ggplot2::ylab("Stability") +
          ggplot2::theme_bw() +
          ggplot2::geom_text(data = slope_text, ggplot2::aes(x = -Inf, y = -Inf, label = slope),
                             x = max(plotdata$Xvariable, na.rm = TRUE),
                             y = min(plotdata$Stability, na.rm = TRUE),
                             size = 5, hjust = 1, vjust = 0.1) +
          ggplot2::theme(
            strip.text   = ggplot2::element_text(size = 13),
            legend.title = ggplot2::element_text(size = 13),
            legend.text  = ggplot2::element_text(size = 12),
            legend.key.size = grid::unit(0.8, "cm"),
            axis.title   = ggplot2::element_text(size = 16)
          ) +
          ggplot2::guides(linetype = ggplot2::guide_legend(keywidth = 2.5))
      }
      
    } else {  # no by_group
      plotout <- ggplot2::ggplot(plotdata, ggplot2::aes(x = Xvariable, y = Stability)) +
        ggplot2::geom_point(size = 2.7) +
        ggplot2::geom_line(ggplot2::aes(x = Xvariable, y = pred, linetype = sign),
                           linewidth = 1.2, color = "black") +
        ggplot2::scale_linetype_manual(values = c("solid", "dashed"), drop = FALSE) +
        ggplot2::facet_wrap(~ Order_q, scales = "fixed",
                            ncol = min(5, length(unique(plotdata$Order_q)))) +
        ggplot2::labs(linetype = "") +
        ggplot2::xlab(x_variable) + ggplot2::ylab("Stability") +
        ggplot2::theme_bw() +
        ggplot2::geom_text(data = slope_text, ggplot2::aes(x = -Inf, y = -Inf, label = slope),
                           x = max(plotdata$Xvariable, na.rm = TRUE),
                           y = min(plotdata$Stability, na.rm = TRUE),
                           size = 5, hjust = 1, vjust = 0.1) +
        ggplot2::theme(
          strip.text   = ggplot2::element_text(size = 13),
          legend.title = ggplot2::element_text(size = 13),
          legend.text  = ggplot2::element_text(size = 12),
          legend.key.size = grid::unit(0.8, "cm"),
          axis.title   = ggplot2::element_text(size = 16)
        ) +
        ggplot2::guides(linetype = ggplot2::guide_legend(keywidth = 2.5))
    } 
  }
  else{ 
    
    # fit model 
    lm_sign <-c() 
    lm_slope <- c() 
    plotdata <- c() 
    plotdata_text_part <- c() 
    part_fit <- c() 
    
    for (qq in unique(output$Order_q)) {
      subdata <- dplyr::filter(output, Order_q == qq)
      
      if (identical(model, "lm")) {
        sign_num <- 4
        MODEL_G <- stats::lm(Gamma ~ Xvariable, subdata);     sum_gamma <- summary(MODEL_G); pred_G <- stats::predict(MODEL_G, newdata = subdata)
        MODEL_A <- stats::lm(Alpha ~ Xvariable, subdata);     sum_alpha <- summary(MODEL_A); pred_A <- stats::predict(MODEL_A, newdata = subdata)
        MODEL_BA<- stats::lm(Beta  ~ Xvariable, subdata);     sum_beta_add <- summary(MODEL_BA); pred_BA <- stats::predict(MODEL_BA, newdata = subdata)
        MODEL_S <- stats::lm(Synchrony ~ Xvariable, subdata); sum_syn <- summary(MODEL_S); pred_S <- stats::predict(MODEL_S, newdata = subdata)
        
      } else {
        sign_num <- 5
        MODEL_G <- lmerTest::lmer(Gamma     ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata); sum_gamma <- summary(MODEL_G); pred_G <- stats::predict(MODEL_G, newdata = subdata, re.form = NA)
        MODEL_A <- lmerTest::lmer(Alpha     ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata); sum_alpha <- summary(MODEL_A); pred_A <- stats::predict(MODEL_A, newdata = subdata, re.form = NA)
        MODEL_BA<- lmerTest::lmer(Beta      ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata); sum_beta_add <- summary(MODEL_BA); pred_BA <- stats::predict(MODEL_BA, newdata = subdata, re.form = NA)
        MODEL_S <- lmerTest::lmer(Synchrony ~ 1 + Xvariable + (1 + Xvariable | Gvariable), subdata); sum_syn <- summary(MODEL_S); pred_S <- stats::predict(MODEL_S, newdata = subdata, re.form = NA)
        
        ran <- c()
        for (bb in rownames(stats::coef(MODEL_G)$Gvariable)) {
          ran <- rbind(ran, range(dplyr::filter(subdata, Gvariable == bb)$Xvariable))
        }
        slope_intercept <- data.frame(
          Order_q  = rep(paste0("q = ", qq), nrow(stats::coef(MODEL_G)$Gvariable) * 4),
          type     = rep(c("Gamma", "Alpha", "Beta", "Synchrony"), each = nrow(stats::coef(MODEL_G)$Gvariable)),
          Gvariable= rep(rownames(stats::coef(MODEL_G)$Gvariable), 4),
          intercept= c(stats::coef(MODEL_G)$Gvariable[,1],
                       stats::coef(MODEL_A)$Gvariable[,1],
                       stats::coef(MODEL_BA)$Gvariable[,1],
                       stats::coef(MODEL_S)$Gvariable[,1]),
          slope    = c(stats::coef(MODEL_G)$Gvariable[,2],
                       stats::coef(MODEL_A)$Gvariable[,2],
                       stats::coef(MODEL_BA)$Gvariable[,2],
                       stats::coef(MODEL_S)$Gvariable[,2]),
          rbind(ran, ran, ran, ran)
        )
        colnames(slope_intercept)[6:7] <- c("x_min", "x_max")
        part_fit <- rbind(part_fit, slope_intercept)
        
        plotdata_text_part <- rbind(
          plotdata_text_part,
          data.frame(
            slope_gamma     = paste0("slope = ", format(round(stats::coef(MODEL_G)$Gvariable[,2], 4), nsmall = 4)),
            slope_alpha     = paste0("slope = ", format(round(stats::coef(MODEL_A)$Gvariable[,2], 4), nsmall = 4)),
            slope_beta      = paste0("slope = ", format(round(stats::coef(MODEL_BA)$Gvariable[,2], 4), nsmall = 4)),
            slope_synchrony = paste0("slope = ", format(round(stats::coef(MODEL_S)$Gvariable[,2], 4), nsmall = 4)),
            Order_q = rep(paste0("q = ", qq), nrow(stats::coef(MODEL_G)$Gvariable)),
            Gvariable = rownames(stats::coef(MODEL_G)$Gvariable)
          )
        )
      }
      
      lm_sign <- rbind(
        lm_sign,
        c(
          gamma     = ifelse(sum_gamma$coefficients[2, sign_num] < 0.05, "significant", "non-significant"),
          alpha     = ifelse(sum_alpha$coefficients[2, sign_num] < 0.05, "significant", "non-significant"),
          beta_add  = ifelse(sum_beta_add$coefficients[2, sign_num] < 0.05, "significant", "non-significant"),
          synchrony = ifelse(sum_syn$coefficients[2, sign_num] < 0.05, "significant", "non-significant")
        )
      )
      lm_slope <- rbind(
        lm_slope,
        c(gamma     = sum_gamma$coefficients[2,1],
          alpha     = sum_alpha$coefficients[2,1],
          beta_add  = sum_beta_add$coefficients[2,1],
          synchrony = sum_syn$coefficients[2,1])
      )
      plotdata <- rbind(plotdata,
                        data.frame(subdata, 
                                   pred_G = pred_G, 
                                   pred_A = pred_A, 
                                   pred_BA = pred_BA, 
                                   pred_S = pred_S)
                        )
    }
    
    
    lm_sign <- as.data.frame(lm_sign) 
    lm_slope <- as.data.frame(lm_slope) 
    rownames(lm_sign) <- paste("q = ", unique(output$Order_q), sep = "") 
    rownames(lm_slope) <- paste("q = ", unique(output$Order_q), sep = "") 
    
    plotdata$Order_q <- paste("q = ", plotdata$Order_q, sep="") 
    plotdata$sign_G <- sapply(plotdata$Order_q, function(yy){ 
      lm_sign[which(rownames(lm_sign)==yy),1] }) 
    plotdata$sign_A <- sapply(plotdata$Order_q, function(yy){ 
      lm_sign[which(rownames(lm_sign)==yy),2] 
    }) 
    plotdata$sign_BA <- sapply(plotdata$Order_q, function(yy){ 
      lm_sign[which(rownames(lm_sign)==yy),3] 
    }) 
    plotdata$sign_S <- sapply(plotdata$Order_q, function(yy){ 
      lm_sign[which(rownames(lm_sign)==yy),4] 
    }) 
    
    plotdata_Stab <- data.frame(Dataset = rep(plotdata$Dataset, 4), 
                                Order_q = rep(plotdata$Order_q, 4), 
                                Stability = c(plotdata$Gamma, 
                                              plotdata$Alpha, 
                                              plotdata$Beta, 
                                              plotdata$Synchrony), 
                                pred = c(plotdata$pred_G, plotdata$pred_A, 
                                         plotdata$pred_BA, plotdata$pred_S), 
                                sign = c(plotdata$sign_G, plotdata$sign_A, 
                                         plotdata$sign_BA, plotdata$sign_S), 
                                type = rep(c("Gamma","Alpha", "Beta","Synchrony"), 
                                           each = nrow(plotdata)), 
                                Xvariable = rep(plotdata$Xvariable, 4)
    ) 
    
    if(is.null(by_group)==FALSE){ 
      
      plotdata_Stab$Gvariable <- rep(plotdata$Gvariable,4) 
      plotdata_Stab$Gvariable <- as.factor(plotdata_Stab$Gvariable) 
      plotdata$Gvariable <- as.factor(plotdata$Gvariable) 
      
    } 
    
    plotdata_Stab$sign <- factor(plotdata_Stab$sign, 
                                 levels=c("significant", "non-significant")
    ) 
    plotdata_Stab$type <- factor(plotdata_Stab$type, 
                                 levels = c("Gamma","Alpha","Beta","Synchrony")
    ) 
    plotdata_Stab$Order_q <- as.factor(plotdata_Stab$Order_q) 
    
    if(!any(plotdata_Stab$sign == "non-significant")) { 
      dummy <- plotdata_Stab[1, ] 
      dummy$Xvariable <- NA 
      dummy$pred <- NA 
      dummy$Stability <- NA 
      dummy$sign <- "non-significant" 
      plotdata_Stab <- rbind(plotdata_Stab, dummy) 
    } 
    if(!any(plotdata_Stab$sign == "significant")) { 
      dummy <- plotdata_Stab[1, ] 
      dummy$Xvariable <- NA 
      dummy$pred <- NA 
      dummy$Stability <- NA 
      dummy$sign <- "significant" 
      plotdata_Stab <- rbind(plotdata_Stab, dummy) 
    } 
    
    slope_text_Stab <- data.frame(slope = paste("slope = ",format(round(as.vector(as.matrix(lm_slope)),4), 
                                                                  nsmall=4),sep=""), 
                                  Order_q = rep(rownames(lm_slope),4), 
                                  type = rep(c("Gamma","Alpha","Beta","Synchrony"), 
                                             each = nrow(lm_slope))
    ) 
    slope_text_Stab$type <- factor(slope_text_Stab$type, 
                                   levels = c("Gamma","Alpha","Beta","Synchrony")
    ) 
    slope_text_Stab$Order_q <- as.factor(slope_text_Stab$Order_q) 
    
    if(model=="LMM"){ 
      plotdata_text_part <- data.frame(slope = c(plotdata_text_part$slope_gamma, 
                                                 plotdata_text_part$slope_alpha, 
                                                 plotdata_text_part$slope_beta, 
                                                 plotdata_text_part$slope_synchrony), 
                                       Order_q = rep(plotdata_text_part$Order_q, 4), 
                                       type = rep(c("Gamma","Alpha","Beta","Synchrony"), 
                                                  each=nrow(plotdata_text_part)
                                       ), 
                                       Gvariable = rep(plotdata_text_part$Gvariable, 4)
                                       )
      plotdata_text_part$type <- as.factor(plotdata_text_part$type) 
      plotdata_text_part$Order_q <- as.factor(plotdata_text_part$Order_q) 
      plotdata_text_part$Gvariable <- as.factor(plotdata_text_part$Gvariable) 
      part_fit$Gvariable <- as.factor(part_fit$Gvariable) 
      part_fit$Order_q <- as.factor(part_fit$Order_q) 
      part_fit$type <- as.factor(part_fit$type) 
    } 
     
    tyG <- min(dplyr::filter(plotdata_Stab, type=="Gamma")$Stability, na.rm = TRUE) 
    tyA <- min(dplyr::filter(plotdata_Stab, type=="Alpha")$Stability, na.rm = TRUE) 
    tyBA <- max(dplyr::filter(plotdata_Stab, type=="Beta")$Stability, na.rm = TRUE) 
    tyS <- min(dplyr::filter(plotdata_Stab, type=="Synchrony")$Stability, na.rm = TRUE) 
    
    if(is.null(by_group)==FALSE){ 
    # Check if the number of unique 'Assemblage' is 4 or less 
      if (length(unique(plotdata_Stab$Gvariable)) <= 4){ 
        cbPalette <- c("#EA0000","#64A600","#0066CC","#fb8500") 
      }
      else{ 
        # If there are more than 4 assemblages, start with the same predefined color palette 
        # Then extend the palette by generating additional colors using the 'ggplotColors' function 
        cbPalette <- c(c("#EA0000","#64A600","#0066CC","#fb8500"), 
                       ggplotColors(length(unique(plotdata_Stab$Gvariable))-4)) 
      }
      if(length(unique(plotdata$Gvariable)) == 4){ 
        gnum <- 4 
        hhjust1 <- 1 
        vvjust1 <- rep(c(-3.5,-3.5,4.1,-3.5), each=2) 
        hhjust2 <- rep(c(2.1,1),16) 
        vvjust2 <- c(rep(c(-1.7,-1.7,0.1,0.1),4), 
                     rep(c(1.1,1.1,2.6,2.6),2), 
                     rep(c(-1.7,-1.7,0.1,0.1),2)) 
      }
      else{ 
        gnum <- length(unique(plotdata$Gvariable)) 
        if(gnum%%2==0){ 
          create <- gnum/2 
        }
        else{ 
          create <- (gnum + 1)/2 
        }
        hhjust1 <- 1 
        vcreate_1 <- 0.1+(-1.8)*c(create:0) 
        vcreate_2 <- 1.1+1.5*c(create:0) 
        vvjust1 <- rep(c(vcreate_1[1], vcreate_1[1], vcreate_2[1], vcreate_1[1]), each=2) 
        if(gnum%%2==1){ 
          hhjust2 <- rep(c(2.1,1),create)[-length(rep(c(2.1,1),create))] 
          hhjust2 <- rep(hhjust2,8) 
          vvjust2_1 <- rep(vcreate_1[-1], each=2)[-length(rep(vcreate_1[-1], each=2))] 
          vvjust2_2 <- rev(rep(vcreate_2[-1], each=2))[-length(rev(rep(vcreate_2[-1], each=2)))] 
          vvjust2 <- c(rep(vvjust2_1, 4), rep(vvjust2_2, 2), rep(vvjust2_1, 2)) 
        }
        else{ 
          hhjust2 <- rep(rep(c(2.1,1),create),8) 
          vvjust2 <- c(rep(rep(vcreate_1[-1], each=2),4), 
                       rep(rep(rev(vcreate_2[-1]), each=2),2), 
                       rep(rep(vcreate_1[-1], each=2),2)) 
        } 
      } 
      if(model=="LMM"){ 
        plotout1 <- ggplot2::ggplot() + 
          ggplot2::geom_point(data = plotdata_Stab, 
                              ggplot2::aes(x=Xvariable, y=Stability, color=Gvariable), size=2.7) + 
          ggplot2::geom_segment(data = part_fit, 
                                ggplot2::aes(x=x_min, xend=x_max, 
                           y=intercept+slope*x_min, yend=intercept+slope*x_max, 
                           color=Gvariable)) + 
          ggplot2::geom_line(data = plotdata_Stab, 
                             ggplot2::aes(x=Xvariable, y=pred, linetype = sign), 
                    linewidth=1.2, color="black") + 
          ggplot2::scale_linetype_manual(values=c("solid","dashed"), drop = FALSE) + 
          ggplot2::scale_color_manual(values = cbPalette) + 
          ggplot2::facet_grid(type~Order_q, scales = "free_y") + 
          ggplot2::labs(linetype="", color=by_group) + 
          ggplot2::xlab(label=x_variable) + 
          ggplot2::ylab(label="Stability and Synchrony") + 
          ggplot2::theme_bw() + 
          ggplot2::geom_text(data=slope_text_Stab, 
                             ggplot2::aes(x = -Inf, y = -Inf, label = slope), 
                    x=max(plotdata_Stab$Xvariable, na.rm = TRUE), 
                    y=rep(c(tyG, tyA, tyBA, tyS),each=2), 
                    size=3.5, hjust=hhjust1, vjust=vvjust1
          )+ 
          ggplot2::geom_text(data=plotdata_text_part, 
                             ggplot2::aes(x = -Inf, y = -Inf, label=slope, color=Gvariable), 
                    x=rep(max(plotdata_Stab$Xvariable, na.rm = TRUE),gnum*8), 
                    y=rep(c(tyG, tyA, tyBA, tyS),each=gnum*2), 
                    size=3.5, hjust=hhjust2, vjust=vvjust2, 
                    key_glyph = ggplot2::draw_key_path) + 
          ggplot2::theme(strip.text = ggplot2::element_text(size=13), 
                legend.title = ggplot2::element_text(size=13), 
                legend.text = ggplot2::element_text(size=12), 
                legend.key.size = grid::unit(0.8, 'cm'), 
                axis.title = ggplot2::element_text(size=16)
          ) + 
          ggplot2::guides(linetype = ggplot2::guide_legend(keywidth = 2.5)) 
      }
      else{ 
        plotout1 <- ggplot2::ggplot() + 
          ggplot2::geom_point(data = plotdata_Stab, 
                              ggplot2::aes(x=Xvariable, y=Stability, color=Gvariable), size=2.7) + 
          ggplot2::geom_line(data = plotdata_Stab, 
                             ggplot2::aes(x=Xvariable, y=pred, linetype = sign), 
                    linewidth=1.2, color="black") + 
          ggplot2::scale_linetype_manual(values=c("solid","dashed"), drop = FALSE) + 
          ggplot2::scale_color_manual(values = cbPalette) + 
          ggplot2::facet_grid(type~Order_q, scales = "free_y") + 
          ggplot2::labs(linetype="", color=by_group) + 
          ggplot2::xlab(label=x_variable) + 
          ggplot2::ylab(label="Stability and Synchrony") + 
          ggplot2::theme_bw() + 
          ggplot2::geom_text(data=slope_text_Stab, 
                             ggplot2::aes(x = -Inf, y = -Inf, label = slope), 
                    x=max(plotdata_Stab$Xvariable, na.rm = TRUE), 
                    y=rep(c(tyG, tyA, tyBA, tyS),each=2), 
                    size=5, hjust=1, vjust=rep(c(0.1,0.1,1.1,0.1),each=2)) + 
          ggplot2::theme(strip.text = ggplot2::element_text(size=13), 
                legend.title = ggplot2::element_text(size=13), 
                legend.text = ggplot2::element_text(size=12), 
                legend.key.size = grid::unit(0.8, 'cm'), 
                axis.title = ggplot2::element_text(size=16)
          ) + 
          ggplot2::guides(linetype = ggplot2::guide_legend(keywidth = 2.5)) 
      } 
      
    }
    else{ 
      plotout1 <- ggplot2::ggplot(plotdata_Stab, 
                                  ggplot2::aes(x=Xvariable, y=Stability)) + 
        ggplot2::geom_point(size=2.7) + 
        ggplot2::geom_line(ggplot2::aes(x=Xvariable, y=pred, linetype = sign), 
                  linewidth=1.2, color="black") + 
        ggplot2::scale_linetype_manual(values=c("solid","dashed"), drop = FALSE) + 
        ggplot2::facet_grid(type~Order_q, scales = "free_y") + 
        ggplot2::labs(linetype="") + 
        ggplot2::xlab(label = x_variable) + 
        ggplot2::ylab(label = "Stability and Synchrony") + 
        ggplot2::theme_bw() + 
        ggplot2::geom_text(data=slope_text_Stab, 
                           ggplot2::aes(x = -Inf, y = -Inf, label = slope), 
                  x = max(plotdata_Stab$Xvariable, na.rm = TRUE), 
                  y = rep(c(tyG, tyA, tyBA, tyS),each=2), 
                  size=5, hjust=1, vjust=rep(c(0.1,0.1,1.1,0.1),each=2)) + 
        ggplot2::theme(strip.text = ggplot2::element_text(size=13), 
              legend.title = ggplot2::element_text(size=13), 
              legend.text = ggplot2::element_text(size=12), 
              legend.key.size = grid::unit(0.8, 'cm'), 
              axis.title = ggplot2::element_text(size=16)) + 
        ggplot2::guides(linetype = ggplot2::guide_legend(keywidth = 2.5)) 
    }
    
    plotout <- plotout1 
    
  } 
  return(plotout) 
}




# Generate Color Palette for ggplot2
#
# This function creates a color palette suitable for ggplot2 visualizations by evenly spacing colors in the HCL color space. The function ensures that the colors are well-distributed and visually distinct, making it ideal for categorical data where each category needs to be represented by a different color.
#
# @param g An integer indicating the number of distinct colors to generate. This value should be a positive integer, with higher values resulting in a broader range of colors.
# @return A vector of color codes in hexadecimal format, suitable for use in ggplot2 charts and plots. The length of the vector will match the input parameter `g`.
# @examples
# # Generate a palette of 5 distinct colors
# ggplotColors(5)
#
# # Use the generated colors in a ggplot2 chart
# library(ggplot2)
# df <- data.frame(x = 1:5, y = rnorm(5), group = factor(1:5))
# ggplot(df, aes(x, y, color = group)) +
#   geom_point() +
#   scale_color_manual(values = ggplotColors(5))
#
ggplotColors <- function(g){
  d <- 360/g # Calculate the distance between colors in HCL color space
  h <- cumsum(c(15, rep(d,g - 1))) # Create cumulative sums to define hue values
  grDevices::hcl(h = h, c = 100, l = 65) # Convert HCL values to hexadecimal color codes
}





## ========== no visible global function definition for R CMD check ========== ##
utils::globalVariables(c("COLOR", "Order_q", "Dataset", "Stability", "Synchrony"
))





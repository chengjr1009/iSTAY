
utils::globalVariables(c(
  "Alpha", "Gamma", "Gvariable", "Hier", "Plot/Community", "Xvariable",
  "coef", "intercept", "lm", "pred", "slope",
  "type", "value", "x_max", "x_min"
))

#' Calculate stability of time-series data for a single assemblage.
#'
#' \code{Stay_Single} is a function that computes the stability of time-series data (e.g., biomass, productivity) for a single assemblage.
#'
#' @param data can be input as a \code{vector} of time series data, or \code{data.frame} (assemblages by times).
#' @param order.q a numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in the data.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in the data.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in the data.
#'
#'
#'
#' @return a dataframe with columns "Plot/Community", "Order_q" and "Stability".
#'
#' @examples
#' # Stability of each single plot
#' data("Jena_plot_biomass_data")
#' single_plot <- do.call(rbind, Jena_plot_biomass_data)
#' output_single_plot <- Stay_Single(data=single_plot, order.q=c(1,2), Alltime=TRUE)
#' output_single_plot
#'
#' # Stability of each single species in each plot
#' data("Jena_species_biomass_data")
#' single_species <- do.call(rbind, Jena_species_biomass_data)
#' output_single_species <- Stay_Single(data=single_species, order.q=c(1,2), Alltime=TRUE)
#' output_single_species
#'
#' @export

Stay_Single <- function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL){
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
  colnames(result)[1] <- c("Plot/Community")
  return(result)
}





#' Calculate stability and synchrony of time-series data for multiple assemblages.
#'
#' \code{Stay_Multiple} is a function that computes gamma, alpha, and beta stability, as well as synchrony, for time-series data (e.g., biomass, productivity) across multiple assemblages.
#'
#' @param data can be input as a \code{data.frame/matrix} (assemblages by times), or a \code{list} of \code{data.frames} with each dataframe representing a assemblages-by-times data.
#' @param order.q a numerical vector specifying the orders of stability and synchrony. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in (every) dataframe.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in (every) dataframe.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in (every) dataframe.
#'
#'
#' @return a dataframe with columns "Site", "Order_q", "Gamma", "Alpha", "Beta" and "Synchrony".
#'
#' @examples
#' # Stability of multiple plots
#' data("Jena_plot_biomass_data")
#' multiple_plot <- Jena_plot_biomass_data
#' output_multi_plot <- Stay_Multiple(data=multiple_plot, order.q=c(1,2), Alltime=TRUE)
#' output_multi_plot
#'
#' # Stability of multiple species in each plot
#' data("Jena_species_biomass_data")
#' multiple_species <- Jena_species_biomass_data
#' output_multi_species <- Stay_Multiple(data=multiple_species, order.q=c(1,2), Alltime=TRUE)
#' output_multi_species
#'
#'
#' @export

Stay_Multiple <- function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL){
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
  colnames(result) <- c("Site", "Order_q", "Gamma", "Alpha", "Beta", "Synchrony")
  return(result)
}

#' Calculate stability of time-series data for hierarchical structure.
#'
#' \code{Stay_Hier} is a function that computes gamma, alpha, and normalized beta stability for time-series data (e.g., biomass, productivity) within hierarchical structures.
#'
#' @param data can be input as \code{data.frame} (assemblages by times).
#' @param mat hierarchical structure of data.
#' @param order.q a numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime \code{TRUE} or \code{FALSE}, to decide whether to use all the times in the data.
#' @param start_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the start column of time in the data.
#' @param end_T (argument only for \code{Alltime = FALSE}) a positive integer specifying the end column of time in the data.
#'
#' @import dplyr
#'
#' @return a dataframe with columns "Hier", "Order_q", and stability "Gamma", "Alpha", "Beta" and "Synchrony".
#'
#' @examples
#'
#' data("Jena_hierarchical_biomass_data")
#' data("Jena_hierarchical_mat")
#' output_hier <- Stay_Hier(data=Jena_hierarchical_biomass_data, mat=Jena_hierarchical_mat,
#'                          order.q=c(1,2), Alltime=TRUE)
#' output_hier
#'
#'
#' @export

Stay_Hier <- function (data, mat, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL){
  if (Alltime == FALSE) {
    data <- data[, start_T:end_T]
  }
  TT <- ncol(data)
  del <- which(apply(data, 1, sum) == 0)
  if (length(del) != 0) {
    data <- data[-del, ]
    mat <- mat[-del, ]
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
  
  nStruc <- ncol(mat) + 1
  
  for (h in 1:(ncol(mat) - 1)) {
    J <- unique(mat[, h])
    NJ <- length(J)
    MJ <- matrix(0, NJ, TT)
    for (j in 1 : NJ) {
      JJ <- which(mat[, h] == J[j])
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
  
  J <- unique(mat[, 1])
  NJ <- length(J)
  MJ <- matrix(0, NJ, TT)
  for (j in 1 : NJ) {
    JJ <- which(mat[, 1] == J[j])
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
  
  for (h in 2 : (ncol(mat))) {
    beta <- out[, (h + 1)] - out[, (h + 2)]
    
    K <- unique(mat[, h - 1])
    NK <- length(K)
    MK <- matrix(0, NK, TT)
    for (k in 1 : NK) {
      KK <- which(mat[, h - 1] == K[k])
      if(length(KK) == 1) {MK[k,] <- data[KK,]}
      else{MK[k,] <- colSums(data[KK,])} 
    }
    
    wk.temp <- rep(0, nrow(data))
    nk.temp <- rep(0, nrow(data))
    for (k in 1:NK) {
      KK <- which(mat[, h - 1] == K[k])
      wk.temp[KK] <- rowSums(MK)[k] / N
      nk.temp[KK] <- rowSums(MK)[k]
    }
    
    J <- unique(mat[, h])
    NJ <- length(J)
    MJ <- matrix(0, NJ, TT)
    wk <- rep(0, NJ)
    nk <- rep(0, NJ)
    for (j in 1 : NJ) {
      JJ <- which(mat[, h] == J[j])
      
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
  
  return(data.frame(Hier = rep(nStruc : 1, each = length(order.q)),
                    Order_q = rep(order.q, nStruc),
                    Gamma = c(out[, 2], out[, 2 : nStruc] |> unlist()),
                    Alpha = c(rep(NA, length(order.q)), out[, 3 : (nStruc + 1)] |> unlist()),
                    Synchrony = c(rep(NA, length(order.q)), out[, (nStruc + 2) : (nStruc * 2)] |> unlist())) |>
           mutate(Beta = Gamma - Alpha, .before = Synchrony))
  
}



#' ggplot2 extension for a Stay_Single, Stay_Multiple or Stay_Hier object with q-profile.
#'
#' \code{ggStay_qprofile} is a graphical function that based on the output from the function \code{Stay_Single}, \code{Stay_Multiple} or \code{Stay_Hier}. It provides to graph the q-profile of stability (and synchrony if is multiple assemblages).
#'
#' @param output the output obtained from \code{Stay_Single}, \code{Stay_Multiple} or \code{Stay_Hier}.
#'
#' @import ggpubr
#'
#'
#' @return For a \code{Stay_Single} object, this function return a figure of q-profile for stability .
#' For a \code{Stay_Multiple} object, this function return a figure that contains q-profile for (Gamma, Alpha, Beta) stability and synchrony.
#' For a \code{Stay_Hier} object, this function return a figure that contains q-profile for gamma stability of highest hierarchical level and alpha stability of other hierarchical level. And it also provides a figure about the relationship between each decomposition of overall stability and the synchrony of each hierarchical level and order.q.
#'
#'
#' @examples
#' data("Jena_plot_biomass_data")
#' data("Jena_species_biomass_data")
#' data("Jena_hierarchical_biomass_data")
#' data("Jena_hierarchical_mat")
#'
#' ## Single assemblage
#' # Stability of each single plot
#' single_plot <- do.call(rbind, Jena_plot_biomass_data)
#' output_single_plot_q <- Stay_Single(data=single_plot[c(12,38),],
#'                                     order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_single_plot_q)
#'
#' # Stability of each single species
#' single_species <- do.call(rbind, Jena_species_biomass_data)
#' output_single_species_q <- Stay_Single(data=single_species[c(40,49),],
#'                                        order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_single_species_q)
#'
#'
#' ## Multiple assemblages
#' # Stability of multiple plots
#' multiple_plot <- Jena_plot_biomass_data
#' output_multi_plot_q <- Stay_Multiple(data=multiple_plot[c(9,11)],
#'                                          order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_multi_plot_q)
#'
#' # Stability of multiple species in plot
#' multiple_species <- Jena_species_biomass_data
#' output_multi_species_q <- Stay_Multiple(data=multiple_species[c(62,70)],
#'                                             order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_multi_species_q)
#'
#'
#' ## Hierarchies
#' output_hier_q <- Stay_Hier(data=Jena_hierarchical_biomass_data,
#'                            mat=Jena_hierarchical_mat,
#'                            order.q=seq(0.1,2,0.1), Alltime=TRUE)
#' ggStay_qprofile(output=output_hier_q)
#'
#' @export

ggStay_qprofile <- function(output){
  
  # local helper: safe color generator (no extra deps)
  ggplotColors <- function(n) {
    grDevices::hcl(h = seq(15, 375, length.out = n + 1), l = 65, c = 100)[1:n]
  }
  
  if (length(which(colnames(output) == "Stability")) != 0) {
    if (length(which(colnames(output) == "Plot/Community")) == 0 ||
        length(which(colnames(output) == "Order_q")) == 0) {
      stop('Please put the complete output of "Stay_Single", "Stay_Multiple" or "Stay_Hier" function.')
    } else {
      outtype <- "single"
    }
  } else if (length(which(colnames(output) == "Hier")) != 0) {
    if (length(which(colnames(output) == "Order_q")) == 0 ||
        length(which(colnames(output) == "Gamma")) == 0 ||
        length(which(colnames(output) == "Alpha")) == 0 ||
        length(which(colnames(output) == "Beta")) == 0 ||
        length(which(colnames(output) == "Synchrony")) == 0) {
      stop('Please put the complete output of "Stay_Single", "Stay_Multiple" or "Stay_Hier" function.')
    } else {
      outtype <- "hier"
    }
  } else {
    if (length(which(colnames(output) == "Site")) == 0 ||
        length(which(colnames(output) == "Order_q")) == 0 ||
        length(which(colnames(output) == "Gamma")) == 0 ||
        length(which(colnames(output) == "Alpha")) == 0 ||
        length(which(colnames(output) == "Beta")) == 0 ||
        length(which(colnames(output) == "Synchrony")) == 0) {
      stop('Please put the complete output of "Stay_Single", "Stay_Multiple" or "Stay_Hier" function.')
    } else {
      outtype <- "multiple"
    }
  }
  
  if (outtype == "single") {
    
    if (length(unique(output$`Plot/Community`)) <= 4) {
      cbPalette <- c("#EA0000","#0066CC","#64A600","#fb8500")
    } else {
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#fb8500"),
                     ggplotColors(length(unique(output$`Plot/Community`)) - 4))
    }
    
    output$`Plot/Community` <- factor(output$`Plot/Community`,
                                      levels = unique(output$`Plot/Community`))
    
    plotout <- ggplot2::ggplot(
      data = output,
      mapping = ggplot2::aes(x = Order_q, y = Stability, color = `Plot/Community`)
    ) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::ylab("Stability") +
      ggplot2::xlab("Order of q") +
      ggplot2::labs(color = "Plot/Community") +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = cbPalette) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(size = 13),
        legend.text  = ggplot2::element_text(size = 12),
        legend.key.size = grid::unit(0.8, "cm"),
        axis.title   = ggplot2::element_text(size = 16)
      )
    
  } else if (outtype == "hier") {
    
    maxhier  <- max(output$Hier)
    hier_num <- unique(output$Hier)
    qq       <- unique(output$Order_q)
    type_name   <- c(paste0("Gamma(", maxhier, ")"),
                     paste0("Alpha(", hier_num[-1], ")"))
    type_name2_2 <- c(paste0("Beta(", hier_num[-1], ")"), "Alpha(1)")
    type_diff    <- paste0("Hier", hier_num[-1])
    
    plotdat1 <- data.frame(
      Order_q  = rep(qq, length(hier_num)),
      Stability = c(dplyr::filter(output, Hier == maxhier)$Gamma,
                    dplyr::filter(output, Hier != maxhier)$Alpha),
      type = rep(type_name, each = length(qq))
    )
    plotdat1$type <- factor(plotdat1$type, levels = type_name)
    
    plotdat2 <- data.frame(
      Order_q  = rep(qq, (length(hier_num) - 1)),
      Synchrony = dplyr::filter(output, Hier != maxhier)$Synchrony,
      type = rep(type_diff, each = length(qq))
    )
    plotdat2$type <- factor(plotdat2$type, levels = type_diff)
    
    plotdat2_2 <- data.frame(
      Order_q  = rep(qq, length(hier_num)),
      Stability = c(dplyr::filter(output, Hier != maxhier)$Beta,
                    dplyr::filter(output, Hier == 1)$Alpha),
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
    
    output$Site <- factor(output$Site, levels = unique(output$Site))
    
    stab_plotdat <- data.frame(
      Site    = rep(output$Site, 4),
      Order_q = rep(output$Order_q, 4),
      value   = c(output$Gamma, output$Alpha, output$Beta, output$Synchrony),
      type    = rep(c("Gamma", "Alpha", "Beta", "Synchrony"), each = nrow(output))
    )
    stab_plotdat$type <- factor(stab_plotdat$type, levels = c("Gamma","Alpha","Beta","Synchrony"))
    
    if (length(unique(stab_plotdat$Site)) <= 4) {
      cbPalette <- c("#EA0000","#0066CC","#64A600","#fb8500")
    } else {
      cbPalette <- c(c("#EA0000","#0066CC","#64A600","#fb8500"),
                     ggplotColors(length(unique(stab_plotdat$Site)) - 4))
    }
    
    plotout <- ggplot2::ggplot(stab_plotdat, ggplot2::aes(x = Order_q, y = value, color = Site)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::facet_wrap(~ type, nrow = 2, scales = "free") +
      ggplot2::ylab("Stability and Synchrony") +
      ggplot2::xlab("Order of q") +
      ggplot2::labs(color = "Site") +
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



#' ggplot2 extension for a Stay_Single or Stay_Multiple object to analysis with an diversity (or other) variable.
#'
#' \code{ggStay_analysis} is a graphical function that based on the output from the function \code{Stay_Single} or \code{Stay_Multiple}. It provides to graph relationships between stability (and synchrony if is multiple assemblages) and an additional diversity (or other) variable .
#'
#' @param output the output obtained from \code{Stay_Single} or \code{Stay_Multiple} and needs to combine with a column that sets as \code{x_variable}. Also, if \code{by_group} is not \code{NULL}, the output also need to combine with the column that sets as \code{by_group}.
#' @param x_variable name of the column of diversity (or other) variable, that will use as the x-axis in the plot.
#' @param by_group name of the column that is a categorical variable for plotting points with different color. And it is required if \code{model = "LMM"}, model uses it as random effect for intercept and slope. Default is \code{NULL}.
#' @param model specifying the fitting model, \code{model = "lm"} for linear model; \code{model = "LMM"} for linear mixed model with random effects for intercept and slope. Default is \code{model = "LMM"}.
#'
#'
#' @import stringr
#' @import lme4
#' @import dplyr
#'
#'
#' @return For an \code{Stay_Single} object, this function return a figure of diversity (or other) variable vs. stability.
#' For an \code{Stay_Multiple} object, this function return a figure that is about diversity (or other) variable vs. (Gamma, Alpha, Beta) stability and synchrony.
#'
#' @examples
#' data("Jena_plot_biomass_data")
#' data("Jena_species_biomass_data")
#' data("Jena_hierarchical_biomass_data")
#' data("Jena_hierarchical_mat")
#'
#' ## Single assemblage
#' # Stability of each single plot
#' single_plot <- do.call(rbind, Jena_plot_biomass_data)
#' output_single_plot_div <- Stay_Single(data=single_plot, order.q=c(1,2), Alltime=TRUE)
#' output_single_plot_div <- data.frame(output_single_plot_div,
#'                                      sowndiv=as.numeric(do.call(rbind,
#'                                        strsplit(output_single_plot_div[,1],"[._]+"))[,2]),
#'                                      block=do.call(rbind,
#'                                        strsplit(output_single_plot_div[,1],"[._]+"))[,1])
#'
#' ggStay_analysis(output=output_single_plot_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' # Stability of each single species
#' single_species <- do.call(rbind, Jena_species_biomass_data)
#' output_single_species_div <- Stay_Single(data=single_species,
#'                                          order.q=c(1,2), Alltime=TRUE)
#' output_single_species_div <- data.frame(output_single_species_div,
#'                               sowndiv=as.numeric(do.call(rbind,
#'                                       strsplit(output_single_species_div[,1],"[._]+"))[,3]),
#'                               block=do.call(rbind,
#'                                     strsplit(output_single_species_div[,1],"[._]+"))[,2])
#'
#' ggStay_analysis(output=output_single_species_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#'
#' ## Multiple assemblages
#' # Stability of multiple plots
#' multiple_plot <- Jena_plot_biomass_data
#' output_multi_plot_div <- Stay_Multiple(data=multiple_plot, order.q=c(1,2), Alltime=TRUE)
#' output_multi_plot_div <- data.frame(output_multi_plot_div, sowndiv=rep(c(16,8,4,2,1),8),
#'                                     block=rep(rep(c("B1","B2","B3","B4"),each=5),2))
#'
#' ggStay_analysis(output=output_multi_plot_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' # Stability of multiple species in plot
#' multiple_species <- Jena_species_biomass_data
#' output_multi_species_div <- Stay_Multiple(data=multiple_species,
#'                                               order.q=c(1,2), Alltime=TRUE)
#' output_multi_species_div <- data.frame(output_multi_species_div,
#'                              sowndiv=as.numeric(do.call(rbind,
#'                                      strsplit(output_multi_species_div[,1],"_"))[,3]),
#'                              block=do.call(rbind,
#'                                    strsplit(output_multi_species_div[,1],"_"))[,2])
#'
#' ggStay_analysis(output=output_multi_species_div, x_variable="sowndiv",
#'                     by_group="block", model="LMM")
#'
#' @export

ggStay_analysis <- function(output, x_variable, by_group = NULL, model = "LMM") {
  
  # local helper: safe color generator
  ggplotColors <- function(n) {
    grDevices::hcl(h = seq(15, 375, length.out = n + 1), l = 65, c = 100)[1:n]
  }
  
  ## --- Basic checks ---------------------------------------------------------
  if (length(which(colnames(output) == "Stability")) != 0) {
    if ((length(which(colnames(output) == "Plot/Community")) == 0 &&
         length(which(colnames(output) == "Plot.Community")) == 0) ||
        length(which(colnames(output) == "Order_q")) == 0) {
      stop('Please put the complete output of "Stay_Single" or "Stay_Multiple" function.')
    }
  } else {
    if (length(which(colnames(output) == "Site")) == 0 ||
        length(which(colnames(output) == "Order_q")) == 0 ||
        length(which(colnames(output) == "Gamma")) == 0 ||
        length(which(colnames(output) == "Alpha")) == 0 ||
        length(which(colnames(output) == "Beta")) == 0  ||
        length(which(colnames(output) == "Synchrony")) == 0) {
      stop('Please put the complete output of "Stay_Single" or "Stay_Multiple" function.')
    }
  }
  
  if (length(which(colnames(output) == x_variable)) == 0) {
    stop('The output data need to combine a column setting as x_variable.')
  } else {
    colnames(output)[which(colnames(output) == x_variable)] <- "Xvariable"
  }
  
  if (!is.null(by_group)) {
    if (length(which(colnames(output) == by_group)) == 0) {
      stop('The output data need to combine a column setting as by_group.')
    } else {
      colnames(output)[which(colnames(output) == by_group)] <- "Gvariable"
    }
  }
  if (identical(model, "LMM") && is.null(by_group)) {
    stop('For linear mixed model, you need to set by_group variable as random effect.')
  }
  
  if (!is.null(by_group)) output$Gvariable <- as.factor(output$Gvariable)
  
  ## ========================================================================
  ## A) Single 
  ## ========================================================================
  if (length(which(colnames(output) == "Stability")) != 0) {
    
    lm_sign  <- c()
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
    
    lm_sign            <- as.data.frame(lm_sign)
    lm_slope           <- as.data.frame(lm_slope)
    rownames(lm_sign)  <- paste0("q = ", unique(output$Order_q))
    rownames(lm_slope) <- paste0("q = ", unique(output$Order_q))
    
    plotdata$Order_q <- paste0("q = ", plotdata$Order_q)
    plotdata$sign <- sapply(plotdata$Order_q, function(yy) lm_sign[which(rownames(lm_sign) == yy), 1])
    plotdata$sign <- factor(plotdata$sign, levels = c("significant", "non-significant"))
    slope_text <- data.frame(slope = paste0("slope = ", round(lm_slope[, 1], 4)),
                             Order_q = rownames(lm_slope))
    
    if (!any(plotdata$sign == "non-significant")) {
      dummy <- plotdata[1, ]; dummy$Xvariable <- NA; dummy$pred <- NA; dummy$Stability <- NA; dummy$sign <- "non-significant"
      plotdata <- rbind(plotdata, dummy)
    }
    if (!any(plotdata$sign == "significant")) {
      dummy <- plotdata[1, ]; dummy$Xvariable <- NA; dummy$pred <- NA; dummy$Stability <- NA; dummy$sign <- "significant"
      plotdata <- rbind(plotdata, dummy)
    }
    
    if (identical(model, "LMM")) {
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
                             size = 3.5, hjust = 1, vjust = -3.5) +
          ggplot2::geom_text(data = plotdata_text_part,
                             ggplot2::aes(x = -Inf, y = -Inf, label = slope, color = Gvariable),
                             x = max(plotdata$Xvariable, na.rm = TRUE),
                             y = min(plotdata$Stability, na.rm = TRUE),
                             size = 3.5, hjust = 2.1, vjust = -1.7,
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
    
    ## ========================================================================
    ## B) Multiple (Site) case
    ## ========================================================================
  } else {
    
    lm_sign  <- c()
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
        c(
          gamma     = sum_gamma$coefficients[2,1],
          alpha     = sum_alpha$coefficients[2,1],
          beta_add  = sum_beta_add$coefficients[2,1],
          synchrony = sum_syn$coefficients[2,1]
        )
      )
      plotdata <- rbind(plotdata,
                        data.frame(subdata, pred_G = pred_G, pred_A = pred_A, pred_BA = pred_BA, pred_S = pred_S))
    }
    
    lm_sign  <- as.data.frame(lm_sign)
    lm_slope <- as.data.frame(lm_slope)
    rownames(lm_sign)  <- paste0("q = ", unique(output$Order_q))
    rownames(lm_slope) <- paste0("q = ", unique(output$Order_q))
    
    plotdata$Order_q <- paste0("q = ", plotdata$Order_q)
    plotdata$sign_G  <- sapply(plotdata$Order_q, function(yy) lm_sign[which(rownames(lm_sign) == yy), 1])
    plotdata$sign_A  <- sapply(plotdata$Order_q, function(yy) lm_sign[which(rownames(lm_sign) == yy), 2])
    plotdata$sign_BA <- sapply(plotdata$Order_q, function(yy) lm_sign[which(rownames(lm_sign) == yy), 3])
    plotdata$sign_S  <- sapply(plotdata$Order_q, function(yy) lm_sign[which(rownames(lm_sign) == yy), 4])
    
    plotdata_Stab <- data.frame(
      Site     = rep(plotdata$Site, 4),
      Order_q  = rep(plotdata$Order_q, 4),
      Stability= c(plotdata$Gamma, plotdata$Alpha, plotdata$Beta, plotdata$Synchrony),
      pred     = c(plotdata$pred_G, plotdata$pred_A, plotdata$pred_BA, plotdata$pred_S),
      sign     = c(plotdata$sign_G, plotdata$sign_A, plotdata$sign_BA, plotdata$sign_S),
      type     = rep(c("Gamma","Alpha","Beta","Synchrony"), each = nrow(plotdata)),
      Xvariable= rep(plotdata$Xvariable, 4)
    )
    
    if (!is.null(by_group)) {
      plotdata_Stab$Gvariable <- as.factor(rep(plotdata$Gvariable, 4))
      plotdata$Gvariable      <- as.factor(plotdata$Gvariable)
    }
    
    plotdata_Stab$sign    <- factor(plotdata_Stab$sign, levels = c("significant", "non-significant"))
    plotdata_Stab$type    <- factor(plotdata_Stab$type, levels = c("Gamma","Alpha","Beta","Synchrony"))
    plotdata_Stab$Order_q <- as.factor(plotdata_Stab$Order_q)
    
    if (!any(plotdata_Stab$sign == "non-significant")) {
      dummy <- plotdata_Stab[1, ]; dummy$Xvariable <- NA; dummy$pred <- NA; dummy$Stability <- NA; dummy$sign <- "non-significant"
      plotdata_Stab <- rbind(plotdata_Stab, dummy)
    }
    if (!any(plotdata_Stab$sign == "significant")) {
      dummy <- plotdata_Stab[1, ]; dummy$Xvariable <- NA; dummy$pred <- NA; dummy$Stability <- NA; dummy$sign <- "significant"
      plotdata_Stab <- rbind(plotdata_Stab, dummy)
    }
    
    slope_text_Stab <- data.frame(
      slope   = paste0("slope = ", format(round(as.vector(as.matrix(lm_slope)), 4), nsmall = 4)),
      Order_q = rep(rownames(lm_slope), 4),
      type    = rep(c("Gamma","Alpha","Beta","Synchrony"), each = nrow(lm_slope))
    )
    slope_text_Stab$type    <- factor(slope_text_Stab$type, levels = c("Gamma","Alpha","Beta","Synchrony"))
    slope_text_Stab$Order_q <- as.factor(slope_text_Stab$Order_q)
    
    if (identical(model, "LMM")) {
      plotdata_text_part <- data.frame(
        slope    = c(plotdata_text_part$slope_gamma,
                     plotdata_text_part$slope_alpha,
                     plotdata_text_part$slope_beta,
                     plotdata_text_part$slope_synchrony),
        Order_q  = rep(plotdata_text_part$Order_q, 4),
        type     = rep(c("Gamma","Alpha","Beta","Synchrony"), each = nrow(plotdata_text_part)),
        Gvariable= rep(plotdata_text_part$Gvariable, 4)
      )
      plotdata_text_part$type     <- as.factor(plotdata_text_part$type)
      plotdata_text_part$Order_q  <- as.factor(plotdata_text_part$Order_q)
      plotdata_text_part$Gvariable<- as.factor(plotdata_text_part$Gvariable)
      part_fit$Gvariable <- as.factor(part_fit$Gvariable)
      part_fit$Order_q   <- as.factor(part_fit$Order_q)
      part_fit$type      <- as.factor(part_fit$type)
    }
    
    tyG <- min(dplyr::filter(plotdata_Stab, type == "Gamma")$Stability,     na.rm = TRUE)
    tyA <- min(dplyr::filter(plotdata_Stab, type == "Alpha")$Stability,     na.rm = TRUE)
    tyBA<- max(dplyr::filter(plotdata_Stab, type == "Beta")$Stability,      na.rm = TRUE)
    tyS <- min(dplyr::filter(plotdata_Stab, type == "Synchrony")$Stability, na.rm = TRUE)
    
    Q    <- length(unique(plotdata_Stab$Order_q))
    gnum <- if (!is.null(by_group)) length(unique(plotdata$Gvariable)) else 0
    
    if (!is.null(by_group)) {
      if (length(unique(plotdata_Stab$Gvariable)) <= 4) {
        cbPalette <- c("#EA0000", "#64A600", "#0066CC", "#fb8500")
      } else {
        cbPalette <- c(c("#EA0000", "#64A600", "#0066CC", "#fb8500"),
                       ggplotColors(length(unique(plotdata_Stab$Gvariable)) - 4))
      }
      
      plotout1 <- ggplot2::ggplot() +
        ggplot2::geom_point(data = plotdata_Stab,
                            ggplot2::aes(x = Xvariable, y = Stability, color = Gvariable), size = 2.7) +
        ggplot2::geom_segment(data = part_fit,
                              ggplot2::aes(x = x_min, xend = x_max,
                                           y = intercept + slope * x_min,
                                           yend = intercept + slope * x_max,
                                           color = Gvariable)) +
        ggplot2::geom_line(data = plotdata_Stab,
                           ggplot2::aes(x = Xvariable, y = pred, linetype = sign),
                           linewidth = 1.2, color = "black") +
        ggplot2::scale_linetype_manual(values = c("solid", "dashed"), drop = FALSE) +
        ggplot2::scale_color_manual(values = cbPalette) +
        ggplot2::facet_grid(type ~ Order_q, scales = "free_y") +
        ggplot2::labs(linetype = "", color = by_group) +
        ggplot2::xlab(x_variable) + ggplot2::ylab("Stability and Synchrony") +
        ggplot2::theme_bw() +
        # slope_text_Stab：y 向量長度 = 4 * Q，與資料列數一致；hjust/vjust 用常數
        ggplot2::geom_text(data = slope_text_Stab,
                           ggplot2::aes(x = -Inf, y = -Inf, label = slope),
                           x = max(plotdata_Stab$Xvariable, na.rm = TRUE),
                           y = rep(c(tyG, tyA, tyBA, tyS), each = Q),
                           size = 3.5, hjust = 1, vjust = -3.5) +
        # plotdata_text_part：y 長度 = gnum * Q * 4；hjust/vjust 用常數
        ggplot2::geom_text(data = plotdata_text_part,
                           ggplot2::aes(x = -Inf, y = -Inf, label = slope, color = Gvariable),
                           x = rep(max(plotdata_Stab$Xvariable, na.rm = TRUE), gnum * 4 * Q),
                           y = rep(c(tyG, tyA, tyBA, tyS), each = gnum * Q),
                           size = 3.5, hjust = 2.1, vjust = -1.7,
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
      
      plotout1 <- ggplot2::ggplot(plotdata_Stab, ggplot2::aes(x = Xvariable, y = Stability)) +
        ggplot2::geom_point(size = 2.7) +
        ggplot2::geom_line(ggplot2::aes(x = Xvariable, y = pred, linetype = sign),
                           linewidth = 1.2, color = "black") +
        ggplot2::scale_linetype_manual(values = c("solid", "dashed"), drop = FALSE) +
        ggplot2::facet_grid(type ~ Order_q, scales = "free_y") +
        ggplot2::labs(linetype = "") +
        ggplot2::xlab(x_variable) + ggplot2::ylab("Stability and Synchrony") +
        ggplot2::theme_bw() +
        ggplot2::geom_text(data = slope_text_Stab,
                           ggplot2::aes(x = -Inf, y = -Inf, label = slope),
                           x = max(plotdata_Stab$Xvariable, na.rm = TRUE),
                           y = rep(c(tyG, tyA, tyBA, tyS), each = Q),
                           size = 5, hjust = 1, vjust = c(0.1, 0.1, 1.1, 0.1)) +
        ggplot2::theme(
          strip.text   = ggplot2::element_text(size = 13),
          legend.title = ggplot2::element_text(size = 13),
          legend.text  = ggplot2::element_text(size = 12),
          legend.key.size = grid::unit(0.8, "cm"),
          axis.title   = ggplot2::element_text(size = 16)
        ) +
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
utils::globalVariables(c("COLOR", "Order_q", "Site", "Stability", "Synchrony"
))





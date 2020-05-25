library(rlist)
library(r.jive)
library(SpatioTemporal)

jive.iter = function (data, rankJ = 1, rankA = rep(1, length(data)), conv = 1e-06, 
                      maxiter = 1000, orthIndiv = TRUE, showProgress = TRUE, vertical = TRUE) 
{
  l <- length(data)
  A <- list()
  for (i in 1:l) {
    A[[i]] <- matrix(0, nrow(data[[i]]), ncol(data[[i]]))
  }
  Xtot <- do.call(rbind, data)
  Jtot <- matrix(0, nrow(Xtot), ncol(Xtot))
  Atot <- matrix(0, nrow(Xtot), ncol(Xtot))
  nrun = 0
  converged = F
  if (orthIndiv) {
    Vind <- list()
  }
  U = list()
  W = list()
  S = list()
  T. = list()
  Tind = list()
  Wind = list()
  # ---------------
  while (nrun < maxiter & !converged) {
    Jlast <- Jtot
    Alast <- Atot
    timeJ <- system.time({
      if (rankJ > 0) {
        temp <- Xtot - Atot
        s <- svdwrapper(temp, nu = rankJ, nv = rankJ)
        Jtot <- s$u[, 1:rankJ] %*% diag(x = s$d[1:rankJ], 
                                        nrow = rankJ) %*% t(s$v[, 1:rankJ])
        U_ <- s$u[, 1:rankJ]
        V <- s$v[, 1:rankJ]
      }
      else {
        Jtot <- matrix(0, nrow(Xtot), ncol(Xtot))
        U_ <- matrix(0, nrow(Xtot), rankJ)
        V <- matrix(0, ncol(Xtot), rankJ)
      }
      temp <- Jtot
      J <- list()
      for (i in 1:l) {
        J[[i]] <- temp[1:nrow(data[[i]]), ]
        temp <- temp[-(1:nrow(data[[i]])), ]
      }
    })
    timeA <- system.time({
      A <- list()
      for (i in 1:l) {
        if (rankA[i] > 0) {
          temp <- (data[[i]] - J[[i]]) %*% (diag(ncol(Xtot)) - V %*% t(V))
          if (orthIndiv & nrun > 0) {
            for (j in (1:l)[-i]) {
              temp <- temp %*% (diag(ncol(Xtot)) - Vind[[j]] %*% 
                                  t(Vind[[j]]))
            }
          }
          s <- svdwrapper(temp, nu = rankA[i], nv = rankA[i])
          if (orthIndiv) {
            Vind[[i]] <- s$v[, 1:rankA[i]]
          }
          A[[i]] <- s$u[, 1:rankA[i]] %*% diag(x = s$d[1:rankA[i]], 
                                               nrow = rankA[i]) %*% t(s$v[, 1:rankA[i]])
          if (vertical){
            Wind[[i]] <- s$v[, 1:rankA[i]]
            Tind[[i]] <- s$u[, 1:rankA[i]]
          }else{
            Wind[[i]] <- s$u[, 1:rankA[i]]
            Tind[[i]] <- s$v[, 1:rankA[i]]
          }
        }
        else {
          A[[i]] <- matrix(0, nrow(data[[i]]), ncol(data[[i]]))
          if (orthIndiv) {
            Vind[[i]] <- matrix(0, ncol(Xtot), rankA[[i]])
          }
        }
      }
    })
    if (orthIndiv & nrun == 0) {
      for (i in 1:l) {
        for (j in (1:l)[-i]) {
          A[[i]] <- A[[i]] %*% (diag(ncol(Xtot)) - Vind[[j]] %*% 
                                  t(Vind[[j]]))
        }
      }
      for (i in 1:l) {
        if (rankA[i] > 0) {
          s <- svdwrapper(A[[i]], nu = rankA[i], nv = rankA[i])
          Vind[[i]] <- s$v[, 1:rankA[i]]
          if (vertical){
            Wind[[i]] <- s$v[, 1:rankA[i]]
            Tind[[i]] <- s$u[, 1:rankA[i]]
          }else{
            Wind[[i]] <- s$u[, 1:rankA[i]]
            Tind[[i]] <- s$v[, 1:rankA[i]]
          }
        }
      }
    }
    Atot <- do.call(rbind, A)
    if (norm(Jtot - Jlast, type = "f") <= conv & norm(Atot - 
                                                      Alast, type = "f") <= conv) {
      converged <- T
    }
    nrun = nrun + 1
    if (vertical){
      U = list.append(U, V)
      S = list.append(S, U_)
    }else{
      U = list.append(U, U_)
      S = list.append(S, V)
    }
    W = list.append(W, Wind)
    T. = list.append(T., Tind)
  }
  if (showProgress) {
    if (converged) {
      cat(paste("JIVE algorithm converged after ", nrun, 
                " iterations.\n"))
    }
    else {
      cat(paste("JIVE algorithm did not converge after ", 
                nrun, " iterations.\n"))
    }
  }
  return(list(data = data, joint = J, individual = A, rankJ, 
              rankA, method = "given", U = U, W = W, T = T., S = S, nrun = nrun))
}


jive.perm = function (data, nperms = 100, alpha = 0.05, est = TRUE, conv = 1e-06, 
          maxiter = 1000, orthIndiv = TRUE, showProgress = TRUE, vertical = TRUE) 
{
  nrun <- 0
  Jperp <- list()
  Aperp <- list()
  for (i in 1:length(data)) {
    Jperp[[i]] <- matrix(0, nrow = nrow(data[[i]]), ncol = ncol(data[[i]]))
    Aperp[[i]] <- matrix(0, nrow = nrow(data[[i]]), ncol = ncol(data[[i]]))
  }
  last <- rep(-2, length(data) + 1)
  current <- rep(-1, length(data) + 1)
  while (!isTRUE(all.equal(last, current)) & nrun < 10) {
    last <- current
    if (showProgress) {
      if (nrun == 0) {
        cat("Estimating  joint and individual ranks via permutation...\n")
      }
      else {
        cat("Re-estimating  joint and individual ranks via permutation...\n")
      }
    }
    full <- list()
    for (i in 1:length(data)) {
      full[[i]] <- data[[i]] - Aperp[[i]]
    }
    n <- ncol(full[[1]])
    actual <- svdwrapper(do.call(rbind, full), nu = 0, nv = 0)$d
    perms <- matrix(NA, nperms, min(n, sum(unlist(lapply(data, 
                                                         nrow)))))
    for (i in 1:nperms) {
      temp <- list()
      for (j in 1:length(data)) {
        temp[[j]] <- full[[j]][, sample(1:n, n, replace = F)]
      }
      perms[i, ] <- svdwrapper(do.call(rbind, temp), nu = 0, 
                               nv = 0)$d
    }
    rankJ <- 0
    for (i in 1:n) {
      if (actual[i] > quantile(perms[, i], 1 - alpha)) {
        rankJ <- rankJ + 1
      }
      else {
        break
      }
    }
    rankJ <- max(rankJ, last[1])
    rankA <- c()
    for (i in 1:length(data)) {
      ind <- data[[i]] - Jperp[[i]]
      actual <- svdwrapper(ind, nu = 0, nv = 0)$d
      perms <- matrix(NA, nperms, min(n, nrow(data[[i]])))
      for (k in 1:nperms) {
        perm <- t(ind)
        pind <- order(c(col(perm)), runif(length(perm)))
        perm <- matrix(perm[pind], nrow = nrow(ind), 
                       ncol = n, byrow = TRUE)
        perms[k, ] <- svdwrapper(perm, nu = 0, nv = 0)$d
      }
      rankA[i] <- 0
      for (j in 1:n) {
        if (actual[j] > quantile(perms[, j], 1 - alpha)) {
          rankA[i] <- rankA[i] + 1
        }
        else {
          break
        }
      }
    }
    current <- c(rankJ, rankA)
    if (!isTRUE(all.equal(last, current))) {
      dataR <- list()
      if (est) {
        u <- list()
        for (i in 1:length(data)) {
          if (nrow(data[[i]]) > ncol(data[[i]])) {
            temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                               nv = ncol(data[[i]]))
            dataR[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                               nrow = ncol(data[[1]])) %*% t(temp$v[, 
                                                                    1:ncol(data[[1]])])
            u[[i]] <- temp$u
          }
          else {
            u[[i]] <- diag(1, nrow(data[[i]]))
            dataR[[i]] <- data[[i]]
          }
        }
      }
      else {
        dataR <- data
      }
      if (showProgress) {
        cat("Running JIVE algorithm for ranks:\njoint rank:", 
            rankJ, ", individual ranks:", rankA, "\n")
      }
      tempjive <- jive.iter(dataR, rankJ, rankA, conv = conv, 
                            maxiter = maxiter, orthIndiv = orthIndiv, showProgress = showProgress, vertical = vertical)
      Jperp = tempjive$joint
      Aperp = tempjive$individual
      if (est) {
        for (i in 1:length(data)) {
          Jperp[[i]] <- u[[i]] %*% Jperp[[i]]
          Aperp[[i]] <- u[[i]] %*% Aperp[[i]]
        }
      }
    }
    nrun <- nrun + 1
  }
  converged <- ifelse(nrun == 10, F, T)
  if (showProgress) {
    cat("Final joint rank:", rankJ, ", final individual ranks:", 
        rankA, "\n")
  }
  return(list(data = data, joint = Jperp, individual = Aperp, 
              rankJ = rankJ, rankA = rankA, converged = converged,
              U = tempjive$U, W = tempjive$W, nrun = tempjive$nrun))
}

bic.jive = function (data, n = unlist(lapply(data, ncol)) * unlist(lapply(data, 
    nrow)), d = unlist(lapply(data, nrow)), conv = 1e-06, maxiter = 1000, 
    orthIndiv = TRUE, showProgress = TRUE, vertical = TRUE) 
{
    l <- length(data)
    nc <- ncol(data[[1]])
    lambda <- log(sum(n))
    sse <- c()
    Xtot <- do.call(rbind, data)
    bic.improve <- T
    rankJ <- 0
    rankA <- rep(0, l)
    if (showProgress) {
        cat("Running JIVE algorithm for ranks:\njoint rank:", 
            rankJ, ", individual ranks:", rankA, "\n")
    }
    current <- jive.iter(data, rankJ, rankA, conv, maxiter, orthIndiv, 
        showProgress = showProgress, vertical = vertical)
    for (k in 1:length(data)) {
        sse[k] <- norm(data[[k]] - current$joint[[k]] - current$individual[[k]], 
            type = "f")^2
    }
    p.jive <- 0
    current.bic <- sum(n * log(sse/n)) + p.jive * lambda
    bic.table <- c(rankJ, rankA, current.bic)
    while (bic.improve) {
        bic.improve <- F
        temp <- list()
        if (showProgress) {
            cat("Running JIVE algorithm for ranks:\njoint rank:", 
                rankJ + 1, ", individual ranks:", rankA, "\n")
        }
        temp[[1]] <- jive.iter(data, rankJ + 1, rankA, conv, 
            maxiter, orthIndiv, showProgress = showProgress, vertical = vertical)
        for (k in 1:length(data)) {
            sse[k] <- norm(data[[k]] - temp[[1]]$joint[[k]] - 
                temp[[1]]$individual[[k]], type = "f")^2
        }
        p.jive <- sum(sum(d):(sum(d) - (rankJ + 1) + 1)) + sum(nc:(nc - 
            (rankJ + 1) + 1)) + pjsum(d, rankA) + pjsum(rep(nc, 
            length(data)) - (rankJ + 1), rankA)
        bic <- sum(n * log(sse/n)) + p.jive * lambda
        bicvec <- bic
        bic.table <- rbind(bic.table, c(rankJ + 1, rankA, bic))
        for (i in 1:l) {
            tempR <- rankA
            tempR[i] <- tempR[i] + 1
            if (tempR[i] < min(n, nrow(data[[i]]))) {
                if (showProgress) {
                  cat("Running JIVE algorithm for ranks:\njoint rank:", 
                    rankJ, ", individual ranks:", tempR, "\n")
                }
                temp[[i + 1]] <- jive.iter(data, rankJ, tempR, 
                  conv, maxiter, orthIndiv, showProgress = showProgress, vertical = vertical)
                for (k in 1:length(data)) {
                  sse[k] <- norm(data[[k]] - temp[[i + 1]]$joint[[k]] - 
                    temp[[i + 1]]$individual[[k]], type = "f")^2
                }
                p.jive <- ifelse(rankJ == 0, 0, sum(sum(d):(sum(d) - 
                  rankJ + 1)) + sum(nc:(nc - rankJ + 1))) + pjsum(d, 
                  tempR) + pjsum(rep(nc, length(data)) - rankJ, 
                  tempR)
                bic <- sum(n * log(sse/n)) + p.jive * lambda
            }
            else {
                bic <- NA
            }
            bicvec <- c(bicvec, bic)
            bic.table <- rbind(bic.table, c(rankJ, tempR, bic))
        }
        lowest.bic <- temp[[which.min(bicvec)]]
        if (min(bicvec, na.rm = T) < current.bic) {
            bic.improve <- T
            current <- lowest.bic
            current.bic <- min(bicvec, na.rm = T)
            if (which.min(bicvec) == 1) {
                rankJ <- rankJ + 1
            }
            else {
                rankA[which.min(bicvec) - 1] <- rankA[which.min(bicvec) - 
                  1] + 1
            }
        }
    }
    return(list(data = data, joint = current$joint, individual = current$individual, 
        rankJ = rankJ, rankA = rankA, bic.table = bic.table, U = current$U, W = current$W, nrun = current$nrun))
}

jive.multigroup = function (data, rankJ = 1, rankA = rep(1, length(data)), method = "perm", 
          dnames = names(data), conv = "default", maxiter = 1000, scale = TRUE, 
          center = TRUE, orthIndiv = TRUE, est = TRUE, showProgress = TRUE) 
{
  l <- length(data)
  d <- c()
  for (i in 1:l) {
    d[i] <- ncol(data[[i]])
  }
  for (i in 1:l) {
    temp <- SVDmiss(data[[i]], ncomp = min(ncol(data[[i]]), 
                                           nrow(data[[i]])))[[1]]
    data[[i]] <- temp$u %*% diag(x = temp$d) %*% t(temp$v)
  }
  centerValues <- list()
  scaleValues <- c()
  for (i in 1:l) {
    if (center) {
      centerValues[[i]] <- apply(data[[i]], 2, mean, na.rm = T)
      data[[i]] <- sweep(data[[i]], 2, centerValues[[i]])
    }
    if (!center) {
      centerValues[[i]] <- rep(0, d[i])
    }
    if (scale) {
      scaleValues[i] <- norm(data[[i]], type = "f")
      data[[i]] <- data[[i]]/scaleValues[i]
    }
    if (!scale) {
      scaleValues[[i]] <- 1
    }
  }
  if (conv == "default") {
    conv = 10^(-6) * norm(do.call(rbind, data), type = "f")
  }
  orig <- data
  if (method == "given") {
    if (est) {
      u <- list()
      for (i in 1:l) {
        if (nrow(data[[i]]) > ncol(data[[i]])) {
          temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                             nv = ncol(data[[i]]))
          data[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                            nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
          u[[i]] <- temp$u
        }
        else {
          u[[i]] <- diag(1, nrow(data[[i]]))
        }
      }
    }
    if (showProgress) {
      cat("Running JIVE algorithm for ranks:\njoint rank:", 
          rankJ, ", individual ranks:", rankA, "\n")
    }
    temp <- jive.iter(data, rankJ, rankA, conv = conv, maxiter = maxiter, 
                      orthIndiv = orthIndiv, showProgress = showProgress, vertical = TRUE)
    joint <- temp$joint
    individual <- temp$individual
    if (est) {
      for (i in 1:l) {
        joint[[i]] <- u[[i]] %*% joint[[i]]
        individual[[i]] <- u[[i]] %*% individual[[i]]
      }
    }
  }
  else if (method == "perm") {
    temp <- jive.perm(data, est = est, conv = conv, maxiter = maxiter, 
                      orthIndiv = orthIndiv, showProgress = showProgress)
    joint <- temp$joint
    individual <- temp$individual
    rankJ <- temp$rankJ
    rankA <- temp$rankA
    converged <- temp$converged
  }
  else if (method == "bic") {
    if (est) {
      u <- list()
      for (i in 1:l) {
        if (nrow(data[[i]]) > ncol(data[[i]])) {
          temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                             nv = ncol(data[[i]]))
          data[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                            nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
          u[[i]] <- temp$u
        }
        else {
          u[[i]] <- diag(1, nrow(data[[i]]))
        }
      }
    }
    temp <- bic.jive(data, n, d, conv = conv, maxiter = maxiter, 
                     orthIndiv = orthIndiv, showProgress = showProgress, vertical = TRUE)
    joint <- temp$joint
    individual <- temp$individual
    rankJ <- temp$rankJ
    rankA <- temp$rankA
    bic.table <- temp$bic.table
    if (est) {
      for (i in 1:l) {
        joint[[i]] <- u[[i]] %*% joint[[i]]
        individual[[i]] <- u[[i]] %*% individual[[i]]
      }
    }
  }
  if (is.null(dnames)) {
    for (i in 1:l) {
      names(orig)[[i]] <- paste("Source", i, sep = "_")
    }
  }
  else {
    names(orig) <- dnames
  }
  result <- list(data = orig, joint = joint, individual = individual, 
                 rankJ = rankJ, rankA = rankA, method = method, U = temp$U, W = temp$W, nrun = temp$nrun)
  if (method == "bic") {
    result$bic.table <- bic.table
  }
  if (method == "perm") {
    result$converged <- converged
  }
  result$scale <- list(center, scale, centerValues, scaleValues)
  names(result$scale) <- c("Center", "Scale", "Center Values", 
                           "Scale Values")
  class(result) <- "jive"
  return(result)
}

jive.multisource = function (data, rankJ = 1, rankA = rep(1, length(data)), method = "perm", 
    dnames = names(data), conv = "default", maxiter = 1000, scale = TRUE, 
    center = TRUE, orthIndiv = TRUE, est = TRUE, showProgress = TRUE) 
{
    l <- length(data)
    d <- c()
    for (i in 1:length(data)) {
        d[i] <- nrow(data[[i]])
    }
    for (i in 1:l) {
        temp <- SVDmiss(data[[i]], ncomp = min(ncol(data[[i]]), 
            nrow(data[[i]])))[[1]]
        data[[i]] <- temp$u %*% diag(x = temp$d) %*% t(temp$v)
    }
    centerValues <- list()
    scaleValues <- c()
    for (i in 1:l) {
        if (center) {
            centerValues[[i]] <- apply(data[[i]], 1, mean, na.rm = T)
            data[[i]] <- sweep(data[[i]], 1, centerValues[[i]])
        }
        if (!center) {
            centerValues[[i]] <- rep(0, d[i])
        }
        if (scale) {
            scaleValues[i] <- norm(data[[i]], type = "f")
            data[[i]] <- data[[i]]/scaleValues[i]
        }
        if (!scale) {
            scaleValues[i] <- 1
        }
    }
    if (conv == "default") {
        conv = 10^(-6) * norm(do.call(rbind, data), type = "f")
    }
    orig <- data
    if (method == "given") {
        if (est) {
            u <- list()
            for (i in 1:l) {
              u[[i]] <- diag(1, nrow(data[[i]]))
                # if (nrow(data[[i]]) > ncol(data[[i]])) {
                #   temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                #     nv = ncol(data[[i]]))
                #   data[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                #     nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
                #   u[[i]] <- temp$u
                # }
                # else {
                #   u[[i]] <- diag(1, nrow(data[[i]]))
                # }
            }
        }
        if (showProgress) {
            cat("Running JIVE algorithm for ranks:\njoint rank:", 
                rankJ, ", individual ranks:", rankA, "\n")
        }
        temp <- jive.iter(data, rankJ, rankA, conv = conv, maxiter = maxiter, 
            orthIndiv = orthIndiv, showProgress = showProgress, vertical = FALSE)
        joint <- temp$joint
        individual <- temp$individual
        if (est) {
            for (i in 1:l) {
                joint[[i]] <- u[[i]] %*% joint[[i]]
                individual[[i]] <- u[[i]] %*% individual[[i]]
            }
        }
    }
    else if (method == "perm") {
        temp <- jive.perm(data, est = est, conv = conv, maxiter = maxiter, 
            orthIndiv = orthIndiv, showProgress = showProgress, vertical = FALSE)
        joint <- temp$joint
        individual <- temp$individual
        rankJ <- temp$rankJ
        rankA <- temp$rankA
        converged <- temp$converged
    }
    else if (method == "bic") {
        if (est) {
            u <- list()
            for (i in 1:l) {
                if (nrow(data[[i]]) > ncol(data[[i]])) {
                  temp <- svdwrapper(data[[i]], nu = ncol(data[[i]]), 
                    nv = ncol(data[[i]]))
                  data[[i]] <- diag(x = temp$d[1:ncol(data[[1]])], 
                    nrow = ncol(data[[1]])) %*% t(temp$v[, 1:ncol(data[[1]])])
                  u[[i]] <- temp$u
                }
                else {
                  u[[i]] <- diag(1, nrow(data[[i]]))
                }
            }
        }
        temp <- bic.jive(data, n, d, conv = conv, maxiter = maxiter, 
            orthIndiv = orthIndiv, showProgress = showProgress, vertical  = FALSE)
        joint <- temp$joint
        individual <- temp$individual
        rankJ <- temp$rankJ
        rankA <- temp$rankA
        bic.table <- temp$bic.table
        if (est) {
            for (i in 1:l) {
                joint[[i]] <- u[[i]] %*% joint[[i]]
                individual[[i]] <- u[[i]] %*% individual[[i]]
            }
        }
    }
    if (is.null(dnames)) {
        for (i in 1:l) {
            names(orig)[[i]] <- paste("Source", i, sep = "_")
        }
    }
    else {
        names(orig) <- dnames
    }
    result <- list(data = orig, joint = joint, individual = individual, 
                   rankJ = rankJ, rankA = rankA, method = method, U = temp$U, W = temp$W, T = temp$T, S = temp$S, nrun = temp$nrun)
    if (method == "bic") {
        result$bic.table <- bic.table
    }
    if (method == "perm") {
        result$converged <- converged
    }
    result$scale <- list(center, scale, centerValues, scaleValues)
    names(result$scale) <- c("Center", "Scale", "Center Values", 
        "Scale Values")
    class(result) <- "jive"
    return(result)
}

C2beta.jive = function(X.list, Y.list, U, W, scale.X = TRUE){
  G = length(X.list)
  centerValues.X <- list()
  scaleValues.X <- list()
  centerValues.Y <- list()
  scaleValues.Y <- list()
  for (g in 1:G){
    # X
    centerValues.X[[g]] = apply(X.list[[g]], 2, mean)
    X.list[[g]] = sweep(X.list[[g]], 2, centerValues.X[[g]])
    if (scale.X){
      scaleValues.X[[g]] = norm(X.list[[g]], type = "f")
    }else{
      scaleValues.X[[g]] = 1
    }
    X.list[[g]] = X.list[[g]]/scaleValues.X[[g]]
    
    # Y
    centerValues.Y[[g]] = mean(Y.list[[g]])
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, centerValues.Y[[g]])
    scaleValues.Y[[g]] = sqrt(mean(Y.list[[g]]^2))
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, scaleValues.Y[[g]], FUN = "/")
  }
  X = do.call(rbind, X.list)
  Y = do.call(rbind, Y.list)
  beta.C = C2beta(X, Y, U, lambda = 0)$beta
  Yhat.homo.list = lapply(1:G, function(g) X.list[[g]]%*%beta.C)
  Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
  beta.Cind = lapply(1:G, function(g) C2beta(X.list[[g]], Y.heter.list[[g]], W[[g]], lambda = 0)$beta)
  
  beta.C0 = beta.C
  beta.Cind0 = beta.Cind
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind))
}


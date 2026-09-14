ESS.Gibbs <- function(b, loglik, sd = 1, max_steps = 300L) {
  slice <- loglik(b) + log(runif(1))
  angle <- runif(1, 0, 2 * pi)
  lower <- angle - 2 * pi
  upper <- angle
  direction <- rnorm(length(b), sd = sd)

  for (step in seq_len(max_steps)) {
    proposal <- b * cos(angle) + direction * sin(angle)
    if (loglik(proposal) > slice) return(proposal)
    if (angle < 0) lower <- angle else upper <- angle
    angle <- runif(1, lower, upper)
  }

  b
}


relabel_mifa <- function(
  z_store, row_store, norm_store, q_store, eig_store,
  psi_store, mu_store, tau_store, sigma_store, mix_store, z_final
) {
  B <- nrow(z_store)
  N <- ncol(z_store)
  P <- dim(row_store)[1]
  G <- dim(row_store)[2]
  K <- dim(eig_store)[1]

  psm <- matrix(0, N, N)
  for (b in seq_len(B)) {
    psm <- psm + outer(z_store[b, ], z_store[b, ], "==")
  }
  psm <- psm / B

  loss <- vapply(seq_len(B), function(b) {
    sum((outer(z_store[b, ], z_store[b, ], "==") - psm)^2)
  }, numeric(1))
  reference <- z_store[which.min(loss), ]
  maps <- matrix(0L, B, G)

  for (b in seq_len(B)) {
    counts <- table(
      factor(z_store[b, ], levels = seq_len(G)),
      factor(reference, levels = seq_len(G))
    )
    map <- as.integer(clue::solve_LSAP(counts, maximum = TRUE))
    permutation <- order(map)
    maps[b, ] <- map
    z_store[b, ] <- map[z_store[b, ]]
    row_store[, , b] <- matrix(row_store[, , b], P, G)[, permutation, drop = FALSE]
    norm_store[, , b] <- matrix(norm_store[, , b], P, G)[, permutation, drop = FALSE]
    q_store[b, ] <- q_store[b, permutation]
    eig_store[, , b] <- matrix(eig_store[, , b], K, G)[, permutation, drop = FALSE]
    psi_store[, , b] <- matrix(psi_store[, , b], P, G)[, permutation, drop = FALSE]
    mu_store[, , b] <- matrix(mu_store[, , b], P, G)[, permutation, drop = FALSE]
    tau_store[, , b] <- matrix(tau_store[, , b], K, G)[, permutation, drop = FALSE]
    sigma_store[b, ] <- sigma_store[b, permutation]
    mix_store[b, ] <- mix_store[b, permutation]
  }

  counts <- table(
    factor(z_final, levels = seq_len(G)),
    factor(reference, levels = seq_len(G))
  )

  list(
    z_store = z_store,
    row_store = row_store,
    norm_store = norm_store,
    q_store = q_store,
    eig_store = eig_store,
    psi_store = psi_store,
    mu_store = mu_store,
    tau_store = tau_store,
    sigma_store = sigma_store,
    mix_store = mix_store,
    z_ref = reference,
    psm = psm,
    maps = maps,
    final_map = as.integer(clue::solve_LSAP(counts, maximum = TRUE))
  )
}


group_sparse_mifa <- function(
  X, G, K = round(3 * log(ncol(X))),
  n_iter = 6000, burn = 2000, thin = 2,
  nu1 = 3, nu2 = 2,
  alpha1 = 2.1, beta1 = 1, alpha2 = 3.1, beta2 = 1,
  rho1 = 3, rho2 = 2, omega1 = 3, omega2 = 2,
  a_psi = 3, bpsi_shape = 4, bpsi_rate = 2,
  mu_sd = 1, mix_alpha = 1, equal_weights = TRUE,
  rank_prop = 0.99, rank_abs = 1.5,
  threshold_mean = 1.28, threshold_sd = 0.5,
  z_init = NULL, update_z = TRUE,
  truth = NULL, pip_cut = 0.5,
  seed = 200, verbose = TRUE
) {
  set.seed(seed)

  X <- as.matrix(X)
  N <- nrow(X)
  P <- ncol(X)
  K <- min(as.integer(K), P - 1L, N - 1L)
  stopifnot(G >= 1L, K >= 1L, n_iter > burn, thin >= 1L)

  feature_names <- colnames(X)
  if (is.null(feature_names)) feature_names <- paste0("V", seq_len(P))

  z <- if (is.null(z_init)) {
    sample(rep(seq_len(G), length.out = N))
  } else {
    as.integer(factor(z_init, levels = sort(unique(z_init))))
  }
  z_initial <- z

  mu <- sapply(seq_len(G), function(g) colMeans(X[z == g, , drop = FALSE]))
  psi <- matrix(1, P, G)
  psi_rate <- rep(bpsi_shape / bpsi_rate, P)
  mix <- rep(1 / G, G)

  threshold_z <- 0
  threshold <- threshold_mean
  alpha <- matrix(threshold + 0.01, P, G)

  omega <- rgamma(1, omega1, omega2)
  sigma <- rgamma(G, rho1, rho2)
  slab <- phi <- delta <- tau <- Lambda <- eta <- vector("list", G)

  log_tail_delta <- pgamma(
    1, alpha2, rate = beta2, lower.tail = FALSE, log.p = TRUE
  )

  for (g in seq_len(G)) {
    delta[[g]] <- c(
      rgamma(1, alpha1, beta1),
      qgamma(
        log(runif(K - 1L)) + log_tail_delta,
        alpha2, rate = beta2, lower.tail = FALSE, log.p = TRUE
      )
    )
    tau[[g]] <- cumprod(delta[[g]])
    phi[[g]] <- matrix(rgamma(P * K, nu1, nu2), P, K)
    precision <- phi[[g]] * matrix(
      rep(omega * sigma[g] * tau[[g]], each = P), P, K
    )
    slab[[g]] <- matrix(rnorm(P * K), P, K) / sqrt(precision)
    Lambda[[g]] <- slab[[g]] * pmax(alpha[, g] - threshold, 0)
    eta[[g]] <- matrix(0, N, K)
  }

  keep <- seq.int(burn + 1L, n_iter, by = thin)
  B <- length(keep)
  z_store <- matrix(0L, B, N)
  row_store <- array(0, c(P, G, B))
  norm_store <- array(0, c(P, G, B))
  q_store <- matrix(0L, B, G)
  eig_store <- array(0, c(K, G, B))
  psi_store <- array(0, c(P, G, B))
  mu_store <- array(0, c(P, G, B))
  tau_store <- array(0, c(K, G, B))
  sigma_store <- matrix(0, B, G)
  mix_store <- matrix(0, B, G)
  omega_store <- threshold_store <- numeric(B)
  save <- 0L

  for (iter in seq_len(n_iter)) {
    if (update_z) {
      logp <- matrix(0, N, G)

      for (g in seq_len(G)) {
        centered <- sweep(X, 2, mu[, g], "-")
        Dinv <- 1 / psi[, g]
        XD <- centered * rep(Dinv, each = N)
        C <- diag(K) + crossprod(Lambda[[g]] * Dinv, Lambda[[g]])
        R <- chol(C)
        Bmat <- XD %*% Lambda[[g]]
        quad <- rowSums(centered * XD) -
          rowSums((Bmat %*% chol2inv(R)) * Bmat)
        logdet <- sum(log(psi[, g])) + 2 * sum(log(diag(R)))

        logp[, g] <- log(mix[g]) -
          0.5 * (P * log(2 * pi) + logdet + quad)
      }

      prob <- exp(logp - apply(logp, 1, max))
      prob <- prob / rowSums(prob)
      z <- apply(prob, 1, function(x) sample.int(G, 1, prob = x))
    }

    if (!equal_weights) {
      mix <- rgamma(G, mix_alpha + tabulate(z, G))
      mix <- mix / sum(mix)
    }

    for (g in seq_len(G)) {
      ind <- which(z == g)
      ng <- length(ind)
      if (!ng) next

      C <- diag(K) + crossprod(Lambda[[g]] / psi[, g], Lambda[[g]])
      V <- chol2inv(chol(C))
      centered <- sweep(X[ind, , drop = FALSE], 2, mu[, g], "-")
      M <- centered %*% (Lambda[[g]] / psi[, g]) %*% V
      eta[[g]][ind, ] <- M + matrix(rnorm(ng * K), ng, K) %*% chol(V)
      H <- eta[[g]][ind, , drop = FALSE]
      HtH <- crossprod(H)

      for (j in seq_len(P)) {
        y <- X[ind, j] - mu[j, g]
        activation <- pmax(alpha[j, g] - threshold, 0)
        prior <- omega * sigma[g] * phi[[g]][j, ] * tau[[g]]
        R <- chol(diag(prior, K) + activation^2 * HtH / psi[j, g])
        mean_b <- chol2inv(R) %*%
          (activation * crossprod(H, y) / psi[j, g])
        slab[[g]][j, ] <- as.numeric(mean_b + backsolve(R, rnorm(K)))

        loglik_alpha <- function(value) {
          activation <- pmax(value - threshold, 0)
          -sum((y - H %*% (activation * slab[[g]][j, ]))^2) /
            (2 * psi[j, g])
        }
        alpha[j, g] <- ESS.Gibbs(
          alpha[j, g], loglik_alpha, max_steps = 200L
        )
        Lambda[[g]][j, ] <-
          pmax(alpha[j, g] - threshold, 0) * slab[[g]][j, ]
      }

      active_rows <- alpha[, g] > threshold
      m <- sum(active_rows)

      phi[[g]] <- matrix(rgamma(P * K, nu1, nu2), P, K)
      if (m > 0L) {
        rate_phi <- nu2 + 0.5 * omega * sigma[g] *
          slab[[g]][active_rows, , drop = FALSE]^2 *
          matrix(rep(tau[[g]], each = m), m, K)
        phi[[g]][active_rows, ] <- matrix(
          rgamma(m * K, nu1 + 0.5, rate_phi), m, K
        )
      }

      for (h in seq_len(K)) {
        tau[[g]] <- cumprod(delta[[g]])
        columns <- h:K
        shape <- if (h == 1L) alpha1 else alpha2
        rate <- if (h == 1L) beta1 else beta2

        if (m > 0L) {
          energy <- colSums(
            phi[[g]][active_rows, columns, drop = FALSE] *
              slab[[g]][active_rows, columns, drop = FALSE]^2
          )
          rate <- rate + 0.5 * omega * sigma[g] * sum(
            tau[[g]][columns] / delta[[g]][h] * energy
          )
        }

        shape <- shape + m * length(columns) / 2
        if (h == 1L) {
          delta[[g]][h] <- rgamma(1, shape, rate)
        } else {
          log_tail <- pgamma(
            1, shape, rate = rate, lower.tail = FALSE, log.p = TRUE
          )
          delta[[g]][h] <- qgamma(
            log(runif(1)) + log_tail,
            shape, rate = rate, lower.tail = FALSE, log.p = TRUE
          )
        }
      }

      tau[[g]] <- cumprod(delta[[g]])
      energy <- if (m > 0L) {
        sum(
          phi[[g]][active_rows, , drop = FALSE] *
            slab[[g]][active_rows, , drop = FALSE]^2 *
            matrix(rep(tau[[g]], each = m), m, K)
        )
      } else 0
      sigma[g] <- rgamma(
        1, rho1 + m * K / 2, rho2 + 0.5 * omega * energy
      )

      if (m < P) {
        precision <- phi[[g]][!active_rows, , drop = FALSE] *
          matrix(rep(omega * sigma[g] * tau[[g]], each = P - m), P - m, K)
        slab[[g]][!active_rows, ] <-
          matrix(rnorm((P - m) * K), P - m, K) / sqrt(precision)
        Lambda[[g]][!active_rows, ] <- 0
      }

      fitted <- H %*% t(Lambda[[g]])
      residual <- X[ind, , drop = FALSE] - fitted
      precision_mu <- ng / psi[, g] + 1 / mu_sd^2
      mean_mu <- colSums(residual) / psi[, g] / precision_mu
      mu[, g] <- rnorm(P, mean_mu, sqrt(1 / precision_mu))

      residual <- sweep(X[ind, , drop = FALSE], 2, mu[, g], "-") - fitted
      psi[, g] <- 1 / rgamma(
        P, a_psi + ng / 2, psi_rate + colSums(residual^2) / 2
      )
    }

    psi_rate <- rgamma(
      P, bpsi_shape + G * a_psi,
      bpsi_rate + rowSums(1 / psi)
    )

    omega_shape <- omega1
    omega_rate <- omega2
    for (g in seq_len(G)) {
      active_rows <- alpha[, g] > threshold
      m <- sum(active_rows)
      if (!m) next

      energy <- sum(
        phi[[g]][active_rows, , drop = FALSE] *
          slab[[g]][active_rows, , drop = FALSE]^2 *
          matrix(rep(tau[[g]], each = m), m, K)
      )
      omega_shape <- omega_shape + m * K / 2
      omega_rate <- omega_rate + 0.5 * sigma[g] * energy
    }
    omega <- rgamma(1, omega_shape, omega_rate)

    loglik_threshold <- function(value) {
      current <- threshold_mean + threshold_sd * value
      out <- 0

      for (g in seq_len(G)) {
        ind <- which(z == g)
        if (!length(ind)) next
        activation <- pmax(alpha[, g] - current, 0)
        L <- slab[[g]] * activation
        fitted <- eta[[g]][ind, , drop = FALSE] %*% t(L)
        residual <- sweep(X[ind, , drop = FALSE], 2, mu[, g], "-") - fitted
        out <- out - sum(sweep(residual^2, 2, psi[, g], "/")) / 2
      }

      out
    }

    threshold_z <- ESS.Gibbs(
      threshold_z, loglik_threshold, max_steps = 300L
    )
    threshold <- threshold_mean + threshold_sd * threshold_z
    for (g in seq_len(G)) {
      Lambda[[g]] <- slab[[g]] * pmax(alpha[, g] - threshold, 0)
    }

    if (iter %in% keep) {
      save <- save + 1L
      z_store[save, ] <- z
      sigma_store[save, ] <- sigma
      mix_store[save, ] <- mix
      omega_store[save] <- omega
      threshold_store[save] <- threshold
      psi_store[, , save] <- psi
      mu_store[, , save] <- mu

      for (g in seq_len(G)) {
        row_store[, g, save] <- alpha[, g] > threshold
        norm_store[, g, save] <- sqrt(rowSums(Lambda[[g]]^2))
        tau_store[, g, save] <- tau[[g]]

        values <- svd(Lambda[[g]] / sqrt(psi[, g]), nu = 0, nv = 0)$d^2
        values <- c(values, rep(0, K - length(values)))
        eig_store[, g, save] <- values
        signal <- values[values > rank_abs]
        q_store[save, g] <- if (!length(signal)) {
          0L
        } else {
          which(cumsum(signal) / sum(signal) >= rank_prop)[1]
        }
      }
    }

    if (verbose && iter %% 100L == 0L) {
      cat(
        "iteration", iter,
        "cluster sizes", tabulate(z, G),
        "active rows", colSums(alpha > threshold),
        "omega", round(omega, 3), "\n"
      )
    }
  }

  post <- relabel_mifa(
    z_store, row_store, norm_store, q_store, eig_store,
    psi_store, mu_store, tau_store, sigma_store, mix_store, z
  )
  order_final <- order(post$final_map)

  z <- post$final_map[z]
  mu <- mu[, order_final, drop = FALSE]
  psi <- psi[, order_final, drop = FALSE]
  alpha <- alpha[, order_final, drop = FALSE]
  mix <- mix[order_final]
  sigma <- sigma[order_final]
  slab <- slab[order_final]
  Lambda <- Lambda[order_final]
  phi <- phi[order_final]
  delta <- delta[order_final]
  tau <- tau[order_final]
  eta <- eta[order_final]

  z_hat_reference <- apply(post$z_store, 2, function(x) {
    as.integer(names(which.max(table(x))))
  })
  truth_map <- seq_len(G)
  z_true <- NULL

  if (!is.null(truth) && !is.null(truth$z_true)) {
    z_true <- as.integer(factor(
      truth$z_true, levels = sort(unique(truth$z_true))
    ))
    counts <- table(
      factor(z_hat_reference, levels = seq_len(G)),
      factor(z_true, levels = seq_len(G))
    )
    truth_map <- as.integer(clue::solve_LSAP(counts, maximum = TRUE))
    truth_order <- order(truth_map)

    post$z_store <- matrix(truth_map[post$z_store], B, N)
    post$row_store <- post$row_store[, truth_order, , drop = FALSE]
    post$norm_store <- post$norm_store[, truth_order, , drop = FALSE]
    post$q_store <- post$q_store[, truth_order, drop = FALSE]
    post$eig_store <- post$eig_store[, truth_order, , drop = FALSE]
    post$psi_store <- post$psi_store[, truth_order, , drop = FALSE]
    post$mu_store <- post$mu_store[, truth_order, , drop = FALSE]
    post$tau_store <- post$tau_store[, truth_order, , drop = FALSE]
    post$sigma_store <- post$sigma_store[, truth_order, drop = FALSE]
    post$mix_store <- post$mix_store[, truth_order, drop = FALSE]
    post$z_ref <- truth_map[post$z_ref]

    z <- truth_map[z]
    mu <- mu[, truth_order, drop = FALSE]
    psi <- psi[, truth_order, drop = FALSE]
    alpha <- alpha[, truth_order, drop = FALSE]
    mix <- mix[truth_order]
    sigma <- sigma[truth_order]
    slab <- slab[truth_order]
    Lambda <- Lambda[truth_order]
    phi <- phi[truth_order]
    delta <- delta[truth_order]
    tau <- tau[truth_order]
    eta <- eta[truth_order]
  }

  row_pip <- apply(post$row_store, c(1, 2), mean)
  row_norm <- apply(post$norm_store, c(1, 2), mean)
  psi_mean <- apply(post$psi_store, c(1, 2), mean)
  mu_mean <- apply(post$mu_store, c(1, 2), mean)
  eigen_mean <- apply(post$eig_store, c(1, 2), mean)
  tau_mean <- apply(post$tau_store, c(1, 2), mean)
  z_hat <- apply(post$z_store, 2, function(x) {
    as.integer(names(which.max(table(x))))
  })
  q_hat <- apply(post$q_store, 2, function(x) {
    as.integer(names(which.max(table(x))))
  })

  cluster_names <- paste0("cluster_", seq_len(G))
  rownames(row_pip) <- rownames(row_norm) <- rownames(psi_mean) <-
    rownames(mu_mean) <- feature_names
  colnames(row_pip) <- colnames(row_norm) <- colnames(psi_mean) <-
    colnames(mu_mean) <- cluster_names
  rownames(eigen_mean) <- rownames(tau_mean) <- paste0("factor_", seq_len(K))
  colnames(eigen_mean) <- colnames(tau_mean) <- cluster_names
  colnames(post$q_store) <- colnames(post$sigma_store) <-
    colnames(post$mix_store) <- cluster_names

  q_summary <- data.frame(
    cluster = seq_len(G),
    q_mean = colMeans(post$q_store),
    q_median = apply(post$q_store, 2, median),
    q_mode = q_hat
  )
  if (!is.null(truth) && !is.null(truth$q)) {
    q_summary$q_true <- truth$q
    q_summary$correct <- q_summary$q_mode == truth$q
  }

  confusion <- NULL
  RI <- ARI <- accuracy <- NA_real_
  if (!is.null(z_true)) {
    confusion <- table(
      truth = factor(z_true, levels = seq_len(G)),
      estimate = factor(z_hat, levels = seq_len(G))
    )
    same_both <- sum(confusion * (confusion - 1) / 2)
    same_truth <- sum(rowSums(confusion) * (rowSums(confusion) - 1) / 2)
    same_estimate <- sum(colSums(confusion) * (colSums(confusion) - 1) / 2)
    total_pairs <- N * (N - 1) / 2
    different_both <- total_pairs - same_truth - same_estimate + same_both
    expected <- same_truth * same_estimate / total_pairs
    maximum <- (same_truth + same_estimate) / 2
    RI <- (same_both + different_both) / total_pairs
    ARI <- (same_both - expected) / (maximum - expected)
    accuracy <- mean(z_hat == z_true)
  }

  gamma_true <- if (!is.null(truth) && !is.null(truth$gamma_true)) {
    as.matrix(truth$gamma_true) != 0
  } else NULL

  selected_rows <- setNames(vector("list", G), cluster_names)
  true_active_rows <- selected_true_rows <- false_positive_rows <-
    false_negative_rows <- setNames(vector("list", G), cluster_names)
  row_summary <- vector("list", G)

  for (g in seq_len(G)) {
    selected <- row_pip[, g] >= pip_cut
    selected_rows[[g]] <- which(selected)
    has_truth <- !is.null(gamma_true)
    active <- if (has_truth) gamma_true[, g] else selected
    noise <- !active
    q_reference <- if (!is.null(truth) && !is.null(truth$q)) {
      truth$q[g]
    } else q_hat[g]
    total_energy <- sum(eigen_mean[, g])
    tail_energy <- if (q_reference >= K || total_energy == 0) {
      0
    } else {
      sum(eigen_mean[(q_reference + 1L):K, g]) / total_energy
    }

    TP <- TN <- FP <- FN <- NA_integer_
    selection_accuracy <- precision <- recall <- specificity <- F1 <- NA_real_
    if (has_truth) {
      TP <- sum(selected & active)
      TN <- sum(!selected & noise)
      FP <- sum(selected & noise)
      FN <- sum(!selected & active)
      selection_accuracy <- (TP + TN) / P
      precision <- if (TP + FP == 0) NA_real_ else TP / (TP + FP)
      recall <- if (TP + FN == 0) NA_real_ else TP / (TP + FN)
      specificity <- if (TN + FP == 0) NA_real_ else TN / (TN + FP)
      F1 <- if (2 * TP + FP + FN == 0) {
        NA_real_
      } else {
        2 * TP / (2 * TP + FP + FN)
      }
      true_active_rows[[g]] <- which(active)
      selected_true_rows[[g]] <- which(selected & active)
      false_positive_rows[[g]] <- which(selected & noise)
      false_negative_rows[[g]] <- which(!selected & active)
    }

    row_summary[[g]] <- data.frame(
      cluster = g,
      selected_n = sum(selected),
      true_active_n = if (has_truth) sum(active) else NA_integer_,
      TP = TP, TN = TN, FP = FP, FN = FN,
      selection_accuracy = selection_accuracy,
      precision = precision,
      recall = recall,
      specificity = specificity,
      F1 = F1,
      active_PIP = if (any(active)) mean(row_pip[active, g]) else NA_real_,
      noise_PIP = if (any(noise)) mean(row_pip[noise, g]) else NA_real_,
      active_norm = if (any(active)) mean(row_norm[active, g]) else NA_real_,
      noise_norm = if (any(noise)) mean(row_norm[noise, g]) else NA_real_,
      active_psi = if (any(active)) mean(psi_mean[active, g]) else NA_real_,
      noise_psi = if (any(noise)) mean(psi_mean[noise, g]) else NA_real_,
      tail_energy = tail_energy
    )
  }
  row_summary <- do.call(rbind, row_summary)

  omega_summary <- data.frame(
    mean = mean(omega_store),
    sd = sd(omega_store),
    q025 = unname(quantile(omega_store, 0.025)),
    median = median(omega_store),
    q975 = unname(quantile(omega_store, 0.975)),
    final = omega
  )
  threshold_summary <- data.frame(
    mean = mean(threshold_store),
    sd = sd(threshold_store),
    q025 = unname(quantile(threshold_store, 0.025)),
    median = median(threshold_store),
    q975 = unname(quantile(threshold_store, 0.975)),
    final = threshold
  )
  sigma_summary <- data.frame(
    cluster = seq_len(G),
    mean = colMeans(post$sigma_store),
    sd = apply(post$sigma_store, 2, sd),
    q025 = apply(post$sigma_store, 2, quantile, 0.025),
    median = apply(post$sigma_store, 2, median),
    q975 = apply(post$sigma_store, 2, quantile, 0.975),
    final = sigma
  )

  tau_final <- delta_final <- matrix(0, K, G)
  for (g in seq_len(G)) {
    tau_final[, g] <- tau[[g]]
    delta_final[, g] <- delta[[g]]
  }
  rownames(tau_final) <- rownames(delta_final) <- paste0("factor_", seq_len(K))
  colnames(tau_final) <- colnames(delta_final) <- cluster_names

  Sigma_final <- column_diagnostics <- vector("list", G)
  covariance_summary <- vector("list", G)
  names(Sigma_final) <- names(column_diagnostics) <- cluster_names

  for (g in seq_len(G)) {
    Sigma_final[[g]] <- tcrossprod(Lambda[[g]]) + diag(psi[, g])
    active <- if (!is.null(gamma_true)) {
      gamma_true[, g]
    } else {
      row_pip[, g] >= pip_cut
    }
    fitted_block <- Sigma_final[[g]][active, active, drop = FALSE]
    fitted_offdiag <- if (sum(active) > 1L) {
      mean(abs(fitted_block[row(fitted_block) != col(fitted_block)]))
    } else NA_real_
    fitted_diagonal <- if (sum(active) > 0L) {
      mean(diag(fitted_block))
    } else NA_real_

    true_eigen <- true_block <- NULL
    true_offdiag <- true_diagonal <- NA_real_
    if (!is.null(truth) && !is.null(truth$Lambda) && !is.null(truth$psi)) {
      psi_true <- if (length(truth$psi) == 1L) {
        rep(truth$psi, P)
      } else if (is.matrix(truth$psi)) {
        truth$psi[, g]
      } else {
        truth$psi
      }
      true_eigen <- svd(
        truth$Lambda[[g]] / sqrt(psi_true), nu = 0, nv = 0
      )$d^2
    }
    if (!is.null(truth) && !is.null(truth$Sigma) && any(active)) {
      true_block <- truth$Sigma[[g]][active, active, drop = FALSE]
      true_offdiag <- if (sum(active) > 1L) {
        mean(abs(true_block[row(true_block) != col(true_block)]))
      } else NA_real_
      true_diagonal <- mean(diag(true_block))
    }

    column_diagnostics[[g]] <- list(
      q_mean = q_summary$q_mean[g],
      q_median = q_summary$q_median[g],
      q_mode = q_summary$q_mode[g],
      q_true = if (!is.null(truth) && !is.null(truth$q)) truth$q[g] else NULL,
      tail_energy = row_summary$tail_energy[g],
      fitted_eigenvalues = eigen_mean[, g],
      true_eigenvalues = true_eigen,
      tau_posterior_mean = tau_mean[, g],
      tau_final = tau[[g]],
      active_covariance_fitted = fitted_block,
      active_covariance_true = true_block
    )
    covariance_summary[[g]] <- data.frame(
      cluster = g,
      mean_abs_offdiag_true = true_offdiag,
      mean_abs_offdiag_fitted = fitted_offdiag,
      mean_diagonal_true = true_diagonal,
      mean_diagonal_fitted = fitted_diagonal
    )
  }
  covariance_summary <- do.call(rbind, covariance_summary)

  truth_info <- NULL
  if (!is.null(truth)) {
    truth_info <- list(
      cluster_sizes = if (!is.null(z_true)) tabulate(z_true, G) else NULL,
      q = truth$q,
      active_rows = if (!is.null(gamma_true)) {
        lapply(seq_len(G), function(g) which(gamma_true[, g]))
      } else NULL,
      loadings = truth$loadings,
      psi = truth$psi,
      covariance_distance = truth$covariance_distance,
      seed = truth$seed
    )
  }

  list(
    model = "Group-ReLU MIFA",
    call = match.call(),
    settings = list(
      G = G, K = K, n_iter = n_iter, burn = burn, thin = thin,
      retained_draws = B, pip_cut = pip_cut,
      rank_prop = rank_prop, rank_abs = rank_abs,
      equal_weights = equal_weights,
      nu1 = nu1, nu2 = nu2,
      alpha1 = alpha1, beta1 = beta1,
      alpha2 = alpha2, beta2 = beta2,
      rho1 = rho1, rho2 = rho2,
      omega1 = omega1, omega2 = omega2,
      a_psi = a_psi, bpsi_shape = bpsi_shape, bpsi_rate = bpsi_rate,
      mu_sd = mu_sd, mix_alpha = mix_alpha,
      threshold_mean = threshold_mean, threshold_sd = threshold_sd,
      seed = seed
    ),
    data_info = list(
      N = N, P = P, G = G, K = K,
      feature_names = feature_names,
      initial_cluster_sizes = tabulate(z_initial, G),
      estimated_cluster_sizes = tabulate(z_hat, G),
      final_state_cluster_sizes = tabulate(z, G)
    ),
    truth_info = truth_info,
    clustering = list(
      RI = RI, ARI = ARI, accuracy = accuracy,
      confusion_matrix = confusion,
      estimated_cluster_sizes = tabulate(z_hat, G),
      true_cluster_sizes = if (!is.null(z_true)) tabulate(z_true, G) else NULL
    ),
    factor_summary = q_summary,
    row_summary = row_summary,
    selected_rows = selected_rows,
    true_active_rows = if (!is.null(gamma_true)) true_active_rows else NULL,
    selected_true_rows = if (!is.null(gamma_true)) selected_true_rows else NULL,
    false_positive_rows = if (!is.null(gamma_true)) false_positive_rows else NULL,
    false_negative_rows = if (!is.null(gamma_true)) false_negative_rows else NULL,
    global_shrinkage = list(
      omega = omega_summary,
      threshold = threshold_summary,
      sigma = sigma_summary,
      tau_posterior_mean = tau_mean,
      tau_final = tau_final,
      delta_final = delta_final
    ),
    column_diagnostics = column_diagnostics,
    covariance_summary = covariance_summary,
    row_pip = row_pip,
    row_norm = row_norm,
    selected = row_pip >= pip_cut,
    mu_mean = mu_mean,
    psi_mean = psi_mean,
    eigen_mean = eigen_mean,
    tau_mean = tau_mean,
    z_store = post$z_store,
    row_store = post$row_store,
    norm_store = post$norm_store,
    q_store = post$q_store,
    eigen_store = post$eig_store,
    psi_store = post$psi_store,
    mu_store = post$mu_store,
    tau_store = post$tau_store,
    sigma_store = post$sigma_store,
    mix_store = post$mix_store,
    omega_store = omega_store,
    threshold_store = threshold_store,
    z_hat = z_hat,
    z = z,
    z_ref = post$z_ref,
    psm = post$psm,
    label_maps = post$maps,
    truth_map = truth_map,
    q_hat = q_hat,
    Sigma = Sigma_final,
    Lambda = Lambda,
    slab = slab,
    phi = phi,
    delta = delta,
    tau = tau,
    sigma = sigma,
    omega = omega,
    psi_rate = psi_rate,
    alpha = alpha,
    activation = pmax(alpha - threshold, 0),
    alpha0 = threshold,
    mu = mu,
    psi = psi,
    mix = mix,
    eta = eta
  )
}

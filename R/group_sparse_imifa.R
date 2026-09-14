# Group-ReLU spectral IMIFA with a Pitman-Yor process over clusters.
# The only package dependency is clue, used for square-assignment relabeling.

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

  G_draw <- apply(z_store, 1, function(x) length(unique(x)))
  G_mode <- as.integer(names(which.max(table(G_draw))))
  conditioned <- which(G_draw == G_mode)

  psm <- matrix(0, N, N)
  for (b in conditioned) {
    psm <- psm + outer(z_store[b, ], z_store[b, ], "==")
  }
  psm <- psm / length(conditioned)

  loss <- vapply(conditioned, function(b) {
    sum((outer(z_store[b, ], z_store[b, ], "==") - psm)^2)
  }, numeric(1))
  reference <- z_store[conditioned[which.min(loss)], ]
  reference <- match(reference, sort(unique(reference)))
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
    G_mode = G_mode,
    conditioned = conditioned,
    psm = psm,
    maps = maps,
    final_map = as.integer(clue::solve_LSAP(counts, maximum = TRUE))
  )
}


# G_init seeds an over-partition; G_max is only the slice/storage ceiling.
# Smaller pyp_concentration favours fewer occupied clusters; positive
# pyp_discount gives the PYP its heavier tail of small clusters.
group_sparse_imifa <- function(
  X, K = round(3 * log(ncol(X))),
  n_iter = 6000, burn = 2000, thin = 2,
  G_init = min(max(ceiling(3 * log(nrow(X))), 25L), nrow(X) - 1L),
  G_max = max(G_init, min(nrow(X) - 1L, 50L)),
  pyp_concentration = 0.5, pyp_discount = 0.1, slice_rho = 0.75,
  nu1 = 3, nu2 = 2,
  alpha1 = 2.1, beta1 = 1, alpha2 = 3.1, beta2 = 1,
  rho1 = 3, rho2 = 2, omega1 = 3, omega2 = 2,
  a_psi = 20, bpsi_shape = 38, bpsi_rate = 2,
  mu_sd = 1,
  rank_prop = 0.99, rank_abs = 1.5,
  threshold_mean = 1.28, threshold_sd = 0.5,
  spectral_init = TRUE,
  z_init = NULL,
  truth = NULL, pip_cut = 0.5,
  seed = 200, verbose = TRUE
) {
  set.seed(seed)

  X <- as.matrix(X)
  N <- nrow(X)
  P <- ncol(X)
  K <- min(as.integer(K), P - 1L, N - 1L)
  G_init <- as.integer(G_init)
  G_max <- as.integer(G_max)
  stopifnot(
    K >= 1L, n_iter > burn, thin >= 1L,
    G_init >= 1L, G_max >= G_init, G_max < N,
    pyp_discount >= 0, pyp_discount < 1,
    pyp_concentration > -pyp_discount,
    slice_rho > 0, slice_rho < 1
  )

  feature_names <- colnames(X)
  if (is.null(feature_names)) feature_names <- paste0("V", seq_len(P))

  z <- if (is.null(z_init)) {
    kmeans(X, centers = G_init, nstart = 10)$cluster
  } else {
    as.integer(factor(z_init, levels = sort(unique(z_init))))
  }
  z_initial <- z
  stopifnot(max(z) <= G_max)

  mu <- matrix(rnorm(P * G_max, 0, mu_sd), P, G_max)
  for (g in sort(unique(z))) {
    mu[, g] <- colMeans(X[z == g, , drop = FALSE])
  }
  psi_rate <- rep(bpsi_shape / bpsi_rate, P)
  psi <- matrix(
    1 / rgamma(P * G_max, a_psi, rep(psi_rate, G_max)), P, G_max
  )
  psi[, sort(unique(z))] <- 1

  threshold_z <- 0
  threshold <- threshold_mean
  alpha <- matrix(rnorm(P * G_max), P, G_max)
  alpha[, sort(unique(z))] <- threshold + 1

  omega <- rgamma(1, omega1, omega2)
  sigma <- rgamma(G_max, rho1, rho2)
  slab <- phi <- delta <- tau <- Lambda <- eta <- vector("list", G_max)
  # Independent slice sequence from Kalli et al. (2011), as used by IMIFA.
  xi <- (1 - slice_rho) * slice_rho^(seq_len(G_max) - 1L)

  log_tail_delta <- pgamma(
    1, alpha2, rate = beta2, lower.tail = FALSE, log.p = TRUE
  )

  for (g in seq_len(G_max)) {
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

  if (spectral_init) {
    for (g in sort(unique(z))) {
      ind <- which(z == g)
      if (length(ind) < 2L) next
      centered <- sweep(X[ind, , drop = FALSE], 2, mu[, g], "-")
      S <- crossprod(centered) / (length(ind) - 1)
      eig <- eigen(S, symmetric = TRUE)
      noise <- median(diag(S))
      values <- pmax(eig$values[seq_len(K)] - noise, 0)
      values[eig$values[seq_len(K)] <=
        noise * (1 + sqrt(P / (length(ind) - 1)))^2] <- 0
      Lambda[[g]] <- sweep(
        eig$vectors[, seq_len(K), drop = FALSE], 2, sqrt(values), "*"
      )
      alpha[, g] <- threshold + 1
      slab[[g]] <- Lambda[[g]]
      psi[, g] <- noise
    }
  }

  keep <- seq.int(burn + 1L, n_iter, by = thin)
  B <- length(keep)
  z_store <- matrix(0L, B, N)
  row_store <- array(NA, c(P, G_max, B))
  norm_store <- array(NA_real_, c(P, G_max, B))
  q_store <- matrix(NA_integer_, B, G_max)
  eig_store <- array(NA_real_, c(K, G_max, B))
  psi_store <- array(NA_real_, c(P, G_max, B))
  mu_store <- array(NA_real_, c(P, G_max, B))
  tau_store <- array(NA_real_, c(K, G_max, B))
  sigma_store <- matrix(NA_real_, B, G_max)
  mix_store <- stick_store <- matrix(0, B, G_max)
  slice_store <- integer(B)
  slice_min_store <- tail_mass_store <- numeric(B)
  omega_store <- threshold_store <- numeric(B)
  save <- 0L

  for (iter in seq_len(n_iter)) {
    # PYP stick posterior: v_g | z ~ Beta(1-d+n_g, a+gd+n_{>g}).
    nn <- tabulate(z, G_max)
    stick <- rbeta(
      G_max,
      1 - pyp_discount + nn,
      pyp_concentration + seq_len(G_max) * pyp_discount + N - cumsum(nn)
    )
    mix <- stick * c(1, cumprod(1 - stick[-G_max]))
    u_slice <- runif(N, 0, xi[z])
    G_active <- min(G_max, sum(min(u_slice) < xi))

    # Occupied components retain the original within-cluster Gibbs updates.
    # Empty slice proposals are redrawn from the same hierarchical base prior.
    for (g in seq_len(G_active)) {
      ind <- which(z == g)
      ng <- length(ind)
      if (!ng) {
        mu[, g] <- rnorm(P, 0, mu_sd)
        psi[, g] <- 1 / rgamma(P, a_psi, psi_rate)
        alpha[, g] <- rnorm(P)
        sigma[g] <- rgamma(1, rho1, rho2)
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
        next
      }

      C <- diag(K) + crossprod(Lambda[[g]] / psi[, g], Lambda[[g]])
      V <- chol2inv(chol(C))
      centered <- sweep(X[ind, , drop = FALSE], 2, mu[, g], "-")
      M <- centered %*% (Lambda[[g]] / psi[, g]) %*% V
      eta[[g]][ind, ] <- M + matrix(rnorm(ng * K), ng, K) %*% chol(V)
      H <- eta[[g]][ind, , drop = FALSE]
      HtH <- crossprod(H)

      for (j in seq_len(P)) {
        y <- X[ind, j] - mu[j, g]
        prior <- omega * sigma[g] * phi[[g]][j, ] * tau[[g]]
        Hty <- crossprod(H, y)

        loglik_alpha <- function(value) {
          activation <- pmax(value - threshold, 0)
          R <- chol(diag(prior, K) + activation^2 * HtH / psi[j, g])
          quad <- sum(y^2) / psi[j, g] - activation^2 / psi[j, g]^2 *
            drop(crossprod(Hty, chol2inv(R) %*% Hty))
          logdet <- ng * log(psi[j, g]) + 2 * sum(log(diag(R))) -
            sum(log(prior))
          -(logdet + quad) / 2
        }
        alpha[j, g] <- ESS.Gibbs(
          alpha[j, g], loglik_alpha, max_steps = 200L
        )
        activation <- pmax(alpha[j, g] - threshold, 0)
        R <- chol(diag(prior, K) + activation^2 * HtH / psi[j, g])
        mean_b <- chol2inv(R) %*% (activation * Hty / psi[j, g])
        slab[[g]][j, ] <- as.numeric(mean_b + backsolve(R, rnorm(K)))
        Lambda[[g]][j, ] <- activation * slab[[g]][j, ]
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

    occupied <- which(nn > 0)
    psi_rate <- rgamma(
      P, bpsi_shape + length(occupied) * a_psi,
      bpsi_rate + rowSums(1 / psi[, occupied, drop = FALSE])
    )

    omega_shape <- omega1
    omega_rate <- omega2
    for (g in occupied) {
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

      for (g in occupied) {
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
    for (g in seq_len(G_active)) {
      Lambda[[g]] <- slab[[g]] * pmax(alpha[, g] - threshold, 0)
    }

    # p(z_i=g | u_i, -) is proportional to pi_g / xi_g times the MFA density.
    logp <- matrix(-Inf, N, G_active)

    for (g in seq_len(G_active)) {
      eligible <- which(u_slice < xi[g])
      centered <- sweep(X[eligible, , drop = FALSE], 2, mu[, g], "-")
      Dinv <- 1 / psi[, g]
      XD <- centered * rep(Dinv, each = length(eligible))
      C <- diag(K) + crossprod(Lambda[[g]] * Dinv, Lambda[[g]])
      R <- chol(C)
      Bmat <- XD %*% Lambda[[g]]
      quad <- rowSums(centered * XD) -
        rowSums((Bmat %*% chol2inv(R)) * Bmat)
      logdet <- sum(log(psi[, g])) + 2 * sum(log(diag(R)))
      logp[eligible, g] <- log(mix[g]) - log(xi[g]) -
        0.5 * (P * log(2 * pi) + logdet + quad)
    }

    z <- max.col(
      logp - log(-log(matrix(runif(N * G_active), N, G_active))),
      ties.method = "first"
    )

    if (iter %in% keep) {
      save <- save + 1L
      z_store[save, ] <- z
      mix_store[save, ] <- mix
      stick_store[save, ] <- stick
      omega_store[save] <- omega
      threshold_store[save] <- threshold
      slice_store[save] <- G_active
      slice_min_store[save] <- min(u_slice)
      tail_mass_store[save] <- 1 - sum(mix)
      occupied <- which(tabulate(z, G_max) > 0)

      for (g in occupied) {
        row_store[, g, save] <- alpha[, g] > threshold
        norm_store[, g, save] <- sqrt(rowSums(Lambda[[g]]^2))
        psi_store[, g, save] <- psi[, g]
        mu_store[, g, save] <- mu[, g]
        tau_store[, g, save] <- tau[[g]]
        sigma_store[save, g] <- sigma[g]

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
      sizes <- tabulate(z, G_max)
      cat(
        "iteration", iter,
        "occupied", sum(sizes > 0),
        "slice proposals", G_active,
        "cluster sizes", sizes[sizes > 0],
        "omega", round(omega, 3), "\n"
      )
    }
  }

  post <- relabel_mifa(
    z_store, row_store, norm_store, q_store, eig_store,
    psi_store, mu_store, tau_store, sigma_store, mix_store, z
  )
  stick_final_raw <- stick
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

  z_hat_reference <- post$z_ref
  truth_map <- seq_len(G_max)
  z_true <- NULL
  G_true <- NULL

  if (!is.null(truth) && !is.null(truth$z_true)) {
    z_true <- as.integer(factor(
      truth$z_true, levels = sort(unique(truth$z_true))
    ))
    G_true <- length(unique(z_true))
    counts <- matrix(0, G_max, G_max)
    counts[, seq_len(G_true)] <- table(
      factor(z_hat_reference, levels = seq_len(G_max)),
      factor(z_true, levels = seq_len(G_true))
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

  # Within-cluster summaries use only draws at modal G where the component exists.
  size_store <- t(vapply(seq_len(B), function(b) {
    tabulate(post$z_store[b, ], G_max)
  }, integer(G_max)))
  occupied_store <- size_store > 0
  G_store <- rowSums(occupied_store)
  G_mode <- post$G_mode
  inferential <- G_store == G_mode

  row_pip <- row_norm <- psi_mean <- mu_mean <- matrix(
    NA_real_, P, G_max
  )
  eigen_mean <- tau_mean <- matrix(NA_real_, K, G_max)
  q_hat_all <- rep(NA_integer_, G_max)

  for (g in seq_len(G_max)) {
    retained <- which(inferential & occupied_store[, g])
    if (!length(retained)) next
    row_pip[, g] <- rowMeans(matrix(
      post$row_store[, g, retained], P, length(retained)
    ))
    row_norm[, g] <- rowMeans(matrix(
      post$norm_store[, g, retained], P, length(retained)
    ))
    psi_mean[, g] <- rowMeans(matrix(
      post$psi_store[, g, retained], P, length(retained)
    ))
    mu_mean[, g] <- rowMeans(matrix(
      post$mu_store[, g, retained], P, length(retained)
    ))
    eigen_mean[, g] <- rowMeans(matrix(
      post$eig_store[, g, retained], K, length(retained)
    ))
    tau_mean[, g] <- rowMeans(matrix(
      post$tau_store[, g, retained], K, length(retained)
    ))
    q_values <- post$q_store[retained, g]
    q_hat_all[g] <- as.integer(names(which.max(table(q_values))))
  }

  z_hat <- post$z_ref
  reported <- sort(unique(z_hat))
  cluster_names <- paste0("cluster_", seq_len(G_max))
  rownames(row_pip) <- rownames(row_norm) <- rownames(psi_mean) <-
    rownames(mu_mean) <- feature_names
  colnames(row_pip) <- colnames(row_norm) <- colnames(psi_mean) <-
    colnames(mu_mean) <- cluster_names
  rownames(eigen_mean) <- rownames(tau_mean) <- paste0("factor_", seq_len(K))
  colnames(eigen_mean) <- colnames(tau_mean) <- cluster_names
  colnames(post$q_store) <- colnames(post$sigma_store) <-
    colnames(post$mix_store) <- cluster_names

  q_summary <- data.frame(
    cluster = reported,
    occupied_probability = colMeans(occupied_store)[reported],
    mean_size = colMeans(size_store)[reported],
    dahl_size = tabulate(z_hat, G_max)[reported],
    q_mean = vapply(reported, function(g) {
      mean(post$q_store[inferential & occupied_store[, g], g])
    }, numeric(1)),
    q_median = vapply(reported, function(g) {
      median(post$q_store[inferential & occupied_store[, g], g])
    }, numeric(1)),
    q_mode = q_hat_all[reported]
  )
  if (!is.null(truth) && !is.null(truth$q)) {
    q_summary$q_true <- ifelse(
      q_summary$cluster <= length(truth$q),
      truth$q[pmin(q_summary$cluster, length(truth$q))], NA
    )
    q_summary$correct <- q_summary$q_mode == q_summary$q_true
  }

  G_summary <- data.frame(
    mean = mean(G_store),
    sd = sd(G_store),
    q025 = unname(quantile(G_store, 0.025)),
    median = median(G_store),
    q975 = unname(quantile(G_store, 0.975)),
    mode = G_mode,
    dahl = length(reported),
    conditioned_draws = sum(inferential),
    truth = if (is.null(G_true)) NA_integer_ else G_true
  )
  G_distribution <- data.frame(
    G = as.integer(names(table(G_store))),
    probability = as.numeric(prop.table(table(G_store)))
  )
  component_summary <- data.frame(
    cluster = seq_len(G_max),
    occupied_probability = colMeans(occupied_store),
    occupied_probability_at_mode = colMeans(occupied_store[inferential, , drop = FALSE]),
    mean_size = colMeans(size_store),
    mean_size_when_occupied = vapply(seq_len(G_max), function(g) {
      if (any(occupied_store[, g])) {
        mean(size_store[occupied_store[, g], g])
      } else NA_real_
    }, numeric(1)),
    dahl_size = tabulate(z_hat, G_max),
    weight_mean = colMeans(post$mix_store)
  )

  confusion <- NULL
  RI <- ARI <- accuracy <- NA_real_
  if (!is.null(z_true)) {
    confusion <- table(
      truth = factor(z_true, levels = seq_len(G_true)),
      estimate = factor(z_hat, levels = reported)
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

  reported_names <- paste0("cluster_", reported)
  selected_rows <- setNames(vector("list", length(reported)), reported_names)
  true_active_rows <- selected_true_rows <- false_positive_rows <-
    false_negative_rows <- setNames(
      vector("list", length(reported)), reported_names
    )
  row_summary <- vector("list", length(reported))

  for (r in seq_along(reported)) {
    g <- reported[r]
    selected <- row_pip[, g] >= pip_cut
    selected_rows[[r]] <- which(selected)
    has_truth <- !is.null(gamma_true) && g <= ncol(gamma_true)
    active <- if (has_truth) gamma_true[, g] else selected
    noise <- !active
    q_reference <- if (!is.null(truth) && !is.null(truth$q) &&
                       g <= length(truth$q)) {
      truth$q[g]
    } else q_hat_all[g]
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
      true_active_rows[[r]] <- which(active)
      selected_true_rows[[r]] <- which(selected & active)
      false_positive_rows[[r]] <- which(selected & noise)
      false_negative_rows[[r]] <- which(!selected & active)
    }

    row_summary[[r]] <- data.frame(
      cluster = g,
      occupied_probability = mean(occupied_store[, g]),
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
  sigma_summary <- do.call(rbind, lapply(reported, function(g) {
    values <- post$sigma_store[inferential & occupied_store[, g], g]
    data.frame(
      cluster = g,
      mean = mean(values), sd = sd(values),
      q025 = unname(quantile(values, 0.025)),
      median = median(values),
      q975 = unname(quantile(values, 0.975)),
      final = sigma[g]
    )
  }))

  tau_final <- delta_final <- matrix(0, K, length(reported))
  for (r in seq_along(reported)) {
    g <- reported[r]
    tau_final[, r] <- tau[[g]]
    delta_final[, r] <- delta[[g]]
  }
  rownames(tau_final) <- rownames(delta_final) <- paste0("factor_", seq_len(K))
  colnames(tau_final) <- colnames(delta_final) <- reported_names

  Sigma_final <- column_diagnostics <- vector("list", length(reported))
  covariance_summary <- vector("list", length(reported))
  names(Sigma_final) <- names(column_diagnostics) <- reported_names

  for (r in seq_along(reported)) {
    g <- reported[r]
    Sigma_final[[r]] <- tcrossprod(Lambda[[g]]) + diag(psi[, g])
    active <- if (!is.null(gamma_true) && g <= ncol(gamma_true)) {
      gamma_true[, g]
    } else {
      row_pip[, g] >= pip_cut
    }
    fitted_block <- Sigma_final[[r]][active, active, drop = FALSE]
    fitted_offdiag <- if (sum(active) > 1L) {
      mean(abs(fitted_block[row(fitted_block) != col(fitted_block)]))
    } else NA_real_
    fitted_diagonal <- if (sum(active) > 0L) {
      mean(diag(fitted_block))
    } else NA_real_

    true_eigen <- true_block <- NULL
    true_offdiag <- true_diagonal <- NA_real_
    if (!is.null(truth) && !is.null(truth$Lambda) &&
        !is.null(truth$psi) && g <= length(truth$Lambda)) {
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
    if (!is.null(truth) && !is.null(truth$Sigma) &&
        g <= length(truth$Sigma) && any(active)) {
      true_block <- truth$Sigma[[g]][active, active, drop = FALSE]
      true_offdiag <- if (sum(active) > 1L) {
        mean(abs(true_block[row(true_block) != col(true_block)]))
      } else NA_real_
      true_diagonal <- mean(diag(true_block))
    }

    column_diagnostics[[r]] <- list(
      q_mean = q_summary$q_mean[r],
      q_median = q_summary$q_median[r],
      q_mode = q_summary$q_mode[r],
      q_true = if (!is.null(truth) && !is.null(truth$q) &&
                   g <= length(truth$q)) truth$q[g] else NULL,
      tail_energy = row_summary$tail_energy[r],
      fitted_eigenvalues = eigen_mean[, g],
      true_eigenvalues = true_eigen,
      tau_posterior_mean = tau_mean[, g],
      tau_final = tau[[g]],
      active_covariance_fitted = fitted_block,
      active_covariance_true = true_block
    )
    covariance_summary[[r]] <- data.frame(
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
      cluster_sizes = if (!is.null(z_true)) tabulate(z_true, G_true) else NULL,
      q = truth$q,
      active_rows = if (!is.null(gamma_true)) {
        lapply(seq_len(ncol(gamma_true)), function(g) which(gamma_true[, g]))
      } else NULL,
      loadings = truth$loadings,
      psi = truth$psi,
      covariance_distance = truth$covariance_distance,
      seed = truth$seed
    )
  }

  list(
    model = "Group-ReLU spectral IMIFA with PYP slice sampling",
    call = match.call(),
    settings = list(
      K = K, G_init = G_init, G_max = G_max,
      n_iter = n_iter, burn = burn, thin = thin,
      retained_draws = B, pip_cut = pip_cut,
      rank_prop = rank_prop, rank_abs = rank_abs,
      pyp_concentration = pyp_concentration,
      pyp_discount = pyp_discount,
      slice_rho = slice_rho,
      nu1 = nu1, nu2 = nu2,
      alpha1 = alpha1, beta1 = beta1,
      alpha2 = alpha2, beta2 = beta2,
      rho1 = rho1, rho2 = rho2,
      omega1 = omega1, omega2 = omega2,
      a_psi = a_psi, bpsi_shape = bpsi_shape, bpsi_rate = bpsi_rate,
      mu_sd = mu_sd,
      threshold_mean = threshold_mean, threshold_sd = threshold_sd,
      spectral_init = spectral_init,
      seed = seed
    ),
    data_info = list(
      N = N, P = P, K = K, G_init = length(unique(z_initial)),
      G_max = G_max, posterior_mode_G = G_mode,
      dahl_G = length(reported),
      feature_names = feature_names,
      initial_cluster_sizes = table(z_initial),
      estimated_cluster_sizes = setNames(
        tabulate(z_hat, G_max)[reported], reported_names
      ),
      final_state_cluster_sizes = table(z)
    ),
    truth_info = truth_info,
    clustering = list(
      G_summary = G_summary,
      G_distribution = G_distribution,
      RI = RI, ARI = ARI, accuracy = accuracy,
      confusion_matrix = confusion,
      estimated_cluster_sizes = setNames(
        tabulate(z_hat, G_max)[reported], reported_names
      ),
      true_cluster_sizes = if (!is.null(z_true)) {
        tabulate(z_true, G_true)
      } else NULL
    ),
    pyp = list(
      concentration = pyp_concentration,
      discount = pyp_discount,
      slice_rho = slice_rho,
      xi = xi,
      G_summary = G_summary,
      G_distribution = G_distribution,
      component_summary = component_summary,
      slice_summary = data.frame(
        mean = mean(slice_store),
        median = median(slice_store),
        max = max(slice_store),
        cap_hit_probability = mean(slice_store == G_max),
        mean_tail_mass = mean(tail_mass_store)
      )
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
      tau_posterior_mean = tau_mean[, reported, drop = FALSE],
      tau_final = tau_final,
      delta_final = delta_final
    ),
    column_diagnostics = column_diagnostics,
    covariance_summary = covariance_summary,
    row_pip = row_pip[, reported, drop = FALSE],
    row_norm = row_norm[, reported, drop = FALSE],
    selected = row_pip[, reported, drop = FALSE] >= pip_cut,
    mu_mean = mu_mean[, reported, drop = FALSE],
    psi_mean = psi_mean[, reported, drop = FALSE],
    eigen_mean = eigen_mean[, reported, drop = FALSE],
    tau_mean = tau_mean[, reported, drop = FALSE],
    z_store = post$z_store,
    G_store = G_store,
    slice_store = slice_store,
    occupied_store = occupied_store,
    inferential_draws = inferential,
    size_store = size_store,
    row_store = post$row_store,
    norm_store = post$norm_store,
    q_store = post$q_store,
    eigen_store = post$eig_store,
    psi_store = post$psi_store,
    mu_store = post$mu_store,
    tau_store = post$tau_store,
    sigma_store = post$sigma_store,
    mix_store = post$mix_store,
    stick_store_raw = stick_store,
    slice_min_store = slice_min_store,
    tail_mass_store = tail_mass_store,
    omega_store = omega_store,
    threshold_store = threshold_store,
    z_hat = z_hat,
    z = z,
    z_ref = post$z_ref,
    psm = post$psm,
    label_maps = post$maps,
    truth_map = truth_map,
    reported_clusters = reported,
    q_hat = setNames(q_hat_all[reported], reported_names),
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
    stick_raw = stick_final_raw,
    eta = eta
  )
}

Rsparse <- function(Xg, q_g, n_iter = 4000, burn = 2000, thin = 2L, v0 = 1e-4, v1 = 1,
                    a_pi = 1, b_pi = 1, a_psi = 2, b_psi = 1,
                    init_Lambda = NULL, init_psi = NULL, init_mu = NULL,
                    verbose = FALSE) {
  Xg <- as.matrix(Xg)
  n <- nrow(Xg)
  p <- ncol(Xg)

  feature_names <- colnames(Xg)
  if (is.null(feature_names)) {
    feature_names <- paste0("V", seq_len(p))
    colnames(Xg) <- feature_names
  }

  q_g <- as.integer(q_g)
  n_save <- floor((n_iter - burn) / thin)

  mu_g <- if (is.null(init_mu)) colMeans(Xg) else init_mu
  Y <- sweep(Xg, 2, mu_g, FUN = "-")

  if (q_g == 0L) {
    psi0 <- if (is.null(init_psi)) {
      pmax(apply(Y, 2, stats::var), 1e-3)
    } else {
      as.numeric(init_psi)
    }

    return(list(
      feature = feature_names,
      mu_g = mu_g,
      q_g = q_g,
      pip = rep(0, p),
      gamma_mode = rep(0L, p),
      mean_row_l2 = rep(0, p),
      mean_row_sqsum = rep(0, p),
      mean_h2 = rep(0, p),
      Lambda_mean = matrix(0, nrow = p, ncol = 0,
                           dimnames = list(feature_names, NULL)),
      psi_mean = psi0,
      pi_mean = 0,
      gamma_samples = matrix(
        0L,
        nrow = p,
        ncol = n_save,
        dimnames = list(feature_names, paste0("draw", seq_len(n_save)))
      )
    ))
  }

  if (is.null(init_Lambda)) {
    sv <- svd(Y)
    r <- min(q_g, length(sv$d))

    Lambda <- matrix(0, nrow = p, ncol = q_g)

    Lambda[, seq_len(r)] <-
      sv$v[, seq_len(r), drop = FALSE] %*%
      diag(
        sqrt(pmax(sv$d[seq_len(r)]^2 / n, 1e-6)),
        nrow = r
      )
  } else {
    Lambda <- as.matrix(init_Lambda)
  }

  rownames(Lambda) <- feature_names

  psi <- if (is.null(init_psi)) {
    pmax(apply(Y, 2, stats::var), 1e-3)
  } else {
    as.numeric(init_psi)
  }

  gamma <- as.integer(sqrt(rowSums(Lambda^2)) > 0.05)
  pi_g <- mean(gamma)
  pi_g <- min(max(pi_g, 1e-4), 1 - 1e-4)

  gamma_draws <- matrix(NA_integer_, nrow = p, ncol = n_save, dimnames = list(feature_names, paste0("draw", seq_len(n_save))))

  Lambda_sum <- matrix(0, nrow = p, ncol = q_g)
  psi_sum <- numeric(p)
  gamma_sum <- numeric(p)
  row_l2_sum <- numeric(p)
  row_sqsum_sum <- numeric(p)
  h2_sum <- numeric(p)
  pi_sum <- 0

  save_id <- 0L

  for (iter in seq_len(n_iter)) {

    prec_eta <- diag(q_g) + crossprod(Lambda, Lambda / psi)

    V_eta <- solve(prec_eta)
    V_eta <- (V_eta + t(V_eta)) / 2

    Eta_mean <- Y %*% (Lambda / psi) %*% V_eta

    Eta <- Eta_mean +
      matrix(rnorm(n * q_g), nrow = n, ncol = q_g) %*% chol(V_eta)

    XtX <- crossprod(Eta)

    for (j in seq_len(p)) {
      vj <- if (gamma[j] == 1L) v1 else v0

      Vj <- solve(XtX / psi[j] + diag(1 / vj, q_g))
      Vj <- (Vj + t(Vj)) / 2

      mj <- Vj %*% (crossprod(Eta, Y[, j, drop = FALSE]) / psi[j])

      Lambda[j, ] <- as.numeric(
        mj + t(chol(Vj)) %*% rnorm(q_g)
      )
    }

    pi_g <- min(max(pi_g, 1e-8), 1 - 1e-8)

    for (j in seq_len(p)) {
      lam_j <- Lambda[j, ]

      log_p1 <- log(pi_g) +
        sum(dnorm(lam_j, mean = 0, sd = sqrt(v1), log = TRUE))

      log_p0 <- log(1 - pi_g) +
        sum(dnorm(lam_j, mean = 0, sd = sqrt(v0), log = TRUE))

      m <- max(log_p0, log_p1)

      prob1 <- exp(log_p1 - m) /
        (exp(log_p0 - m) + exp(log_p1 - m))

      gamma[j] <- rbinom(1, size = 1, prob = prob1)
    }

    pi_g <- rbeta(
      1,
      shape1 = a_pi + sum(gamma),
      shape2 = b_pi + p - sum(gamma)
    )

    resid <- Y - Eta %*% t(Lambda)
    rss <- colSums(resid^2)

    psi <- 1 / rgamma(
      p,
      shape = a_psi + n / 2,
      rate = b_psi + rss / 2
    )

    psi <- pmax(psi, 1e-8)

    if (iter > burn && ((iter - burn) %% thin == 0L)) {
      save_id <- save_id + 1L

      gamma_draws[, save_id] <- gamma

      Lambda_sum <- Lambda_sum + Lambda
      psi_sum <- psi_sum + psi
      gamma_sum <- gamma_sum + gamma
      pi_sum <- pi_sum + pi_g

      row_sqsum_now <- rowSums(Lambda^2)
      row_l2_now <- sqrt(row_sqsum_now)
      h2_now <- row_sqsum_now / (row_sqsum_now + psi)

      row_sqsum_sum <- row_sqsum_sum + row_sqsum_now
      row_l2_sum <- row_l2_sum + row_l2_now
      h2_sum <- h2_sum + h2_now
    }

    if (verbose && iter %% 500 == 0L) {
      message(sprintf(
        "[Rsparse] iter = %d / %d, mean(gamma) = %.3f, pi_g = %.3f",
        iter, n_iter, mean(gamma), pi_g
      ))
    }
  }

  Lambda_mean <- Lambda_sum / n_save
  rownames(Lambda_mean) <- feature_names

  list(
    feature = feature_names,
    mu_g = mu_g,
    q_g = q_g,
    pip = gamma_sum / n_save,
    gamma_mode = as.integer((gamma_sum / n_save) >= 0.5),
    mean_row_l2 = row_l2_sum / n_save,
    mean_row_sqsum = row_sqsum_sum / n_save,
    mean_h2 = h2_sum / n_save,
    Lambda_mean = Lambda_mean,
    psi_mean = psi_sum / n_save,
    pi_mean = pi_sum / n_save,
    gamma_samples = gamma_draws
  )
}

Rsparse_imifa <- function(dat, res, n_iter = 4000, burn = 2000, thin = 2L,
                          v0 = 1e-4, v1 = 1, a_pi = 1, b_pi = 1, a_psi = 2, b_psi = 1,
                          verbose = FALSE) {
  X <- as.matrix(dat)

  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  
  zMAP <- as.integer(res$Clust$MAP)
  q_g <- as.integer(res$GQ.results$Q)
  G <- length(q_g)

  fits <- vector("list", G)
  names(fits) <- paste0("Cluster", seq_len(G))

  for (g in seq_len(G)) {
    idx <- which(zMAP == g)
    Xg <- X[idx, , drop = FALSE]

    if (verbose) {
      message(sprintf(
        "[Rsparse_imifa] Cluster%d: n_g = %d, q_g = %d",
        g, nrow(Xg), q_g[g]
      ))
    }

    fits[[g]] <- Rsparse(Xg = Xg, q_g = q_g[g], n_iter = n_iter, burn = burn, thin = thin,
                        v0 = v0, v1 = v1, a_pi = a_pi, b_pi = b_pi, a_psi = a_psi, b_psi = b_psi,
                        init_Lambda = as.matrix(res$Loadings$post.load[[g]]), init_psi = as.numeric(res$Uniquenesses$post.psi[, g]), init_mu = colMeans(Xg),
                        verbose = verbose)
  }

  list(fits = fits, zMAP = zMAP, q_g = q_g, cluster_sizes = table(factor(zMAP, levels = seq_len(G))))
}

summarise_RS <- function(row_sparse_fit, p_cut = 0.5) {
  fits <- row_sparse_fit$fits

  tables <- lapply(seq_along(fits), function(g) {
    fit_g <- fits[[g]]

    df <- data.frame(
      cluster = names(fits)[g],
      feature = fit_g$feature,
      p_inclusion = fit_g$pip,
      selected = fit_g$pip >= p_cut,
      mean_row_l2 = fit_g$mean_row_l2
    )

    df[order(df$p_inclusion, decreasing = TRUE), ]
  })

  names(tables) <- names(fits)

  summary <- do.call(
    rbind,
    lapply(tables, function(df) {
      data.frame(
        cluster = unique(df$cluster),
        n_features = nrow(df),
        n_selected = sum(df$selected),
        mean_pip_selected = mean(df$p_inclusion[df$selected])
      )
    })
  )

  rownames(summary) <- NULL

  list(
    tables = tables,
    summary = summary
  )
}
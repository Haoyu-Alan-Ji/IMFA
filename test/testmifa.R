source(here::here('R', 'group_sparse_mifa.R'))


# This is the simfun from the supplied QMD.
simfun <- function(n = 50, p = 60, q = c(2, 2, 3), active = 5,
                   loadings = c(2, 3, 4), psi = 1, seed = 200) {
  set.seed(seed)
  G <- length(q)
  stopifnot(length(loadings) == G)

  Xg <- vector("list", G)
  eta <- vector("list", G)
  eps <- vector("list", G)
  Lambda <- vector("list", G)
  Sigma <- vector("list", G)
  allocation <- vector("list", G)

  for (g in seq_len(G)) {
    L <- matrix(0, p, q[g])
    A <- matrix(0, active, q[g])

    for (j in seq_len(active)) {
      main <- 1L + (j - 1L) %% q[g]
      rest <- setdiff(seq_len(q[g]), main)
      rest <- rest[sample.int(length(rest))]
      order <- c(main, rest)
      left <- loadings[g]^2
      share <- numeric(q[g])

      if (q[g] == 1L) {
        share[main] <- left
      } else {
        for (k in seq_len(q[g] - 1L)) {
          part <- sample(60:70, 1) / 100
          share[order[k]] <- left * part
          left <- left - share[order[k]]
        }
        share[order[q[g]]] <- left
      }

      A[j, ] <- sqrt(share)
      if (q[g] > 1L) {
        A[j, rest] <- A[j, rest] *
          sample(c(-1, 1), length(rest), replace = TRUE)
      }
    }

    L[seq_len(active), ] <- A
    eta[[g]] <- matrix(rnorm(n * q[g]), n, q[g])
    eps[[g]] <- matrix(rnorm(n * p, sd = sqrt(psi)), n, p)
    Xg[[g]] <- eta[[g]] %*% t(L) + eps[[g]]
    Lambda[[g]] <- L
    Sigma[[g]] <- tcrossprod(L) + diag(psi, p)
    allocation[[g]] <- A^2 / rowSums(A^2)
  }

  X <- do.call(rbind, Xg)
  z <- rep(seq_len(G), each = n)
  ind <- sample.int(G * n)
  X <- X[ind, , drop = FALSE]
  z <- z[ind]
  colnames(X) <- paste0("V", seq_len(p))

  distance <- matrix(0, G, G)
  for (g in seq_len(G)) {
    for (h in seq_len(G)) {
      distance[g, h] <- sqrt(sum((Sigma[[g]] - Sigma[[h]])^2))
    }
  }

  list(
    X = X, z_true = z, X_cluster = Xg, eta = eta, noise = eps,
    Lambda = Lambda, Sigma = Sigma, allocation = allocation,
    gamma_true = matrix(rep(c(rep(1L, active), rep(0L, p - active)), G),
                        nrow = p, ncol = G),
    q = q, n = n, p = p, G = G, active = seq_len(active),
    loadings = loadings, psi = psi, covariance_distance = distance, seed = seed
  )
}


selection_summary <- function(fit, truth) {
  out <- vector("list", ncol(truth))

  for (g in seq_len(ncol(truth))) {
    selected <- fit$row_pip[, g] >= 0.5
    actual <- truth[, g] == 1L
    TP <- sum(selected & actual)
    FP <- sum(selected & !actual)
    FN <- sum(!selected & actual)

    out[[g]] <- data.frame(
      cluster = g,
      selected = sum(selected),
      TP = TP,
      FP = FP,
      FN = FN,
      precision = ifelse(TP + FP == 0, 0, TP / (TP + FP)),
      recall = TP / (TP + FN),
      F1 = ifelse(2 * TP + FP + FN == 0, 0,
                  2 * TP / (2 * TP + FP + FN))
    )
  }
  do.call(rbind, out)
}


test_group_sparse_mifa <- function() {
  # Verify the supplied simulation before fitting anything.
  sim0 <- simfun()
  stopifnot(dim(sim0$X)[1] == 150L, dim(sim0$X)[2] == 60L)
  stopifnot(all(tabulate(sim0$z_true) == 50L))
  stopifnot(all(vapply(sim0$Lambda, function(x) qr(x)$rank,
                       integer(1)) == sim0$q))
  stopifnot(all(vapply(seq_along(sim0$Lambda), function(g) {
    max(abs(sqrt(rowSums(sim0$Lambda[[g]][1:5, , drop = FALSE]^2)) -
              sim0$loadings[g]))
  }, numeric(1)) < 1e-10))
  stopifnot(all(sim0$covariance_distance[upper.tri(
    sim0$covariance_distance
  )] > 0))

  # A smaller draw from exactly the same simfun keeps the smoke tests quick.
  sim <- simfun(
    n = 12, p = 12, q = c(2, 2, 3), active = 5,
    loadings = c(2, 2.5, 3), seed = 301
  )

  common <- list(
    X = sim$X,
    G = sim$G,
    q_start = 3,
    q_cap = 4,
    adapt = FALSE,
    n_iter = 50,
    burn = 25,
    thin = 5,
    z_init = sim$z_true,
    equal_weights = TRUE,
    verbose = FALSE
  )

  fit_indicator <- do.call(group_sparse_mifa, c(common, list(
    gate = "indicator", update_z = TRUE, seed = 401
  )))
  fit_relu <- do.call(group_sparse_mifa, c(common, list(
    gate = "relu", update_z = FALSE, seed = 402
  )))

  for (fit in list(fit_indicator, fit_relu)) {
    stopifnot(
      all(dim(fit$row_pip) == c(sim$p, sim$G)),
      all(dim(fit$z_store) == c(5, sim$n * sim$G)),
      all(fit$row_pip >= 0 & fit$row_pip <= 1),
      all(is.finite(fit$row_norm)),
      all(is.finite(fit$threshold)),
      all(fit$q_store == 3L),
      all(vapply(seq_len(sim$G), function(g) {
        isTRUE(all.equal(fit$tau[[g]], cumprod(fit$delta[[g]])))
      }, logical(1)))
    )

    for (g in seq_len(sim$G)) {
      activation <- if (fit$gate == "indicator") {
        as.numeric(fit$alpha[, g] > fit$alpha0)
      } else {
        pmax(fit$alpha[, g] - fit$alpha0, 0)
      }
      stopifnot(isTRUE(all.equal(
        fit$Lambda[[g]], fit$slab[[g]] * activation,
        tolerance = 1e-10
      )))
    }
  }

  # Directly exercise adaptive addition/removal and the q_g = 0 boundary.
  fit_adapt <- group_sparse_mifa(
    sim$X, sim$G, gate = "indicator",
    q_start = 3, q_cap = 4,
    adapt = TRUE, factor_eps = 1e6,
    n_iter = 30, burn = 15, thin = 3,
    z_init = sim$z_true, update_z = FALSE,
    equal_weights = TRUE, seed = 403, verbose = FALSE
  )
  stopifnot(
    all(fit_adapt$q_store >= 0L),
    all(fit_adapt$q_store <= 4L),
    any(fit_adapt$q_store == 0L)
  )

  # Known label permutation: postprocessing must recover one common labeling.
  z_fake <- rbind(c(1L, 1L, 2L, 2L), c(2L, 2L, 1L, 1L))
  row_fake <- norm_fake <- array(0, c(2, 2, 2))
  row_fake[, , 1] <- norm_fake[, , 1] <- matrix(c(1, 1, 2, 2), 2, 2)
  row_fake[, , 2] <- norm_fake[, , 2] <- matrix(c(2, 2, 1, 1), 2, 2)
  q_fake <- rbind(c(1L, 2L), c(2L, 1L))
  sigma_fake <- rbind(c(1, 2), c(2, 1))

  post <- relabel_mifa(
    z_fake, row_fake, norm_fake, q_fake, sigma_fake, z_fake[2, ]
  )
  stopifnot(
    identical(as.integer(post$z_store[1, ]),
              as.integer(post$z_store[2, ])),
    identical(as.integer(post$q_store[1, ]),
              as.integer(post$q_store[2, ])),
    identical(as.numeric(post$row_store[, , 1]),
              as.numeric(post$row_store[, , 2]))
  )

  message("simfun, group gates, MGP updates, adaptive q, and relabeling: passed")
  invisible(list(
    sim = sim,
    indicator = fit_indicator,
    relu = fit_relu,
    adaptive = fit_adapt,
    indicator_selection = selection_summary(fit_indicator, sim$gamma_true),
    relu_selection = selection_summary(fit_relu, sim$gamma_true)
  ))
}

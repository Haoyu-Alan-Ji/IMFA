simfun <- function(n = 50, p = 60, q = c(2, 2, 3), active = 5,
                   loadings = c(2, 3, 4), psi = 1, seed = 200) {
  set.seed(seed)
  G <- length(q)
  stopifnot(length(loadings) == G)

  Xg <- eta <- eps <- Lambda <- Sigma <- allocation <- vector("list", G)

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
  order <- sample.int(G * n)
  X <- X[order, , drop = FALSE]
  z <- z[order]
  colnames(X) <- paste0("V", seq_len(p))

  covariance_distance <- matrix(0, G, G)
  for (g in seq_len(G)) {
    for (h in seq_len(G)) {
      covariance_distance[g, h] <-
        sqrt(sum((Sigma[[g]] - Sigma[[h]])^2))
    }
  }

  list(
    X = X,
    z_true = z,
    X_cluster = Xg,
    eta = eta,
    noise = eps,
    Lambda = Lambda,
    Sigma = Sigma,
    allocation = allocation,
    gamma_true = matrix(
      rep(c(rep(1L, active), rep(0L, p - active)), G), p, G
    ),
    q = q,
    n = n,
    p = p,
    G = G,
    active = seq_len(active),
    loadings = loadings,
    psi = psi,
    covariance_distance = covariance_distance,
    seed = seed
  )
}

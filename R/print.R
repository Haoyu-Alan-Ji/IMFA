q_summary <- function(fit) {

  apply(
    fit$q_store,
    2,
    function(x) {
      c(
        mean = mean(x),
        median = median(x),
        mode = as.integer(
          names(
            which.max(
              table(x)
            )
          )
        )
      )
    }
  )
}

print_spectrum <- function(
  fit,
  sim
) {

  cat("\nSVD\n")
  cat("-----------------------------\n")

  for (g in seq_len(sim$G)) {

    L <- fit$Lambda[[g]][
      sim$active,
      ,
      drop = FALSE
    ]

    Ltrue <- sim$Lambda[[g]][
      sim$active,
      ,
      drop = FALSE
    ]

    cat("\ncluster", g, "\n")

    cat(
      "stored q =",
      ncol(L),
      "\n"
    )

    cat(
      "estimated:",
      paste(
        round(
          svd(L)$d,
          3
        ),
        collapse = " "
      ),
      "\n"
    )

    cat(
      "true:     ",
      paste(
        round(
          svd(Ltrue)$d,
          3
        ),
        collapse = " "
      ),
      "\n"
    )
  }
}

print_tau <- function(
  fit,
  sim
) {

  cat("\nMGP tau\n")
  cat("-----------------------------\n")

  for (g in seq_len(sim$G)) {

    cat(
      "cluster", g, ":",
      paste(
        round(
          fit$tau[[g]],
          3
        ),
        collapse = " "
      ),
      "\n"
    )
  }
}

print_clustering <- function(
  fit,
  sim
) {

  cat("\nClustering\n")
  cat("-----------------------------\n")

  ari <- mclust::adjustedRandIndex(
    fit$z_hat,
    sim$z_true
  )

  cat(
    "ARI =",
    round(
      ari,
      4
    ),
    "\n\n"
  )

  print(
    table(
      truth = sim$z_true,
      estimate = fit$z_hat
    )
  )
}

print_support <- function(
  fit,
  sim,
  cutoff = 0.5
) {

  cat("\nRow selection\n")
  cat("-----------------------------\n")

  truth <- sim$gamma_true == 1

  for (g in seq_len(sim$G)) {

    selected <-
      fit$row_pip[, g] >
      cutoff

    TP <- sum(
      selected &
        truth[, g]
    )

    FP <- sum(
      selected &
        !truth[, g]
    )

    FN <- sum(
      !selected &
        truth[, g]
    )

    precision <-
      TP /
      (TP + FP)

    recall <-
      TP /
      (TP + FN)

    f1 <-
      2 *
      precision *
      recall /
      (
        precision +
          recall
      )

    cat(
      "cluster", g,
      " selected =", sum(selected),
      " TP =", TP,
      " FP =", FP,
      " FN =", FN,
      " F1 =", round(f1, 3),
      "\n"
    )
  }

  cat("\nfirst 10 row PIPs:\n")

  print(
    round(
      fit$row_pip[
        seq_len(
          min(
            10,
            nrow(fit$row_pip)
          )
        ),
        ,
        drop = FALSE
      ],
      3
    )
  )
}

print_fit <- function(
  name,
  fit,
  sim,
  clustering = FALSE,
  support = FALSE,
  tau = TRUE
) {

  cat("\n")
  cat("========================================\n")
  cat(name, "\n")
  cat("========================================\n")

  cat("\nq summary\n")
  cat("-----------------------------\n")

  print(
    round(
      q_summary(fit),
      3
    )
  )

  if (clustering) {
    print_clustering(
      fit,
      sim
    )
  }

  if (support) {
    print_support(
      fit,
      sim
    )
  }

  if (tau) {
    print_tau(
      fit,
      sim
    )
  }

  print_spectrum(
    fit,
    sim
  )
}
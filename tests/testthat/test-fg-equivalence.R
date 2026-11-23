test_that("U2 components: C++ kernel equals the R reference", {
  set.seed(42)
  
  n <- 40
  z1 <- rbinom(n, 1, 0.5)
  z2 <- rbinom(n, 1, 0.4)
  z  <- cbind(z1, z2)
  colnames(z) <- c("trt", "sex")
  
  # Survival with a treatment effect; random censoring
  t      <- rexp(n, rate = exp(0.3 * z1))
  status <- rbinom(n, 1, 0.8)
  qt     <- pmin(t, t * runif(n, 0.4, 1))
  
  lt <- log(t)
  lq <- log(qt)
  
  beta <- c(0.3, -0.2)
  L    <- stats::quantile(t, 0.9)
  hL   <- log(L)
  
  km_surv <- QAST:::prep_km_by_column(z, t, status)
  
  # Grid of theta for the point-by-point comparison
  theta_grid <- rbind(
    c(0.0, 0.0),
    c(0.2, -0.1),
    c(-0.4, 0.5),
    c(0.7, 0.7)
  )
  
  for (phi in c("identity", "piecewise")) {
    for (norm_min in c("L1", "L2")) {
      
      for (i in seq_len(nrow(theta_grid))) {
        theta <- theta_grid[i, ]
        
        parts_cpp <- QAST:::U2_parts_fg_cpp(theta, lt, lq, z, status,
                                            beta, hL, km_surv, phi)
        parts_ref <- U2_parts_fg_ref(theta, lt, lq, z, status,
                                            beta, hL, km_surv, phi)
        expect_equal(parts_cpp, parts_ref, tolerance = 1e-10)
        
        obj_cpp <- QAST:::U2_obj_fg_cpp(theta, lt, lq, z, status,
                                        beta, hL, km_surv, phi, norm_min)
        
        # Replicates the aggregation of the norm in the reference
        if (norm_min == "L1") {
          obj_ref <- sum(abs(parts_ref)) / sqrt(length(parts_ref))
        } else {
          obj_ref <- sqrt(sum(parts_ref^2)) / sqrt(length(parts_ref))
        }
        expect_equal(obj_cpp, obj_ref, tolerance = 1e-10)
      }
    }
  }
})
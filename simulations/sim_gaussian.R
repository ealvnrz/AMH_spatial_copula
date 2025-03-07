# Libraries
rm(list = ls())
require(GeoModels)
require(dotCall64)
require(parallel)

################################################################################
# Objective function for copula-based pairwise likelihood
################################################################################
complogcopulik2 <- function(param, npairs, fixed, namesfixed, namescorr, namesnuis,
                            namesparam, corrmodel, colidx, rowidx, lags, data1, data2, X, model) {
  names(param) <- namesparam
  param <- c(param, fixed)
  
  # Check for NAs in paramcorr
  paramcorr <- as.numeric(param[namescorr])
  if (any(is.na(paramcorr))) {
    stop("paramcorr contains NA values, check your parameters.")
  }
  
  nuisance <- param[namesnuis]
  sel <- substr(names(nuisance), 1, 4) == "mean"
  mm <- as.numeric(nuisance[sel])  # linear mean
  Mean <- c(X %*% mm)
  other_nuis <- as.numeric(nuisance[!sel])
  nu <- other_nuis[3]  # adjust in C code if needed
  
  result <- dotCall64::.C64(
    "Comp_Pair_cop2",
    SIGNATURE = c(
      "integer", "integer", "double", "double", "double", "double",
      "double", "double", "double", "double", "double", "integer"
    ),
    corrmodel, npairs, lags, data1, data2, paramcorr,
    res = dotCall64::numeric_dc(1), Mean[colidx], Mean[rowidx],
    other_nuis, nu, model,
    INTENT = c("r", "r", "r", "r", "r", "r", "rw", "r", "r", "r", "r", "r")
  )$res
  
  return(-result)
}

################################################################################
# Main script
################################################################################

# Load dynamic library
dyn.load("parte.so")

# Simulation parameters
N <- 400
coords <- cbind(runif(N), runif(N))

nu <- 1
mean <- 0
sill <- 1
scale <- 0.2 / 3
smooth <- 0.5
nugget <- 0

pcorr <- list(scale = scale, smooth = smooth)
param_marg <- list(mean = mean, nugget = nugget, sill = sill, nu = nu)
param_simu <- c(param_marg, pcorr)

model <- "Gaussian"
corrmodel <- 14  # Matern

X <- matrix(rep(1, N), ncol = 1)

I <- 3.5
lower <- c(-I, 0, -I, 0)
upper <- c(I, I, I, I)

param <- unlist(param_simu)
start <- param[c(1, 3, 4, 5)]
fixed <- param[c(2, 6)]

set.seed(7)

# Function to run one iteration
run_iteration <- function(iter) {
  # Simulate data
  data <- GeoSimCopula(
    coordx = coords,
    corrmodel = "Matern",
    model = model,
    param = param_simu,
    copula = "AMH"
  )$data
  
  # Neighbors selection
  sel <- GeoNeighIndex(coordx = coords, neighb = 2)
  data1 <- data[sel$colidx]
  data2 <- data[sel$rowidx]  
  
  # Check for NA in simulated data
  if (any(is.na(data1)) || any(is.na(data2))) {
    stop("Simulated data contains NA values.")
  }
  
  # Minimize the objective function
  b <- nlminb(
    start = start,
    objective = complogcopulik2,
    lower = lower,
    upper = upper,
    npairs = length(sel$lags),
    fixed = fixed,
    data1 = data1,
    data2 = data2,
    X = X,
    corrmodel = corrmodel,
    colidx = sel$colidx,
    rowidx = sel$rowidx,
    lags = sel$lags,
    model = CkModel(model),
    namescorr = names(unlist(pcorr)),
    namesnuis = names(unlist(param_marg)),
    namesfixed = names(fixed),
    namesparam = names(start)
  )
  
  # Save the result to the file after each iteration
  write.table(
    t(b$par),
    file = "results_gaussian_400.csv",
    append = TRUE,
    sep = ",",
    col.names = !file.exists("results_gaussian_400.csv"),
    row.names = FALSE
  )
  
  # Print progress
  cat(paste("Iteration", iter, "completed at", Sys.time(), "\n"))
  
  return(b$par)
}

# Parallelization of iterations using mclapply
num_iterations <- 100
num_cores <- detectCores() - 4

# Run iterations in parallel
results <- mclapply(1:num_iterations, run_iteration, mc.cores = num_cores)

# Convert the results to a matrix
res <- do.call(rbind, results)

# Save the final results to a file
write.csv(res, "final_results_gaussian_400.csv", row.names = FALSE)

# Unload dynamic library
dyn.unload("parte.so")

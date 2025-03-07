# Clear the environment
rm(list = ls())

# Load required packages
require(GeoModels)
require(dotCall64)

# Auxiliary function
complogcopulik2 <- function(param, npairs, fixed, namesfixed, namescorr, namesnuis, namesparam,
                            corrmodel, colidx, rowidx, lags, data1, data2, X, model) {
  names(param) <- namesparam
  param <- c(param, fixed)
  paramcorr <- as.numeric(param[namescorr])
  nuisance <- param[namesnuis]
  sel <- substr(names(nuisance), 1, 4) == "mean"
  mm <- as.numeric(nuisance[sel])  # linear mean
  Mean <- c(X %*% mm)
  other_nuis <- as.numeric(nuisance[!sel])
  nu <- other_nuis[2]  # note: adjust if needed in corresponding C code
  
  result <- dotCall64::.C64(
    "Comp_Pair_cop2",
    SIGNATURE = c(
      "integer", "integer", "double", "double", "double", 
      "double", "double", "double", "double", "double",
      "double", "integer"
    ),
    corrmodel, npairs, lags, data1, data2, paramcorr,
    res = dotCall64::numeric_dc(1), Mean[colidx], Mean[rowidx],
    other_nuis, nu, model,
    INTENT = c("r", "r", "r", "r", "r", "r", "rw", "r", "r", "r", "r", "r")
  )$res
  
  return(-result)
}

# Load compiled shared object
dyn.load("parte.so")

# Data
data <- read.csv("data/data_z.csv")[, 2]
coords <- read.csv("data/data_coords.csv")[, 2:3]
maxdist <- max(dist(coords))

# Plots
quilt.plot(coords, data)

# Empirical Semivariogram
a <- GeoVariogram(coordx = coords, data = data, maxdist = maxdist / 3)
plot(a, pch = 20, ylim = c(0, 0.0025), xlab = "Distance", ylab = "Semivariogram")

# Check Beta model as a marginal model (assuming independence)
corrmodel <- "Matern"
scale <- 200
smooth <- 0
power2 <- 4
nugget <- 0
model <- "Beta2"
min_1 <- -1
max_1 <- 1
shape <- 1
mean <- 0
nugget <- 0
optimizer <- "nlminb"

I <- 5
fixed <- list(nugget = nugget, sill = 1, scale = scale, smooth = smooth, min = -1, max = 1)
start <- list(shape = shape, mean = mean)
lower <- list(shape = 0, mean = -I)
upper <- list(shape = 1000, mean = I)

# Maximum independence likelihood
fit0 <- GeoFit(
  data = data, coordx = coords, corrmodel = corrmodel, model = model,
  likelihood = "Marginal", type = "Independence",
  optimizer = optimizer, lower = lower, upper = upper,
  start = start, fixed = fixed
)

# Histogram of data with fitted Beta density
hist(data, xlim = c(0.1, 0.75), xlab = "NDVI", main = "Histogram of NDVI",
     ylim = c(0, 10), freq = FALSE, breaks = 20)

pmin <- -1
pmax <- 1
ll <- seq(0.1, 0.7, 0.001)
sh <- fit0$param$shape
mm <- 1 / (1 + exp(-fit0$param$mean))
ds <- dbeta((ll - pmin) / (pmax - pmin), shape1 = mm * sh, shape2 = (1 - mm) * sh) / (pmax - pmin)
lines(ll, ds, lwd = 2)

# Normal score transformation
bb <- GeoPit(fit0, type = "Gaussian")
hist(bb$data, freq = FALSE)
GeoScatterplot(bb$data, coords, neighb = c(1, 2, 3, 4))

# Gaussian Copula Beta Model
copula <- "Gaussian"
model <- "Beta2"
corrmodel <- "Matern"
smooth <- 0.5

I <- 7
lower <- list(shape = 0, scale = 0, mean = -I)
upper <- list(shape = 5000, scale = 1000000, mean = I)
fixed <- list(sill = 1, smooth = smooth, min = -1, max = 1, nugget = 0)
start <- list(mean = 0, shape = 6, scale = 1500)

fit1 <- GeoFit2(
  data = data, coordx = coords, corrmodel = corrmodel, model = model,
  neighb = 2, likelihood = "Marginal", type = "Pairwise",
  optimizer = "nlminb", lower = lower, upper = upper,
  copula = copula, start = start, fixed = fixed, sensitivity = TRUE
)

# AMH Copula
nu <- 1
min <- 0
max <- 1
N <- length(data)
X <- matrix(rep(1, N), ncol = 1)

pcorr <- list(scale = scale, smooth = smooth)
param_marg <- list(mean = mean, nugget = nugget, nu = nu, shape = shape, min = min, max = max)
param_simu <- c(param_marg, pcorr)

model <- "Beta2"
corrmodel <- 14  # Matern
param <- unlist(param_simu)

start <- param[c(1, 3, 4, 7)]
fixed <- param[c(2, 5, 6, 8)]

sel <- GeoNeighIndex(coordx = coords, neighb = 2)
data1 <- data[sel$colidx]
data2 <- data[sel$rowidx]

b <- nlminb(
  start = start,
  objective = complogcopulik2,
  lower = lower, upper = upper,
  npairs = length(sel$lags),
  fixed = fixed,
  data1 = data1, data2 = data2, X = X,
  corrmodel = corrmodel, colidx = sel$colidx, rowidx = sel$rowidx,
  lags = sel$lags, model = CkModel(model),
  namescorr = names(unlist(pcorr)),
  namesnuis = names(unlist(param_marg)),
  namesfixed = names(fixed),
  namesparam = names(start)
)

print(b$par)
print(b$objective)

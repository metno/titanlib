titanlib_path <- "../build/extras"
dyn.load( file.path( titanlib_path, 
                     paste("SWIG/R/titanlib", .Platform$dynlib.ext, sep="")))
source(   file.path( titanlib_path,"SWIG/R/titanlib.R"))

#---------------------------------------------------------------
# Test with small vectors
P <- 11
lats <- rep( 60, P)
lons <- seq(10, length=P, by=0.005)
elevs <- rep(0, length(lats))
values <- round( 10*sin( lons * 2 * base::pi / ( max(lons) - min(lons))), 2)
# this is a bad observation
values[2]<-20
# check all obs
obs_to_check <- rep(1, length(lats))
# -- buddy check parameters --
radius <- 20000 # km
num_iterations <- 10
num_min <- 3 
threshold <- 2
max_elev_diff <- 100000 # m, if negative then elevation is not considered
elev_gradient <- -0.0065 # used only when max_elev_diff is positive
min_std <- 1
# call titanlib
points <- Points(lats, lons, elevs)
flags <- buddy_check(points, values, radius, num_min, threshold, max_elev_diff, elev_gradient, min_std, num_iterations)
print(flags)

#------------------------------------------
# Test with bigger vectors
P <- 50000
pGE <- 0.3 # probability of gross error 0.3 = 30%
lats <- runif(P, min = 55, max = 60)
lons <- runif(P, min = 5, max = 10)
elevs <- runif(P, min = 0, max = 2500)
# simple vertical profile
values <- 30 - 0.0065 * elevs
idx <- sample(1:P,ceiling(P*pGE))
true_flags <- values; true_flags[] <- 0
true_flags[idx] <- 1
values[idx]<-runif(ceiling(P*pGE), min = -50, max = 50)
obs_to_check = rep(1,P)
# buddy check parameters
radius <- 20000
num_iterations <- 10
num_min <- 3
threshold <- 2
max_elev_diff <- 100000
elev_gradient <- -0.0065
min_std <- 1
# call titanlib
points <- Points(lats, lons, elevs)
t0<-Sys.time()
flags <- buddy_check(points, values, radius, num_min, threshold, max_elev_diff, elev_gradient, min_std, num_iterations)
t1<-Sys.time()
print(t1-t0)
a <- length( which( true_flags == 1 & flags == 1)) # hit (correct bad)
c <- length( which( true_flags == 1 & flags == 0)) # miss (false good)
b <- length( which( true_flags == 0 & flags == 1)) # false bad
d <- length( which( true_flags == 0 & flags == 0)) # correct good
rand <- (a+c) * (a+b) / (a+b+c+d)
ets <- (a-rand) / (a+b+c-rand)
acc <- (a+d)/(a+b+c+d)
pod <- a/(a+c)
pofa <- b/(b+d)
cat( paste0( "TOT / corr.bad(true_bad) false_bad false_good  corr.good: ",a+b+c+d," / ",a," (",length(which( true_flags == 1)),") ",b," ",c," ",d,"\n"))
cat( paste("scores / acc pod pofa ets:", round(acc,2), round(pod,2), round(pofa,2), round(ets,2),"\n"))


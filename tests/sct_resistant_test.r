titanlib_path <- "../build/extras"
dyn.load( file.path( titanlib_path, 
                     paste("SWIG/R/titanlib", .Platform$dynlib.ext, sep="")))
source(   file.path( titanlib_path,"SWIG/R/titanlib.R"))

#---------------------------------------------------------------
# Test with small vectors
lats <- c(60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60)
lons <- 10 + 0.005 * 0:(length(lats)-1)
elevs <- rep(0, length(lats))
values <- round( 10*sin(lons*2*base::pi/(max(lons)-min(lons))), 2)
#values = c(0, 0, 0, 0, -100)
obs_to_check = rep(1, length(lats))
background_values = rep(0, length(lats))
background_elab_type = "MedianOuterCircle"
N = length(lats)
num_min_outer = 3
num_max_outer = 10
inner_radius = 20000
outer_radius = 50000
num_iterations = 10
num_min_prof = 1
min_elev_diff = 100
min_horizontal_scale = 250 
max_horizontal_scale = 100000
kth_closest_obs_horizontal_scale = 2
vertical_scale = 200
tpos = rep(1,N) * 16
tneg = rep(1,N) * 16
t_sod = rep(1,N) * 4
eps2 = rep(1,N) * 0.5
values_mina = values - 20
values_maxa = values + 20
values_minv = values - 1
values_maxv = values + 1
debug = T
points = Points(lats, lons, elevs)
res <- sct_resistant(points, values, obs_to_check, background_values, background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, values_mina, values_maxa, values_minv, values_maxv, eps2, tpos, tneg, debug)
# check the results of sct woth the following OI
#  first create the data_inner.txt file from sct output with debug=T
#d<-read.table(file="data_inner.txt",header=F,stringsAsFactors=F,strip.white=T)
#lats<-d$V4
#lons<-d$V3
#x<-(lons-10)/0.005*278.329 #278.329*0:(length(lats)-1)
#disth2<-outer(x,x,FUN="-")**2
#dh2<-278.329**2
#S <- exp( -0.5*( disth2 / dh2))
#SRinv <- chol2inv( chol( ( S + diag( x=0.5, length(lats)))))
#yb <- rep( median(d$V5), length(lons))
#yo<-d$V5
#inno <- yo - yb
#SRinv_di <- crossprod( SRinv, inno) 
#ya  <- yb + S %*% SRinv_di
#yav <- yo - 1 / diag( SRinv) * SRinv %*% inno
#plot(d$V3,d$V5)
#lines(d$V3,yb)
#points(d$V3,ya,col="blue")
#points(d$V3,d$V7,col="cyan")
#points(d$V3,yav,col="red")
#points(d$V3,d$V8,col="red")
#points(d$V3,d$V8,col="pink")
#mu <- median(chi)
#sigma <- as.numeric( diff( quantile( chi, probs=c( 0.25, 0.75))))
#sigma_mu <- sigma / length(lats)
#--------------------------------------------------------
# Test with larger vectors
P = 5000
pGE = 0.3 # probability of gross error 0.3 = 30%
lats = runif(P, min = 55, max = 70)
lons = runif(P, min = 5, max = 30)
elevs = runif(P, min = 0, max = 2500)
# simple vertical profile
values <- 30 - 0.0065 * elevs
#values <- (103-1.333*lats) - 0.0065 * elevs
idx <- sample(1:P,ceiling(P*pGE))
true_flags <- values; true_flags[] <- 0
true_flags[idx] <- 1
values[idx]<-runif(ceiling(P*pGE), min = -50, max = 50)
obs_to_check = rep(1,P)
background_values = 0
background_elab_type = "VerticalProfileTheilSen"
tpos = rep(3,P)
tneg = rep(3,P)
eps2 = rep(0.5,P)
values_mina = values - 20
values_maxa = values + 20
values_minv = values - 1
values_maxv = values + 1
debug = F
num_min_outer = 3
num_max_outer = 50
inner_radius = 30000
outer_radius = 50000
num_iterations = 20
num_min_prof = 10
min_elev_diff = 500
min_horizontal_scale = 500
max_horizontal_scale = 10000
kth_closest_obs_horizontal_scale = 3
vertical_scale = 600
points = Points(lats, lons, elevs)
t0<-Sys.time()
res<-sct_resistant(points, values, obs_to_check, background_values, background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, values_mina, values_maxa, values_minv, values_maxv, eps2, tpos, tneg, debug)
t1<-Sys.time()
print(t1-t0)
flags<-res[[1]]
score<-res[[2]]
a <- length( which( true_flags == 1 & flags == 1))
c <- length( which( true_flags == 1 & flags == 0))
b <- length( which( true_flags == 0 & flags == 1))
d <- length( which( true_flags == 0 & flags == 0))
rand <- (a+c) * (a+b) / (a+b+c+d)
ets <- (a-rand) / (a+b+c-rand)
acc <- (a+d)/(a+b+c+d)
pod <- a/(a+c)
pofa <- b/(b+d)
print( paste("a(bad) b c d", a,"(",length(which( true_flags == 1)),")", b, c, d, a+b+c+d))
print( paste("acc pod pofa ets", round(acc,2), round(pod,2), round(pofa,2), round(ets,2)))

#save.image("test_sct.rdata")

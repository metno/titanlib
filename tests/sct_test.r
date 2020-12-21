#titanlib_path <- "/home/cristianl/projects/titanlib/build/extras"
titanlib_path <- "/home/cristianl/projects/titanlib/build/extras"
#titanlib_path <- "/home/lineb/projects/titanlib/titanlib/build/extras"
dyn.load( file.path( titanlib_path, 
                     paste("SWIG/R/titanlib", .Platform$dynlib.ext, sep="")))
source(   file.path( titanlib_path,"SWIG/R/titanlib.R"))

#---------------------------------------------------------------
# Test with small vectors
lats = c(60, 60, 60, 60, 60)
lons = c(10, 10.005, 10.01, 10.015, 10.02)
elevs = c(0, 0, 0, 0, 0)
values = c(0, 1, 10, -10, -100)
#values = c(0, 0, 0, 0, -100)
obs_to_check = c(1, 1, 1, 1, 1)
background_values = c(0, 0, 0, 0, 0)
background_elab_type = "vertical_profile"
N = length(lats)
num_min_outer = 3
num_max_outer = 10
inner_radius = 20000
outer_radius = 50000
num_iterations = 10
num_min_prof = 0
min_elev_diff = 100
min_horizontal_scale = 10000
max_horizontal_scale = 100000
kth_closest_obs_horizontal_scale = 2
vertical_scale = 200
tpos = rep(1,N) * 16
tneg = rep(1,N) * 16
t_sod = rep(1,N) * 4
eps2 = rep(1,N) * 0.5
value_min = -50
value_max = 50
values_mina = values - 20
values_maxa = values + 20
values_minv = values - 1
values_maxv = values + 1
debug = T
res<-sct(lats, lons, elevs, values, obs_to_check, background_values, background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, value_min, value_max, values_mina, values_maxa, values_minv, values_maxv, eps2, tpos, tneg, debug)
print(res[[1]])

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
background_elab_type = "vertical_profile_Theil_Sen"
tpos = rep(3,P)
tneg = rep(3,P)
eps2 = rep(0.5,P)
value_min = -50
value_max = 50
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

t0<-Sys.time()
res<-sct(lats, lons, elevs, values, obs_to_check, background_values, background_elab_type, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, num_min_prof, min_elev_diff, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, value_min, value_max, values_mina, values_maxa, values_minv, values_maxv, eps2, tpos, tneg, debug)
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

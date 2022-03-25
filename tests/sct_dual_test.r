titanlib_path <- "../build"
dyn.load( file.path( titanlib_path, 
                     paste("swig/R/titanlib", .Platform$dynlib.ext, sep="")))
source(   file.path( titanlib_path,"swig/R/titanlib.R"))

#---------------------------------------------------------------
# Test with small vectors
lats <- c(60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60)
lons <- 10 + 0.005 * 0:(length(lats)-1)
elevs <- rep(0, length(lats))
#values <- round( 10*sin(lons*2*base::pi/(max(lons)-min(lons))), 2)
values <- rep( 10, length( lats))
values[5] <- 0
values[10] <- 0
#values = c(0, 0, 0, 0, -100)
obs_to_check = rep(1, length(lats))
event_thresholds = 0.1
test_thresholds = 0.8
condition = "Gt"
N = length(lats)
num_min_outer = 3
num_max_outer = 10
inner_radius = 20000
outer_radius = 50000
num_iterations = 10
min_horizontal_scale = 250 
max_horizontal_scale = 100000
kth_closest_obs_horizontal_scale = 2
vertical_scale = 200
debug = T
points = Points(lats, lons, elevs)
res <- sct_dual(points, values, obs_to_check, event_thresholds, condition, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, test_thresholds, debug)
print(res)
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
P = 100000
pGE = 0.01 # probability of gross error 0.3 = 30%
lat_mn <- 55
lat_mx <- 70
lon_mn <- 5
lon_mx <- 30
lats = runif( P, min = lat_mn, max = lat_mx)
lons = runif( P, min = lon_mn, max = lon_mx)
elevs = runif( P, min = 0, max = 2500)
# 
values <- rep(10,P)
values[which(lons>=(lon_mn+(lon_mx-lon_mn)/2))] <- 0
#
true_flags <- values; true_flags[] <- 0
values_or <- values
idx <- sample( which( values_or == 10), ceiling(P*pGE/2))
values[idx] <- 0
true_flags[idx] <- 1
idx <- sample( which( values_or == 0), ceiling(P*pGE/2))
values[idx] <- 10
true_flags[idx] <- 1
#
obs_to_check = rep(1,P)
event_thresholds = 0.1
test_thresholds = 0.5
condition = "Gt"
debug = F
num_min_outer = 3
num_max_outer = 50
inner_radius = 30000
outer_radius = 50000
num_iterations = 20
min_horizontal_scale = 500
max_horizontal_scale = 10000
kth_closest_obs_horizontal_scale = 3
vertical_scale = 600
points = Points(lats, lons, elevs)
t0<-Sys.time()
res <- sct_dual(points, values, obs_to_check, event_thresholds, condition, num_min_outer, num_max_outer, inner_radius, outer_radius, num_iterations, min_horizontal_scale, max_horizontal_scale, kth_closest_obs_horizontal_scale, vertical_scale, test_thresholds, debug)
t1<-Sys.time()
print(t1-t0)
flags<-res
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

q()
#save.image("test_dual_sct.rdata")
# plots from "test_dual_sct.rdata"
load("test_dual_sct.rdata")
png(file="test_dual_sct.png",width=800,height=800)
par(mar=c(2,2,1,1))
plot( lons, lats, col="white",xlab="",ylab="")
ix<-which(values==10)
points(lons[ix],lats[ix],cex=1.5,col="darkblue",pch=21,bg="blue")
ix<-which(values==0)
points(lons[ix],lats[ix],cex=1.5,col="darkred",pch=21,bg="red")
ix<-which(values==0 & true_flags==1)
points(lons[ix],lats[ix],cex=1.5,col="darkred",pch=21,bg="red")
ix<-which(values==10 & true_flags==1)
points(lons[ix],lats[ix],cex=1.5,col="darkblue",pch=21,bg="blue")
ix<-which(flags==1 & true_flags==1)
points(lons[ix],lats[ix],cex=1.5,col="gold",pch=4)
ix<-which(flags==1 & true_flags==0)
points(lons[ix],lats[ix],cex=1.5,col="white",pch=4)
ix<-which(flags==0 & true_flags==1)
points(lons[ix],lats[ix],cex=1.5,col="white",pch=0)
dev.off()

png(file="sct_dual_If1wrtf2.png",width=800,height=800)
par(mar=c(5,5,1,1))
f1<-0.9
f2<-seq(0.01,1,by=0.001)
m<-(-1)
plot( f2, f1*log(f1/f2), lwd=3, main="", xlab="f2", ylab="f1=0.9, f1 * ln( f1 / f2)",cex.lab=2,cex.axis=1.5,col="white")
lines(c(0,f2),c(f1,-f2+f1))
abline(h=seq(-10,10,by=0.1),v=seq(-10,10,by=0.1),lty=2,col="gray")
abline(h=0,v=0)
lines( f2, f1*log(f1/f2), col="black",lwd=5)
dev.off()

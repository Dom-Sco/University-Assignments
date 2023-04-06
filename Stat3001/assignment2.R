#data
x = c(8.23, 7.58, 7.39, 9.02, 6.69, 8.05, 8.38, 8.03, 9.54, 12.10)
n = length(x)
#calculate Fn and F
xgrid = seq(floor(min(x)), ceiling(max(x)), by=0.01)
Fn = c()
for(i in 1:length(xgrid)) Fn[i] = mean(x <= xgrid[i]) #proportion of X’s less than x
plot(xgrid, Fn, type="n", xlab="x", ylab="F_N")
lines(xgrid, Fn)

#add the true cdf
cdf <- pnorm(sort(xgrid), mean=9, sd=1)
lines(xgrid, cdf, col="red")

#simulate Dn
dn = NULL; K = 100000;
for(k in 1:K) {
  u = runif(n, min=0, max=1);  #random uniform sample
  i = 1:n;                     #index
  u.sorted = sort(u);          #sort u from smallest to largest
  dn[k] = max( max(abs(u.sorted-i/n)), max(abs(u.sorted-(i-1)/n)) ) 
}
#test statistic
dn.max = max(abs(Fn-cdf)); dn.max
#p-value
p.value= mean(dn >= dn.max); p.value


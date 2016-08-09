
# plot the data from quadtree::testr() 
# in order to verify results 

x <- c(0.1, 0.15, 0.4, 0.6, 0.7, 0.71, 0.8, 0.74, 0.48, 0.41)
y <- c(0.6, 0.85, 0.65, 0.3, 0.6, 0.1, 0.4, 0.33, 0.42, 0.57)

xcom <- 0.62
ycom <- 0.423214

ptfar <- 2

dist <- sqrt( (x[ptfar]-xcom)^2 + (y[ptfar]-ycom)^2 )
print(dist)

x11()
plot(x,y)
points(xcom, ycom, pch=23)

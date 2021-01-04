y <- function(x,b1,b2){
  y <- b1*x + b2*I(x^2) - 2 
}


x <- 1:50
xrange <- c(0,30)
yrange <- c(0,60)
plot(x, y(x,5,-0.2), ylim=yrange, xlim=xrange)
points(x, y(x,6,-0.25), ylim=yrange, xlim=xrange )
points(x, y(x,5,-0.3), ylim=yrange, xlim=xrange)

plot(x, y(x,6,-0.2), ylim=yrange, xlim=xrange)
points(x, y(x,5.5,-0.4), ylim=yrange, xlim=xrange )
points(x, y(x,5,-0.2), ylim=yrange, xlim=xrange)

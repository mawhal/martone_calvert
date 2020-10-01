
Evar <- function( x ){
  S = length( x[x>0] )
  1 - 2/pi*atan( sum((log(x[ x>0 ]) - sum(log(x[ x>0 ]))/S)^2)/S ) 
}



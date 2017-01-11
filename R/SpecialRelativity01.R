#############
# variables #
#############
#v <- 0.9 # the speed of the moving frame of reference with respect to the frame of reference at rest expresses in light years per year (= ratio with the speed of light)

############
# packages #
############
# library(ggplot2)

#############
# constants #
#############
c <- 1 # speed of light expressed in light years per year

#############
# functions #
#############
gamma <- function(v = O) {
  # Calculates the Lorentz time dilation / length contraction factor gamma
  # Input: v --> the relative velocity compared to the frame of reference at rest in light years per year
  return (1/sqrt(1 - v^2/c^2))
}

gammaInverse <- function(gamma = 1) {
  # Calculates the velocity needed to obtain the Lorentz time dilation / length contraction factor gamma
  # Input: gamma --> the Lorentz time dilation / length contraction factor gamma
  return (c * sqrt(1 - 1/gamma^2))
}

xa <- function(xb = O, tb = 0, v = 0) {
  # Calculates the spatial coordinate in the frame of reference at rest of an event at a certain time and location in the moving frame of reference
  # Input: xb --> spatial coordinate (light years) in the frame of reference moving at speed v with respect to the rest frame of reference
  #        tb --> temporal coordinate (years) in the frame of reference moving at speed v with respect to the rest frame of reference
  #        v  --> the speed fo the moving frame of reference with respect to the rest frame of reference in light years per year
  return (gamma(v) * (xb + v*tb))
}

ta <- function(xb = O, tb = 0, v = 0) {
  # Calculates the temporal coordinate in the frame of reference at rest of an event at a certain time and location in the moving frame of reference
  # Input: xb --> spatial coordinate (light years) in the frame of reference moving at speed v with respect to the rest frame of reference
  #        tb --> temporal coordinate (years) in the frame of reference moving at speed v with respect to the rest frame of reference
  #        v  --> the speed fo the moving frame of reference with respect to the rest frame of reference in light years per year
  return (gamma(v) * (tb + v*xb/c^2))
}

plotSpaceTimeDiagram <- function(v = 0, t = -10:10, x = -10:10, color = 'red') {
  # Plots the spacetime diagram of the rest and the moving frame of reference over a certain time range (specified for the rest frame of reference)
  # Input: v     --> the speed fo the moving frame of reference with respect to the rest frame of reference
  #        t     --> the array of points in time in the rest frame of reference which are to be shown on the time axis (in years)
  #        x     --> the array of points in space in the rest frame of reference which are to be shown on the space (x) axis (in light years)
  #        color --> the color in which the moving frame of reference is to be displayed
  stopifnot(v < 1 && -1 < v) # v must be between -1 and +1
  lx <- length(x)
  lt <- length(t)
  print(paste("Gamma factor =",gamma(v = v)))
  # Create the plot -- just plotting the end points of the supplied x and t vector
  plot(x = c(x[1],x[lx],rep(0,2)), y = c(rep(0,2),t[1],t[lt])
       , xlab = 'spacial coordinate rest frame (in lightyears)'
       , ylab = 'temporal coordinate rest frame (in years)')
  title(main = paste('space time diagram velocity =',v,'times c'), sub = paste('gamma =',gamma(v)))
  # add the rest frame axis
  lines(x = c(x[1],x[lx]), y = rep(0,2), lwd = 2)
  lines(x = rep(0,2), y = c(t[1],t[lt]), lwd = 2, lty = 3)
  # add the light cone
  lines(x = c(-100,100),y = c(-100,100), lwd = 0.5, col = 'green')
  lines(x = c(-100,100),y = c(100,-100), lwd = 0.5, col = 'green')
  # add the rest frame lines of identical time and location
  for (xi in x) {
    lines(x = rep(xi,2), y = c(t[1],t[lt]), lwd = 1, lty = 3)
  }
  for (ti in t) {
    lines(x = c(x[1],x[lx]), y = rep(ti,2), lwd = 1)
  }
  # add the moving frame axis
  lines(x = c(xa(x[1],0,v),xa(x[lx],0,v)), y = c(ta(x[1],0,v),ta(x[lx],0,v)), lwd = 2, col = color)
  lines(x = c(xa(0,t[1],v),xa(0,t[lt],v)), y = c(ta(0,t[1],v),ta(0,t[lt],v)), lwd = 2, col = color, lty = 3)
  # add the moving frame lines of identical time and location
  # extend the moving coordinates to make sure there are sufficient equal distance and time lines
  x_alt <-seq.default(from = -abs(max(x)), to = max(x), by = x[2] - x[1])
  t_alt <-seq.default(from = -abs(max(t)), to = max(t), by = t[2] - t[1])
  lx_alt <- length(x_alt)
  lt_alt <- length(t_alt)
  for (xi in x_alt) {
    lines(x = c(xa(xi,t_alt[1],v),xa(xi,t_alt[lt_alt],v)), y = c(ta(xi,t_alt[1],v),ta(xi,t_alt[lt_alt],v)), lwd = 1, col = color, lty = 3)
  }
  for (ti in t_alt) {
    lines(x = c(xa(x_alt[1],ti,v),xa(x_alt[lx_alt],ti,v)), y = c(ta(x_alt[1],ti,v),ta(x_alt[lx_alt],ti,v)), lwd = 1, col = color)
  }
}

################
# working code #
################
plotSpaceTimeDiagram(v = 0.1, x = 0:5, t = 0:5)
plotSpaceTimeDiagram(v = 0.8, x = 0:10, t = c(0,1,2,3,4,5))
plotSpaceTimeDiagram(v = 0.8, x = c(0,2,4,6,8,10), t = c(0,1,2,3,4,5))
plotSpaceTimeDiagram(v = 0.8, x = c(0,2,4,6,8,10), t = seq.default(from = 0, to = 10, by = 0.5))
plotSpaceTimeDiagram(v = 0.8, x = seq.default(from = 0, to = 2, by = 0.5), seq.default(from = 0, to = 2, by = 0.5))
plotSpaceTimeDiagram(v = 0.942809, x = seq.default(from = -2, to = 6, by = 1), seq.default(from = -2, to = 6, by = 1))

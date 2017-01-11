#############
# variables #
#############

#############
# packages  #
#############

#############
# constants #
#############
c <- 1 # speed of light expressed in light years per year

#############
# classes   #
#############
SpaceTimeDiagram <- setRefClass("SpaceTimeDiagram",
                                fields = list(v = "numeric",
                                              x = "numeric",
                                              t = "numeric",
                                              edge = "numeric"
                                             ),
                                methods = list(initialize = function (vi,xi,ti) {
                                                 # initialise an instance of this class
                                                 stopifnot(vi < 1 && -1 < vi) # v must be between -1 and +1
                                                 v <<- vi
                                                 if (0 < xi[1]) {
                                                   xi <- seq.default(from = 0, to = xi[length(xi)], by = xi[2] - xi[1])
                                                 }
                                                 if (xi[length(xi)] < 0) {
                                                   xi <- seq.default(from = xi[1], to = 0, by = xi[2] - xi[1])
                                                 }
                                                 x <<- xi
                                                 if (0 < ti[1]) {
                                                   ti <- seq.default(from = 0, to = ti[length(ti)], by = ti[2] - ti[1])
                                                 }
                                                 if (ti[length(ti)] < 0) {
                                                   ti <- seq.default(from = ti[1], to = 0, by = ti[2] - ti[1])
                                                 }
                                                 t <<- ti
                                                 edge <<- c(c(min(xi),min(ti)),c(max(xi),min(ti)),c(max(xi),max(ti)),c(min(xi),max(ti)))
                                               },
                                               summary = function () {
                                                 return (paste('Space Time Diagram with velocity:', v))},
                                               gamma = function (){
                                                 # Calculates the Lorentz time dilation / length contraction factor gamma
                                                 return (1/sqrt(1 - v^2/c^2))
                                               },
                                               lorentz_transform_x = function(x, t) {
                                                 # Calculates the spatial coordinate in the frame of reference at rest of an event at a certain time and location in the moving frame of reference
                                                 # Input: x --> spatial coordinate (light years) in the frame of reference moving at speed v with respect to the rest frame of reference
                                                 #        t --> temporal coordinate (years) in the frame of reference moving at speed v with respect to the rest frame of reference
                                                 return (gamma() * (x + v*t))
                                               },
                                               lorentz_transform_t = function(x, t) {
                                                 # Calculates the temporal coordinate in the frame of reference at rest of an event at a certain time and location in the moving frame of reference
                                                 # Input: x --> spatial coordinate (light years) in the frame of reference moving at speed v with respect to the rest frame of reference
                                                 #        t --> temporal coordinate (years) in the frame of reference moving at speed v with respect to the rest frame of reference
                                                 return (gamma() * (t + v*x/c^2))
                                               },
                                               plotSpaceTimeDiagram = function(color = 'red') {
                                                 # Plots the spacetime diagram of the rest and the moving frame of reference over a certain time range (specified for the rest frame of reference)
                                                 # Input: color --> the color in which the moving frame of reference is to be displayed
                                                 lx <- length(x)
                                                 lt <- length(t)
                                                 # Create the plot -- just plotting the edges and the origin
                                                 plot(x = c(0,edge[c(1,3,5,7)]), y = c(0,edge[c(2,4,6,8)])
                                                      ,xlab = 'spacial coordinate rest frame (in lightyears)'
                                                      ,ylab = 'temporal coordinate rest frame (in years)')
                                                 title(main = paste('space time diagram velocity =',v,'c'), sub = paste('gamma =',gamma()))
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
                                                 lines(x = c(lorentz_transform_x(x[1],0),lorentz_transform_x(x[lx],0)), y = c(lorentz_transform_t(x[1],0),lorentz_transform_t(x[lx],0)), lwd = 2, col = color)
                                                 lines(x = c(lorentz_transform_x(0,t[1]),lorentz_transform_x(0,t[lt])), y = c(lorentz_transform_t(0,t[1]),lorentz_transform_t(0,t[lt])), lwd = 2, col = color, lty = 3)
                                                 # add the moving frame lines of identical time and location
                                                 # extend the moving coordinates to make sure there are sufficient equal distance and time lines
                                                 x_alt <-seq.default(from = -abs(max(x)), to = max(x), by = x[2] - x[1])
                                                 t_alt <-seq.default(from = -abs(max(t)), to = max(t), by = t[2] - t[1])
                                                 lx_alt <- length(x_alt)
                                                 lt_alt <- length(t_alt)
                                                 for (xi in x_alt) {
                                                   lines(x = c(lorentz_transform_x(xi,t_alt[1]),lorentz_transform_x(xi,t_alt[lt_alt])), y = c(lorentz_transform_t(xi,t_alt[1]),lorentz_transform_t(xi,t_alt[lt_alt])), lwd = 1, col = color, lty = 3)
                                                 }
                                                 for (ti in t_alt) {
                                                   lines(x = c(lorentz_transform_x(x_alt[1],ti),lorentz_transform_x(x_alt[lx_alt],ti)), y = c(lorentz_transform_t(x_alt[1],ti),lorentz_transform_t(x_alt[lx_alt],ti)), lwd = 1, col = color)
                                                 }
                                               }
                                              )
                               )

#############
# functions #
#############

################
# working code #
################
spaceTimeDiagram30k <- SpaceTimeDiagram$new(v = 0.1, x = 2:5, t = -5:-2)
spaceTimeDiagram240k <- SpaceTimeDiagram$new(v = 0.8, x = seq.default(from = 0, to = 2, by = 0.5), t = seq.default(from = 0, to = 2, by = 0.5))

spaceTimeDiagram240k$v
spaceTimeDiagram30k$x
spaceTimeDiagram30k$t
spaceTimeDiagram240k$summary()
spaceTimeDiagram240k$gamma()
spaceTimeDiagram240k$lorentz_transform_x(rep(0,6),c(0,1,2,3,4,5))
spaceTimeDiagram240k$lorentz_transform_t(rep(0,6),c(0,1,2,3,4,5))
spaceTimeDiagram240k$lorentz_transform_x(c(0,1,2,3,4,5),rep(0,6))
spaceTimeDiagram240k$lorentz_transform_t(c(0,1,2,3,4,5),rep(0,6))
spaceTimeDiagram240k$edge
spaceTimeDiagram240k$plotSpaceTimeDiagram(color = 'purple')
spaceTimeDiagram30k$edge
spaceTimeDiagram30k$plotSpaceTimeDiagram()
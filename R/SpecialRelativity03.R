#' @export SpaceTimeDiagram
 
##############
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
                                fields = list(v           = "numeric",
                                              x           = "numeric",
                                              t           = "numeric",
                                              voyage_dist = "numeric",
                                              voyage_time = "numeric",   
                                              edge        = "numeric"
                                ),
                                methods = list(
                                  initialize = function (v,x_to,t_to,voyage=t_to,x_step=1,t_step=1,x_negative=0) {
                                    # initialise an instance of this class
                                    # initial verifications
                                    stopifnot(all(v < 1) && all(-1 < v)) # v must be between -1 and +1
                                    stopifnot(0 < x_to || 0 < t_to || 0 < x_step || 0 < t_step) # additional requirements
                                    stopifnot(length(v) == length(voyage)) # v must have the same length as voyage(age)
                                    # assign speed vector
                                    v <<- v
                                    # assign spatial coordinate vector
                                    if (x_negative == 0) {
                                      x <<- seq.default(from = 0, to = x_to, by = x_step)
                                    }
                                    else {
                                      x <<- seq.default(from = - x_to, to = x_to, by = x_step)
                                    }
                                    # assign temporal coordiante vector
                                    t <<- seq.default(from = 0, to = t_to, by = t_step)
                                    # assign edge of the coordinates
                                    edge <<- c(min(x),min(t),max(x),min(t),max(x),max(t),min(x),max(t))
                                    # assign the spatial and temporal coordinates where the voyage is at the moment when the speed changes in the rest frame of reference
                                    vl <- length(voyage)
                                    distance <- v * voyage # the distance per speed section
                                    distance_total <- rep(0,vl)
                                    voyage_total <- rep(0,vl)
                                    for (i in 1:vl) {
                                      distance_total[i] <- sum(distance[1:i]) # the cumulative distance for each speed section
                                      voyage_total[i] <- sum(voyage[1:i]) # the cumulative time for each speed section
                                    }
                                    voyage_dist <<- c(0,distance_total)
                                    voyage_time <<- c(0,voyage_total)
                                  },
                                  summary = function () {
                                    return (paste('Space Time Diagram with velocity:', v))},
                                  gamma = function (){
                                    # Calculates the Lorentz time dilation / length contraction factor gamma
                                    return (1/sqrt(1 - v^2/c^2))
                                  },
                                  lorentz_transform_x = function(x, t, i = 0) {
                                    # Calculates the spatial coordinate in the frame of reference at rest of an event at a certain time and location in the moving frame of reference
                                    # Input: x --> spatial coordinate (light years) in the frame of reference moving at speed v with respect to the rest frame of reference
                                    #        t --> temporal coordinate (years) in the frame of reference moving at speed v with respect to the rest frame of reference 
                                    if (i == 0) {
                                      stopifnot(length(x) == length(t) && length(x) == length(v))
                                      return (gamma() * (x + v*t))
                                    }
                                    stopifnot(1 <= i && i <= length(v) && round(i) == i)
                                    return (gamma()[i] * (x + v[i]*t))
                                  },
                                  lorentz_transform_t = function(x, t, i = 0) {
                                    # Calculates the temporal coordinate in the frame of reference at rest of an event at a certain time and location in the moving frame of reference
                                    # Input: x --> spatial coordinate (light years) in the frame of reference moving at speed v with respect to the rest frame of reference
                                    #        t --> temporal coordinate (years) in the frame of reference moving at speed v with respect to the rest frame of reference
                                    if (i == 0) {
                                      stopifnot(length(x) == length(t) && length(x) == length(v))
                                      return (gamma() * (t + v*x/c^2))
                                    }
                                    stopifnot(1<= i && i <= length(v) && round(i) == i)
                                    return (gamma()[i] * (t + v[i]*x/c^2))
                                  },
                                  voyage_add = function() {
                                    # Plots the voyage and the light cone in the frame of reference at rest
                                    lines(x = voyage_dist, y = voyage_time, lwd = 3, col = "blue")
                                    # add the light cone
                                    vl <- length(voyage_dist) - 1
                                    for (i in 1:vl) {
                                      lines(x = c(voyage_dist[i],voyage_dist[i+1]), y = c(voyage_time[i],voyage_time[i] + abs(voyage_dist[i+1] - voyage_dist[i])), lwd = 3, col = "green")
                                      #lines(x = c(c(0,distance_total)[i],distance_total[i]), y = c(c(0,voyage_total)[i],c(0,voyage_total)[i] + abs(distance[i])), lwd = 3, col = "green")
                                    }
                                  },
                                  moving_frame = function(color) {
                                    vl <- length(voyage_dist) - 1
                                    spatial_axis_x_coordinate <- 1/gamma() * lorentz_transform_x(voyage_dist[2:(vl+1)] - voyage_dist[1:vl],rep(0,vl))
                                    spatial_axis_t_coordinate <- 1/gamma() * lorentz_transform_t(voyage_dist[2:(vl+1)] - voyage_dist[1:vl],rep(0,vl))
                                    temporal_axis_x_coordinate <- 1/gamma() * lorentz_transform_x(rep(0,vl),voyage_time[2:(vl+1)] - voyage_time[1:vl])
                                    temporal_axis_t_coordinate <- 1/gamma() * lorentz_transform_t(rep(0,vl),voyage_time[2:(vl+1)] - voyage_time[1:vl])
                                    for (i in 1:vl) {
                                      # spatial axis (t is 0 in moving reference frame)
                                      lines(x = c(voyage_dist[i],voyage_dist[i] + spatial_axis_x_coordinate[i]), y = c(voyage_time[i],voyage_time[i] + spatial_axis_t_coordinate[i]), lwd = 2, col = color)
                                      # equal distance to spatial axis
                                      lines(x = c(voyage_dist[i+1],voyage_dist[i+1] - spatial_axis_x_coordinate[i]), y = c(voyage_time[i+1],voyage_time[i+1] - spatial_axis_t_coordinate[i]), lwd = 1, col = color)
                                      #t_diff <- voyage_time[i+1] - voyage_time[i]
                                      #t_step <- t[2] - t[1]
                                      #ni <- floor(t_diff / t_step)
                                      #for (step_ti in 1:ni) {
                                        #lines(x = sp_x_coor, y = sp_y_coor + step_ti * t_step, lwd = 1, col = color)
                                      #}
                                      # temporal axis (overlaps with the voyage = x is 0 in moving reference frame)
                                      lines(x = c(voyage_dist[i],voyage_dist[i] + temporal_axis_x_coordinate[i]), y = c(voyage_time[i],voyage_time[i] + temporal_axis_t_coordinate[i]), lwd = 2, col = color, lty = 3)
                                      # equal distance to temporal axis
                                      # temporal axis (overlaps with the voyage = x is 0 in moving reference frame)
                                      lines(x = c(voyage_dist[i] + spatial_axis_x_coordinate[i],voyage_dist[i] + spatial_axis_x_coordinate[i] - temporal_axis_x_coordinate[i]), y = c(voyage_time[i] + spatial_axis_t_coordinate[i],voyage_time[i] + spatial_axis_t_coordinate[i] - temporal_axis_t_coordinate[i]), lwd = 1, col = color, lty = 3)
                                      #x_diff <- voyage_time[i+1] - temp_y_coor
                                      #x_step <- x[2] - x[1]
                                      #mi <- floor(x_diff / x_step)
                                      #for (step_xi in 1:ni) {
                                      #  lines(x = temp_x_coor + step_xi * x_step, y = temp_y_coor, lwd = 1, col = color, lty = 3)
                                      #}
                                    }
                                  },
                                  plotSpaceTimeDiagram = function(color = 'red') {
                                    # Plots the spacetime diagram of the rest and the moving frame of reference over a certain time range (specified for the rest frame of reference)
                                    # Input: color --> the color in which the moving frame of reference is to be displayed
                                    lx <- length(x)
                                    lt <- length(t)
                                    # Create the plot -- just plotting the edges and the origin
                                    plot(x = edge[c(1,3,5,7)]
                                         ,y = edge[c(2,4,6,8)]
                                         ,xlab = 'spacial coordinate rest frame (in lightyears)'
                                         ,ylab = 'temporal coordinate rest frame (in years)')
                                    title(main = paste('space time diagram velocity =',v,'c'), sub = paste('gamma =',paste(gamma())))
                                    # add the rest frame axis
                                    lines(x = c(x[1],x[lx]), y = rep(0,2), lwd = 2)
                                    lines(x = rep(0,2), y = c(t[1],t[lt]), lwd = 2, lty = 3)
                                    # add the rest frame lines of identical time and location
                                    for (xi in x) {
                                      lines(x = rep(xi,2), y = c(t[1],t[lt]), lwd = 1, lty = 3)
                                    }
                                    for (ti in t) {
                                      lines(x = c(x[1],x[lx]), y = rep(ti,2), lwd = 1)
                                    }
                                    # add the voyage and the light cone
                                    voyage_add()
                                    # add the moving frame of reference
                                    moving_frame(color = color)
                                  }
                                )
)

#############
# functions #
#############

################
# working code #
################
mySTD <- SpaceTimeDiagram$new(v = c(0.8,-0.75,0,2/3), x_to = 5, t_to = 10, voyage = c(5,2,1.5,1.5), x_step = 0.5, t_step = 0.5)
mySTD <- SpaceTimeDiagram$new(v = c(0.8,-0.5), x_to = 5, t_to = 9, voyage = c(5,4), x_step = 0.5, t_step = 0.5)
mySTD <- SpaceTimeDiagram$new(v = 0.8, x_to = 5, t_to = 5, t_step = 0.5)
mySTD$plotSpaceTimeDiagram()

mySTD$voyage_dist
mySTD$voyage_time
mySTD$gamma()

mySTD$lorentz_transform_x(c(0,1,2,3,4,5),rep(0,6), i = 2)
mySTD$lorentz_transform_t(c(0,1,2,3,4,5),rep(0,6), i = 2)

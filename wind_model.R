# See Alessandrini et. al. 2012

sigma_wind_speed <- 3.7
sigma_bird_speed <- 3.5 

# TO MIKE: Define with beta
mean_bird_speed <- 40 


# To Mike: We need to add the variance, not the standard deviation
sigma_ground_speed <- sqrt(sigma_bird_speed^2 + sigma_wind_speed^2)

wind_model <- function(x, z, dt, bird_vec_angle, sigma_ground_speed) {
  
  #' The function computes the probability of a migratory route, 
  #' given light measurements, the bird's activity and the wind conditions along the route
  #' @param x List of x positions along the route (for definition, 
  #' see Bayesian Estimation of Animal Movement from Archival and Satellite Tag by Sumner et al., 2009)
  #' @param z List of z positions between the x positions (for definition see Sumner et al. 2009)
  #' @param bird_vec_angle List of angles in radians, indicating the bird's flight direction
  #' @param dt Time span during which the bird was actively flying, given by the activity sensor
  #' @param sigma_ground_speed expected standard deviation of the ground speed [m/s]
  #' @return prob_route List of log-probabilities for each point along the migratory route


  spd <- pmax.int(trackDist2(x, z), 1e-06)/dt
  spd <- ifelse(is.nan(spd),0, spd)

  list[wind_speeds, wind_direction] <- get_wind(z, wind_data, wind_positions)
  
  # To Mike: alpha hatten wir schon für die Gammaverteilung, nimm lieber theta
  #compute theta, the angel between the flight direction and the wind direction
  theta <- ifelse(bird_vec_angle > wind_directions, 
                  bird_vec_angle - wind_directions, 
				  wind_directions - bird_vec_angle)
  
  # Compute compensation
  mag_comp <- sin(alpha) * wind_speeds
  bird_speed_compensated <- sqrt(mean_bird_speed^2 - mag_comp^2)
  bird_speed_compensated <- ifelse(is.nan(compensation), 0, compensation)
  
  # Compute ground_speed 
  ground_speed <- cos(theta) * wind_speeds + bird_speed_compensated
		
  #Beta and alpha for gamma distribution
  beta <- ground_speed/(sigma_ground_speed^2)
  alpha  <- ground_speed * beta
  
  #compute probabilities of each point along the migration route
  prob_list <- dgamma(spd, alpha, beta, 
                      log = T)
					  
  prob_list <- ifelse(is.na(prob_list), 
					  -1000, prob_list)
  return(prob_list)
}

# 
# To MIKE: pass geo_indexed wind_data to function 

get_wind <- function(z, wind_data, wind_positions, sRad=2, k=9) {
  #' The function computes both wind speed and direction along the migratory route
  #' @param z List of z positions along the migratory route (for definition see Sumner et al. 2009)
  #' @param wind_data List of wind measurements 
  #' @param wind_positions List of positions where wind_data were measured
  #' @param sRad search radius around bird position in degrees
  #' @param k number of nearest wind measurements to interpolate wind at each z 
  #' @param geoIndex I have no clue. 
  
  #' @return prob_route List of log-probabilities for each point along the migratory route

  position_idx <- seq(1:(length(z[,1])))
  
  # MIKE: Remove  
  # Compute subset of wind positions: only wind measurements near the z positions are considered
  pos_subset <- lapply(position_idx, FUN = function(x){
		which(
			wind_positions[,1] <= z[x,] + sRad & 
			wind_positions[,2] <= z[x,] + sRad & 
			wind_positions[,1] >= z[x,] - sRad & 
			wind_positions[,2] >= z[x,] - sRad)})

  # Compute spherical distance from each point in z to each wind position in the subset
  dist_to_wind <- lapply(position_idx, FUN = function(x){
		from <- wind_positions[pos_subset[[x]],]
		to <- z[x,]
		suppressMessages(geodist(from, to, measure='cheap'))
  })
  
  # Use the distance to find the k-nearest wind positions 
  k_nn_data <- lapply(position_idx, function(x){
 
    # Find the value of the k smallest distances 
    # (partial = k implies  that the sorting stops once the kth nearest neighbor has been found)
    k_th <- sort(dist[[x]], partial = k)[k]
    # Find the indices of all elements smaller than k_nn (not ordered!!)
    k_nn_idx <- which(dist[[x]] <= k_th)
    #Filter the nine closest distances
    k_nn_dist<- dist[[x]][k_nn_idx]
    
	cbind(k_nn_idx, k_nn_dist)})
  
  # Wind measurments are taken at different altitude levels
  # filter the wind speeds at the correct height at each position
  # return wind in x and y direction
 
  wind_vector <- lapply(position_idx, FUN= function(x){
  
	nn_positions <- pos_subset[[x]][k_nn_data[[x]][,1]]
    
	wind_x <- lapply(nn_positions, FUN= function(p){
		wind_data[geoIndex[[p]][x],2]})
    
	wind_y <- lapply(nn_positions,FUN=function(p){
		wind_data[geoIndex[[p]][x],3]})
    
	#bind the x and y column to one data frame
    cbind(unlist(wind_x), unlist(wind_y))})
  
  # Wind measurements are weighted inverse to their distance to z
  dist_min <- lapply(position_idx, FUN=function(x){
    min(k_nn_data[[x]][,2])})
  
  dist_max <- lapply(position_idx, FUN=function(x){
    max(k_nn_data[[x]][,2])})
  
  weights <- lapply(position_idx, FUN=function(x){
		(k_nn_data[[x]][,2] - dist_max[[x]]) / 
		(dist_min[[x]] - dist_max[[x]])})
  
  #Sum the weights at each position along the track
  sum_weights <- lapply(weights, sum)
  
  # To Mike: I made this more compact. Check.
  # Compute the weighted wind on each position along the chain
  wind_vector_weighted <- lapply(position_idx, FUN=function(x){
	apply(wind_vector[[x]] * weights[[x]], 2,sum) / sum_weights[[x]]
  })
    
  
  #Extract wind speed [m/s] and wind direction at each position
  wind_speeds <- unlist(lapply(position_idx, FUN=function(x){
		sqrt(wind_vector_weighted[[x]][[1]]^2 + 
		     wind_vector_weighted[[x]][[2]]^2) * 3.6}),recursive = T, use.names = T)
  
  wind_directions <- unlist(lapply(position_idx,FUN = function(x){
		atan2(wind_vector_weighted[[x]][[2]], 
		      wind_vector_weighted[[x]][[1]])}), recursive = T, use.names = T)
  
  return(list(wind_speeds = wind_speeds, 
              wind_directions = wind_directions))
 }
 
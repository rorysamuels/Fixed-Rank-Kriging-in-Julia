# Functions to for the construction, placement, and evaluation of basis
# functions and basis function matrix.

# Common choices for basis functions. Here, r is the scale parameter.
bisquare(d, r) = d <= (r) ? (1 - ( d / (r) )^2 )^2 : 0
gaussian(d,r) = exp( -(d^2) / (2*r^2) )
exponential(d,r) = exp( -d/r )
matern32(d,r) = ( 1 + (sqrt(3)*d)/r )*exp( -(sqrt(3)*d)/r )


# Create a set n points in cartesian 3d space, evenly distributed on the sphere,
# and convert to degrees lat/lon.
# Taken from http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
function Create_basis_centers(n)
	goldenRatio = (1 + 5^0.5)/2;
	i = 1:n;
	phi = acos.(round.(collect(1 .- 2 .* i ./ n), digits=7));
	theta = collect(2 * pi * i ./ goldenRatio);
	x = cos.(theta) .* sin.(phi);
	y = sin.(theta) .* sin.(phi);
	z = cos.(phi);
	# Unable to convert to lat/lon when z is directly on the boundary.
	r = 6377.9999;
	Cartesian3d_to_latlon(hcat(r * x, r * y, r * z))
end


# Create a rectangluar grid of (n_lat x n_lon) evenly distributed points.
# Returns (n_lat x n_lon) x 2 matrix of lat/lon coordinates.
function Create_basis_centers(n_lat, n_lon, lat, lon)

	lat_points = collect(range(lat..., length = n_lat))
	lon_points = collect(range(lon..., length = n_lon))

	hcat(repeat(lat_points, inner = length(lon_points)),
		 repeat(lon_points, outer = length(lat_points)))
end



# Evaluate basis functions of a certain class 'type', given scale parameters,
# centers, and data locations. For 'type', any named function can be referenced,
# even user defined functions that are not specified in the above 4.
function Eval_basis(U, C, scales; basis_type = "bisquare")

	Ut = U'
	Ct = C'
	n_obs = size(U,1)
	n_centers = size(C, 1)

	distance_matrix = spzeros(n_obs, n_centers)
	S = spzeros(n_obs, n_centers)

	# Define f() to be the requested class of basis functions. Must exist in the
	# global environment.
	f(d,r) = getfield(Main, Symbol(basis_type))(d,r)

	for j in 1:n_centers
		for i in 1:n_obs
			distance_matrix[i,j] = hsdist(Ut[:,i], Ct[:,j])
		end
		S[:,j] = f.(distance_matrix[:,j],scales[j])
	end

	return S
end


# Manually construct basis function matrix; i.e. number of centers at each
# resolution and the scale parameter value at each resolution are both
# specified by the user. Returns both a matrix and a DataFrame containing
# information about funcitons.

function Create_basis(U, n_centers::AbstractVector, scales::AbstractVector; basis_type = "bisquare", prune0=true)

	Ut = U'
	n_obs = size(U,1)
	n_res = size(n_centers,1)

	basis_centers = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [2])
	min_dists = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [1])
	S = BlockArray(spzeros(n_obs,sum(n_centers)), [n_obs], n_centers)

	for k in 1:n_res
		basis_centers[Block(k,1)] = Create_basis_centers(n_centers[k])
		S[Block(1,k)] = Eval_basis(U,basis_centers[Block(k,1)],rep([scales[k]],n_centers[k]), basis_type =basis_type)
	end

	S, basis_centers = Array(S), Array(basis_centers)

	# If prune0, remove basis functions which evaluate to zero at every location.
	if prune0
	    ind = vec(mapslices(col -> sum(col) .> .001, S, dims = 1))
	else
		ind = : ;
	end

	# Return S matrix and data frame containing basis function centers, the
	# value of their scale parameters, and their resolution group.
	return S[:, ind],
	       DataFrame(basis_lat = basis_centers[:,1],
			   		 basis_lon = basis_centers[:,2],
	           		 scale = rep(scales,n_centers),
	            	 res = rep(1:n_res, n_centers))[ind,:]
end



# Manually construct basis function matrix; i.e. number of centers at each
# resolution and the scale parameter value at each resolution are both
# specified by the user. This method of Create_basis() allows for specification
# of a lat/lon boundary in which to place the centers of the basis functions.
# Returns both a matrix and a DataFrame containing information about funcitons.

function Create_basis(U, n_centers::AbstractVector, scales::AbstractVector, lat_range, lon_range,; basis_type = "bisquare", prune0=true)

	Ut = U'
	n_obs = size(U,1)
	n_res = size(n_centers,1)

	n_latlon = split_int.(n_centers)
	n_centers = Int.(prod.(n_latlon))

	basis_centers = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [2])
	min_dists = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [1])
	S = BlockArray(spzeros(n_obs,sum(n_centers)), [n_obs], n_centers)

	for k in 1:n_res
		basis_centers[Block(k,1)] = Create_basis_centers(Int(n_latlon[k][1]),Int(n_latlon[k][2]),lat_range,lon_range)
		S[Block(1,k)] = Eval_basis(U,basis_centers[Block(k,1)],rep([scales[k]],n_centers[k]), basis_type =basis_type)
	end

	S, basis_centers = Array(S), Array(basis_centers)

	# If prune0, remove basis functions which evaluate to zero at every location.
	if prune0
	    ind = vec(mapslices(col -> sum(col) .> .001, S, dims = 1))
	else
		ind = : ;
	end

	# Return S matrix and data frame containing basis function centers, the
	# value of their scale parameters, and their resolution group.
	return S[:, ind],
	       DataFrame(basis_lat = basis_centers[:,1],
			   		 basis_lon = basis_centers[:,2],
	           		 scale = rep(scales,n_centers),
	            	 res = rep(1:n_res, n_centers))[ind,:]
end




# Automatically construct basis funtion matrix; i.e. number of centers at each
# resolution and the scale parameter value at each resolution are both
# decided automatically. The user must specify the number of basis functions
# at the coarsest resolution and thenumber of resolutions.
# Returns both a matrix and a DataFrame containing information about funcitons.

function Create_basis(U, lowres_centers::Integer, n_res::Integer; basis_type = "bisquare", prune0=true)

	Ut = U'
	n_obs = size(U,1)
	n_centers = lowres_centers .* (2 .^ (0:(n_res-1)))

	basis_centers = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [2])
	min_dists = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [1])
	scales = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [1])
	S = BlockArray(spzeros(n_obs,sum(n_centers)), [n_obs], n_centers)

	for k in 1:n_res
		basis_centers[Block(k,1)] = Create_basis_centers(n_centers[k])
		min_dists[Block(k,1)] = Minimum_distances(basis_centers[Block(k,1)])
		if basis_type == "bisquare"
			scales[Block(k,1)] = 1.5*min_dists[Block(k,1)]
		else
			scales[Block(k,1)]  = min_dists[Block(k,1)]
		end
		S[Block(1,k)] = Eval_basis(U,basis_centers[Block(k,1)],scales[Block(k,1)], basis_type =basis_type)
	end

	S, basis_centers, scales = Array(S), Array(basis_centers), vec(Array(scales))

	# If prune0, remove basis functions which evaluate to zero at every location.
	if prune0
	    ind = vec(mapslices(col -> sum(col) .> .001, S, dims = 1))
	else
		ind = : ;
	end

	# Return S matrix and data frame containing basis function centers, the
	# value of their scale parameters, and their resolution group.
	return S[:, ind],
	       DataFrame(basis_lat = basis_centers[:,1],
			   		 basis_lon = basis_centers[:,2],
	           		 scale = scales,
	            	 res = rep(1:n_res, n_centers))[ind,:]
end


# Automatically construct basis function matrix; i.e. number of centers at each
# resolution and the scale parameter value at each resolution are both
# decided automatically. The user must specify the number of basis functions
# at the coarsest resolution and thenumber of resolutions.
# # This method of Auto_basis() allows for specification of a lat/lon boundary
# in which to place the centers of the basis functions.
# Returns both a matrix and a DataFrame containing information about funcitons.

function Create_basis(U, approx_lowres_centers::Integer, n_res::Integer,  lat_range, lon_range; basis_type = "bisquare", prune0=true)


	Ut = U'
	n_obs = size(U,1)

	lowres_centers = prod(split_int(approx_lowres_centers))
	approx_n_centers = lowres_centers .* (2 .^ (0:(n_res-1)))
	n_latlon = split_int.(approx_n_centers)
	n_centers = Int.(prod.(n_latlon))

	basis_centers = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [2])
	min_dists = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [1])
	scales = BlockArray(undef_blocks, Matrix{Float32}, n_centers, [1])
	S = BlockArray(spzeros(n_obs,sum(n_centers)), [n_obs], n_centers)

	for k in 1:n_res
		basis_centers[Block(k,1)] =  Create_basis_centers(Int(n_latlon[k][1]),Int(n_latlon[k][2]),lat_range,lon_range)
		min_dists[Block(k,1)] = Minimum_distances(basis_centers[Block(k,1)])
		if basis_type == "bisquare"
			scales[Block(k,1)] = 1.5*min_dists[Block(k,1)]
		else
			scales[Block(k,1)]  = min_dists[Block(k,1)]
		end
		S[Block(1,k)] = Eval_basis(U,basis_centers[Block(k,1)],scales[Block(k,1)], basis_type =basis_type)
	end

	S, basis_centers, scales = Array(S), Array(basis_centers), vec(Array(scales))

	# If prune0, remove basis functions which evaluate to zero at every location.
	if prune0
	    ind = vec(mapslices(col -> sum(col) .> .001, S, dims = 1))
	else
		ind = : ;
	end

	# Return S matrix and data frame containing basis function centers, the
	# value of their scale parameters, and their resolution group.
	return S[:, ind],
	       DataFrame(basis_lat = basis_centers[:,1],
			   		 basis_lon = basis_centers[:,2],
	           		 scale = scales,
	            	 res = rep(1:n_res, n_centers))[ind,:]
end

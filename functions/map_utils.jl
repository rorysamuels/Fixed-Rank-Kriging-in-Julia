# Utility functions useful for computing distances, transforming coordinate
# systems, etc.


# Compute great circle distance between two points using haversine formula.
# Points must be given as [lat, lon]. Final distance given in km.
function hsdist(x, y; r = 6378)
    λ₁, ϕ₁ = x
    λ₂, ϕ₂ = y
	Δϕ = ϕ₂ - ϕ₁ # Difference latitude
    Δλ = λ₂ - λ₁ # Difference longitude

    h = sind(Δλ/2)^2 + cosd(λ₁)*cosd(λ₂)*sind(Δϕ/2)^2
    2 * r * asin(sqrt(h))
end


# Change latitudes and longitudes to a 3d coordinate system so that points nearby on
# the sphere are nearby in the new coordinate system. Without this, locations near the
# boundaries are far apart when in reality they are not. V is a matrix with lats
# in the first column and lons in the second column.
# All given in degrees.
# Formulas from http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node42.html
function LatLon_to_cartesian3d(V)
	R = 6378 # radius of the Earth in km
	x = R .* cosd.(V[:,2]) .* sind.(90. .- V[:,1])
   	y = R .* sind.(V[:,2]) .* sind.(90. .- V[:,1])
	z = R .* cosd.(90. .- V[:,1])
    hcat(x, y, z)
end

# Convert cartesian coordinates to lat-lon. xyz is a matrix with x coordinate
# in the first column, y coordinate in the second column, and z coordinate in
# the third column. Returns a matrix of the corresponding lats and lons.
# Formulas derived from http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node42.html:
function Cartesian3d_to_latlon(xyz)
	x = xyz[:,1]; y = xyz[:,2]; z = xyz[:,3];
	R = 6378; # radius of the Earth in km
   	lat = 90. .- acosd.(z ./ R);
	a = sind.(90. .- lat);

	# Find indices in k for which a != 0; computation below can't divide by zero.
	f(w) = w != 0; k = findall(f,a);
	bx = x[k] ./ (R .* a[k]);
	lat = lat[k];
	# Find elements b with disallowed values.
	g1(w) = w > 1; g2(w) = w < -1;
	j1bx = findall(g1,bx);  # Unusable indices of *bx*.
	j2bx = findall(g2,bx);  # Unusable indices of *bx*.
	bx[j1bx] = floor.(bx[j1bx]); bx[j2bx] = ceil.(bx[j2bx]);

	# Function for the longitude calculations.
	h1(w) = w >= 0; h2(w) = w < 0;

	# 4 combinations of +/- z and y correspond to 4 hemispheres.
	k11 = intersect(findall(h1,z[k]), findall(h1,y[k])); # z pos, y pos NE
	k12 = intersect(findall(h1,z[k]), findall(h2,y[k])); # z pos, y neg NW
	k21 = intersect(findall(h2,z[k]), findall(h1,y[k])); # z neg, y pos SE
	k22 = intersect(findall(h2,z[k]), findall(h2,y[k])); # z neg, y neg SW

	lon = zeros(length(lat));
	lon[k11] = acosd.(bx[k11]); # ok NE
	lon[k12] = -acosd.(bx[k12]); # NW
	lon[k21] = acosd.(bx[k21]); # ok SE
	lon[k22] = -acosd.(bx[k22]); # SW

	# Combine results into single lat and lon vectors.
	hcat(lat, lon)
end

# Return minimum (nonzero) distance for each location (latitude/longitude)
# from neighboring locations.
function Minimum_distances(A)
	At = transpose(A);
	n = size(At,2);
	D = zeros(n,n);
	for j in 1:n
		s = At[:,j];
		for i in (j+1):n
			u = At[:,i]; D[j,i] = hsdist(s, u); D[i,j] = D[j,i];
		end
	end
	mapslices(col -> minimum(filter(x -> x != 0, col)), D, dims = 1)'
end

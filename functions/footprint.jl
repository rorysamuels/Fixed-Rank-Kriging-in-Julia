# Fcuntion for computing approximate footprint and evauating if there is overlap
# with BAUs.



# Construct approximate footprint as union of BAUs. Assumes BAU resolution, L_P,
# is finer than the resolution of the footprint, L_Z.
# Returns vector of string-valued hexadecimal H3 hexagon IDs that make up each
# observations footprint. All resolutions are defined by H3 grid system.
function approx_footprint(U, L_P, L_Z)
	U_Geo = mapslices(u -> GeoCoord(deg2rad.(u)...), U, dims = 2)
	index = geoToH3.(U_Geo,L_P)
	footprints = hex.(reduce(vcat,kRing.(index,(L_P - L_Z))))
end


# Return vector indicating which BAUs overlap with observation footprints.
# Both BAU_hex_ids and footprints are assumed to be vectors of string-valued
# hexadecimal H3 hexagon IDs.
function BAU_overlap(BAU_hex_ids, footprints; parallel = false)

	n_bau = size(BAU_hex_ids,1)
	c = zeros(n_bau)
	if parallel
		Threads.@threads for i in 1:n_bau
			c[i] = Int(in(BAU_hex_ids[i], footprints))
		end
	else
		for i in 1:n_bau
			c[i] = Int(in(BAU_hex_ids[i], footprints))
		end
	end
	return c
end

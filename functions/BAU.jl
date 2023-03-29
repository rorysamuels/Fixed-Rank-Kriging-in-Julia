


# Get latitude and longitude for geoocord type.
get_lat(geocoord) = rad2deg(geocoord.lat)
get_lon(geocoord) = rad2deg(geocoord.lon)

# Force preservation of hexadecimal hexago IDs (output is string).
hex(s) = string(s, base = 16)



# Return all valid hexagon IDs at resolution L as a vector.
function Get_hex_id(L)
	res0indexes = getRes0Indexes()
	childrenL = unique(reduce(vcat,h3ToChildren.(res0indexes,L)))
	return childrenL[findall((h3IsValid.(childrenL) .> 0))]
end


# Return DataFrame of all valid hexagon IDs at resolution L along with their
# latitudes and longitudes.
function BAU_grid(L; parallel = false)
	hex_id = Get_hex_id(L)
	h3_lons = get_lon.(h3ToGeo.(hex_id))
	h3_lats = get_lat.(h3ToGeo.(hex_id))

	@pipe DataFrame(lat = h3_lats,
	 				lon = h3_lons,
					hex_id = hex.(hex_id))
end


# Return DataFrame of valid hexagon IDs within a lat/lon boundaty at resolution
# L along with their latitudes and longitudes.
function BAU_grid(L, lat_range, lon_range; parallel = false)
	n0 = numHexagons(0)
	parent = getRes0Indexes()


	if parallel
		df = [DataFrame() for i in 1:n0]
		Threads.@threads for i in 1:n0
			child = @pipe unique(h3ToChildren(parent[i],L)) |>
							filter(x -> (h3IsValid.(x)),_)
			child_geo = h3ToGeo.(child)

			child_df = @pipe DataFrame(lat = get_lat.(child_geo),
									   lon = get_lon.(child_geo),
									   hex_id = hex.(child)) |>
									   	filter(:lon => lon -> is_between(lon,lon_range), _) |>
									   	filter(:lat => lat -> is_between(lat,lat_range), _)

			df[i] = child_df
		end
		return reduce(vcat,df)
	else
		df = DataFrame()
		for id in parent
			child = @pipe unique(h3ToChildren(id,L)) |>
							filter(x -> (h3IsValid.(x)),_)
			child_geo = h3ToGeo.(child)

			child_df = @pipe DataFrame(lat = get_lat.(child_geo),
									   lon = get_lon.(child_geo),
									   hex_id = hex.(child)) |>
									   	filter(:lon => lon -> is_between(lon,lon_range), _) |>
									   	filter(:lat => lat -> is_between(lat,lat_range), _)

			df = vcat(df, child_df)
		end
		return df
	end

end

# Functions for making predictions and calculating uncertainties.


# For each of the below functions:
# Compute predictions and/or mspe given new locations, data, locations, and
# estimated parameters. Return value is a DataFrame containing predictions
# and/or mspe along with the associated latitude and longitude of new locations.


## Prediction considering overlap of footprint and BAU.
# The scalar/vector, c, indicates whether or not BAU overlaps with approx
# footprint.


# Method for prediction at a single location, with single covariate.
function Calculate_pred_mspe(s0::AbstractVector, x0::AbstractFloat, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type, c; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = 1


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(reshape(s0,1,2),hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type =  basis_type)
	k_s0 = mapslices(x -> x .+ c*σ_ξ, S*K*S_s0', dims = 2)

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) and (n x n_bau).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 * α .+ k_s0'*Σ_inv*(Z-(X .* α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[1,:], lons = s0[2,:])
	elseif pred return  DataFrame(preds = vec(x0 * α .+ k_s0'*Σ_inv*(Z-(X .* α))), lons = s0[1,:], lats = s0[2,:])
	elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[1,:], lats = s0[2,:])
	end
end



# Method for prediction at a single location, with multiple covariates.
function Calculate_pred_mspe(s0::AbstractVector, x0::AbstractVector, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type, c; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = size(s0,1)


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(s0,hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type =  basis_type)
	k_s0 = mapslices(x -> x .+ c*σ_ξ, S*K*S_s0', dims = 2)

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) or (n x n).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[:,1], lons = s0[:,2])
			elseif pred return  DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), lons = s0[:,1], lats = s0[:,2])
			elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[:,1], lats = s0[:,2])
	end
end


# Method for prediction at multiple locations, with a single covariate.
function Calculate_pred_mspe(s0::AbstractMatrix, x0::AbstractVector, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type, c; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = size(s0,1)


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(s0,hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type =  basis_type)
	k_s0 = mapslices(x -> x + vec(c).*σ_ξ, S*K*S_s0', dims = 2)

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) or (n x n).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 .* α + k_s0'*Σ_inv*(Z-(X .* α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[:,1], lons = s0[:,2])
			elseif pred return  DataFrame(preds = vec(x0 .* α + k_s0'*Σ_inv*(Z-(X .* α))), lons = s0[:,1], lats = s0[:,2])
			elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[:,1], lats = s0[:,2])
	end
end



# Method for prediction at multiple locations, with a multiple covariates.
function Calculate_pred_mspe(s0::AbstractMatrix, x0::AbstractMatrix, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type, c; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = size(s0,1)


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(s0,hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type =  basis_type)
	k_s0 = mapslices(x -> x + vec(c).*σ_ξ, S*K*S_s0', dims = 2)

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) or (n x n).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[:,1], lons = s0[:,2])
			elseif pred return  DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), lons = s0[:,1], lats = s0[:,2])
			elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[:,1], lats = s0[:,2])
	end
end



## Prediction without considering overlap of footprint and BAU.


# Method for prediction at a single location, with single covariate.
function Calculate_pred_mspe(s0::AbstractVector, x0::AbstractFloat, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = 1


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(reshape(s0,1,2),hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type =  basis_type)
	k_s0 = S*K*S_s0'

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) or (n x n).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 * α .+ k_s0'*Σ_inv*(Z-(X .* α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[1,:], lons = s0[2,:])
	elseif pred return  DataFrame(preds = vec(x0 * α .+ k_s0'*Σ_inv*(Z-(X .* α))), lons = s0[1,:], lats = s0[2,:])
	elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[1,:], lats = s0[2,:])
	end
end



# Method for prediction at a single location, with multiple covariates.
function Calculate_pred_mspe(s0::AbstractVector, x0::AbstractVector, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = size(s0,1)


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(s0,hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type = basis_type)
	k_s0 = S*K*S_s0'

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) or (n x n).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[:,1], lons = s0[:,2])
			elseif pred return  DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), lons = s0[:,1], lats = s0[:,2])
			elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[:,1], lats = s0[:,2])
	end
end


# Method for prediction at multiple locations, with a single covariate.
function Calculate_pred_mspe(s0::AbstractMatrix, x0::AbstractVector, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = size(s0,1)


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(s0,hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type = basis_type)
	k_s0 = S*K*S_s0'

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) or (n x n).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 .* α + k_s0'*Σ_inv*(Z-(X .* α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[:,1], lons = s0[:,2])
			elseif pred return  DataFrame(preds = vec(x0 .* α + k_s0'*Σ_inv*(Z-(X .* α))), lons = s0[:,1], lats = s0[:,2])
			elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[:,1], lats = s0[:,2])
	end
end



# Method for prediction at multiple locations, with a multiple covariates.
function Calculate_pred_mspe(s0::AbstractMatrix, x0::AbstractMatrix, U, Z, X, S, K, α, σ_ξ, Σ_inv, BasisDF::DataFrame, basis_type; pred = true, mspe = false)

	if !(mspe || pred)
		error("Nothing will be returned. Specify 'pred = true' and/or 'mspe = true.'")
	end

	n = size(U,1)
	n_bau = size(s0,1)


	# Evaluate basis function matrix at BAU centroids
	S_s0 = Eval_basis(s0,hcat(BasisDF[:,1],BasisDF[:,2]),BasisDF[:,3],basis_type = basis_type)
	k_s0 = S*K*S_s0'

	# Compute diagonal of triple matrix product using dot-products to conserve
	# memory - this way the largest matrix stored is (n_bau x n) instead of
	# (n_bau x n_bau) or (n x n).
	if mspe
		ABt, CDt = (S_s0*K)', (k_s0'Σ_inv)'
		At, Ct = S_s0', k_s0
		aba, cdc = zeros(n_bau), zeros(n_bau)
		for i in 1:n_bau
	 		aba[i] = ABt[:,i]'*At[:,i]
			cdc[i] = CDt[:,i]'*Ct[:,i]
		end
	end

	if (mspe && pred)
		return DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), mspe = vec(aba - cdc .+ σ_ξ), lats = s0[:,1], lons = s0[:,2])
			elseif pred return  DataFrame(preds = vec(x0 * α + k_s0'*Σ_inv*(Z-(X * α))), lons = s0[:,1], lats = s0[:,2])
			elseif mspe return DataFrame(mspe = vec(aba - cdc .+ σ_ξ), lons = s0[:,1], lats = s0[:,2])
	end
end

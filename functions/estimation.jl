# Functions for estimating unkown via Expectation-Maximization algorithm.
# This estiamtion procedure assumes homoskedastic variance for both the
# measurement error and fine-scale variation.


# Compute trace terms faster.
tr_fast(A,B) = sum(A .* B)
trQ_fast(S, Σ_η, μ_η, Z, T, α) = (tr_fast(S * Σ_η, S) + sum((S * μ_η)' * (S * μ_η))
    - 2 * sum((Z - (T * α))' * (S * μ_η)) + sum( (Z - (T * α))' * (Z - (T * α))))




# Method of EM_FRK() for case of single covariate stored as a vector.
function EM_FRK(Z, X::Vector, S, V_ξ, V_ϵ, σ_ϵ, K_0, α_0, σ_ξ_0; max_iter = 1000, tol = .0005, monitor = false)

	n, r = size(S);
	K_last = K = K_0;
	α_last = α = α_0;
	σ_ξ_last = σ_ξ = σ_ξ_0;
	invΣz = I


	i = 1; ΔK = 1; Δα = 1; Δσ_ξ = 1;
	while ( (i <= max_iter) && ( (ΔK > tol) || (Δα > tol) || (Δσ_ξ > tol) ) )
		if monitor
			println("i = ",i)
		end
		# E-step
		Dz = (σ_ξ * V_ξ) + (σ_ϵ * V_ϵ);
		invDz = Diagonal(1 ./ Vector(diag(Dz)));
		invK = inv(K);
		a = S' * invDz; b = a * S;
		Σ_η = inv( b + invK );
		μ_η = Σ_η * a * (Z - (X .* α));
		invΣz = invDz - invDz*S*Σ_η*S'invDz;

		# M-step
		K = Σ_η + (μ_η * μ_η');
		a = invDz * (Z - (S * μ_η));
		b = X' * invDz;
		α = inv( b * X) * (X' * a);

		σ_ξ = (trQ_fast(S, Σ_η, μ_η, Z, X, α) / n) - σ_ϵ;

		# Check convergence
		ΔK = norm(K .- K_last) / norm(K_last);
		Δα = norm(α .- α_last) / norm(α_last);
		Δσ_ξ = abs(σ_ξ - σ_ξ_last) / σ_ξ_last;
		if monitor && i > 1
			println("ΔK = ",  ΔK, " Δα = ",  Δα, " Δσ_ξ = ",  Δσ_ξ)
		end

		K_last = K;
		α_last = α;
		σ_ξ_last = σ_ξ;
		i = i + 1;
	end

	return (K+K')/2, α, maximum([σ_ξ,0]), invΣz
end


# Method of EM_FRK() for case of one or more covariates stored as a matrix.
function EM_FRK(Z, X::Matrix, S, V_ξ, V_ϵ, σ_ϵ, K_0, α_0, σ_ξ_0; max_iter = 1000, tol = .0005, monitor = false)

	n, r = size(S);
	K_last = K = K_0;
	α_last = α = α_0;
	σ_ξ_last = σ_ξ = σ_ξ_0;
	invΣz = I


	i = 1; ΔK = 1; Δα = 1; Δσ_ξ = 1;
	while ( (i <= max_iter) && ( (ΔK > tol) || (Δα > tol) || (Δσ_ξ > tol) ) )
		if monitor
			println("i = ",i)
		end
		# E-step
		Dz = (σ_ξ * V_ξ) + (σ_ϵ * V_ϵ);
		invDz = Diagonal(1 ./ Vector(diag(Dz)));
		invK = inv(K);
		a = S' * invDz; b = a * S;
		Σ_η = inv( b + invK );
		μ_η = Σ_η * a * (Z - (X * α));
		invΣz = invDz - invDz*S*Σ_η*S'invDz;

		# M-step
		K = Σ_η + (μ_η * μ_η');
		a = invDz * (Z - (S * μ_η));
		b = X' * invDz;
		α = inv( b * X) * (X' * a);

		σ_ξ = (trQ_fast(S, Σ_η, μ_η, Z, X, α) / n) - σ_ϵ;

		# Check convergence
		ΔK = norm(K .- K_last) / norm(K_last);
		Δα = norm(α .- α_last) / norm(α_last);
		Δσ_ξ = abs(σ_ξ - σ_ξ_last) / σ_ξ_last;
		if monitor && i > 1
			println("ΔK = ",  ΔK, " Δα = ",  Δα, " Δσ_ξ = ",  Δσ_ξ)
		end

		K_last = K;
		α_last = α;
		σ_ξ_last = σ_ξ;
		i = i + 1;
	end

	return (K+K')/2, α, maximum([σ_ξ,0]), invΣz
end

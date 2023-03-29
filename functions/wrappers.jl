# Wrapper functions for the high-level tasks of fitting and predicting.

function FRK_fit(U, Z, X; basis_type = "bisquare", basis_args = [], inits = [],
	 known_params = [], max_iter = 1000, tol = .005, monitor = false, msg = false)

	 if basis_args == [] error("Specify basis function arguments!") end
	 if msg println("\n Constructing basis function matrix...\n") end
	 S, Basis_DF = Create_basis(U,basis_args..., basis_type = basis_type)
	 n, r = size(S)

	 if known_params == []
	 known_params = [Matrix(I,n,n),Matrix(I,n,n),1]
	 end
	 # Set resonable initial values for EM if nothing specified.
	  if inits == []
		   V2 = Z'Z/n
		   inits = [.9*V2*Matrix(I,r,r), zeros(size(X,2)), .1*V2]
	  end

	  if msg println("\n Estimating unkown parameters...\n") end
	  K_hat, α_hat, σ_ξ_hat, Σ_inv_hat = EM_FRK(Z,X,S,known_params...,inits...,
											   max_iter = max_iter, tol = tol,
											   monitor = monitor)

	 return [U,Z,X,S,K_hat, α_hat, σ_ξ_hat, Σ_inv_hat, Basis_DF, basis_type]

end







function FRK_inner_predict(model, BAU_args; pred, mspe, footprint, FP_LVL,
        max_store, parallel, msg)

    if msg println("\n Generating Lvl $(BAU_args[1]) BAU grid for prediction...\n") end

    pred_BAU_df = BAU_grid(BAU_args..., parallel = parallel)
    pred_centers = Array(pred_BAU_df[:,[:lat,:lon]])
    X_pred = Array(pred_BAU_df[:,:lat])
    n_bau = nrow(pred_BAU_df)


    n_blocks = Int(round(n_bau/max_store))
    close_div = Int(floor(n_bau/n_blocks))
    div_remain = n_bau - sum(rep([close_div],n_blocks - 1))
    pred_blocks = BlockArray(pred_centers,rep([close_div, div_remain],[n_blocks-1,1]),[2])
    X_blocks = BlockArray(reshape(X_pred,n_bau,1) ,rep([close_div, div_remain],[n_blocks-1,1]),[1])

    if msg println("\n Making predictions/calculating uncertaintites...\n") end

    if footprint
            pred_hex_ids = Array(pred_BAU_df[:,:hex_id])
            fps = approx_footprint(U,BAU_args[1],FP_LVL)
            c = BAU_overlap(pred_hex_ids,fps, parallel = parallel)
            c_blocks = BlockArray(reshape(c,n_bau,1) ,rep([close_div, div_remain],[n_blocks-1,1]),[1])
            if parallel
                    if !(pred && mspe)
                            res_df = BlockArray(zeros(n_bau,3),rep([close_div, div_remain],[n_blocks-1,1]),[3])
                            Threads.@threads for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model...,c_blocks[Block(k,1)], pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:Value,:lat,:lon])
                    elseif (pred && mspe)
                            res_df = BlockArray(zeros(n_bau,4),rep([close_div, div_remain],[n_blocks-1,1]),[4])
                            Threads.@threads for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model...,c_blocks[Block(k,1)], pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:preds,:mspe,:lat,:lon])
                    end
            else
                    if !(pred && mspe)
                            res_df = BlockArray(zeros(n_bau,3),rep([close_div, div_remain],[n_blocks-1,1]),[3])
                            for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model...,c_blocks[Block(k,1)], pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:Value,:lat,:lon])
                    elseif (pred && mspe)
                            res_df = BlockArray(zeros(n_bau,4),rep([close_div, div_remain],[n_blocks-1,1]),[4])
                            for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model...,c_blocks[Block(k,1)], pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:preds,:mspe,:lat,:lon])
                    end
            end
    else
            if parallel
                    if !(pred && mspe)
                            res_df = BlockArray(zeros(n_bau,3),rep([close_div, div_remain],[n_blocks-1,1]),[3])
                            Threads.@threads for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model..., pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:Value,:lat,:lon])
                    elseif (pred && mspe)
                            res_df = BlockArray(zeros(n_bau,4),rep([close_div, div_remain],[n_blocks-1,1]),[4])
                            Threads.@threads for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model..., pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:preds,:mspe,:lat,:lon])
                    end
            else
                    if !(pred && mspe)
                            res_df = BlockArray(zeros(n_bau,3),rep([close_div, div_remain],[n_blocks-1,1]),[3])
                            for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model..., pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:Value,:lat,:lon])
                    elseif (pred && mspe)
                            res_df = BlockArray(zeros(n_bau,4),rep([close_div, div_remain],[n_blocks-1,1]),[4])
                            for k in 1:n_blocks
                                res_df[Block(k,1)] = Calculate_pred_mspe(pred_blocks[Block(k,1)],X_blocks[Block(k,1)], model..., pred = pred, mspe = mspe)
                            end
                            res_df = DataFrame(res_df,[:preds,:mspe,:lat,:lon])
                    end
            end
    end

    if msg println("\n Done! \n") end

    return innerjoin(res_df,pred_BAU_df, on = [:lat, :lon])
end


function FRK_predict(model::AbstractVector; BAU_args = [], pred = true, mspe=false,
	footprint=false, FP_LVL =  BAU_args[1]-1, max_store = 5000, parallel = false,
	show_basis = false, print = true, save = false, path = "", dataname = "", msg = false)

	if BAU_args == [] error("Specify BAU arguments!") end
	res_df = FRK_inner_predict(model, BAU_args, pred = pred, mspe = mspe,
	 footprint = footprint, FP_LVL = FP_LVL, max_store = max_store,
	 parallel = parallel, msg = msg)

	 if save
 		if show_basis
 			CSV.write("$path/$(dataname)_LVL_$(BAU_args[1]).csv", res_df)
 			CSV.write("$path/$(dataname)_basis_info.csv", Basis_DF)
 			println("Saved results at $path/$(dataname)_LVL_$(BAU_args[1]).csv")
 			println("Saved basis info at $path/$(dataname)_basis_info.csv")
 		else
 			CSV.write("$path/$(dataname)_LVL_$(BAU_args[1]).csv", res_df)
 			println("Saved results at $path/$(dataname)_LVL_$(BAU_args[1]).csv")
 		end
 	end

 	if print
 		if show_basis
 			return res_df, model[end-1]
 		else return res_df
 		end
 	end
end



function FRK(U, Z, X; basis_type = "bisquare", basis_args = [], inits = [],
	 known_params = [], max_iter = 1000, tol = .005, monitor = false,
	 pred = true, mspe=false, BAU_args = [], footprint=false, FP_LVL =  BAU_args[1]-1,
	 max_store = 5000, parallel = false, show_basis = false, print = true,
	 save = false, path = "", dataname = "", msg = true)

	 if BAU_args == [] error("Specify BAU arguments!") end

	 model = FRK_fit(U,Z,X, basis_type = basis_type, basis_args = basis_args, inits = inits,
	 	known_params = known_params, max_iter = max_iter, tol = tol, monitor = monitor,
	 	msg = msg)

	return FRK_predict(model, BAU_args = BAU_args, pred = pred, mspe = mspe,
	 footprint = footprint, FP_LVL = FP_LVL, max_store = max_store,
	 parallel = parallel, print = print, save = save,
	 path = path, dataname = dataname, msg = msg)

 end

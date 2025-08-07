""" This module implements the possibility to have an outer action acting on cells in given regions of a grid.
"""

struct OuterActionParams
    regions
    action :: Union{Function, Nothing}
end

function inside(low, hi, point)
    n = length(low)

    for i in 1:n
        if point[i] < low[i] || point[i] > hi[i]
            return false
        end
    end

    return true
end

function outer_action_post_mcs_listeners(model :: GridModel{T, N, CPMState{Tdc, Trnd, N}, CPMCfg}) where {T, Tdc, N, Trnd}    
    params = model.cfg.params[:outer_action]
    action = params.action

    if action === nothing
        return
    end

    centroids::Dict{T,Vector{Float64}} = get_stat(model, :centroid)

    for r in params.regions
        # application of action to cells in region
        for (cid, centroid) in centroids
            if inside(r[1], r[2], centroid)
                action(model, r, cid, centroid)
            end
        end

        # application of action to region
        action(model, r, nothing, nothing)
    end
end

const OuterActionCfg = mk_cfg(
    post_mcs_listeners = [outer_action_post_mcs_listeners]
)

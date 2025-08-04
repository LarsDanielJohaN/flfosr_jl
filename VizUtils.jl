
using Statistics,  Base.Threads  , Plots


#=
postf, compute bounds for (1-alpha)% credible bands per points in design grid and mean function. 

F: Matrix{Float64} of dimensions Tn x Nf where Tn is the number of points in the design grid and Nf is the number of functions. 
alpha: Float64 or Nothing, significance level, defaults to 0.05. 
=#

function postf(F, alpha::Union{Float64, Nothing} = 0.05)
    @assert (alpha >= 0.0) & (alpha <= 1.0) "Fatal error on postf, alpha should satisfy 0.0 <= alpha <= 1.0" 

    Tn, Nf = size(F) #Obtains dimensions of F
    q_alpha_half= zeros(Tn) #Declares vector of dimension Tn to store empiricial quantiles at alpha/2 
    q_one_m_alpha_half = zeros(Tn) #Declares vector of dimension Tn to store empirical quantiles at 1 - alpha/2
    mu = mean(F, dims = 2) #Computes point wise means. 

    Threads.@threads for tn in 1:Tn #Computes quantiles per design point in parallel. 
        q_alpha_half[tn] = quantile( F[tn, :], alpha/2)
        q_one_m_alpha_half[tn] = quantile( F[tn, :], 1- alpha/2)
    end

    return Dict("mu"=>mu, "q_alpha_half"=>q_alpha_half, "q_one_m_alpha_half"=>q_one_m_alpha_half)

end



function plotf_mean_q(F::Array{Float64}, T_grid::Vector{Float64} = range(0, 1, length = size(F)[1]) |> collect,  alpha::Float64 = 0.05, alpha_viz::Float64 = 0.5, add_area::Bool = true , title = "")
    @assert length(T_grid) == size(F)[1] "Fatal error on plotf! T_grid should be of size Tn. "

    dict_postf = postf(F, alpha) #Obtains pointwise quantile values and mean function. 
    pq = plot(T_grid, dict_postf["q_one_m_alpha_half"], fillrange = dict_postf["q_alpha_half"], fillalpha = alpha_viz, color = :red, label = false, linealpha = 0, title = title) #Adds the filling between confidence bands.
    pq =plot!(pq, T_grid, [dict_postf["mu"], dict_postf["q_alpha_half"], dict_postf["q_one_m_alpha_half"]], label = ["mu" "q: $(alpha/2)" "q: $(1-alpha/2)"], color = [:red :blue :blue]) #Adds confidence bands. 
    return pq 


end
#=
plotf returns a plot object with a given collection of functions in it. 
F: Matrix{Float64} of dimensions Tn x Nf where Tn is the number of design points and Nf the total number of functions to plot. 
T_grid: Union{Vector{Float64}, Nothing}, vector of length Tn with design points. Defaults to an equally spaced grid of points in the unit interval of size Tn. 
add_postf: Union{Bool, Nothing}, boolean value which tells whether to add confidence bands and mean value, defaults to true. If add_postf = true, the functions in F 
will be ploted with a lower intensity than the mean and confidence bands. 
alpha: Union{Float64, Nothing} significance level for the confidence bands to use in case add_postf = true, defaults to 0.05. 
alpha_viz: Union{Float64, Nothing}, intensity of curve values. e.g. the alpha value specified to plot. 
=#

function plotf(F::Matrix{Float64}, T_grid::Union{Vector{Float64}, Nothing} = range(0, 1, length = size(F)[1]) |> collect, add_postf::Union{Bool, Nothing} = true, alpha::Union{Float64, Nothing} = 0.05, alpha_viz::Union{Float64, Nothing} = 0.2)
    @assert length(T_grid) == size(F)[1] "Fatal error on plotf! T_grid should be of size Tn. "
    Tn, Nf = size(F)
    p = plot(T_grid, F, label = false, alpha = alpha_viz) #Plots functions on p
    if add_postf #If add_postf = true, then obtain confidence bands, mean function and plot them. 
        dict_postf = postf(F, alpha)
        p =plot!(p, T_grid, [dict_postf["mu"], dict_postf["q_alpha_half"], dict_postf["q_one_m_alpha_half"]], label = ["mu" "q: $(alpha/2)" "q: $(1-alpha/2)"], color = [:red :blue :blue])
    end

    return p

end
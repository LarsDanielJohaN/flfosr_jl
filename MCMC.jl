using LinearAlgebra, Distributions ,  Base.Threads, StatsBase, RCall

#=
sample_MVN_canonical. Obtains sample from MVN with mean inv(Q)b, and precision matrix Q. 
This is performed using the procedure described on Algorithm 2.5 from Rue H. and  Held L.'s 
Gaussian Markov Random Fields book. 
=#
function sample_MVN_canonical(;Q::Matrix{Float64}, b::Vector{Float64} ) #, L::Matrix{Float64})
    L = cholesky( (Q + Q')./2 ).L
    return (L' \ (L\b)) + (L'\randn(size(Q)[1]))
end

function get_Z_mat(M_rep)
    M = sum(M_rep)
    N = length(M_rep)
    Z = zeros(M, N) #Creates M x N zeros matrix where M = M1 + M2 + ... + MN
    start_idx = 1
    for n in 1:N
        len = M_rep[n]
        Z[ start_idx:(start_idx + len -1), n] .= 1 #Sets the corresponding block diagonal. 
        start_idx += len
    end
    return Z

end


#=
Y: Matrix{Float64} of observations  of dimension T x M with M = M1 + M2 + ... + MN
X: Matrix{Float64} design matrix of dimension M x L+1 with M = M1 + M2 + ... + MN
M_rep: Vector{Int64} of size N whose ith entry corresponds to Mi, e.g. M_rep[i] = Mi
K: Int64, number of basis functions. 
S: Int64, number of MCMC iterations. 
center_base: Bool, boolean value that indicates whether reparametrized base should consider a centered design matrix (for functions). 
tol_sm: Float64, Zero tolernce for sm. e.g. in factorization all values lower than tol_sm are considered to be zero. 

Notes: 
-It is assumed that the order in Y's columns goes in accordance to the entries in M_rep. This in the sense
that, knowing that M_rep[1] = M1, then Y[:, 1:M1] are the observations of the 1st subject, Y[:, M1+1: M1+M2] correspond to 
the obsevations of subject 2, etc. 

-It is assumed that X's rows are ordered in accordance to the columns of Y. That is, if 
Y[:,l] corresponds to the j-th visit of the i-th subject, then X[l,:] has its corresponding covariates. 

-For the basis re-paremetrization, we replicate the sm(.) function from R's spikeSlabGAM library. Specifically, we 
replicate the function with the "ortho" decomposition option, spline.degree = 3, diff.ord = 2, centerx = x,
and rankZ = 0.999. 
=#

function flfosr(; Y::Matrix{Float64}, X::Matrix{Float64}, M_rep::Vector{Int64}, K::Int64 = 10, S::Int64=2000 , S_burn::Int64 = 1000 , a_alph::Float64 = 0.1, b_alph::Float64 = 0.1, a_gamm::Float64 = 0.1, b_gamm::Float64 = 0.1, a_omeg::Float64 = 0.1,  b_omeg::Float64 = 0.1, tol_sm::Float64 = (1/10^10))
    #Makes sure that appropiate inputs are recieved. 
    T,M_Y = size(Y)
    M_X, L_p_one = size(X)
    M_M_rep = sum(M_rep)
    @assert M_Y == M_X "Fatal error on flfosr!! Y should be of dimensions (T x M) and X of (M x L+1) but the M's dont coincide. "
    @assert M_M_rep == M_X "Fatal error on flfosr!! X should be of dimensions (M x L+1) and the entries of M_rep (M1 + M2 + ... + MN) should sum to M. They dont. "

    Tn, M_Y = size(Y) #Obtains number of evaluation points with Tn and total number of functions M_Y
    N = length(M_rep) #Gets total number of subjects.
    tau = range(0, 1, length = Tn) |> collect  #Obtains the observation points.  
    idx = inverse_rle( 1:N |> collect, M_rep )

    
    R"library(spikeSlabGAM) "

    R"""
    B <- cbind(1/sqrt($Tn), poly($tau, 1), sm($tau, K = $K, rankZ = .99999999,  spline.degree = 3, diff.ord = 2, centerBase = T))
    B <- B/sqrt(sum(diag(crossprod(B))))
    Dk <- diag(crossprod(B))
    Bk <- B%*%diag(1/sqrt(Dk))
        """
    B = rcopy(Matrix{Float64}, R"B"   )
    B_proj = rcopy(Matrix{Float64}, R"Bk")'
    B = B_proj'

    Y_proj = B_proj*Y #Obtains projected data, now an K x M matrix. 
    Z = get_Z_mat(M_rep)

    Alpha = zeros(K, L+1, S +1 ) #Creates K x L+1 x S tensor, each slice K x L+1 x s for s = 1,2,..., S represents the iter's sample of fixed effect coeffients. 
    Gamma = zeros(K, N, S+1) #Creates K x N x S tensor, each slice K x N x s represents s'ths values for the subject specific effects. 
    Omega = zeros(K, M_Y, S+1) #Creates a K x M x S tensor, each slice K x M_Y x s represents the s´ths vales for the specific specific coefficients. 

    Gamma[:, :, 1] =  ((rowsum( Matrix( Y'),   idx,    1:N |>collect    ) ./M_rep )*B)'
    Omega[:, :, 1] =  B'*(Y - B*Gamma[:, idx, 1])
    
    Sig_Eps = zeros(S+1) #Creates a length S vector to store the realizations for the noise variance. 
    Sig_Eps[1] = 0.1
    Sig_Alpha = zeros(L+1, S+1) #Creates an L x S matrix to store the realizations of the fixed effect coefficient variances. 
    Sig_Alpha[:, 1] =2.0 * ones(L+1)
    Sig_Gamma = zeros(S+1) #Creates a length S vector to store the realization for the subject efffect coefficient variances. 
    Sig_Gamma[1] = 1.0
    Sig_Omega = zeros(N, S+1) #Creates a N x S matrix to store the realizations of the visit effect coefficient variances. 
    Sig_Omega[:, 1] = 0.5 .* ones(N)

    ##---------------------------------------------------------------------------------------------------------------------------------
    #Reserve memory space or calculate terms which will be used repeatedly during the MCMC scheme (i.e. to reduce allocation time and calculation time per iteration)

    ell_alpha_k = zeros(L+1, M_Y) #Creates fixed memory for storing the terms dkX'(Vk - Wk).
    Q_inv_gamma_k = zeros(N) #Creates fixed memory for storing the terms dk Z'Vk(yk - Xalphak).
    Q_inv_omega_k = zeros(M_Y) #Creates fixed memory for storing terms diag({ dk/sigma_eps^2 + 1/sigma_omegai^2  })
    su = 1
    ##---------------------------------------------------------------------------------------------------------------------------------

    Alphaf = zeros(Tn, L+1, S - S_burn) #Creates a Tn x L+1 x S- S_burn tensor to store Alpha realizations in their functional form. 
    @views begin
        for s in 2:(S+1) #Do MCMC iterations. 
            if mod(s, 200) == 0 
                println("On iteration $(s) from $(S)")
            end

            Threads.@threads for k in 1:K #Sample coefficients for each individual function in parallel. 

                ##---------------------------------------------------------------------------------------------------------------------------------
                ell_alpha_k = X' *(  Diagonal(  inverse_rle(1 ./ (Sig_Eps[s-1] .+ Sig_Omega[:, s-1] .+ M_rep .*Sig_Gamma[s-1] ), M_rep) )) #Good to go
                Alpha[k, :, s] =  sample_MVN_canonical( Q = Diagonal(1 ./ Sig_Alpha[:, s-1]) + ell_alpha_k*X ,  b = ell_alpha_k*Y_proj[k, :] ) #Good to go. 

                ##---------------------------------------------------------------------------------------------------------------------------------
                Q_inv_gamma_k = 1 ./ ((1/Sig_Gamma[s-1]) .+ M_rep .* (   1 ./(    Sig_Eps[s-1] .+ Sig_Omega[:, s-1]  ) ) ) #Good to go. 
                #Vk = inverse_rle(1 ./ (sig_eps[s-1] .+ d[k].*Sig_Omega[:, s-1]), M_rep) 
                Gamma[k, :, s] = rand(MvNormal(  Diagonal(Q_inv_gamma_k)* (Z')*(  Diagonal(inverse_rle(1 ./ (Sig_Eps[s-1] .+ Sig_Omega[:, s-1]), M_rep) )*(Y_proj[k, :] - X*Alpha[k, :, s]) ), Diagonal(Q_inv_gamma_k)                          )         ,1)

                ##---------------------------------------------------------------------------------------------------------------------------------
                Q_inv_omega_k = inverse_rle(   1 ./ ( (1/Sig_Eps[s-1]) .+ (1 ./ Sig_Omega[:, s-1]) ), M_rep)
                #ell_omegak = (d[k]/Sig_Eps[s-1]).*(Y_proj[k, :] - X* Alpha[k, :, s] - Z*Gamma[k, :, s] )
                Omega[k , :, s] = rand(MvNormal(Diagonal(Q_inv_omega_k) * (  (1/Sig_Eps[s-1]).*(Y_proj[k, :] - X* Alpha[k, :, s] - Gamma[k, idx, s] )   ),  Diagonal(Q_inv_omega_k)    ), 1)
            end 

            #=
            Note that Y - B*( Alpha[:, :, s]*X' + Omega[:, :, s] + Gamma[:,idx, s] constitutes the residuals of estimiating Y with the current values for the coefficients. 
            therefore, taking the norm of the previous consitutues the sum squared residuals. 
            =#
            Sig_Eps[s] = 1/rand(Distributions.Gamma(Tn*M_Y/2, 1/( sum(    ((Y - B*( Alpha[:, :, s]*X' + Omega[:, :, s] + Gamma[:,idx, s])).^2 ) ./ 2 )   )), 1)[1] #Obtain new realization for the observation error variance. 

            ##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            #lam_sig_alpha = 1 ./ vec( b_alph .+ sum(Alpha[:, 1:(L+1), s].^2, dims=1)/2  )
            Threads.@threads for l in 1:(L+1) #Obtains new samples for fixed effect coefficient variances. 
                Sig_Alpha[l, s] =  1.0/rand(Distributions.Gamma(  a_alph+ K/2,     1/(b_alph + sum( (Alpha[:, l,s].^2) ./ 2))  ) ,1)[1] 
            end
            ##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

            #=
            Note that sum(Gamma[:, :, s] .^2 ) takes all coefficients for subject effect functions, squares them and smbs them. 
            =#
            Sig_Gamma[s] =  1/rand(Distributions.Gamma(a_gamm + N*K/2,     1/(b_gamm + sum(  (Gamma[:, :, s].^2) ./ 2 ))    ),     1)[1]   #Obtain new relizations for the subject effect coefficient variance/smoothing parameter.

            ##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            su = 1 
            Sig_Omega[1, s] = 1.0/ rand(Distributions.Gamma(a_omeg + M_rep[1]*K/2,  1/(    b_omeg + sum( (Omega[ : ,1:M_rep[1] , s].^2) ./2    )       )   ), 1)[1]
            for n in 2:N 
                Sig_Omega[n, s] = 1.0/ rand(Distributions.Gamma(a_omeg + M_rep[n]*K/2,  1/(    b_omeg + sum(  (Omega[ : ,su + 1:su + M_rep[n] , s].^2) ./ 2    )      )   ), 1)[1]
                su+= M_rep[n]
            end 
            ##----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        end 
    end
    for s in 1:(S - S_burn)  
        Alphaf[: ,:, s] = B*Alpha[: ,:, s+S_burn-1]        
    end

    return Dict("X"=>X, "B"=>B_proj',  "w_post"=>Omega[:, :, S_burn:end], "ga_post"=>Gamma[:, :, S_burn:end], "alpha_post"=>Alpha[:, :, S_burn:end], "alpha_postf"=>Alphaf, "sig_eps_post"=>Sig_Eps[S_burn:end], "sig_alpha_post"=>Sig_Alpha[:, S_burn:end], "sig_gamma_post"=>Sig_Gamma[S_burn:end], "sig_omega_post"=>Sig_Omega[:, S_burn:end] )

end 

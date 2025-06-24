
using LinearAlgebra, Statistics, StatsBase 

#=
m: value such that m+1 is the degree of B^m(x)
i: element in the basis
T: vector of knots
x: value to evaluate. 
=#
function eval_B_Spline(x::Float64, T::Vector{Float64}, i::Int64, m::Int64)
    #T = range(0, stop = 10, length = 20) |> collect
    #@assert T isa Vector{Float64} "Fatal error on eval_B_Spline!! T should be of type Vector{Float64}"
    #@assert x isa Float64 "Fatal error on eval_B_Spline! x should be of type Float64"
    #@assert (i isa Int64) & (m isa Int64) "Fatal error on eval_B_Spline!! i and m should be of type Int64"
    if m == -1 
        return Float64(  (x < T[i+1]) & (x >= T[i])  )
    else 
        z0 = (x - T[i])/(T[i+m+1] - T[i])
        z1 = (T[i+m+2] - x)/(T[i+m+2] - T[i+1])
        return z0*eval_B_Spline(x,T, i, m-1) + z1*eval_B_Spline(x, T, i+1, m-1)
    end  
end

#=
get_des mat: given a collection of evaluation points x_e, obtains the design matrix of a B-Spline basis of order m+1, k basis functions
the B-Spline basis defined on knots on T. 
T: design knots. 
x_e: points where to evaluate functions.
k: number of elements in basis.
m: number such that m+1 is the degree of the resulting polinomial.

=#

function get_des_mat_B_Spline(x_e::Vector{Float64}, K::Int64, m::Int64, T::Vector{Float64})

    #@assert x_e isa Vector{Float64} "Fatal error on get_des_mat_B_Spline!! x_e should be of type Vector{Float64}"
    #@assert (k isa Int64) & (m isa Int64) "Fatal error on get_des_mat_B_Spline!! Both k and m should be of type Int64"
    
    B1(x) = eval_B_Spline(x, T, 1, m) #Creates function that evaluates function 1 on x
    X = B1.(x_e) #Evaluates x_e with function B1. 
    for i in 2:K 
        Bi(x) = eval_B_Spline(x,T, i, m)#Creates function that evaluates function i on x
        X = hcat( X, Bi.(x_e))  #Evaluates x_e with function Bi. 
    end 
    return X

end


#=
x: vector of knots to evaluate basis
K: number of elements on the basis
spline_deg: degree of the resulting spline
=#
function psBasis(x::Vector{Float64}, K::Int64 = length(x), spline_degree::Int64 = 3, diff_ord::Int64 = 2 )
    @assert K > spline_degree "Fatal error on psBasis!! K should be different from spline_degree"
    knots_no = K - spline_degree + 1 
    xl = min(x...)
    xr = max(x...)
    xmin = xl - (xr - xl)/100 
    xmax = xr + (xr - xl)/100
    dx = (xmax - xmin)/(knots_no - 1)
    T = sort( range(xmin - spline_degree * dx, stop = xmax + spline_degree * dx, step = dx) |> collect ) #Obtains vector of knots. 
    X= get_des_mat_B_Spline(x, K, spline_degree-1, T)
    return X
end
#=
get_scnd_pen_mat creates penalty matrix on second squared diferences on the coeficeints. 
=#

function get_pen_mat(k::Int64)
    D = zeros(k-1, k)

    for i in 1:(k-1)
        D[i,i] = -1 
        D[i, i+1] = 1
    end
    return D'*D
end

#=
sm 
Warning: sm, in our current implementation, only performs part of the functionalities in spikeSlabGAM's sm.
=#

function sm(; x::Vector{Float64}, K::Union{Int64, Nothing} = min(length(x), 20), spline_deg::Union{Int64, Nothing} = 3, diff_ord::Union{Int64, Nothing} = 2,  center_base::Union{Bool, Nothing} = true,  tol::Float64 = (1/10^10))
    B = psBasis(x, K, spline_deg, diff_ord) 
    P = get_pen_mat(K) #Checar si el diff_ord no causa problemas. 
    if center_base
        B = centerBase(B, x, diff_ord - 1 )
    end
    B = orthoDesign(P = P, X = B, tol = tol)
    B = scaleMat(B)
    return B 

end


#=
orthoDesign gets reparametrized centered basis with procedure described by Scheipl F., Fahrmeir L. & Kneib T. (2012)
D: Penalization matrix. 
X: Design matrix. 
tol: To be considered zero tolserance. 
=#
function orthoDesign(; P::Matrix{Float64}, X::Matrix{Float64}, tol::Float64 = (1/10^10))
    @assert (tol >= 0.0) "Fatal error on sm!!, tol should satisfy  tol >= 0.0"
    Pi = pinv(P) #Gets generalized inverse (Moore-Penrose) of penalization matrix. 
    M = X*Pi*X' #Computes covariance matrix of XB, with B the coeficients. 
    U, V, Ut = svd(M) #Computes SVD of X*Pi*X', V i san array. 
    keep_me_indx =  V  .> tol #Gets which values are greater than the tolerance. 

    Vp = diagm(V[keep_me_indx]) #Gets singular values that are greater than zero (in accordance to the tolerance)
    Up = U[:, keep_me_indx] #Gets vectors corresponding to the 

    Xm = Up*sqrt.(Vp) #Computes new basis. 
    return Xm

end


#=
center_base (intended to replicate the one from R's spikeSlabGAM)

center a basis matrix B s.t. its column space does not
include polynomials of a vector x up to degree 'degree', or, if x is a matrix,
s.t. B's new column space is orthogonal to x's column space
=#

function centerBase(B::Matrix{Float64}, x::Union{Matrix{Float64}, Vector{Float64}}, degree::Int64)
    if x isa Vector{Float64}
        Xx = hcat(ones(length(x), 1 ), poly(x, degree))
    else 
        Xx = x
    end
    qrX = qr(Xx)
    Bc = qr_resid(qrX, Xx, B)

    return Bc     
end

#=
poly, intended to replicate functionalities of R's poly. It only replicates the use of poly for the case where coefficients 
arent provided and raw evaluation inst specified. Additionally, poly only takes one dimensional vectors as input. 
x: Vector{Float64}, 
=#

function poly( x::Vector{Float64}, degree::Int64 )
    ##@assert x isa Vector{Float64} "Fatal error on poly!! x must be a Vector{Float64}"
    @assert (degree isa Int64) & (degree < length(unique(x)) ) "Fatal error on poly! Either degree isnt an Int64, or degree >= length(unique(x)). "
    @assert sum(  isnan.(x)  ) == 0 "Fatal error on poly!! x should not contain nan's."
    x = x .- mean(x)
    X = reduce(hcat, (x.^k for k in 0:degree))
    QR = qr(X)
    return QR.Q[:, 2:(degree + 1)]
end

function qr_resid(F,X, Y)
    b = F \ Y #Computes OLS solutions columnwise. 
    r = Y - X*b #Computes residuals. 
    return r 
end

function scaleMat(X::Matrix{Float64}, factor::Int64 = 2)
    r,c = size(X)
    return X ./ (factor*  sqrt(tr(X*X')/r)    )
end


function rowsum(A::Matrix{Float64}, idx::Vector{Int64}, id::Vector{Int64})
    r,c = size(A)
    T = Matrix{Float64}(undef, length(id), c)
    for i in id
        T[i, :] = sum( A[idx .== i, :],  dims = 1  )
    end
    return T 
end
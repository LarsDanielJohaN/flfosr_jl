using RCall, StatsBase, Plots
include("VizUtils.jl")


R"""
    #----------------------------------------------------------------------------
    # Helper functions
    #----------------------------------------------------------------------------


    psBasis <- function(x, K = length(x), spline.degree = 3, diff.ord = 2,
                        knots = NULL) {
    if (is.null(knots)) {
        knots.no <- K - spline.degree + 1
        xl <- min(x)
        xr <- max(x)
        xmin <- xl - (xr - xl) / 100
        xmax <- xr + (xr - xl) / 100
        dx <- (xmax - xmin) / (knots.no - 1)
        knots <- seq(xmin - spline.degree * dx, xmax + spline.degree * dx, by = dx)
    }

    X <- splines::spline.des(knots, x, spline.degree + 1, outer.ok = TRUE)$design
    P <- diag(K) # precision
    if (diff.ord > 0) {
        for (d in 1:diff.ord) P <- diff(P)
        P <- crossprod(P)
    }
    return(list(
        X = X, P = P, knots = knots, K = K, spline.degree = spline.degree,
        diff.ord = diff.ord
    ))
    }


    makeZ <- function(mi){
    z <- matrix(0, nrow = sum(mi), ncol = length(mi))
    k <- 1
    for(i in 1:length(mi)){
        z[k:(k+mi[i]-1),i] <- 1
        k <- k+mi[i]
    }
    z
    }


    rowrep <- function(X, ntimes){
    #as.matrix(as.data.frame(lapply(as.data.frame(X), rep, ntimes)))
    
    X[rep(seq_along(ntimes), ntimes), ]
    }


    makeYall <- function(){
    
    X2 <- as.data.frame(rowrep(X, Mi))
    X2$subject <- rep(1:N, Mi)
    X2$subject <- factor(X2$subject)
    Yall <- X2
    Yall$Y <- t(unname(Y))
    class(Yall$Y) = class(Yall$Y)[-1]
    Yall$curve <- 1:nrow(Yall)
    Yall$curve <- factor(Yall$curve)
    Yall
    }

    fosrcoef <- function(alphaf_post){
    
    S <- ncol(alphaf_post)
        
    alphaf_pci <- apply(alphaf_post[,], 1, function(x) quantile(x, c(.025,.975)))
    alphaf_mean <- rowMeans(alphaf_post[,])

    output <- list(alphaf_pci = alphaf_pci,
                    alphaf_mean = alphaf_mean)
    output
    }

    fmse <- function(alpha_true, alpha_hat){
    
        sum((alpha_true - alpha_hat)^2)/length(alpha_true)
    
    }

    mciw <- function(u, l){
    mean(u - l)
    
    }

    ecp <- function(u, l, alpha){
    
    mean(alpha < u & alpha > l)
    
    
    }

    fosryhat <- function(m1){
    
    yhat_draws <- array(NA, dim = c(nrow(m1$B), nrow(m1$w_post[[1]]), length(m1$w_post[(S/2 + 1):S])))
    
    for(i in (S/2 + 1):S){
        yhat_draws[,,i - (S/2)] <- m1$B%*%t(rowrep(m1$X,Mi)%*%(m1$alpha_post[[i]][,]) + rowrep(m1$ga_post[[i]], Mi) + m1$w_post[[i]])
    }
    output <- rowMeans( yhat_draws , dims = 2 )
    output
    }

    #' Compute Simultaneous Credible Bands
    #'
    #' Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
    #'
    #' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
    #' @param alpha confidence level
    #'
    #' @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
    #'
    #' @note The input needs not be curves: the simultaneous credible "bands" may be computed
    #' for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
    #' level across all components of the vector.
    #'
    #' @export
    credBands = function(sampFuns, alpha = .05){
    
    N = nrow(sampFuns); m = ncol(sampFuns)
    
    # Compute pointwise mean and SD of f(x):
    Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)
    
    # Compute standardized absolute deviation:
    Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)
    
    # And the maximum:
    Maxfx = apply(Standfx, 1, max)
    
    # Compute the (1-alpha) sample quantile:
    Malpha = quantile(Maxfx, 1-alpha)
    
    # Finally, store the bands in a (m x 2) matrix of (lower, upper)
    t(cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx))
    }
 """



R"""
    library(fda)
    library(mgcv)
    library(spikeSlabGAM)

    #----------------------------------------------------------------------------
    # Fast longitudinal function-on-scalar regression (FLFOSR)
    #----------------------------------------------------------------------------

    ### Y: Tn x M Matrix of functional observations,
    # Each column corresponds to one curve from a subject.
    # Each row is a measurement at a single timepoint on a common grid.
    ### X: M x L+1 or N x L+1 design matrix of fixed effects
    # Should include column of 1s for the intercept
    ### z: Length M vector specifying group memberships of each curve
    # e.g. if there are 20 subjects with 5 repeated measurements each, z=rep(1:20, each = 5)
    ### k: number of basis functions
    ### S: number of total MCMC iterations
    ### S_burn: burn first # of MCMC iterations for warm-up
    ### a_a, b_a, a_g, b_g, a_w, b_w: gamma hyperparameters for variance priors of alpha, gamma and omega

flfosr1 <- function(Y, X, z, k = 10, S = 2000, S_burn = S/2,
                       a_a = .1, b_a = .1, a_g = .1, b_g = .1, a_w = .1, b_w = .1){
  

  Mi <- table(z)
  MM <- sum(Mi)
  N <- length(Mi)
  Tn <- nrow(Y)
  
  
  # tau is the obs points (length Tn)
  tau <- seq(0, 1, by = 1/(Tn-1))
  Bmat = psBasis(1:Tn, K = k)$X
  Kbmat <- ncol(Bmat)
  P <- psBasis(1:Tn, K = k)$P

  B = cbind(1/sqrt(Tn), poly(tau, 1), sm(tau, K = k, rankZ = .99999999,  spline.degree = 3, diff.ord = 2, centerBase = T))
  B = B/sqrt(sum(diag(crossprod(B))))
  K <- ncol(B)
  L <- ifelse(is.matrix(X) == T, ncol(X) - 1, 0)
  Dk <- diag(crossprod(B))
  D <-  rep(Dk, rep.int(MM, ncol(B)))
  Bk <- B%*%diag(1/sqrt(Dk))
  Yk <- crossprod(Y, Bk)
  group1 <- z
  
  Nx <- nrow(X)

  if(Nx == N){ #If rows of X is the number of subjects (e.g. potentially rows(X) != M1 + M2 + ... + MN) update design matrix. 
    ZX <- rowrep(X, Mi)
  }else{ #In other case, keep design matrix as it is. 
    ZX <- X
  }

  ## The following section performs the MCMC
  ##....................................................................................................................................................
  
  
  ##----------------------------------------------------------------------------------------------------------------------------------------------------
  ##Sets the hyperparameters for the variance components. 
  a_alph <- a_a
  b_alph <- b_a
  a_omega <- a_w
  b_omega <- b_w
  a_ga <- a_g
  b_ga <- b_g
  
  e <- 1
  alpha <- matrix(0, nrow = L+1 , K) #dim(alpha) =  (L+1, K), creates matrix to store alpha values. 
  ga <- (rowsum(t(Y - Bk%*%t(ZX%*%alpha)), group1)/c(Mi))%*%Bk #dim(ga) = (N, K), creates matrix to store gamma coefficent values. 
  w <- t(Y - Bk%*%t(ZX%*%alpha + rowrep(ga,Mi)))%*%Bk #dim(w) = (M,K) creates matrix to store omega coefficient values. 
  
  sig_e <- .1
  sig_alpha <- c(matrix(.0001, nrow = L+1)) #
  sig_ga <- rowSums(ga^2)/K #length(sig_ga) = N ? Sets initial values for sig_ga
  sig_w <- rowSums(w^2)/K #length(sig_w) = M ? Sets initial values for sig_w
  
  S <- S
  
  w_post <- list() #List to store posterior draws from omega coefficients. 
  ga_post <- list() #List to store posterior draws from gamma coefficients. 
  alpha_post <- list() #List to store posterior draws from alpha coefficients. 
  sig_e_post <- rep(NA, S) #Create vector to store posterior 
  sig_alpha_post <- matrix(NA, nrow = S, ncol = L+1)
  sig_ga_post <- matrix(NA, nrow = S, ncol = N)
  sig_w_post <- matrix(NA, nrow = S, ncol = MM)
  
  nog <- which(Mi == 1)
  progress <- floor(seq(1, S, length.out= 11))
  


  for(s in 1:S){ #For loop that performs MCMC iterations. 
    
    ##----------------------------------------------------------------------------------------------------------------------------------------------------
    ##Sample alpha coefficients. 

      eG <- matrix(rep(1/(sig_e + sig_w), K), ncol = K)
      sumeG <- rowsum(eG, group1)
      eGh <- eG - (rowrep(sumeG/((1/sig_ga)+sumeG), Mi) )*eG
      
      if(Nx == N){
      Q_alpha <- matrix(rep(c(diag(1/sig_alpha) + crossprod(X*((sumeG - (sumeG^2)/(1/sig_ga+sumeG))[,1]), X)), K), ncol = K)
      }else{
      Q_alpha <- matrix(rep(c(diag(1/sig_alpha) + crossprod(ZX, ZX*eGh[,1])), K), ncol = K)
      }

      l_alpha <- crossprod((ZX), Yk*eGh)
      
      alpha <- sapply(1:K, function(x) { #Open function. 

        if(L <= Nx){
        
          ch_Q <- chol(matrix(Q_alpha[,x], nrow=L+1, ncol=L+1))
          alpha <- backsolve(ch_Q,
                             forwardsolve(t(ch_Q), l_alpha[,x]) +
                               rnorm((L+1)))

          alpha
          
        }else{
          H <- Map('*', 1/((1/sig_ga)+sumeG[,x]), tapply(eG[,x],  group1, function(v)
            if(length(v) > 1) {outer(v, v)}
            else{v}))
          
          u <- rnorm(L+1, 0, sqrt(sig_alpha))
          delta <- rnorm(N)

          Phi <- X*(sqrt(sumeG[,x] -  sapply(H, sum)))
          v = Phi%*%u + delta

          pw <- solve(tcrossprod(Phi*rep(sqrt(sig_alpha), each = N)) + diag(N),
                      (rowsum(Yk[,x],group1)/c(Mi))*(sqrt(sumeG[,x] -  sapply(H, sum))) - v)
          alpha <- u + sig_alpha*t(Phi)%*%pw
          

          alpha
        }
      })#Close function and sapply
      
      
      ##----------------------------------------------------------------------------------------------------------------------------------------------------
      ##Sample gamma coefficients. 
      

      Q_gaij <- rep(1/sig_ga, times = K) + rowsum(eG, group1) #Isnt the dk missing here?
      l_gaij <- rowsum(eG*(Yk - ZX%*%alpha), group1)
      ga <- matrix(rnorm(N*K, (1/Q_gaij)*l_gaij, sqrt(1/Q_gaij)), nrow = N , K)
      ga[nog,] <- 0
      
      
      ##----------------------------------------------------------------------------------------------------------------------------------------------------
      ##Sample omega coefficients. 
      
      
      
      

      Q_wij <- 1/sig_w + 1/sig_e

        if(Nx == N){
        w <- matrix(rnorm(MM*K, (1/(sig_e))*(Yk -   rowrep(X%*%alpha + ga, Mi))*c(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
        }else{
        w <- matrix(rnorm(MM*K, (1/(sig_e))*(Yk -   ZX%*%alpha - rowrep(ga, Mi))*c(1/Q_wij), sd = sqrt(1/Q_wij)  ), nrow = MM, ncol = K)
        }

      ##----------------------------------------------------------------------------------------------------------------------------------------------------
      ##Sample variance components. 
      
      
      sig_w <- rep(1/rgamma(N, a_omega + Mi*K/2, rate = b_omega + rowSums(rowsum(w^2, group1))/2), Mi)
      sig_ga <- rep(1/rgamma(1, a_ga + N*K/2, rate = b_ga + sum(((ga))^2)/2), N)
      sig_alpha <- 1/rgamma((L+1), a_alph + K/2, rate = b_alph + rowSums((alpha)^2)/2)
      
      
      if(MM > 50000){
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(B, (w + rowrep((X%*%alpha + ga), Mi)) ))^2)/2)
      }else{
        if(Nx == N){
        sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(Bk, (w +  rowrep(X%*%alpha + ga, Mi)) ))^2)/2)
        }else{
          sig_e <- 1/rgamma(1, Tn*MM/2, rate = sum((Y - tcrossprod(Bk, (w +  ZX%*%alpha + rowrep(ga, Mi)) ))^2)/2)
        }
      }

      w_post[[s]] <- w
      ga_post[[s]] <- ga
      alpha_post[[s]] <- alpha
      sig_e_post[s] <- sig_e
      sig_alpha_post[s,] <- sig_alpha
      sig_ga_post[s,] <- sig_ga
      sig_w_post[s,] <- sig_w

      if(s %in% progress){
        print(paste0("MCMC draws: [", s, "/", S, "]"))
      }
   
  } #Closes loop for iterations. 
  ##Section that performs MCMC ends. 
  ##....................................................................................................................................................
  
  
  
  
  ##----------------------------------------------------------------------------------------------------------------------------------------------------
  #store MCMC draws of fixed effects functions
  alphaf_post <- list()
  for(i in 1:(L+1)){
    alphaf_post[[i]] <- (Bk)%*%t(do.call(rbind, lapply(alpha_post, function(x) x[i,]))[(S_burn + 1):S,])
  }
  
  m1 <- list(X = X,
             B = Bk,
             w_post = w_post,
             ga_post = ga_post,
             alpha_post = alpha_post,
             sig_e_post = sig_e_post,
             sig_alpha_post = sig_alpha_post,
             sig_ga_post = sig_ga_post,
             sig_w_post = sig_w_post,
             alphaf_post = alphaf_post)
  
  return(m1)
  
}
"""



function test_psBasis(; tau::Vector{Float64} = range(0, 1, length = 100) |> collect, k::Int64 = 10)
  return rcopy(Matrix{Float64},  R"psBasis(x = $tau, K = $k)$X")
end

function test_sm(; tau::Vector{Float64} = range(0, 1, length = 100) |> collect, k::Int64 = 10)
  return rcopy(Matrix{Float64}, R"sm($tau, K = $k, rankZ = .99999999,  spline.degree = 3, diff.ord = 2, centerBase = T)")
end 


function flfosr1_Sun_Kowal(;Y::Matrix{Float64}, X::Matrix{Float64}, z::Vector{Int64})

    M, Lp1 = size(X)
    Tn, M = size(Y)
    R"""
    start_time <- proc.time()
    run_flfosr <- flfosr1(Y = $Y, X = $X, z = $z)
    end_time <- proc.time()
    print(end_time-start_time)
    uz <- unique($z)
    """
    omega_post = Array{Float64}(undef,  M, 10, 1000 )
    gamma_post = Array{Float64}(undef, length(rcopy(R"uz")),10, 1000)
    alphaf_post = Array{Float64}(undef, Lp1, Tn, 1000 )

    for s in 1000:1999

        omega_post[:, :,s -1000 + 1] = rcopy(R"run_flfosr$w_post[[$s]]")
        gamma_post[:, :,s -1000 + 1] = rcopy(R"run_flfosr$ga_post[[$s]]" )

    end

    for l in 1:Lp1
        alphaf_post[l, :, :] = rcopy(R"run_flfosr$alphaf_post[[$l]]")
    end 
    return Dict("X"=> rcopy(Matrix{Float64}, R"run_flfosr$X"), "B"=> rcopy(Matrix{Float64}, R"run_flfosr$B"), "w_post"=>omega_post, "ga_post"=> gamma_post, "alphaf_post"=> alphaf_post,  "sig_eps_post"=> rcopy(R"run_flfosr$sig_e_post"), "sig_alpha_post"=> rcopy(Matrix{Float64},R"run_flfosr$sig_alpha_post")  ,"run_time" =>rcopy(R"start_time-end_time") )

end


function compare_bands(;Y1::Matrix{Float64}, Y2::Matrix{Float64}, alpha::Float64 = 0.05)
  @assert size(Y1) == size(Y2)
  Tn, Nf = size(Y1)
  d1 = postf(Y1, alpha)
  d2 = postf(Y2, alpha)

  d1_bandwidth = abs.(  d1["q_alpha_half"]  -  d1["q_one_m_alpha_half"]  )
  d2_bandwidth = abs.(  d2["q_alpha_half"]  -  d2["q_one_m_alpha_half"]  )

  return plot(plot(range(0,1 ,length = Tn) |>collect, d1_bandwidth, color = "red"), plot(range(0,1 ,length = Tn) |>collect, d2_bandwidth, color = "purple"), ylims = (min(min(d1_bandwidth ...) ,min(d2_bandwidth ...) ),   max(max(d1_bandwidth ...) ,max(d2_bandwidth ...) )   )     )

end




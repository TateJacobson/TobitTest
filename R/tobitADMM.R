tobitADMM = function(x, y, left, C, M, t, a = 3.7, 
                     nlambda = 100, lambda.factor = 0.01, lambda = NULL, 
                     rho = 1, abs_tol = 1e-4, rel_tol = 1e-2, nr_tol = 1e-6, theta_tol = 1e-6, 
                     full.model = FALSE, maxit = 1000){ 
    p = ncol(x)
    n = nrow(x)
    
    if( !isTRUE(all.equal( M, sort(M) )) ) stop("Elements of M must be ordered from low to high.")

    #Define M complement
    Mc = setdiff(0:p, M)
    #Check that x is standardized
    x_colsds = apply(x, 2, sd)
    if( max(x_colsds)/min(x_colsds) > 5 ) warning("The columns of x are on very different scales. You may want to standardize the columns of x before using this function.")
    
    d = (y > left)
    
    #Setup for computing lambda max
    x_M = x[ , setdiff(M,0), drop = F]
    mod_M = tobitnet_innerC(xin = x_M, yin = y, cin = left, lambda1 = 0, lambda2 = 0, 
                            pf1 = rep(1,ncol(x_M)), pf2 = rep(1,ncol(x_M)),delta_init = rep(0,ncol(x_M)), standardize = F)
    r_M = rep(mod_M$d0, n) + x_M%*%mod_M$delta
    gamma_M = mod_M$gamma
    
    mod_null = tobitnet_innerC(xin = matrix(0, 0, 0), yin = y, cin = left, lambda1 = 0, lambda2 = 0, 
                               pf1 = rep(1,p), pf2 = rep(1,p), delta_init = rep(0,p), standardize = F)
    r_null = rep(mod_null$d0, n)
    
    gamma_null = mod_null$gamma
    
    #Set lambda solution path if not provided
    if(is.null(lambda)){
        #Compute lambda max
        lmax_null = max( vapply(1:p, function(j) abs( elln_prime(x0j = x[!d,j], x1j = x[d,j], y1 = y[d] - left, r0 = r_null[!d], r1 = r_null[d], gamma = gamma_null) ), FUN.VALUE = numeric(1) ) )
        lmax_M = max( vapply(setdiff(Mc,0), function(j) abs( elln_prime(x0j = x[!d,j], x1j = x[d,j], y1 = y[d] - left, r0 = r_M[!d], r1 = r_M[d], gamma = gamma_M) ), FUN.VALUE = numeric(1) ) )
        lmax = max(lmax_null, lmax_M)
        
        lmin = lmax*lambda.factor
        
        llseq = seq(log(lmin), log(lmax), length.out = nlambda)
        lambda = exp(llseq)
    } else {
        nlambda = length(lambda)
    }
    
    #Set inits
    delta_init = rep(0,p+1)
    gamma_init = 1
    eta_init = rep(0, length(Mc))
    nu1_init = rep(0, nrow(C))
    nu2_init = rep(0, length(Mc))
    
    sigma_vec = rep(0, nlambda)
    b0_vec = rep(0, nlambda)
    beta_mat = matrix(0, nrow = p, ncol = nlambda)
    
    iter_vec = rep(0, nlambda)
    eta_mat = matrix(0, nrow = length(Mc), ncol = nlambda)
    nu_mat = matrix(0, nrow = nrow(C) + length(Mc), ncol = nlambda) 
    
    #Add an intercept column to x
    x = cbind(rep(1,n), x)
    
    for(l in 1:nlambda){
        tadmm = tobitADMM_C(x = x, yin = y, left = left, Cin = C, M = M, Mc = Mc, tin = t, a = a, lambda = lambda[l], rho = rho, 
                            delta_init = delta_init, gamma_init = gamma_init, eta_init = eta_init, nu1_init = nu1_init, nu2_init = nu2_init,
                            abs_tol = abs_tol, rel_tol = rel_tol, nr_tol = nr_tol, theta_tol = theta_tol,
                            full_model = full.model, maxit = maxit)
         
        #Save output
        sigma_vec[l] = tadmm$sigma
        b0_vec[l] = tadmm$b0
        beta_mat[,l] = tadmm$beta
        
        iter_vec[l] = tadmm$admm_iter
        eta_mat[,l] = tadmm$eta
        nu_mat[,l] = tadmm$nu
        
        #Update inits for warm starts
        delta_init = tadmm$delta_temp
        gamma_init = tadmm$gamma
        eta_init = tadmm$eta
        nu1_init = tadmm$nu[1:nrow(C)]
        nu2_init = tadmm$nu[(nrow(C) + 1):(nrow(C) + length(Mc))]
    }
    
    return(
        structure(
            list(b0 = b0_vec,
                 beta = beta_mat,
                 sigma = sigma_vec,
                 lambda = lambda,
                 left = left, C = C, M = M, t = t,
                 full.model = full.model
                 ),
            class = "tobitADMM"
            )
        )
}


#Function to pick optimal lambda for tobitADMM object based on GIC
pick_lambda_tobit = function(tADMM, x, y, ic = c("max", "GIC", "BIC")){
    ic = match.arg(ic)
    stopifnot(class(tADMM) == "tobitADMM")
    d = (y > tADMM$left)
    
    n = nrow(x)
    p = ncol(x)
    
    if(ic == "max"){
        cn = max(log(n), log(log(n))*log(p))
    } else if (ic =="GIC"){
        cn = log(log(n))*log(p)
    } else if (ic == "BIC"){
        cn = log(n)
    }
    
    nlambda = length(tADMM$lambda)
    gic = numeric( nlambda )
    neg.loglik = numeric( nlambda )
    beta.zero.norm = numeric( nlambda )
    
    for(l in 1:nlambda){
        #get fitted values (in terms of delta)
        tADMM.fitted = tADMM$b0[l]/tADMM$sigma[l] +  x%*%tADMM$beta[,l]/tADMM$sigma[l]
        
        #compute -loglik
        neg.loglik[l] = n*elln(y1 = y[d], r0 = tADMM.fitted[!d], r1 = tADMM.fitted[d], gamma = 1/tADMM$sigma[l], left = tADMM$left)
        beta.zero.norm[l] = sum(tADMM$beta[,l] != 0)
        
        gic[l] = neg.loglik[l] + cn*beta.zero.norm[l]
    }
    
    ind = which.min(gic)
    
    return(
        list(
            gic = gic,
            cn = cn,
            neg.loglik = neg.loglik,
            beta.zero.norm = beta.zero.norm,
            ind = ind, 
            lambda = tADMM$lambda[ind], 
            b0 = tADMM$b0[ind],
            beta = tADMM$beta[,ind],
            sigma = tADMM$sigma[ind]
            )
        )
}

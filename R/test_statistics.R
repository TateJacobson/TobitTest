#Compute inverse of Hessian
compute_hessian_inv = function(x, y, left, delta, gamma, M, S){
    d = (y > left)
    y = y - left
    
    n1 = sum(d)
    y1 = y[d]
    x0 = cbind(1, x[!d,])
    x1 = cbind(1, x[d,])
    
    r0 = x0%*%delta - left*gamma
    
    hessian_inv = 
        qr.solve( 
            rbind( 
                cbind(
                    t(x0[, union(M, S) ])%*%diag( as.vector(h(-r0)) )%*%x0[, union(M, S) ] + crossprod(x1[, union(M, S) ]), 
                    crossprod(x1[, union(M, S) ], -y1) 
                    ), 
                cbind(
                    crossprod(-y1, x1[, union(M, S) ]), 
                    sum(y1^2) +  n1*gamma^(-2)
                    )
                )
            )
    
    return(hessian_inv)
}

#Compute Wald test statistic
compute_wald = function(tADMM, x, y, ic = c("max", "GIC", "BIC")){
    ic = match.arg(ic)
    stopifnot(class(tADMM) == "tobitADMM")
    stopifnot(tADMM$full.model)
    
    #pull parameters from tADMM
    left = tADMM$left; C = tADMM$C; M = tADMM$M; t = tADMM$t
    
    #select model with GIC
    tADMM.selection = pick_lambda_tobit(tADMM, x, y, ic = ic)
    delta = c(tADMM.selection$b0, tADMM.selection$beta)/tADMM.selection$sigma
    gamma = 1/tADMM.selection$sigma
    
    Sa = setdiff(which(tADMM.selection$beta != 0), M)
    
    #Add 1 to all indices in M and S to account for intercept column.
    Mnew = M + 1
    Snew = Sa + 1
    
    #Want 0 in either M or Sa
    if(!(0 %in% M)) Snew = c(1, Snew)
    
    omega_a = compute_hessian_inv(x, y, left = left, delta = delta, gamma = gamma, M = Mnew, S = Snew)
    
    gamma_ind = nrow( omega_a )
    
    #Account for 0 in M
    m = length(Mnew)
    inds = 1:m
    
    tob_wald = 
        t( C%*%delta[Mnew] - gamma*t  )%*%
        qr.solve( cbind(C, -t)%*%omega_a[c(inds, gamma_ind), c(inds, gamma_ind)]%*%t( cbind(C, -t)) )%*%
        ( C%*%delta[Mnew] - gamma*t  )
    
    return(
        list(TW = tob_wald,
             selected_model = tADMM.selection,
             tADMM = tADMM
        )
    )
}

#Compute score test statistic
compute_score = function(tADMM, x, y, ic = c("max", "GIC", "BIC")){
    ic = match.arg(ic)
    stopifnot(class(tADMM) == "tobitADMM")
    stopifnot(!tADMM$full.model)
    
    #pull parameters from tADMM
    left = tADMM$left; C = tADMM$C; M = tADMM$M; t = tADMM$t
    
    #select model with GIC
    tADMM.selection = pick_lambda_tobit(tADMM, x, y, ic = ic)
    delta = c(tADMM.selection$b0, tADMM.selection$beta)/tADMM.selection$sigma
    gamma = 1/tADMM.selection$sigma
    
    S0 = setdiff(which(tADMM.selection$beta != 0), M)
    
    #Add 1 to all indices in M and S to account for intercept column.
    Mnew = M + 1
    Snew = S0 + 1
    
    #Want 0 in either M or Sa
    if(!(0 %in% M)) Snew = c(1, Snew)
    
    omega_0 = compute_hessian_inv(x, y, left, delta = delta, gamma = gamma, M = Mnew, S = Snew)
    
    d = (y > left)
    y = y - left
    
    n1 = sum(d)
    y1 = y[d]
    x0 = cbind(1, x[!d,])
    x1 = cbind(1, x[d,])
    
    r0 = x0%*%delta - left*gamma
    r1 = x1%*%delta - left*gamma
    
    score = c(crossprod(x1[,union(Mnew, Snew)], (gamma*y1 - r1)) - crossprod(x0[,union(Mnew, Snew)],  g(-r0)),
              n1*gamma^(-1) - crossprod(y1, (gamma*y1 - r1) ) )
    tob_score = t(score)%*%omega_0%*%score
    
    return(
        list(TS = tob_score,
             selected_model = tADMM.selection,
             tADMM = tADMM
        )
    )
}

#Compute LRT statistic
compute_lrt = function(tADMM.full, tADMM.null, x, y, ic = c("max", "GIC", "BIC")){
    ic = match.arg(ic)
    stopifnot(class(tADMM.full) == "tobitADMM")
    stopifnot(tADMM.full$full.model)
    stopifnot(class(tADMM.null) == "tobitADMM")
    stopifnot(!tADMM.null$full.model)
    stopifnot(all.equal(list(tADMM.full$left, tADMM.full$M), list(tADMM.null$left, tADMM.null$M) ))
    
    #Select models with GIC
    tob.null.selection = pick_lambda_tobit(tADMM.null, x, y, ic = ic)
    tob.full.selection = pick_lambda_tobit(tADMM.full, x, y, ic = ic)
    
    full_delta = c(tob.full.selection$b0, tob.full.selection$beta)/tob.full.selection$sigma
    full_gamma = 1/tob.full.selection$sigma
    
    null_delta = c(tob.null.selection$b0, tob.null.selection$beta)/tob.null.selection$sigma
    null_gamma = 1/tob.null.selection$sigma
    
    n = length(y)
    left = tADMM.full$left
    d = (y > left)
    
    x = cbind(1, x)
    
    full_r = x%*%full_delta
    null_r = x%*%null_delta
    
    tob_lrt = -2*n*( elln(y1 = y[d], r0 = full_r[!d], r1 = full_r[d], gamma = full_gamma, left = left) -
                         elln(y1 = y[d], r0 = null_r[!d], r1 = null_r[d], gamma = null_gamma, left = left) )

    return(
        list(TL = tob_lrt,
             selected_full_model = tob.full.selection,
             selected_null_model = tob.null.selection,
             tADMM.full = tADMM.full,
             tADMM.null = tADMM.null 
        )
    )    
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//Helper functions
//[[Rcpp::export]]
arma::vec g(const arma::vec& x){
	int n = x.n_elem;
	arma::vec gout(n);
	for(int i = 0; i < n; i++){
		if(x(i) < -37){
			gout(i) = -x(i);
		} else {
			gout(i) = arma::normpdf(x(i))/arma::normcdf(x(i));
		}
	}
	return gout;
}

//[[Rcpp::export]]
arma::vec h(const arma::vec& x){
	arma::vec hout = g(x)%(x + g(x));
	return hout;
}

double update_etaj(const double& z, const double& lambda, const double& rho, const double& a){
	double etaj;
	double sgn = (z > 0) - (z < 0);
	if( std::abs(z) <= lambda*(rho+1)/rho ){
		etaj = sgn*std::max(std::abs(z) - lambda/rho, 0.0);
	} else if ((lambda*(rho+1)/rho < std::abs(z)) && (std::abs(z) <= a*lambda) ){
		etaj = rho/(rho*(a-1) - 1)*((a-1)*z - sgn*a*lambda/rho);
	} else {
		etaj = z;
	}
	return etaj;
}

arma::vec scad_penalty(const arma::vec& t, const double& lambda, const double& a){
	int n = t.n_elem;
	arma::vec pt(n);
	for(int i = 0; i < n; i++){
		if( std::abs(t(i)) <= lambda ){
			 pt(i) = lambda*std::abs(t(i)) ;
		} else if (lambda < std::abs(t(i)) && std::abs(t(i)) <= a*lambda){
			pt(i) = ( - t(i)*t(i) + 2*a*lambda*std::abs(t(i)) - lambda*lambda )/(2*(a-1));
		} else {
			pt(i) = (a+1)*lambda*lambda/2;
		}
	}
	return pt;
}

//[[Rcpp::export]]
double elln(const arma::vec& y1, const arma::vec& r0, const arma::vec& r1, const double& gamma, const double& left = 0){
	double n0 = r0.n_elem; 
	double n1 = r1.n_elem;
	arma::vec logL0 = arma::log( arma::normcdf(-(r0 - gamma*left)) );
	arma::vec logL1 = arma::square(gamma*y1 - r1);
	return (0.5*n1*std::log(2*M_PI) - n1*std::log(gamma) + 0.5*arma::sum( logL1 ) - arma::sum( logL0 ))/(n0 + n1) ;
}

//[[Rcpp::export]]
double elln_prime(const arma::vec& x0j, const arma::vec& x1j, const arma::vec& y1, const arma::vec& r0, const arma::vec& r1, const double& gamma, const double& left = 0){
	double n0 = r0.n_elem; 
	double n1 = r1.n_elem;
	arma::vec logL0_prime = g(-(r0 - gamma*left) )%x0j;
	arma::vec logL1_prime = -(gamma*y1 - r1)%x1j;
	return ( arma::sum(logL0_prime) + arma::sum(logL1_prime) )/(n0 + n1);
}

// [[Rcpp::export]]
Rcpp::List tobitADMM_C (const arma::mat& x, const arma::colvec& yin, const double& left, 
						const arma::mat& Cin, arma::uvec& M, arma::uvec& Mc, const arma::colvec& tin, 
						const double& a, const double& lambda, const double& rho, 
						const arma::vec& delta_init, const double& gamma_init, const arma::vec& eta_init, const arma::vec& nu1_init, const arma::vec& nu2_init, 
						const double& abs_tol, const double& rel_tol, const double& nr_tol, const double& theta_tol, const bool& full_model, const int& maxit){

	double n = x.n_rows;
    int p = x.n_cols - 1;
	
	int m = M.n_elem;
	int mc = Mc.n_elem;
	
	//Modify t to account for shift in the intercept if 0 in M
	arma::vec t(tin.size());
	arma::mat C( size(Cin) );

	int nconst = C.n_rows;
	if(full_model){
		t.zeros(); 
		C.zeros();
	} else {
		if( 0 == M(0) ){
			t = tin - left*Cin.col(0);
		} else {
			t = tin;
		}
		C = Cin;
	}

	//left-shift y
	arma::vec y = yin - left;
	
	arma::uvec d = arma::find(y > 0);
	arma::uvec dc = arma::find(y <= 0);
	
	arma::mat xtemp(n, p+1, arma::fill::none);
    //Reorder columns so those in M come first
    xtemp.cols(0,m-1) = x.cols(M);
    xtemp.cols(m,p) = x.cols(Mc);
    
    //Separate censored and uncensored observations
    arma::colvec y1 = y.elem(d); 
    
    arma::mat x0 = xtemp.rows(dc);
    arma::mat x1 = xtemp.rows(d);
	
    double n1 = d.n_elem;
    
	//Create theta = (delta, gamma)
	arma::colvec theta_current(p+2), theta_new(p+2), theta_step(p+2), theta_prev(p+2);
	
	//Warm start theta
	//reshuffle inits to reflect x column reordering
	theta_current.subvec(0,m-1) = delta_init.elem(M); 
	theta_current.subvec(m,p) = delta_init.elem(Mc);
	theta_current(p+1) = gamma_init;
	
	//Warm start eta and nu
	arma::colvec eta_current = eta_init;
	arma::colvec eta_new(mc), eta_step(mc);
	
	arma::colvec nu1_current = nu1_init;
    arma::colvec nu1_step(nconst);
	
	if(full_model){
		nu1_current.zeros();
		nu1_step.zeros();
	} 
	
    arma::colvec nu2_current = nu2_init;
	arma::colvec nu2_step(mc);
    
    //r
	arma::colvec r0 = x0*theta_current.subvec(0,p);
    arma::colvec r1 = x1*theta_current.subvec(0,p);
	
	//Initialize tolerances and norms
	double nr_norm = 9999;
	double pri_tol = 0;
	double dual_tol = 0;
	double pri_norm = 9999;
	double dual_norm = 9999;
	double theta_step_norm = 9999;
    
	arma::mat K = arma::join_cols(arma::join_rows(C, arma::zeros(nconst, mc), -t), 
								  arma::join_rows(arma::zeros(mc, m), arma::eye(mc, mc), arma::zeros(mc, 1)) );
	
	if( 0 == Mc(0) ){
		K(nconst,m) = 0;
    }
	
	int admm_iter = 0;
	
    //outer loop: ADMM
    while(pri_norm > pri_tol || dual_norm > dual_tol || theta_step_norm > theta_tol ){
		theta_prev = theta_current;
		
		double objective_current = elln(y1, r0, r1, theta_current(p+1)) + rho*0.5*arma::sum( arma::square(C*theta_current.subvec(0,m-1) - theta_current(p+1)*t + nu1_current/rho) );
		if( 0 == Mc(0) ){
			objective_current += rho*0.5*arma::sum( arma::square(theta_current.subvec(m+1, p) - eta_current.tail(mc-1) + nu2_current.tail(mc-1)/rho) );
		} else {
			objective_current += rho*0.5*arma::sum( arma::square(theta_current.subvec(m, p) - eta_current + nu2_current/rho) );
		}

		//NR update for theta
        do{
			if( 0 == Mc(0) ){
				nu2_current(0) = 0;
				eta_current(0) = 0;
			}
			
			//Setup for backtracking
			double nr_step_size = 1;
			double objective_new;
			
			//NR calculations
			arma::mat U = arma::join_cols(arma::join_rows( x0.t()*diagmat(h(-r0))*x0 + x1.t()*x1, -x1.t()*y1 ),
										  arma::join_rows( -y1.t()*x1, y1.t()*y1 + n1*std::pow(theta_new(p+1), -2.0) ) )/n
						  + rho*K.t()*K;
			
			arma::vec S = arma::join_cols(- x1.t()*(theta_current(p+1)*y1 - r1) + x0.t()*g(-r0), 
										  - n1/theta_current(p+1) + y1.t()*(theta_current(p+1)*y1 - r1) )/n 
						  + rho*K.t()*K*theta_current
					      - rho*arma::join_cols(arma::zeros(m, 1), eta_current, arma::zeros(1,1)) 
						  + arma::join_cols(C.t()*nu1_current, nu2_current, -t.t()*nu1_current);
			
			arma::vec nr_step = - arma::inv_sympd(U)*S;
			
			do{
				theta_new = theta_current + nr_step_size*nr_step;
				
				objective_new = elln(y1, x0*theta_new.subvec(0,p), x1*theta_new.subvec(0,p), theta_new(p+1)) + rho*0.5*arma::sum(arma::square(C*theta_new.subvec(0,m-1) - theta_new(p+1)*t + nu1_current/rho));
				if(0 == Mc(0)){
					objective_new += rho*0.5*arma::sum(arma::square(theta_new.subvec(m+1, p) - eta_current.tail(mc-1) + nu2_current.tail(mc-1)/rho));
				} else {
					objective_new += rho*0.5*arma::sum(arma::square(theta_new.subvec(m, p) - eta_current + nu2_current/rho));
				}

				//If the objective decreases with the step, then take the step
				if(objective_new < objective_current){
					break;
				//if even a tiny step does not decrease the objective, then do not take a step	
				} else if (arma::sum( arma::square(theta_new - theta_current) ) < nr_tol/10){
					theta_new = theta_current;
					objective_new = objective_current;
					break;
				}
				nr_step_size = 0.5*nr_step_size;
			} while(true);
			
			//Update objective_current, theta_current, r0, and r1
			objective_current = objective_new;
			theta_step = theta_new - theta_current;
			theta_current = theta_new;
			r0 = x0*theta_current.subvec(0,p);
			r1 = x1*theta_current.subvec(0,p);

			nr_norm = arma::sum( arma::square(theta_step) );
								
        } while( nr_norm > nr_tol );
		
		theta_step_norm = arma::sum( arma::square(theta_current - theta_prev) );
		
        //eta update
        if( 0 == Mc(0) ){ 
			for(int j = 1; j < mc; j++ ){
				eta_new(j) = update_etaj(theta_current(m+j) + nu2_current(j)/rho, lambda, rho, a); 
			}
        } else {
			for(int j = 0; j < mc; j++ ){
				eta_new(j) = update_etaj(theta_current(m+j) + nu2_current(j)/rho, lambda, rho, a);
			}
        }
        
        eta_step = eta_new - eta_current;
        eta_current = eta_new; 
		
        //Dual (nu) update
        if(!full_model) nu1_step = rho*(C*theta_current.subvec(0, m-1) - theta_current(p+1)*t);
        nu2_step = rho*( theta_current.subvec(m, p) - eta_current );
        if( 0 == Mc(0) ){
            nu2_step(0) = 0;
        }
        
        if(!full_model) nu1_current = nu1_current + nu1_step; 
        nu2_current = nu2_current + nu2_step; 
		
		//Compute primal and dual residual norms
		pri_norm = arma::norm(arma::join_cols(nu1_step, nu2_step), 2)/rho;
		dual_norm = rho*arma::norm(eta_step, 2);
		
		//Compute primal and dual tolerances
		pri_tol = std::sqrt(nconst + mc)*abs_tol + 
				  rel_tol*std::max({ arma::norm(arma::join_cols(C*theta_current.subvec(0, m-1), theta_current.subvec(m, p)), 2), arma::norm(eta_current, 2), arma::norm(theta_current(p+1)*t, 2) }) ;
		dual_tol = std::sqrt(p+1)*abs_tol + rel_tol*arma::norm( arma::join_cols( C.t()*nu1_current, nu2_current ), 2);
		
		//Increment ADMM iteration count
		admm_iter++;
		
		if(admm_iter > maxit - 1){
		    break;
		}
    }
    
    //clean up output
    arma::vec delta(p+1), delta_temp(p+1);
    delta.elem(M) = theta_current.subvec(0,m-1);
	delta_temp.elem(M) = theta_current.subvec(0,m-1);
	if(lambda == 0){
		delta.elem(Mc) = theta_current.subvec(m, p);
	} else {
		//alternative: return eta instead of delta_Mc
		delta.elem(Mc) = eta_current; 
	}
	delta_temp.elem(Mc) = theta_current.subvec(m, p);
	
	if( 0 == Mc(0) ){
		delta(0) = theta_current(m);
		delta_temp(0) = theta_current(m);
	} 
	
	double gamma = theta_current(p+1);
    double b0 = delta(0)/gamma;
    arma::vec beta = delta.subvec(1,p)/gamma;
    double sigma = 1/gamma;
    
    b0 = b0 + left;
	
	return Rcpp::List::create(
			Rcpp::Named("delta") = delta,
			Rcpp::Named("delta_temp") = delta_temp,
			Rcpp::Named("gamma") = gamma,
			Rcpp::Named("eta") = eta_current,
			Rcpp::Named("nu") = arma::join_cols(nu1_current, nu2_current),
			Rcpp::Named("admm_iter") = admm_iter,
			Rcpp::Named("b0") = b0,
            Rcpp::Named("beta") = beta,
            Rcpp::Named("sigma") = sigma);
}

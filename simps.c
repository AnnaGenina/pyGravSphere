#include <math.h>  // include file for fabs and sin
#include <stdio.h> // include file for /////////////////printf and perror
#include <errno.h>
#include <time.h>
#include <float.h>

#define PI                         (3.14159265358979312)
#define tol_beta                   (0.001)
#define G                          (0.0000043009125)



const double half_pi = PI/2.;



double calc_sin(double x){
    int i = 1;
    double cur = x;
    double acc = 1;
    double fact= 1;
    double power = x;
    while (fabs(acc) > .00001 &&   i < 100){
        fact *= ((2*i)*(2*i+1));
        power *= -1 * x*x;
        acc =  power / fact;
        cur += acc;
        i++;
    }
    return cur;

}

double calc_cos(double x){
    int i = 1;
    double cur = 1;
    double acc = 1;
    double fact= 1;
    double power = 1;
    while (fabs(acc) > .00001 &&   i < 100){
        fact *= ((2*i -1)*(2*i));
        power *= -1 * x*x;
        acc =  power / fact;
        cur += acc;
        i++;
    }
    return cur;

}

double table_lookup(double check, int len_table, double *table){

	double interval = half_pi/(len_table - 1);

	int counter = check/interval;

  double val;

	if (counter >= (len_table - 1)){
		counter = len_table - 2;}


  val = table[counter] + (table[counter+1] - table[counter])/interval * (check - counter*interval);

	return val;

}

double table_lookup_dist(double check, int len_table, double table[len_table][2]){
	double max_rad = table[len_table-1][0];
	double min_rad = table[0][0];
	double interval = (max_rad - min_rad)/(len_table - 1);
  int counter = (check-min_rad)/interval;
  //printf("Counter %d \n", counter);
  double val;

	if (counter >= (len_table - 1)){
		counter = len_table - 2;}


  val = table[counter][1] + (table[counter+1][1] - table[counter][1])/interval * (check - table[counter][0]);


  //double val = table[counter];
	return pow(10,val);

}







double adaptiveSimpsonsAux(double (*f)(double, double, double*, double*, double*, int, double[*][2]), double a, double fa, double b, double fb,  double eps,
                          double whole, double m, double fm, int rec, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){
       	double lm = (a+m)/2.; double flm=f(lm,R,rho_params,beta_params,plum_params,rows,mass_table); double left = fabs(m-a)/ 6. * (fa + 4. * flm + fm);
	double rm = (b+m)/2.; double frm = f(rm,R,rho_params,beta_params,plum_params,rows,mass_table); double right = fabs(b-m)/ 6. * (fb + 4. * frm + fm);


	if ( a == lm ) { errno = EDOM; return whole; }

        double delta = left + right - whole;


	if (rec <= 0 && errno != EDOM) errno = ERANGE;
        if (rec <= 0 || (fabs(delta) <= eps * fabs(whole)))
              return left + right;

        return adaptiveSimpsonsAux(f, a, fa, m, fm, eps, left, lm, flm , rec-1, R,rho_params,beta_params,plum_params,rows,mass_table) +
                adaptiveSimpsonsAux(f, m,fm,b,fb,eps,right,rm,frm,rec-1,R,rho_params,beta_params,plum_params,rows,mass_table);

}

double adaptiveSimpsons(double (*f)(double,double,double*,double*,double*, int, double[*][2]),     // function ptr to integrate
                       double a, double b,      // interval [a,b]
                       double epsilon,         // error tolerance
                       int maxRecDepth,double R, double *rho_params, double *beta_params, double *plum_params,  int rows, double mass_table[rows][2]) {     // recursion cap


    errno = 0;
    double h = fabs(a-b);
    if (h == 0) return 0;
    double fa = f(a,R,rho_params,beta_params,plum_params,rows, mass_table); double fb = f(b,R,rho_params,beta_params,plum_params, rows,mass_table);
    double m = (a+b)/2.; double fm = f(m,R,rho_params,beta_params,plum_params,rows,mass_table);
    double whole = fabs(a-b)/6.*(fa + 4.*fm + fb);
    return adaptiveSimpsonsAux(f,  a, fa,  b,  fb,  epsilon,
                          whole,m, fm, maxRecDepth,   R,rho_params,beta_params,plum_params, rows, mass_table);
}


double Simpsons(double (*f)(double,double,double*,double*,double*, int, double[*][2]),     // function ptr to integrate
                       double a, double b,      // interval [a,b]
                       int bins,         // error tolerance
                       double R, double *rho_params, double *beta_params, double *plum_params,  int rows, double mass_table[rows][2]){

int n = bins -1;
double h = (b - a)/n;
double y0 = f(a+0*h, R,rho_params,beta_params,plum_params,rows,mass_table);
double yn = f(a+n*h,  R,rho_params,beta_params,plum_params,rows,mass_table);
double sum1 = 0;
double sum2 = 0;
double tot = 0;

for(int i=1; i<n; i++){
		if(i%2 == 0){
			sum1 = sum1 + f(a+i*h,  R,rho_params,beta_params,plum_params,rows,mass_table);}
		else{
			sum2 = sum2 + f(a+i*h,  R,rho_params,beta_params,plum_params,rows,mass_table); }}

	tot =(h/3)*(y0 + yn + 2*sum1 +4*sum2);

  return tot;
}


double adaptiveSimpsonsAuxSingle(double (*f)(double, double*), double a, double fa, double b, double fb,  double eps,
                          double whole, double m, double fm, int rec, double *args){
       	double lm = (a+m)/2.; double flm=f(lm,args); double left = fabs(m-a)/ 6. * (fa + 4. * flm + fm);
	double rm = (b+m)/2.; double frm = f(rm,args); double right = fabs(b-m)/ 6. * (fb + 4. * frm + fm);


	if ( a == lm ) { errno = EDOM; return whole; }

        double delta = left + right - whole;


	if (rec <= 0 && errno != EDOM) errno = ERANGE;
        if (rec <= 0 || (fabs(delta) <= eps * fabs(whole)))
              return left + right;

        return adaptiveSimpsonsAuxSingle(f, a, fa, m, fm, eps, left, lm, flm , rec-1, args) +
                adaptiveSimpsonsAuxSingle(f, m,fm,b,fb,eps,right,rm,frm,rec-1,args);

}

double adaptiveSimpsonsSingle(double (*f)(double,double*),     // function ptr to integrate
                       double a, double b,      // interval [a,b]
                       double epsilon,         // error tolerance
                       int maxRecDepth,double *args) {     // recursion cap

    errno = 0;
    double h = fabs(a-b);
    if (h == 0) return 0;
    double fa = f(a,args); double fb = f(b,args);
    double m = (a+b)/2.; double fm = f(m,args);
    double whole = fabs(a-b)/6.*(fa + 4.*fm + fb);
    return adaptiveSimpsonsAuxSingle(f,  a, fa,  b,  fb,  epsilon,
                          whole,m, fm, maxRecDepth,   args);
}

double mass(double x, double* args){
	// Note log
	x = pow(10,x);
	double result = x * log(10)* 4. * PI * x * x * pow(10,args[1]) * pow(x/pow(10,args[2]), -args[5]) * pow(1. + pow(x/pow(10,args[2]), args[3]), -(args[4] - args[5])/args[3] );
	return result ;
}




double BrokenPowerLaw(double x, double *args){

	double bins[]={0.25*args[0],0.5*args[0],1.0*args[0],2.0*args[0],4*args[0]};

	double gammas[] = {args[2],args[3],args[4],args[5],args[6]};

	double rho = pow(10, args[1])	;



    	if (x < bins[0]){
        return rho * pow(x/bins[0],-gammas[0]); }


	else if ((bins[0] <= x ) &  ( x < bins[1])){
            return rho * pow(x/bins[1], -gammas[1]) * pow(bins[1]/bins[0], -gammas[1] ) ; }
        else if  ((bins[1] <= x ) & (x < bins[2])) {
            return rho* pow(x/bins[2],-gammas[2]) * pow(bins[2]/bins[1], -gammas[2]) * pow(bins[1]/bins[0], -gammas[1] ) ; }
        else if ((bins[2] <= x) & ( x < bins[3])) {
            return  rho * pow(x/bins[3],-gammas[3])* pow(bins[3]/bins[2], -gammas[3]) * pow(bins[2]/bins[1], -gammas[2]) * pow(bins[1]/bins[0],-gammas[1]) ; }

        else if (( bins[3] <= x )& (x < bins[4])){
            return rho * pow(x/bins[4],-gammas[4]) * pow(bins[4]/bins[3],-gammas[4])* pow(bins[3]/bins[2],-gammas[3]) * pow(bins[2]/bins[1],-gammas[2]) * pow(bins[1]/bins[0],-gammas[1]) ; }

       else{
           return rho * pow(x/bins[4],-gammas[4]) * pow(bins[4]/bins[3],-gammas[4])* pow(bins[3]/bins[2],-gammas[3]) * pow(bins[2]/bins[1],-gammas[2]) * pow(bins[1]/bins[0],-gammas[1]) ; }

}



double BrokenPowerLawMass(double x, double *args){


	double bins[]={0.25*args[0],0.5*args[0],1.0*args[0],2.0*args[0],4*args[0]};

	double gammas[] = {args[2],args[3],args[4],args[5],args[6]};

	



	double rho = pow(10, args[1])* 4 * PI	;



    if (x < bins[0]){
        return rho * pow(bins[0], gammas[0]) * pow(x, 3 - gammas[0])/(3 - gammas[0]) ;}

    else if ((bins[0] <= x ) &  ( x < bins[1])){
        return rho * pow(bins[0], gammas[0]) * pow(bins[0], 3 - gammas[0])/(3 - gammas[0]) + rho *pow(bins[1], gammas[1]) * pow(bins[1]/bins[0], - gammas[1]) * ( pow(x, 3-gammas[1]) -  pow(bins[0], 3-gammas[1])  )/(3. -gammas[1]) ; }

    else if ((bins[1] <= x ) &  ( x < bins[2])){
        return rho * pow(bins[0], gammas[0]) * pow(bins[0], 3 - gammas[0])/(3 - gammas[0]) + rho *pow(bins[1], gammas[1]) * pow(bins[1]/bins[0], - gammas[1]) * ( pow(bins[1], 3-gammas[1]) -  pow(bins[0], 3-gammas[1])  )/(3. -gammas[1]) + rho * pow(bins[2], gammas[2]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2]) * ( pow(x, 3-gammas[2]) -  pow(bins[1], 3-gammas[2])  )/(3. -gammas[2]) ;}

    else if ((bins[2] <= x ) &  ( x < bins[3])){
        return rho * pow(bins[0], gammas[0]) * pow(bins[0], 3 - gammas[0])/(3 - gammas[0]) + rho *pow(bins[1], gammas[1]) * pow(bins[1]/bins[0], - gammas[1]) * ( pow(bins[1], 3-gammas[1]) -  pow(bins[0], 3-gammas[1])  )/(3. -gammas[1]) + rho * pow(bins[2], gammas[2]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2]) * ( pow(bins[2], 3-gammas[2]) -  pow(bins[1], 3-gammas[2])  )/(3. -gammas[2]) + rho *pow(bins[3], gammas[3]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * ( pow(x, 3-gammas[3]) -  pow(bins[2], 3-gammas[3])  )/(3. -gammas[3])  ;}

    else if ((bins[3] <= x ) &  ( x < bins[4])){
        return rho * pow(bins[0], gammas[0]) * pow(bins[0], 3 - gammas[0])/(3 - gammas[0]) + rho *pow(bins[1], gammas[1]) * pow(bins[1]/bins[0], - gammas[1]) * ( pow(bins[1], 3-gammas[1]) -  pow(bins[0], 3-gammas[1])  )/(3. -gammas[1]) + rho * pow(bins[2], gammas[2]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2]) * ( pow(bins[2], 3-gammas[2]) -  pow(bins[1], 3-gammas[2])  )/(3. -gammas[2]) + rho *pow(bins[3], gammas[3]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * ( pow(bins[3], 3-gammas[3]) -  pow(bins[2], 3-gammas[3])  )/(3. -gammas[3]) + rho *pow(bins[4], gammas[4]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * pow(bins[4]/bins[3], - gammas[4]) * (pow( x, 3-gammas[4]) -  pow(bins[3], 3-gammas[4])  )/(3. -gammas[4]) ; }

    else{
        return rho * pow(bins[0], gammas[0]) * pow(bins[0], 3 - gammas[0])/(3 - gammas[0]) + rho *pow(bins[1], gammas[1]) * pow(bins[1]/bins[0], - gammas[1]) * ( pow(bins[1], 3-gammas[1]) -  pow(bins[0], 3-gammas[1])  )/(3. -gammas[1]) + rho * pow(bins[2], gammas[2]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2]) * ( pow(bins[2], 3-gammas[2]) -  pow(bins[1], 3-gammas[2])  )/(3. -gammas[2]) + rho *pow(bins[3], gammas[3]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * ( pow(bins[3], 3-gammas[3]) -  pow(bins[2], 3-gammas[3])  )/(3. -gammas[3]) + rho *pow(bins[4], gammas[4]) * pow(bins[1]/bins[0], - gammas[1]) * pow(bins[2]/bins[1], - gammas[2])* pow(bins[3]/bins[2], - gammas[3]) * pow(bins[4]/bins[3], - gammas[4]) * (pow( x, 3-gammas[4]) -  pow(bins[3], 3-gammas[4])  )/(3. -gammas[4]) ;}

    }






double Plummer2d(double r, double* args){
	return (1./(PI * args[1]*args[1])) * pow(1. + (r/args[1])*(r/args[1]),-2) ;}

double Plummer3d(double r, double* args){
	return (3./(4.*PI * args[1]*args[1]*args[1])) * pow(1. + (r/args[1])*(r/args[1]),-2.5); }

double Plummer2dM(double r, double* args){
	return (args[0]/(PI * args[1]*args[1])) * pow(1. + (r/args[1])*(r/args[1]),-2) ;}

double Plummer3dM(double r, double* args){
	return ((3.*args[0])/(4.*PI * args[1]*args[1]*args[1])) * pow(1. + (r/args[1])*(r/args[1]),-2.5); }


double nu_mass(double r, double* args){
	 return args[0]*r*r*r/(args[1]*args[1]*args[1]) * pow((1. + r*r/(args[1]*args[1])),-1.5);
}

double mass_star(double x, double* args){
	double sum0 = 0;
  int count = 0;

	for(int i=0; i < 3; i++){
		double args2[] = {args[count], args[count + 1]} ;
		
		sum0 = sum0 + nu_mass(x,args2);
		count = count + 2;
			}



	return sum0;
}



double Plummer2dSph(double r, double* args){
	double sum0 = 0;
	int count = 0;

	for(int i=0; i < 3; i++){
		double args2[] = {args[count], args[count + 1]} ;
		
		sum0 = sum0 + Plummer2dM(r,args2);
		count = count + 2;
			}

	return sum0;
}

double Plummer3dSph(double r, double* args){
	double sum0 = 0;
	int count = 0;

	for(int i=0; i < 3; i++){
		double args2[] = {args[count], args[count + 1]} ;
		
		sum0 = sum0 + Plummer3dM(r,args2);
		count = count + 2;
			}

	return sum0;
}


double Plummer3dSphConst(double r, double* args){
	double sum0 = 0;
	int count = 0;

	for(int i=0; i < 3; i++){
		double args2[] = {args[count], args[count + 1]} ;
		
		sum0 = sum0 + Plummer3dM(r,args2);
		count = count + 2;
			}

	return sum0;
}





double Baes(double r, double* args){

	//double u = pow((r/args[2]),args[3]);

	return  args[0] + (args[1]-args[0])*(1.0/(1.0 + pow((args[2]/r),args[3]))) ; 
	//return (args[0] +  args[1]*u)/(1+ u);


}


double IntBaes(double r, double* args){

    	//return  pow( r, 2*args[0] ) * pow((1 + pow((r/args[2]),args[3])), (2. * (args[1] - args[0])/args[3]) )  ;

	return pow(r,(2.0*args[1])) *pow( (pow((args[2]/r),args[3])+1.0),(2.0/args[3]*(args[1]-args[0]))) ;

}


double IntConstant(double r, double* args){

	return pow(r,(2*args[0]));

}

double nusigintSph(double theta, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){


	double sin_calc = sin(theta);
	double cos_calc = cos(theta);


	double r = pow(10, R)/cos_calc;
	double check = log10(r);
        double mass_r = table_lookup_dist(check, rows, mass_table);
	return IntBaes(r,beta_params) * Plummer3dSph(r,plum_params) * mass_r  * sin_calc;

}


double nusigintSphConstant(double theta, double R,  double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){


	double sin_calc = sin(theta);
	double cos_calc = cos(theta);

	double r = pow(10, R)/cos_calc;

	double check = log10(r);

	double mass_r = table_lookup_dist(check, rows, mass_table);

	return IntConstant(r,beta_params) * Plummer3dSph(r,plum_params) * mass_r * sin_calc;

}

double nusigintSphPowerLaw(double theta, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){


	double r = pow(10, R)/cos(theta);
	double check = log10(r);
  double mass_r = table_lookup_dist(check, rows, mass_table);
	return IntBaes(r,beta_params) * Plummer3dSph(r,plum_params) * (BrokenPowerLawMass(r, rho_params) + mass_r)  * sin(theta) ;

}


double nusigintSphPowerLawConstant(double theta, double R, double* rho_params, double* beta_params, double* plum_params, int rows, double mass_table[rows][2]){

	double r = pow(10, R)/cos(theta);
  double check = log10(r);
	double mass_r = table_lookup_dist(check, rows, mass_table);

	return IntConstant(r,beta_params) * Plummer3dSph(r,plum_params) *  (BrokenPowerLawMass(r, rho_params) + mass_r) * sin(theta);

}



double vsp1int(double r, double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){
	r = pow(10, r);

	double check = log10(r);

	double mass_r = table_lookup_dist(check, rows, mass_table);
	return r * log(10) * mass_r * (5. - 2*Baes(r, beta_params)) * r ;


}

double vsp1intConstant(double r, double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){
	r = pow(10, r);

	double check = log10(r);

	double mass_r = table_lookup_dist(check, rows, mass_table);
	
	return r * log(10) * mass_r * (5. - 2*beta_params[0]) * r ;


}




double vsp2int(double r, double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){

	r = pow(10, r);

	double check = log10(r);

	double mass_r = table_lookup_dist(check, rows, mass_table);

	return r * log(10) * mass_r * (7. - 6*Baes(r, beta_params)) * r*r*r;


}





double vsp2intConstant(double r,double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){

	r = pow(10, r);

	double check = log10(r);

	double mass_r = table_lookup_dist(check, rows, mass_table);

	return r * log(10) * mass_r * (7. - 6*beta_params[0]) * r*r*r ;


}



double resBaes(double theta,  double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){


	//double cth = table_lookup(theta, st, costhetas);

	double cth = cos(theta);
	double r = R/cth;

	double check = log10(r);

  double mass_r = table_lookup_dist(check, rows, mass_table);


	return 2.  * mass_r * (1. - Baes(r,beta_params) * cth * cth) / (cth * cth) ;

}


double resConstant(double theta, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]){



	double cth = calc_cos(theta);

	double r = R/cth;

	double check = log10(r);

	double mass_r = table_lookup_dist(check, rows, mass_table);
	return 2.  * mass_r * (1. - beta_params[0]  * cth * cth)  / (cth * cth) ;

}








double round2(double val, double dig){
	double fac = pow(10, dig);
	
	return round(fac * val)/fac;

}




double likelihood(double *sig2, double err2,  int tracers, double *vel_rel){

	double like =  0.0;
	for (int i = 0; i < tracers; i++){
		like = like + log(  1./sqrt(2.*PI * (sig2[i] + err2)))  + (-0.5 * (vel_rel[i])*(vel_rel[i])/(sig2[i] + err2)) ;  }

	return like;
}


double likelihoodBin(double *sig2, double *sigmas, double *sigerrs,  int tracers){

	double like =  0.0;
	for (int i = 0; i < tracers; i++){
		like = like + log(  1./sqrt(2.*PI * (sigerrs[i]*sigerrs[i])))  + (-0.5 * (   (sqrt(sig2[i]) - sigmas[i]) * (sqrt(sig2[i]) - sigmas[i])   )/(sigerrs[i] * sigerrs[i])   ) ;  }
	

	return like;
}

double likelihoodBinVSP(double *sig2, double *sigmas, double *sigerrs,  int tracers, double vsp1, double vsp2, double *vir_shape){
	double like =  0.0;
	for (int i = 0; i < tracers; i++){
		like = like  + (-0.5 * (   (sqrt(sig2[i]) - sigmas[i]) * (sqrt(sig2[i]) - sigmas[i])   )/(sigerrs[i] * sigerrs[i])   ) ;  }


		like = like - 0.5 * (vir_shape[0] - vsp1) * (vir_shape[0] - vsp1)/(vir_shape[1] * vir_shape[1]) ; 
		like = like - 0.5 * (vir_shape[2] - vsp2) * (vir_shape[2] - vsp2)/(vir_shape[3]*vir_shape[3]) ;   
	
	return like;
}






double GetLikeBinAnisSphVSP(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r,  double rh, int rows, double min_rad, double max_rad){

	    double like = 0;
	    double upper_limit = half_pi - 1e-5;



	    if (beta_params[1] != beta_params[0]){




         int bins = rows;
		     double mass_table[rows][2];

		     double spacing = (max_rad - min_rad)/rows;
         double bot;
		     double top;
         bot = min_rad;
         top = min_rad + spacing;
		     double mdm = 0.	;
		     double mstar;
		     double tot_mass_star = mass_star(50, plum_params);
		     double u;
		     //double star_extent = log10(50);
		     mass_table[0][0] = bot;
		     mdm = mdm + adaptiveSimpsonsSingle(mass, -15, bot, 1e-05, 25, rho_params);
		     mass_table[0][1] = log10(mdm + mass_star(pow(10,bot), plum_params) /tot_mass_star  * mass_star_3r );
      
		     for (int i = 1; i < rows; i++){
			      u = pow(10, top);
			      mass_table[i][0] = top;
			      mdm = mdm + adaptiveSimpsonsSingle(mass, bot,top, 1e-05, 25, rho_params);
			      
			    mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
				
			      mass_table[i][1] = log10(mdm + mstar);

			      bot = top;
			      top = top + spacing;
					}


		     double nusig_table[rows][2];
		      spacing =  (max_rad - min_rad)/rows;
		      bot = min_rad;
		      top = min_rad + spacing;
		      double ns;



		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;
			      ns = Simpsons(nusigintSph, 0, upper_limit, bins,bot, rho_params, beta_params, plum_params, rows, mass_table);

		        nusig_table[i][1] = log10(ns/IntBaes(pow(10, nusig_table[i][0]),beta_params )/ pow(10,bot));

			      bot = top;
			      top = top + spacing;

					}





		


		    double nusig_m_table[rows][2];
		     for (int i = 0; i < rows; i++){

			nusig_m_table[i][0] = nusig_table[i][0];

			nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1]));



					}
	


		     double results[tracers];

		     bot = min_rad;

		    for (int i = 0; i < tracers; i++){


			         results[i] = G *1./Plummer2dSph(R_list[i],plum_params) * R_list[i] * Simpsons(resBaes, 0, upper_limit, bins, R_list[i], rho_params, beta_params, plum_params, rows,nusig_table);

					}

		    bot = min_rad;

	double vsp1 = 0.;
	double vsp2 = 0.;

	if (vir_shape[0] != 0.0){

	      
		vsp1 = vsp1 + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp1 = vsp1 * (2./5.) * G * G;}

	if (vir_shape[2] != 0.0){
	      
		vsp2 = vsp2 + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp2 = vsp2 * (4./35.) * G * G; }


	

	like = likelihoodBinVSP(results, sigmas, sigerrs, tracers,vsp1, vsp2,vir_shape);


					   }

	else{



         int bins = rows;
		     double mass_table[rows][2];
		     double spacing = (max_rad - min_rad)/(rows-1);
		     double bot;
		     double top;
         bot = min_rad;
         top = min_rad + spacing;
		     double mdm = 0.	;
		     double mstar;
		     double tot_mass_star = mass_star(50, plum_params);
		     double u;
		    

		     mass_table[0][0] = bot;
		     mdm = mdm + adaptiveSimpsonsSingle(mass, -15,bot, 1e-05, 25, rho_params);
		     mass_table[0][1] = log10(mdm + mass_star(pow(10,bot), plum_params) /tot_mass_star  * mass_star_3r );

		     for (int i = 1; i < rows; i++){
			      u = pow(10, top);
			      mass_table[i][0] = top;
			      mdm = mdm + adaptiveSimpsonsSingle(mass, bot,top, 1e-05, 25, rho_params);

			    
				    mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
					
			      mass_table[i][1] = log10(mdm + mstar);

			    
			      bot = top;
			      top = top + spacing;
					}


		     double nusig_table[rows][2];
		     spacing = (max_rad - min_rad)/(rows-1);
		      bot = min_rad;
		      top = min_rad + spacing;
		      double ns = 0.	;



		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;
			      ns = Simpsons(nusigintSphConstant,0,upper_limit, bins,bot,  rho_params, beta_params, plum_params, rows, mass_table);
	          nusig_table[i][1] = log10(ns/IntConstant(pow(10, nusig_table[i][0]),beta_params )/ pow(10,bot));
			      bot = top;
			      top = top + spacing;
				

					}


	


		    double nusig_m_table[rows][2];
		     for (int i = 0; i < rows; i++){

			nusig_m_table[i][0] = nusig_table[i][0];

			nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1]));



					}
		

		     double results[tracers];

		     bot = min_rad;





			double res = 0;
		    for (int i = 0; i < tracers; i++){

			res = 0;
				res = res + Simpsons(resConstant,0, upper_limit, bins,R_list[i], rho_params, beta_params, plum_params, rows,nusig_table);


				res = res * G * 1./Plummer2dSph(R_list[i],plum_params) * R_list[i] ;
				results[i] = res;


					}


		     bot = min_rad;

	double vsp1 = 0.;
	double vsp2 = 0.;

	 if (vir_shape[0] != 0.0){

	      
		vsp1 = vsp1 + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp1 = vsp1 * (2./5.) * G * G;}

	if (vir_shape[2] != 0.0){
	      
		vsp2 = vsp2 + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp2 = vsp2 * (4./35.) * G * G; }


	

	like = likelihoodBinVSP(results, sigmas, sigerrs, tracers,vsp1, vsp2,vir_shape);

			}

      return like;
}

void ZhaoFit(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r,  double rh, int rows, double *results, double *vsp1, double *vsp2, double min_rad, double max_rad){

	    
	    double upper_limit = half_pi - 1e-5;



	    if (beta_params[1] != beta_params[0]){




         int bins = rows;
		     double mass_table[rows][2];

		     double spacing = (max_rad - min_rad)/rows;
         double bot;
		     double top;
         bot = min_rad;
         top = min_rad + spacing;
		     double mdm = 0.	;
		     double mstar;
		     double tot_mass_star = mass_star(50, plum_params);
		     double u;
		     
		     mass_table[0][0] = bot;
		     mdm = mdm + adaptiveSimpsonsSingle(mass, -15, bot, 1e-05, 25, rho_params);
		     mass_table[0][1] = log10(mdm + mass_star(pow(10,bot), plum_params) /tot_mass_star  * mass_star_3r );
      
		     for (int i = 1; i < rows; i++){
			      u = pow(10, top);
			      mass_table[i][0] = top;
			      mdm = mdm + adaptiveSimpsonsSingle(mass, bot,top, 1e-05, 25, rho_params);
			     
				    mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
				
			      mass_table[i][1] = log10(mdm + mstar);

			      bot = top;
			      top = top + spacing;
					}


		     double nusig_table[rows][2];
		      spacing =  (max_rad - min_rad)/rows;
		      bot = min_rad;
		      top = min_rad + spacing;
		      double ns;



		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;
			      ns = Simpsons(nusigintSph, 0, upper_limit, bins,bot, rho_params, beta_params, plum_params, rows, mass_table);

		        nusig_table[i][1] = log10(ns/IntBaes(pow(10, nusig_table[i][0]),beta_params )/ pow(10,bot));

			      bot = top;
			      top = top + spacing;

					}





		
		    double nusig_m_table[rows][2];
		     for (int i = 0; i < rows; i++){

			nusig_m_table[i][0] = nusig_table[i][0];

			nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1]));



					}
		

		     

		     bot = min_rad;

		    for (int i = 0; i < tracers; i++){


			         results[i] = sqrt(G *1./Plummer2dSph(R_list[i],plum_params) * R_list[i] * Simpsons(resBaes, 0, upper_limit, bins, R_list[i], rho_params, beta_params, plum_params, rows,nusig_table));

					}

		    bot = min_rad;

	    if (vir_shape[0] != 0.0){

			    

	    	vsp1[0] = vsp1[0] + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
	    	   	

		vsp1[0] = vsp1[0] * (2./5.) * G * G;
			}

	if (vir_shape[2] != 0){

		vsp2[0] = vsp2[0] + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		   
		vsp2[0] = vsp2[0] * (4./35.) * G * G;
			}

		    


					   }

	else{



         int bins = rows;
		     double mass_table[rows][2];
		     double spacing = (max_rad - min_rad)/(rows-1);
		     double bot;
		     double top;
         bot = min_rad;
         top = min_rad + spacing;
		     double mdm = 0.	;
		     double mstar;
		     double tot_mass_star = mass_star(50, plum_params);
		     double u;
		     //double star_extent = log10(50);

		     mass_table[0][0] = bot;
		     mdm = mdm + adaptiveSimpsonsSingle(mass, -15,bot, 1e-05, 25, rho_params);
		     mass_table[0][1] = log10(mdm + mass_star(pow(10,bot), plum_params) /tot_mass_star  * mass_star_3r );

		     for (int i = 1; i < rows; i++){
			      u = pow(10, top);
			      mass_table[i][0] = top;
			      mdm = mdm + adaptiveSimpsonsSingle(mass, bot,top, 1e-05, 25, rho_params);

			     
				    mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
					
			      mass_table[i][1] = log10(mdm + mstar);

			      //////printf("Mass table %f %f \n", mass_table[i][0], mass_table[i][1]);

			      bot = top;
			      top = top + spacing;
					}


		     double nusig_table[rows][2];
		     spacing = (max_rad - min_rad)/(rows-1);
		      bot = min_rad;
		      top = min_rad + spacing;
		      double ns = 0.	;



		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;
			      ns = Simpsons(nusigintSphConstant,0,upper_limit, bins,bot,  rho_params, beta_params, plum_params, rows, mass_table);
	          nusig_table[i][1] = log10(ns/IntConstant(pow(10, nusig_table[i][0]),beta_params )/ pow(10,bot));
			      bot = top;
			      top = top + spacing;
				

					}


	


		    double nusig_m_table[rows][2];
		     for (int i = 0; i < rows; i++){

			nusig_m_table[i][0] = nusig_table[i][0];

			nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1]));



					}
		

		    

		     bot = min_rad;





			
		    for (int i = 0; i < tracers; i++){

			
				results[i] =  sqrt(Simpsons(resConstant,0, upper_limit, bins,R_list[i], rho_params, beta_params, plum_params, rows,nusig_table) *G * 1./Plummer2dSph(R_list[i],plum_params) * R_list[i]);


					}


		     bot = min_rad;



	if (vir_shape[0] != 0.0){

			    

	    	vsp1[0] = vsp1[0] + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
	    	   	

		vsp1[0] = vsp1[0] * (2./5.) * G * G;
			}

	if (vir_shape[2] != 0){

		vsp2[0] = vsp2[0] + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		   
		vsp2[0] = vsp2[0] * (4./35.) * G * G;
			}

		    

			}

     
}



double GetLikeBinAnisSphVSPPowerLaw(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r, double rh, int rows, double min_rad, double max_rad){

	double like = 0;
	double upper_limit = half_pi - 1e-5;

	if (beta_params[1] != beta_params[0]){


		 for (int i=2;i < 7;i++){
			   if (fabs(3-rho_params[i]) < 1e-5){
				if (rho_params[i] > 3.){
					rho_params[i] = rho_params[i]+1e-5;	}
				else{
					rho_params[i] = rho_params[i]-1e-5;     }
				     
								}


					}

		 double mass_table[rows][2];

		 double mass_table_star[rows][2];
		 double spacing = (max_rad - min_rad)/(rows-1);  //spacing is the same
    		 int bins = rows;


	     double bot;
	     double top;
	     bot = min_rad;
	     top = min_rad + spacing;
		 double mdm = 0.	;
		 double tot_mass_star = 0;
		 tot_mass_star = mass_star(50, plum_params);
		 double mstar = 0;
		 double u;
		 //double star_extent = log10(50);


		 for (int i = 0; i < rows; i++){
			    u = pow(10,bot);
			    mass_table[i][0] = bot;
			    mass_table_star[i][0] = bot;
			    mdm = BrokenPowerLawMass(u, rho_params);
			   mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
				                        

			    

			        mass_table[i][1] = log10(mdm + mstar);
			      	mass_table_star[i][1] = log10(mstar);
		  		    bot = top;
			        top = top + spacing;  }






		    double nusig_table[rows][2];
		    spacing = (max_rad - min_rad)/(rows-1);
		    bot = min_rad;
		    top = min_rad + spacing;
		    double ns ;

		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;

			      ns = Simpsons(nusigintSphPowerLaw, 0, upper_limit, bins, bot,rho_params, beta_params,plum_params, rows, mass_table_star);

			      nusig_table[i][1] = log10(ns/IntBaes(pow(10, nusig_table[i][0]),beta_params ) / pow(10,bot) );

			      bot = top;
			      top = top + spacing;

					}



		    double nusig_m_table[rows][2];



		    for (int i = 0; i < rows; i++){

			      nusig_m_table[i][0] = nusig_table[i][0];
			      nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1]));   //Dm+ stars


                                      }


		     double results[tracers];
         bot = min_rad;


		    for (int i = 0; i < tracers; i++){

			      results[i] = G * 1./Plummer2dSph(R_list[i],plum_params) * R_list[i] *Simpsons(resBaes, 0, upper_limit, bins, R_list[i], rho_params, beta_params,plum_params, rows,nusig_table);

				                               	}


		    bot = min_rad;


	double vsp1 = 0.;
	
	double vsp2 = 0.;
	

        if (vir_shape[0] != 0.0){

	        
		vsp1 = vsp1 + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp1 = vsp1 * (2./5.) * G * G;}

	if (vir_shape[2] != 0.0){
	        
		vsp2 = vsp2 + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp2 = vsp2 * (4./35.) * G * G; 
		}


	
	
	like = likelihoodBinVSP(results, sigmas, sigerrs, tracers,vsp1, vsp2,vir_shape);




	}

		   
			


	else {         // Anisotropy is constant



		for (int i=2;i < 7;i++){
			   if (fabs(3-rho_params[i]) < 1e-5){
				if (rho_params[i] > 3.){
					rho_params[i] = rho_params[i]+1e-5;	}
				else{
					rho_params[i] = rho_params[i]-1e-5;     }
				     
								}


					}




   int bins = rows;
   double mass_table[rows][2];
   double mass_table_star[rows][2];
   double spacing = (max_rad - min_rad)/(rows-1);
   double bot;
   double top;
   bot = min_rad;
   top = min_rad + spacing;
   double mdm = 0.	;
   double tot_mass_star = 0;
   tot_mass_star = mass_star(50, plum_params);
   double mstar = 0;
   double u;
   //double star_extent = log10(50);

	 for (int i = 0; i < rows; i++){
			 u = pow(10,bot);
			 mass_table[i][0] = bot;
			 mass_table_star[i][0] = bot;
			 mdm = BrokenPowerLawMass(u, rho_params);
			
			 mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
				      
			     	   

			      mass_table[i][1] = log10(mdm + mstar);
			      mass_table_star[i][1] = log10(mstar);
			      bot = top;
			      top = top + spacing;

					}



     double nusig_table[rows][2];
     spacing = (max_rad - min_rad)/(rows-1);
     bot = min_rad;
     top = min_rad + spacing;
     double ns	;


		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;
			      ns = Simpsons(nusigintSphPowerLawConstant, 0, upper_limit, bins, bot,rho_params, beta_params,plum_params, rows, mass_table_star);

			      nusig_table[i][1] = log10(ns/IntConstant(pow(10, nusig_table[i][0]),beta_params )/ pow(10,bot));
			      bot = top;
			      top = top + spacing;

					}


		    double nusig_m_table[rows][2];
		     for (int i = 0; i < rows; i++){

			nusig_m_table[i][0] = nusig_table[i][0];

			nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1])); }

			
		     double results[tracers];

		     bot = min_rad;



		    for (int i = 0; i < tracers; i++){

			       results[i] = Simpsons(resConstant,0, upper_limit,bins,R_list[i], rho_params, beta_params,plum_params, rows,nusig_table)* G * 1./Plummer2dSph(R_list[i],plum_params) * R_list[i] ;


					}


		    bot = min_rad;



	double vsp1 = 0.;
	double vsp2 = 0.;

	

        if (vir_shape[0] != 0.0){

	      
		vsp1 = vsp1 + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp1 = vsp1 * (2./5.) * G * G;}

	if (vir_shape[2] != 0.0){
	      
		vsp2 = vsp2 + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		vsp2 = vsp2 * (4./35.) * G * G; }


	

	like = likelihoodBinVSP(results, sigmas, sigerrs, tracers,vsp1, vsp2,vir_shape);


		}
      return like;
}



void PowerLawFit(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r, double rh, int rows, double* results, double* vsp1, double* vsp2, double min_rad, double max_rad){

	
	double upper_limit = half_pi - 1e-5;

	if (beta_params[1] != beta_params[0]){


		 for (int i=2;i < 7;i++){
			   if (fabs(3-rho_params[i]) < 1e-5){
				if (rho_params[i] > 3.){
					rho_params[i] = rho_params[i]+1e-5;	}
				else{
					rho_params[i] = rho_params[i]-1e-5;     }
				     
								}


					}

		 double mass_table[rows][2];

		 double mass_table_star[rows][2];
		 double spacing = (max_rad - min_rad)/(rows-1);  
     int bins = rows;


     double bot;
     double top;
     bot = min_rad;
     top = min_rad + spacing;
		 double mdm = 0.	;
		 double tot_mass_star = 0;
		 tot_mass_star = mass_star(50, plum_params);
		 double mstar = 0;
		 double u;
		 


		 for (int i = 0; i < rows; i++){
			    u = pow(10,bot);
			    mass_table[i][0] = bot;
			    mass_table_star[i][0] = bot;
			    mdm = BrokenPowerLawMass(u, rho_params);
			    
	                    mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
				                        

			    
			mass_table[i][1] = log10(mdm + mstar);
		      	mass_table_star[i][1] = log10(mstar);
	  		    bot = top;
			top = top + spacing;  }






		    double nusig_table[rows][2];
		    spacing = (max_rad - min_rad)/(rows-1);
		    bot = min_rad;
		    top = min_rad + spacing;
		    double ns ;

		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;

			      ns = Simpsons(nusigintSphPowerLaw, 0, upper_limit, bins, bot,rho_params, beta_params,plum_params, rows, mass_table_star);

			      nusig_table[i][1] = log10(ns/IntBaes(pow(10, nusig_table[i][0]),beta_params ) / pow(10,bot) );

			      bot = top;
			      top = top + spacing;

					}



		    double nusig_m_table[rows][2];



		    for (int i = 0; i < rows; i++){

			      nusig_m_table[i][0] = nusig_table[i][0];
			      nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1]));   //Dm+ stars


                                      }


		    
         		bot = min_rad;


		    for (int i = 0; i < tracers; i++){

			      results[i] = sqrt(G * 1./Plummer2dSph(R_list[i],plum_params) * R_list[i] *Simpsons(resBaes, 0, upper_limit, bins, R_list[i], rho_params, beta_params,plum_params, rows,nusig_table));

				                               	}


		    bot = min_rad;

	

        if (vir_shape[0] != 0.0){

			    

	    	vsp1[0] = vsp1[0] + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
	    	   	

		vsp1[0] = vsp1[0] * (2./5.) * G * G;
			}

	if (vir_shape[2] != 0.0){

		vsp2[0] = vsp2[0] + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		   
		vsp2[0] = vsp2[0] * (4./35.) * G * G;
			}
			

}


	else {         // Anisotropy is constant



		for (int i=2;i < 7;i++){
			   if (fabs(3-rho_params[i]) < 1e-5){
				if (rho_params[i] > 3.){
					rho_params[i] = rho_params[i]+1e-5;	}
				else{
					rho_params[i] = rho_params[i]-1e-5;     }
				     
								}


					}




   int bins = rows;
   double mass_table[rows][2];
   double mass_table_star[rows][2];
   double spacing = (max_rad - min_rad)/(rows-1);
   double bot;
   double top;
   bot = min_rad;
   top = min_rad + spacing;
   double mdm = 0.	;
   double tot_mass_star = 0;
   tot_mass_star = mass_star(50, plum_params);
   double mstar = 0;
   double u;
   //double star_extent = log10(50);

	 for (int i = 0; i < rows; i++){
			 u = pow(10,bot);
			 mass_table[i][0] = bot;
			 mass_table_star[i][0] = bot;
			 mdm = BrokenPowerLawMass(u, rho_params);
			mstar = mass_star(u, plum_params) /tot_mass_star  * mass_star_3r   ;
				       

			mass_table[i][1] = log10(mdm + mstar);
			mass_table_star[i][1] = log10(mstar);
			bot = top;
			top = top + spacing;

					}



     double nusig_table[rows][2];
     spacing = (max_rad - min_rad)/(rows-1);
     bot = min_rad;
     top = min_rad + spacing;
     double ns	;


		     for (int i = 0; i < rows; i++){

			      nusig_table[i][0] = bot;
			      ns = Simpsons(nusigintSphPowerLawConstant, 0, upper_limit, bins, bot,rho_params, beta_params,plum_params, rows, mass_table_star);

			      nusig_table[i][1] = log10(ns/IntConstant(pow(10, nusig_table[i][0]),beta_params )/ pow(10,bot));
			      bot = top;
			      top = top + spacing;

					}


		    double nusig_m_table[rows][2];
		     for (int i = 0; i < rows; i++){

			nusig_m_table[i][0] = nusig_table[i][0];

			nusig_m_table[i][1] = log10(pow(10,nusig_table[i][1]) * pow(10, mass_table[i][1])); }

			
		     bot = min_rad;

		    for (int i = 0; i < tracers; i++){

			       results[i] = sqrt(Simpsons(resConstant,0, upper_limit,bins,R_list[i], rho_params, beta_params,plum_params, rows,nusig_table)* G * 1./Plummer2dSph(R_list[i],plum_params) * R_list[i] );


					}


		    bot = min_rad;

			
	      if (vir_shape[0] != 0.0){

			    

	    	vsp1[0] = vsp1[0] + Simpsons(vsp1int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
	    	   	

		vsp1[0] = vsp1[0] * (2./5.) * G * G;
			}

	if (vir_shape[2] != 0.0){

		vsp2[0] = vsp2[0] + Simpsons(vsp2int, min_rad, max_rad, bins,0, rho_params, beta_params,plum_params, rows,nusig_m_table);
		   
		vsp2[0] = vsp2[0] * (4./35.) * G * G;
			}
			
		}
     
}





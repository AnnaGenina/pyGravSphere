/* simps.i */
%module gravsphere

%{
 #define SWIG_FILE_WITH_INIT
 #include "gravsphere.h"
%}



%include "numpy.i"

%init %{
import_array();
%}


%apply (double* IN_ARRAY1, int DIM1) {(double* R_list, int tracers1)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigmas, int tracers2)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigerrs, int tracers3)};

%apply (double* IN_ARRAY1, int DIM1){(double* rho_params, int l1)}
%apply (double* IN_ARRAY1, int DIM1){(double* beta_params, int l2)}
%apply (double* IN_ARRAY1, int DIM1){(double* plum_params, int l3)}
%apply (double* IN_ARRAY1, int DIM1){(double* vir_shape, int l4)}


%inline %{
   double GetLikeBinAnisSphVSPfunc(double *R_list, int tracers1,double *sigmas,int tracers2, double *sigerrs, int tracers3, double* rho_params, int l1, double* beta_params, int l2, double* plum_params, int l3,  double* vir_shape, int l4, double stellar_mass_3r, double rh, int rows, double min_rad, double max_rad){

	double like;
	like =  GetLikeBinAnisSphVSP(R_list, sigmas, sigerrs, tracers1, rho_params,beta_params,plum_params,vir_shape, stellar_mass_3r,  rh, rows, min_rad, max_rad);
	return like;


}


%}



%apply (double* IN_ARRAY1, int DIM1) {(double* R_list, int tracers1)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigmas, int tracers2)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigerrs, int tracers3)};

%apply (double* IN_ARRAY1, int DIM1){(double* rho_params, int l1)}
%apply (double* IN_ARRAY1, int DIM1){(double* beta_params, int l2)}
%apply (double* IN_ARRAY1, int DIM1){(double* plum_params, int l3)}
%apply (double* IN_ARRAY1, int DIM1){(double* vir_shape, int l4)}

%inline %{
   double GetLikeBinAnisSphVSPPowerLawfunc(double *R_list, int tracers1,double *sigmas,int tracers2, double *sigerrs, int tracers3, double* rho_params, int l1,  double* beta_params, int l2,  double* plum_params, int l3, double* vir_shape,int l4, double stellar_mass_3r, double rh, int rows, double min_rad, double max_rad){

	double like;
	like =  GetLikeBinAnisSphVSPPowerLaw(R_list, sigmas, sigerrs, tracers1,  rho_params,beta_params,plum_params,vir_shape, stellar_mass_3r,  rh, rows, min_rad, max_rad);
	return like;


}


%}


%apply (double* IN_ARRAY1, int DIM1) {(double* R_list, int tracers1)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigmas, int tracers2)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigerrs, int tracers3)};

%apply (double* IN_ARRAY1, int DIM1){(double* rho_params, int l1)}
%apply (double* IN_ARRAY1, int DIM1){(double* beta_params, int l2)}
%apply (double* IN_ARRAY1, int DIM1){(double* plum_params, int l3)}
%apply (double* IN_ARRAY1, int DIM1){(double* vir_shape, int l4)}
%apply (double* INPLACE_ARRAY1, int DIM1){(double* results, int tracers)}
%apply (double *INPLACE_ARRAY1, int DIM1) { (double *vsp1, int lv1 )};
%apply (double *INPLACE_ARRAY1, int DIM1) { (double *vsp2, int lv2 )};

%inline %{
   extern void PowerLawFitfunc(double *R_list, int tracers1,double *sigmas,int tracers2, double *sigerrs, int tracers3, double* rho_params, int l1,  double* beta_params, int l2,  double* plum_params, int l3, double* vir_shape,int l4, double stellar_mass_3r, double rh, int rows, double* results, int tracers, double* vsp1,int lv1, double* vsp2, int lv2, double min_rad, double max_rad){
	
	
	PowerLawFit(R_list, sigmas, sigerrs, tracers1,  rho_params,beta_params,plum_params,vir_shape, stellar_mass_3r,  rh, rows, results, vsp1,vsp2, min_rad, max_rad);

}



%}


%apply (double* IN_ARRAY1, int DIM1) {(double* R_list, int tracers1)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigmas, int tracers2)};
%apply (double* IN_ARRAY1, int DIM1) {(double* sigerrs, int tracers3)};

%apply (double* IN_ARRAY1, int DIM1){(double* rho_params, int l1)}
%apply (double* IN_ARRAY1, int DIM1){(double* beta_params, int l2)}
%apply (double* IN_ARRAY1, int DIM1){(double* plum_params, int l3)}
%apply (double* IN_ARRAY1, int DIM1){(double* vir_shape, int l4)}
%apply (double* INPLACE_ARRAY1, int DIM1){(double* results, int tracers)}
%apply (double *INPLACE_ARRAY1, int DIM1) { (double *vsp1, int lv1 )};
%apply (double *INPLACE_ARRAY1, int DIM1) { (double *vsp2, int lv2 )};

%inline %{
   extern void ZhaoFitfunc(double *R_list, int tracers1,double *sigmas,int tracers2, double *sigerrs, int tracers3, double* rho_params, int l1,  double* beta_params, int l2,  double* plum_params, int l3, double* vir_shape,int l4, double stellar_mass_3r, double rh, int rows, double* results,int tracers,  double* vsp1, int lv1, double* vsp2, int lv2, double min_rad, double max_rad){

	
	ZhaoFit(R_list, sigmas, sigerrs, tracers1,  rho_params,beta_params,plum_params,vir_shape, stellar_mass_3r,  rh, rows, results, vsp1, vsp2, min_rad, max_rad);
	

}
%}


#include <math.h> 
#include <stdio.h> 
#include <errno.h>
#include <time.h>
#include <float.h>

#define PI                         (3.14159265358979312)
#define tol_beta                   (0.001)
#define G                          (0.0000043009125)


double calc_sin(double x);

double calc_cos(double x);

double table_lookup(double check, int len_table, double *table);

double adaptiveSimpsonsAux(double (*f)(double, double, double*, double*, double*, int, double[*][2], double*, double*, int), double a, double fa, double b, double fb,  double eps,
                          double whole, double m, double fm, int rec, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double adaptiveSimpsons(double (*f)(double,double,double*,double*,double*, int, double[*][2], double*, double*, int),     // function ptr to integrate
                       double a, double b,      // interval [a,b]
                       double epsilon,         // error tolerance
                       int maxRecDepth,double R, double *rho_params, double *beta_params, double *plum_params,  int rows, double mass_table[rows][2]) ;


double adaptiveSimpsonsAuxSingle(double (*f)(double, double*), double a, double fa, double b, double fb,  double eps,
                          double whole, double m, double fm, int rec, double *args);

double adaptiveSimpsonsSingle(double (*f)(double,double*),     // function ptr to integrate
                       double a, double b,      // interval [a,b]
                       double epsilon,         // error tolerance
                       int maxRecDepth,double *args) ;

double mass(double x, double* args);

double BrokenPowerLaw(double x, double *args);

double BrokenPowerLawMass(double x, double *args);

double Plummer2d(double r, double* args);

double Plummer3d(double r, double* args);

double Plummer2dM(double r, double* args);

double Plummer3dM(double r, double* args);

double nu_mass(double r, double* args);

double mass_star(double x, double* args);

double Plummer2dSph(double r, double* args);

double Plummer3dSph(double r, double* args);

double Plummer3dSphConst(double r, double* args);

double Baes(double r, double* args);

double IntBaes(double r, double* args);

double IntConstant(double r, double* args);

double nusigintSph(double theta, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double nusigintSphConstant(double theta, double R,  double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double nusigintSphPowerLaw(double theta, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double nusigintSphPowerLawConstant(double theta, double R, double* rho_params, double* beta_params, double* plum_params, int rows, double mass_table[rows][2]);

double vsp1int(double r, double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double vsp1intConstant(double r, double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);
double vsp2int(double r, double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double vsp2intConstant(double r,double filler, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);
double resBaes(double theta,  double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double resConstant(double theta, double R, double *rho_params, double *beta_params, double *plum_params, int rows, double mass_table[rows][2]);

double integrand(double r, double *args,  int rows, double mass_table[rows][2]);

double integrandSph(double r, double *args,  int rows, double mass_table[rows][2]);

double round2(double val, double dig);

double likelihood(double *sig2, double err2,  int tracers, double *vel_rel);

double likelihoodBin(double *sig2, double *sigmas, double *sigerrs,  int tracers);

double likelihoodBinVSP(double *sig2, double *sigmas, double *sigerrs,  int tracers, double vsp1, double vsp2, double *vir_shape);

double GetLikeBinAnisSphVSP(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r, double rh, int rows, double min_rad, double max_rad);

void ZhaoFit(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r,  double rh, int rows, double *results, double* vsp1, double* vsp2, double min_rad, double max_rad);

double GetLikeBinAnisSphVSPPowerLaw(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r, double rh, int rows, double min_rad, double max_rad);

void PowerLawFit(double *R_list, double *sigmas, double *sigerrs, int tracers, double *rho_params,  double *beta_params, double *plum_params, double *vir_shape, double mass_star_3r, double rh, int rows, double* results, double* vsp1, double* vsp2, double min_rad, double max_rad);

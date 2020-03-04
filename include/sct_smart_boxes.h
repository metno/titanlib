#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_blas.h>

/* Runs the SCT by splitting into boxes first

Arguments:
   n: Number of stations, specifies lengths of x, y, z, t arrays:
   x: Array of x-position [m]
   y: Array of y-position [m]
   z: Array of station altitudes [m]
   t: Array of temperatures (in any units)
   nmax: Pointer to one int. Split a box into two if it contains more than this number of stations. Suggested value: 1000.
   nmin: Don't allow boxes to have fewer than this. Suggested value: 100.
   nminprof: Revert to basic atmospheric profile if fewer than this number of stations. Suggested value: 50.
   dzmin: minimum elevation range to fit a vertical profile [m]. Suggested value: 30
   dhmin: minimum value for OI horizontal decorellation length [m]. Suggested value: 10000
   dz: OI vertical decorellation length [m]. Suggested value: 200
   t2pos: Flag an observation if it exceeds this SCT value above the background [^oC]
     One value for each station. Suggested values: 4
   t2neg: Flag an observation if it is below this SCT value below the background [^oC]
     One value for each station. Suggested values: 4
   eps2: Epsilon squared. One value for each station. Suggested values: 0.5

Outputs:
   flags: Output array of flags, one for each station. 0 means passed the SCT, 1 means fails.
   sct: Output array of sct-values (cross validation deviation), one for each station.
   rep: Output array of station coefficient of representativity, one for each station.

Compiling to library:
   gcc SCT_wrapper.c -fPIC -lgslcblas -lgsl -lblas -L/usr/lib -lm -shared -o SCT_wrapper.so
*/
void sct_smart_boxes(int *n,
                 double *x,
                 double *y,
                 double *z,
                 double *t,
                 int *nmax,
                 int *nmin,
                 int *nminprof,
                 double* dzmin,
                 double* dhmin,
                 double* dz,
                 double *t2pos,
                 double *t2neg,
                 double *eps2,
                 int *flags,
                 double *sct,
                 double *rep,
                 int *boxids);

// Structure to contain station information for a box
struct Box {
  int  n; // Number of stations in box
  double *x;
  double *y;
  double *z;
  double *t;
  int *i;  // Index into global array
};

// Structure to contain a list of boxes
struct BoxList {
   int n; // Number of boxes
   struct Box* boxes;
};

// Run the SCT on a particular box
void spatial_consistency_test(struct Box *currentBox, int *nminprof, double* dzmin, double* dhmin, double* dz, double *t2pos, double *t2neg, double *eps2, int *flags, double* sct_out, double* rep_out);
void spatial_consistency_test_mod(int *N, int *obs_to_test, double *x, double *y, double *z, double *elevs, double *t, int *nminprof, double* dzmin, double* dhmin, double* dz, double *t2pos, double *t2neg, double *eps2, int *flags, double* sct_out, double* rep_out);

int compute_vertical_profile(struct Box *box, double meanT, double gamma, double a, double exact_p10, double exact_p90, int nminprof, double dzmin, double *vp);
// optimizer functions
double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data);
double vertical_profile_optimizer_function(const gsl_vector *v, void *data); // GSL format
// vp profile calculations
void basic_vertical_profile(int nz, double *z, double t0, double *t_out);
void vertical_profile(int nz, double *z, double t0, double gamma, double a, double h0, double h1i, double *t_out);

// box division controller (other functions used by this controller not forward declared)
struct BoxList control_box_division(int maxNumStationsInBox, int minNumStationsInBox, struct Box inputBox);

// Helper functions
double compute_quantile(double quantile, double *array, int sizeArray);
double mean(const double *array, int sizeArray);
double max(const double *array, int sizeArray);
double min(const double *array, int sizeArray);
void print_vector(double *vector, int size);
void print_gsl_vector(gsl_vector *vector, int size);
void print_matrix(double **matrix, int rows, int columns);
void print_gsl_matrix(gsl_matrix *matrix, int rows, int columns);
void print_sub_gsl_matrix(gsl_matrix *matrix, int start, int stop);
gsl_matrix* inverse_matrix(const gsl_matrix *matrix);

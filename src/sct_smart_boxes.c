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
#include "sct_smart_boxes.h"

void sct_smart_boxes(int *n, double *x, double *y, double *z, double *t, int *nmax, int *nmin,
  int *nminprof, double* dzmin, double* dhmin, double* dz, double *t2pos, double *t2neg, double *eps2,
  int *flags, double *sct, double *rep, int *boxids) {
   if(n[0] == 0)
      return;

  // put the input box in the struct
  struct Box inputBox;
  inputBox.n = n[0];
  inputBox.x = x;
  inputBox.y = y;
  inputBox.z = z;
  inputBox.t = t;
  inputBox.i = malloc(sizeof(int) * n[0]);
  for(int i = 0; i < n[0]; i++)
     inputBox.i[i] = i;

  // split the box if needed
  struct BoxList box_list = control_box_division(nmax[0], nmin[0], inputBox);

  int nB = box_list.n;
  // printf("SCT wrapper - splitting into %i boxes\n", nB);

  // loop over the boxes to call SCT
  for(int i=0; i<nB; i++) {
    int box_n = box_list.boxes[i].n;
    int * box_i = box_list.boxes[i].i;
    // printf("box: %i \n", i);

    // Run the SCT on the current box
    clock_t start = clock(), diff;
    int* local_flags = malloc(sizeof(int) * box_n);
    double* local_t2pos = malloc(sizeof(double) * box_n);
    double* local_t2neg = malloc(sizeof(double) * box_n);
    double* local_eps2 = malloc(sizeof(double) * box_n);
    double* local_rep = malloc(sizeof(double) * box_n);
    double* local_sct = malloc(sizeof(double) * box_n);
    // printf("box_n %i \n", box_n);

    // Extract values for the current box into contiguous arrays
    for(int r = 0; r < box_n; r++) {
       assert(box_i[r] < n[0]);
       local_t2pos[r] = t2pos[box_i[r]];
       local_t2neg[r] = t2neg[box_i[r]];
       local_eps2[r] = eps2[box_i[r]];
       local_rep[r] = rep[box_i[r]];
    }
    spatial_consistency_test(&box_list.boxes[i], nminprof, dzmin, dhmin, dz, local_t2pos, local_t2neg, local_eps2, local_flags, local_sct, local_rep);

    // Merge flags back into global array
    for(int r = 0; r < box_n; r++) {
       assert(box_i[r] < n[0]);
       flags[box_i[r]] = local_flags[r];
       rep[box_i[r]] = local_rep[r];
       sct[box_i[r]] = local_sct[r];
       boxids[box_i[r]] = i;
    }
    free(local_t2pos);
    free(local_t2neg);
    free(local_eps2);
    free(local_flags);
    free(local_rep);
    free(local_sct);

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    // printf("SCT wrapper - Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);

  }
  free(inputBox.i);
  return;
}


//----------------------------------------------------------------------------//
int compute_vertical_profile(struct Box *box, double meanT, double gamma, double a, double exact_p10, double exact_p90, int nminprof, double dzmin, double *vp)
{
  // optimize inputs for VP (using Nelder-Mead Simplex algorithm)
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss;
  gsl_multimin_function vp_optim;

  int iter = 0;
  int status;
  double size;

  /* Set initial step sizes to 1 */

  // data (params) that needs to be passed into vp
  double nd = (double) box[0].n; // cast + resize to double
  double *z = box[0].z;
  double *t = box[0].t;
  double * data[4] = {&nd, z, t, vp};


  // Check if terrain is too flat
  double z05 = compute_quantile(0.05, z, nd);
  double z95 = compute_quantile(0.95, z, nd);


  bool use_basic = box[0].n < nminprof || (z95 - z05) < dzmin;
  // printf("use_basic=%d\n", use_basic);

  /* Initialize method and iterate */
  gsl_vector* input;
  if(use_basic) {
    vp_optim.n = 1;
    input = gsl_vector_alloc(vp_optim.n);
    gsl_vector_set(input, 0, meanT);
    vp_optim.f = basic_vertical_profile_optimizer_function;
  }
  else {
    vp_optim.n = 5;
    input = gsl_vector_alloc(vp_optim.n);
    gsl_vector_set(input, 0, meanT);
    gsl_vector_set(input, 1, gamma);
    gsl_vector_set(input, 2, a);
    gsl_vector_set(input, 3, exact_p10);
    gsl_vector_set(input, 4, exact_p90);
    vp_optim.f = vertical_profile_optimizer_function;
  }
  ss = gsl_vector_alloc (vp_optim.n);

  gsl_vector_set_all (ss, 1.0);
  vp_optim.params = data;

  s = gsl_multimin_fminimizer_alloc (T, vp_optim.n);
  gsl_multimin_fminimizer_set (s, &vp_optim, input, ss);
  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1e-2);

    if (status == GSL_SUCCESS) {
     // printf ("converged to minimum at\n");
    }
  }
  while (status == GSL_CONTINUE && iter < 100);


  // now actually calculate the vertical profile using these optimized variables
  if(use_basic) { // basic vp
    basic_vertical_profile(box[0].n, box[0].z, gsl_vector_get(s->x, 0), vp);

    // printf("Basic profile: t0: %.4f\n", gsl_vector_get(s->x, 0));
  }
  else { // more complicated vp
    vertical_profile(box[0].n, box[0].z, gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1),
            gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4), vp);

    // printf ("Full profile: t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n", gsl_vector_get(s->x, 0), gsl_vector_get(s->x, 1), gsl_vector_get(s->x, 2), gsl_vector_get(s->x, 3), gsl_vector_get(s->x, 4));
    // printf("VP: %f %f\n", z[0], vp[0]);

  }

  gsl_vector_free(input);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}

//----------------------------------------------------------------------------//
/*
#+ cost function used for optimization of tvertprof parameter
tvertprofbasic2opt<-function(par) {
  te<-tvertprof_basic(z=zopt,t0=par[1],gamma=argv$gamma.standard)
  return(log((mean((te-topt)**2))**0.5))
}
*/
double basic_vertical_profile_optimizer_function(const gsl_vector *v, void *data)
{
  double **p = (double **)data;
  int n = (int) *p[0]; // is of type double but should be an int
  double *z = p[1];
  double *t = p[2];
  double *t_out = p[3];

  // the parameters to mess with
  double t0 = gsl_vector_get(v,0);

  // give everything to vp to compute t_out
  basic_vertical_profile(n, z, t0, t_out);
  // RMS
  double total = 0;
  for(int i=0; i<n; i++) {
    total += pow((t_out[i]-t[i]),2);
  }
  double value = log(pow((total / n),0.5));

  return value;
}

/*
#+ vertical profile of temperature (linear)
tvertprof_basic<-function(z,t0,gamma) {
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
# Output
#  t= array. temperature [K or degC]
#------------------------------------------------------------------------------
  return(t0+gamma*z)
}
*/
void basic_vertical_profile(int nz, double *z, double t0, double *t_out)
{
  double gamma = -0.0065;
  for(int i=0; i<nz; i++) {
    t_out[i] = t0 + gamma*z[i];
  }
}


//----------------------------------------------------------------------------//
/*
#+ cost function used for optimization of tvertprof parameter
tvertprof2opt<-function(par) {
te<-tvertprof(z=zopt,t0=par[1],gamma=par[2],a=par[3],h0=par[4],h1i=par[5])
return(log((mean((te-topt)**2))**0.5))
}
*/
// vector (double t0, double gamma, double a, double h0, double h1i)
// data (int n, double *z, double *t, double *t_out)
double vertical_profile_optimizer_function(const gsl_vector *v, void *data)
{
  // my input dat:a double * data[4] = {&nd, z, t, t_out};
  double **p = (double **)data;
  int n = (int) *p[0]; // is of type double but should be an int
  double *z = p[1];
  double *t = p[2];
  double *t_out = p[3];

  // the parameters to mess with
  double t0 = gsl_vector_get(v,0);
  double gamma = gsl_vector_get(v,1);
  double a = gsl_vector_get(v,2);
  double h0 = gsl_vector_get(v,3);
  double h1i = gsl_vector_get(v,4);

  //printf("t0: %f gamma: %f a: %f h0: %f h1i: %f\n", t0, gamma, a, h0, h1i);

  // give everything to vp to compute t_out
  vertical_profile(n, z, t0, gamma, a, h0, h1i, t_out);
  // RMS
  double total = 0;
  for(int i=0; i<n; i++) {
    total += pow((t_out[i]-t[i]),2);
  }
  double value = log(pow((total / n),0.5));

  //printf("first values in t_out %f %f %f\n", t_out[0], t_out[1], t_out[2]);
  //printf("first values in t %f %f %f\n", t[0], t[1], t[2]);
  //printf("first values in z %f %f %f\n", z[0], z[1], z[2]);
  //printf("optimizer value: %f\n", value);
  return value;
}


/*
#+ vertical profile of temperature (Frei, 2014)
tvertprof<-function(z,t0,gamma,a,h0,h1i) {
# ref:
# Frei, C. (2014). Interpolation of temperature in a mountainous region
#  using nonlinear profiles and non‐Euclidean distances.
#  International Journal of Climatology, 34(5), 1585-1605.
# input
#  z= array. elevations [m amsl]
#  t0= numeric. temperature at z=0 [K or degC]
#  gamma=numeric. temperature lapse rate [K/m]
#  a= numeric. inversion (spatial) length
#  h0= numeric. z where inversion starts [m]
#  h1i= numeric. h0+h1i is z where inversion stops [m]
#       (Frei uses h1 directly, I use an increment to h0 so to avoid ending
#        up with h1<=h0 during the optimization)
# Output
#  t= array. temperature [K or degC]
*/
void vertical_profile(int nz, double *z,
    double t0, double gamma, double a, double h0, double h1i, double *t_out)
{
  // printf("VERTICAL PROFILE %d %f %f %f %f %f\n", nz, t0, gamma, a, h0, h1i);
  double h1 = h0 + fabs(h1i); // h1<-h0+abs(h1i)
  // loop over the array of elevations (z)
  // t_out is an empty array of length nz (e.g. the same length as z)
  for(int i=0; i<nz; i++) {
    // define some bools
    bool z_le_h0 = z[i] <= h0; // z.le.h0<-which(z<=h0)
    bool z_ge_h1 = z[i] >= h1; // z.ge.h1<-which(z>=h1)
    bool z_in = (z[i]>h0 && z[i]<h1); // z.in<-which(z>h0 & z<h1)
    if(z_le_h0) {
      t_out[i] = t0-gamma*z[i]-a;
    }
    if(z_ge_h1) {
      t_out[i] = t0-gamma*z[i];
    }
    if(z_in) {
      t_out[i] = t0-gamma*z[i]-a/2*(1+cos(M_PI*(z[i]-h0)/(h1-h0)));
    }
  }
}

/*
Lussana, C., Uboldi, F., & Salvati, M. R. (2010). A spatial consistency
test for surface observations from mesoscale meteorological networks.
Quarterly Journal of the Royal Meteorological Society, 136(649), 1075-1088.
*/
void spatial_consistency_test(struct Box *currentBox, int *nminprof, double* dzmin, double* dhmin, double* dz, double *t2pos, double *t2neg, double *eps2, int *flags, double* sct_out, double* rep_out)
{
  // break out the box for simplicity
  int n = currentBox[0].n;
  double *x = currentBox[0].x;
  double *y = currentBox[0].y;
  double *z = currentBox[0].z;
  double *t = currentBox[0].t;

  // fill with 0 to keep track of stations that are flagged
  for(int i=0; i<n; i++) {
    flags[i] = 0;
  }

  /*
  Stuff for VP
  */
  double gamma = -0.0065;
  double a = 5.0;
  double meanT = mean(t,n);
  // allocate for output
  double *vp = malloc(sizeof(double) * n);
  for(int i=0; i<n; i++) {
    vp[i] = -999;
  }
  double exact_p10 = compute_quantile(0.10, z, n);
  double exact_p90 = compute_quantile(0.90, z, n);

  // vector (double t0, double gamma, double a, double h0, double h1i)
  // Starting point for optimization
  // printf ("VP input vector set = t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n", meanT, gamma, a, exact_p10, exact_p90);

  // calculate background
  int status = compute_vertical_profile(currentBox, meanT, gamma, a, exact_p10, exact_p90, nminprof[0], dzmin[0], vp);

  // printf("t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n", meanT, gamma, a, exact_p10, exact_p90);

  // printf("status optimizer: %d\n", status);
  // now have temperature profile (vp)
  for(int i=0; i<n; i++) {
    assert(vp[i] !=-999);
  }
  // done calculating VP for the first time
  int sizeWhenProfileCalculated = currentBox[0].n;

  // distance matrices
  double** disth = malloc(sizeof(double*)*n);
  double** distz = malloc(sizeof(double*)*n);
  double *Dh = malloc(sizeof(double)*n);

  // no need to select j since already only have those for a particular box
  // outer product of the matrices
  for(int i=0; i<n; i++) {
    disth[i] = malloc(sizeof(double)*n);
    distz[i] = malloc(sizeof(double)*n);
    double *Dh_vector = malloc(sizeof(double)*(n-1)); // need to remove one since not considering the diagonal
    for(int j=0; j<n; j++) {
      disth[i][j] = pow((pow((x[i]-x[j]),2)+pow((y[i]-y[j]),2)),0.5);
      distz[i][j] = abs(z[i]-z[j]);
      if(i != j) { // do not want to consider the diagonal
        if(i < j) {
          Dh_vector[j-1] = disth[i][j];
        }
        else if(i > j) {
          Dh_vector[j] = disth[i][j];
        }
      }
    }
    Dh[i] = compute_quantile(0.10, Dh_vector, n-1);
    free(Dh_vector);
  }

  // set to optimal Dh for the SCT
  // either Dhmin or the average 10-percentile of the set of distances between a
  // station and all the others
  double Dh_mean = mean(Dh,n);
  free(Dh);
  // printf("Dh: %f\n", Dh_mean);
  if(Dh_mean < dhmin[0]) {
    Dh_mean = dhmin[0];
  }
  // printf("Dh_mean: %f\n", Dh_mean);

  // background error correlation matrix
  float Dz = 200; // good value (use this default)
  gsl_matrix *S, *Sinv;
  S = gsl_matrix_alloc(n,n);
  double **S_global = malloc(sizeof(double*)*n);
  for(int i=0; i<n; i++) {
    S_global[i] = malloc(sizeof(double)*n);
  }
  Sinv = gsl_matrix_alloc(n,n);
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      double value = exp(-.5*pow((disth[i][j]/Dh_mean),2)-.5*pow((distz[i][j]/(dz[0])),2));
      S_global[i][j] = value;
      if(i == j) { // weight the diagonal?? (0.5 default)
        value = value + eps2[i];
      }
      gsl_matrix_set(S,i,j,value);
      //gsl_matrix_set(Sinv,i,j,value); // not yet inverted, but need a copy of S
    }
  }
  // printf("created the S matrix - size1 %lu size2 %lu \n", S->size1, S->size2);
  //print_gsl_matrix(S,n,n);

  // d<-topt-tb
  gsl_vector *d;
  d = gsl_vector_alloc(n);
  double *d_global = malloc(sizeof(double)*n);
  for(int i=0; i<n; i++) {
    gsl_vector_set(d,i,(t[i]-vp[i])); // difference between actual temp and temperature from vertical profile
    d_global[i] = (t[i]-vp[i]);
  }

  /* ---------------------------------------------------
  Beginning of real SCT looping
  ------------------------------------------------------*/
  bool first = true;
  int current_n = n;
  int throwOut = 0;
  // loop for SCT
  // note that current_n should be used inside this loop!!!
  while(1) {
    if(first) {
      // if first time then invert matrix
      //gsl_matrix_memcpy (Sinv, S);
      clock_t start = clock(), diff;
      Sinv = inverse_matrix(S);
      // if use this, then make Sinv a copy of S
      //gsl_linalg_cholesky_decomp1(Sinv);
      //gsl_linalg_cholesky_invert(Sinv); // in place
      diff = clock() - start;
      int msec = diff * 1000 / CLOCKS_PER_SEC;
      // printf("Time taken to invert matrix %d seconds %d milliseconds \n", msec/1000, msec%1000);

      // UN-weight the diagonal of S
      for(int i=0; i<current_n; i++) {
        double value = gsl_matrix_get(S,i,i) - eps2[i];
        gsl_matrix_set(S,i,i,value);
      }
      //printf("S first\n");
      //print_gsl_matrix(S, current_n, current_n); //(int rows, int columns, gsl_matrix *matrix)
      //printf("Sinv first\n");
      //print_gsl_matrix(Sinv, current_n, current_n);

      // no longer the first iteration
      first = false;
    }
    else { // not first time
      if (throwOut > 0) { // change to "else if ( >1 ) if implement the shortcut
        /*
        S<-S[-indx,-indx]
        eps2.vec<-eps2.vec[-indx]
        diag(S)<-diag(S)+eps2.vec
        SRinv<-chol2inv(chol(S))
        # from S+R go back to S
        diag(S)<-diag(S)-eps2.vec
        */
        /* --------------------------------------------------
        have 3 sizes:
        n = original
        current_n = current size of vectors / matrices
        newSize = size that will be current n after throwing things out
        ---------------------------------------------------- */
        int newSize = current_n - throwOut;
        // we now have an empty box!!!
        if(newSize == 0) {
          break;
        }
        assert(newSize > 0);
        gsl_vector *d_temp = gsl_vector_alloc(newSize);
        gsl_matrix *s_temp = gsl_matrix_alloc(newSize, newSize);

        //printf("d before: ");
        //for(int vec=0; vec<current_n; vec++) {
          //printf(" %f", gsl_vector_get(d, vec));
        //}
        //printf("\n");
        // loop over the original length
        //for(int original=0; original<n; original++) {}

        int li = 0;
        for(int i=0; i<n; i++) { // this needs to loop over current_n
          if(flags[i] == 2) { // stations flagged for removal
            // printf("Removing column - li: %i, i: %i \n", li, i);
            flags[i] = 1; // newly removed station now set to 1
          }
          else if(flags[i] == 0){ // all rows and columns that we want to keep
            // update d
            gsl_vector_set(d_temp,li, d_global[i]);
            int lj = 0;
            for(int j=0; j<n; j++) { // this needs to loop over current_n
              if(flags[j] == 0) {
                // update S
                gsl_matrix_set(s_temp, li, lj, S_global[i][j]);
                lj++;
              }
              else if (flags[j] == 1){
                //printf("Removing row - lj: %i, j: %i \n", lj, j);
              }
            }
            assert(lj == newSize);
            li++;
          }
        }
        assert(li == newSize);
        current_n = newSize;
        gsl_vector_free(d);
        d = d_temp;
        gsl_matrix_free(S);
        S = s_temp;
        //printf("d after: ");
        //for(int vec=0; vec<current_n; vec++) {
          //printf(" %f", gsl_vector_get(d, vec));
        //}
        //printf("\n");
        assert(d->size == current_n);
        assert(S->size1 == current_n);
        assert(S->size2 == current_n);

        // chech size has not changed too much, if it has then recompute VP
        float percentageSizeChange = (float)(sizeWhenProfileCalculated - current_n) / sizeWhenProfileCalculated;
        if(0 && percentageSizeChange > 0.1) {
          printf("size of box has changed significantly (recalculate vp): %f \n", percentageSizeChange);
          // gsl_vector *vp_input, double *vp
          int status = compute_vertical_profile(currentBox, meanT, gamma, a, exact_p10, exact_p90, nminprof[0], dzmin[0], vp);
          printf("status optimizer: %d\n", status);
          sizeWhenProfileCalculated = current_n;
        }

        // weight the diagonal again
        li = 0;
        for(int i=0; i<n; i++) {
          if(flags[i] == 0){
             double value = gsl_matrix_get(S,li,li) + eps2[i];
             gsl_matrix_set(S,li,li,value);
             li++;
          }
        }
        gsl_matrix_free(Sinv); // free the old Sinv
        Sinv = gsl_matrix_alloc(current_n,current_n);
        // invert the matrix again
        //gsl_matrix_memcpy (Sinv, S);
        clock_t start = clock(), diff;
        Sinv = inverse_matrix(S);
        // if use this, then make Sinv a copy of S
        //gsl_linalg_cholesky_decomp1(Sinv);
        //gsl_linalg_cholesky_invert(Sinv); // in place
        diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        // printf("Time taken to invert matrix %d seconds %d milliseconds \n", msec/1000, msec%1000);

        // UN-weight the diagonal of S
        li = 0;
        for(int i=0; i<n; i++) {
          if(flags[i] == 0){
             double value = gsl_matrix_get(S,li,li) - eps2[i];
             gsl_matrix_set(S,li,li,value);
             li++;
          }
        }

      }
      //printf("S\n");
      //print_gsl_matrix(S, current_n, current_n); //(int rows, int columns, gsl_matrix *matrix)
      //printf("Sinv\n");
      //print_gsl_matrix(Sinv, current_n, current_n);

    } // end else
    // printf("Current n (end of matrix and vector updates) %i \n", current_n);
    // printf("d size %lu \n", d->size);
    // printf("S size1 %lu size2 %lu \n", S->size1, S->size2);
    assert(d->size == current_n);
    assert(S->size1 == current_n);
    assert(S->size2 == current_n);
    assert(Sinv->size1 == current_n);
    assert(Sinv->size2 == current_n);
    assert(current_n > 0); // Should hopefully not throw out all the stations...

    gsl_vector *Zinv, *Sinv_d, *ares_temp, *ares;
    Zinv = gsl_vector_alloc(current_n);
    Sinv_d = gsl_vector_alloc(current_n);
    ares_temp = gsl_vector_alloc(current_n);
    ares = gsl_vector_alloc(current_n);

    // (SRinv.d<-crossprod(SRinv,d[sel]))
    // this function does not appear to work properly!!!
    //gsl_blas_dgemv(CblasNoTrans, 1, Sinv, d, 1, Sinv_d); // crossprod Sinv & d to create Sinv_d
    //gsl_blas_dgemv(CblasNoTrans, 1, S, Sinv_d, 1, ares_temp); // crossprod S and Sinv_d to create ares_temp
    int li = 0;
    for(int i=0; i<n; i++) {
      double acc = 0;
      int lj = 0;
      if(flags[i] == 0) {
        for(int j=0; j<n; j++) {
          if(flags[j] == 0) {
            acc += gsl_matrix_get(Sinv,li,lj)*gsl_vector_get(d,lj);
            lj++;
          }
        }
        gsl_vector_set(Sinv_d, li, acc);
        li++;
      }
    }
    li = 0;
    for(int i=0; i<n; i++) {
      double acc = 0;
      int lj = 0;
      if(flags[i] == 0) {
        for(int j=0; j<n; j++) {
          if(flags[j] == 0) {
            acc += gsl_matrix_get(S,li,lj)*gsl_vector_get(Sinv_d,lj);
            lj++;
          }
        }
        gsl_vector_set(ares_temp, li, acc);
        li++;
      }
    }
    li = 0;
    for(int i=0; i<n; i++) {
      if(flags[i] == 0) {
        gsl_vector_set(Zinv,li,(1/gsl_matrix_get(Sinv,li,li))); //Zinv<-1/diag(SRinv)
        gsl_vector_set(ares,li,(gsl_vector_get(ares_temp,li)-gsl_vector_get(d,li))); // ares<-crossprod(S,SRinv.d)-d[sel]
        li++;
      }
    }
    gsl_vector_free(ares_temp);
    //printf("Zinv: ");
    //print_gsl_vector(Zinv,current_n);
    //printf("Sinv_d: ");
    //print_gsl_vector(Sinv_d,current_n);

    // cvres<--Zinv*SRinv.d
    gsl_vector *cvres;
    cvres = gsl_vector_alloc(current_n);
    for(int li=0; li<current_n; li++) {
       gsl_vector_set(cvres,li,-1*gsl_vector_get(Zinv,li));
    }
    gsl_vector_mul(cvres,Sinv_d); // multiplies -Zinv(initial cvres) by Sinv_d (result stored in cvres)
    // print_gsl_vector(cvres, current_n);

    // sig2o<-mean(d[sel]*(-ares))
    gsl_vector *sig2o_temp, *negAres_temp;
    sig2o_temp = gsl_vector_alloc(current_n);
    negAres_temp = gsl_vector_alloc(current_n);
    for(int li=0; li<current_n; li++) {
       gsl_vector_set(negAres_temp,li,-1*gsl_vector_get(ares,li));
    }
    gsl_vector_memcpy(sig2o_temp,d); // copies d into sig2o_temp
    gsl_vector_mul(sig2o_temp,negAres_temp); // multiplies d by -ares

    double sig2o = 0;
    for(int li=0; li<current_n; li++) {
        sig2o = sig2o + gsl_vector_get(sig2o_temp,li);
    }
    sig2o = sig2o/current_n;
    // printf("sig2o: %f\n", sig2o);

    gsl_vector_free(sig2o_temp);
    gsl_vector_free(negAres_temp);
    if(sig2o < 0.01) {
      sig2o = 0.01;
    }

    // pog[sel]<-(ares*cvres)/sig2o
    gsl_vector *pog, *pog_temp;
    pog = gsl_vector_alloc(current_n);
    pog_temp = gsl_vector_alloc(current_n);
    gsl_vector_memcpy(pog_temp,ares); // copies ares into pog_temp
    gsl_vector_mul(pog_temp,cvres); // multiplies ares by cvres

    for(int li=0; li<current_n; li++) {
      gsl_vector_set(pog,li,(gsl_vector_get(pog_temp,li)/sig2o));
    }
    gsl_vector_free(pog_temp);

    // figure out if we should flag a station
    throwOut = 0; // reset this number
    li = 0;
    for(int i=0; i<n; i++) {
      if(flags[i] == 0) {
         assert(sig2o > 0);
         rep_out[i] = gsl_vector_get(d, li)*gsl_vector_get(ares, li) * -1 /sig2o;
         sct_out[i] = gsl_vector_get(pog, li);
        // does it fail the test
        double curr_cvres = gsl_vector_get(cvres,li);
        if((curr_cvres < 0 && gsl_vector_get(pog,li) > t2pos[i]) ||
          (curr_cvres >= 0 && gsl_vector_get(pog,li) > t2neg[i])) {
          // printf("throw out this piece of data: %f cvres=%f sct=%f rep=%f\n", t[i], curr_cvres, gsl_vector_get(pog,li), rep_out[i]);
          throwOut = throwOut + 1;
          flags[i] = 2; // temporarily set to 2 so we know its a newly flagged station
        }
        li++;
      }
    }
    // printf("throw out: %i \n",throwOut);

    // FREE ALL THE MEMORY !!! (but not S or d)
    gsl_vector_free(Zinv);
    gsl_vector_free(Sinv_d);
    gsl_vector_free(ares);
    gsl_vector_free(cvres);
    gsl_vector_free(pog);

    if(throwOut == 0) { // no stations to remove, done SCT iterations
      break;
    }
  } // end of while SCT loop
  // FREE ALL THE MEMORY !!!
  free(vp);
  for(int i=0; i<n; i++) {
     free(disth[i]);
     free(distz[i]);
  }
  free(disth);
  free(distz);
  gsl_vector_free(d);
  gsl_matrix_free(S);
  gsl_matrix_free(Sinv);
  free(d_global);
}


struct Box merge_boxes(struct Box box1, struct Box box2) {

  struct Box mergedBox;
  mergedBox.n = (box1.n+box2.n); // set initial n
  mergedBox.x = malloc(sizeof(double) * mergedBox.n);
  mergedBox.y = malloc(sizeof(double) * mergedBox.n);
  mergedBox.z = malloc(sizeof(double) * mergedBox.n);
  mergedBox.t = malloc(sizeof(double) * mergedBox.n);
  mergedBox.i = malloc(sizeof(double) * mergedBox.n);

  for(int i=0; i<box1.n; i++) {
    mergedBox.x[i] = box1.x[i];
    mergedBox.y[i] = box1.y[i];
    mergedBox.z[i] = box1.z[i];
    mergedBox.t[i] = box1.t[i];
    mergedBox.i[i] = box1.i[i];
  }
  for(int j=box1.n; j<mergedBox.n; j++) {
    mergedBox.x[j] = box2.x[j-box1.n];
    mergedBox.y[j] = box2.y[j-box1.n];
    mergedBox.z[j] = box2.z[j-box1.n];
    mergedBox.t[j] = box2.t[j-box1.n];
    mergedBox.i[j] = box2.i[j-box1.n];
  }

  return mergedBox;
}

// recursive function that will keep splitting boxes until they are the right size
// returns a box at the "leaf"
void split_box(int maxNumStationsInBox, int minNumStationsInBox, struct Box inputBox, struct BoxList* box_list) {

  struct Box * boxes;
  // allocate memory, currently for 4 boxes
  boxes = malloc(sizeof(struct Box) * 4);
  for(int i=0; i<4; i++) {
    //boxes[i].n = malloc(sizeof(int));
    boxes[i].n = 0; // set initial n
    boxes[i].x = malloc(sizeof(double) * inputBox.n);
    boxes[i].y = malloc(sizeof(double) * inputBox.n);
    boxes[i].z = malloc(sizeof(double) * inputBox.n);
    boxes[i].t = malloc(sizeof(double) * inputBox.n);
    boxes[i].i = malloc(sizeof(double) * inputBox.n);
  }
  double maxX = max(inputBox.x,inputBox.n);
  double maxY = max(inputBox.y,inputBox.n);
  double minX = min(inputBox.x,inputBox.n);
  double minY = min(inputBox.y,inputBox.n);
  // halfway between min and max
  double halfwayX = minX + abs(abs(maxX)-abs(minX))/2;
  double halfwayY = minY + abs(abs(maxY)-abs(minY))/2;
  //printf("halfway x: %f y: %f \n", halfwayX, halfwayY);

  // new boxes
  for(int i=0; i<inputBox.n; i++) {
    // (0,0)
    if(inputBox.x[i] < halfwayX && inputBox.y[i] < halfwayY) {
      boxes[0].x[boxes[0].n] = inputBox.x[i];
      boxes[0].y[boxes[0].n] = inputBox.y[i];
      boxes[0].z[boxes[0].n] = inputBox.z[i];
      boxes[0].t[boxes[0].n] = inputBox.t[i];
      boxes[0].i[boxes[0].n] = inputBox.i[i];
      boxes[0].n++;
    }
    // (0,1)
    if(inputBox.x[i] >= halfwayX && inputBox.y[i] < halfwayY) {
      boxes[1].x[boxes[1].n] = inputBox.x[i];
      boxes[1].y[boxes[1].n] = inputBox.y[i];
      boxes[1].z[boxes[1].n] = inputBox.z[i];
      boxes[1].t[boxes[1].n] = inputBox.t[i];
      boxes[1].i[boxes[1].n] = inputBox.i[i];
      boxes[1].n++;
    }
    // (1,0)
    if(inputBox.x[i] < halfwayX && inputBox.y[i] >= halfwayY) {
      boxes[2].x[boxes[2].n] = inputBox.x[i];
      boxes[2].y[boxes[2].n] = inputBox.y[i];
      boxes[2].z[boxes[2].n] = inputBox.z[i];
      boxes[2].t[boxes[2].n] = inputBox.t[i];
      boxes[2].i[boxes[2].n] = inputBox.i[i];
      boxes[2].n++;
    }
    // (1,1)
    if(inputBox.x[i] >= halfwayX && inputBox.y[i] >= halfwayY) {
      boxes[3].x[boxes[3].n] = inputBox.x[i];
      boxes[3].y[boxes[3].n] = inputBox.y[i];
      boxes[3].z[boxes[3].n] = inputBox.z[i];
      boxes[3].t[boxes[3].n] = inputBox.t[i];
      boxes[3].i[boxes[3].n] = inputBox.i[i];
      boxes[3].n++;
    }
  }
  // printf("4 way split - 0: %i 1: %i 2: %i 3: %i \n", boxes[0].n, boxes[1].n, boxes[2].n, boxes[3].n);
  int numBoxes = 4;

  // what kind of aspect ratio does input box have
  double diffX = abs(maxX - minX);
  double diffY = abs(maxY - minY);

  // first check which way it makes more sense to merge
  // ok to merge the way that keeps the boxes squarer?
  if(diffY < diffX) {
    // wide boxes, so preferable to merge 0,2 + 1,3 (vertically)
    int n1 = boxes[0].n + boxes[2].n;
    int n2 = boxes[1].n + boxes[3].n;
    if(n1 > minNumStationsInBox && n2 > minNumStationsInBox) {
      // merge vertically
      // printf("best to merge vertically \n");
      boxes[0] = merge_boxes(boxes[0],boxes[2]);
      boxes[1] = merge_boxes(boxes[1],boxes[3]);
    }
    else {
      // choose the way that keeps the best symmetry
      int n3 = boxes[0].n + boxes[1].n;
      int n4 = boxes[2].n + boxes[3].n;
      int hori = abs(n3-n4);
      int vert = abs(n1-n2);
      if(vert > hori) {
        // merge horizontally
        // printf("had to merge horizontally \n");
        boxes[0] = merge_boxes(boxes[0],boxes[1]);
        boxes[1] = merge_boxes(boxes[2],boxes[3]);
      }
      else {
        // merge vertically
        boxes[0] = merge_boxes(boxes[0],boxes[2]);
        boxes[1] = merge_boxes(boxes[1],boxes[3]);
      }
    }
  }
  else { // diffY > diffX
    // tall boxes, so preferable to merge 0,1 + 2,3 (horizontally)
    int n1 = boxes[0].n + boxes[1].n;
    int n2 = boxes[2].n + boxes[3].n;
    if(n1 > minNumStationsInBox && n2 > minNumStationsInBox) {
      // merge horizontally
      // printf("best to merge horizontally \n");
      boxes[0] = merge_boxes(boxes[0],boxes[1]);
      boxes[1] = merge_boxes(boxes[2],boxes[3]);
    }
    else {
      // choose the way that keeps the best symmetry
      int n3 = boxes[0].n + boxes[2].n;
      int n4 = boxes[1].n + boxes[3].n;
      int vert = abs(n3-n4);
      int hori = abs(n1-n2);
      if(hori > vert) {
        // merge vertically
        // printf("had to merge vertically \n");
        boxes[0] = merge_boxes(boxes[0],boxes[2]);
        boxes[1] = merge_boxes(boxes[1],boxes[3]);
      }
      else {
        // merge horizontally
        boxes[0] = merge_boxes(boxes[0],boxes[1]);
        boxes[1] = merge_boxes(boxes[2],boxes[3]);
      }
    }
  }
  // don't need these boxes anymore
  for(int i=2; i<numBoxes; i++) {
    free(boxes[i].x);
    free(boxes[i].y);
    free(boxes[i].z);
    free(boxes[i].t);
    free(boxes[i].i);
  }
  // printf("2 way split - 0: %i 1: %i \n", boxes[0].n, boxes[1].n);
  numBoxes = 2;

  // loop over the boxes
  for(int i=0; i<numBoxes; i++) {
    int n_temp = boxes[i].n;
    if(n_temp > maxNumStationsInBox) {
      // printf("box still too big %i (being further recursively split)\n", n_temp);
      // split the box further
      int before = box_list->n;
      split_box(maxNumStationsInBox, minNumStationsInBox, boxes[i], box_list);
      int after = box_list->n;
    }
    else {
      // printf("box size: %i, being returned \n", n_temp);
      int current_n = box_list->n;
      // add to list of boxes (finalBoxes)
      box_list->boxes[current_n].n = boxes[i].n;
      box_list->boxes[current_n].x = boxes[i].x;
      box_list->boxes[current_n].y = boxes[i].y;
      box_list->boxes[current_n].z = boxes[i].z;
      box_list->boxes[current_n].t = boxes[i].t;
      box_list->boxes[current_n].i = boxes[i].i;
      // increment the number of boxes
      box_list->n++;
    }
  }
}

struct BoxList control_box_division(int maxNumStationsInBox, int minNumStationsInBox, struct Box inputBox) {

  // check is isn't already smaller than the max
  if(inputBox.n < maxNumStationsInBox) {
    struct BoxList box_list;
    // just return
    box_list.n = 1;
    box_list.boxes = malloc(sizeof(struct Box));
    box_list.boxes[0].n = inputBox.n;
    box_list.boxes[0].x = inputBox.x;
    box_list.boxes[0].y = inputBox.y;
    box_list.boxes[0].z = inputBox.z;
    box_list.boxes[0].t = inputBox.t;
    box_list.boxes[0].i = inputBox.i;

    return box_list;
  }

  struct BoxList box_list;
  int maxNumBoxes = floor(inputBox.n/minNumStationsInBox);
  // printf("allocating memory for potential max number of boxes: %i \n",  maxNumBoxes);
  box_list.boxes = malloc(sizeof(struct Box) * maxNumBoxes);
  box_list.n = 0;

  // pass the outputs in by reference (so the recursive function can add to them)
  split_box(maxNumStationsInBox, minNumStationsInBox, inputBox, &box_list);

  // printf("total number of boxes: %i \n", box_list.n);
  return box_list;
}

//----------------------------------------------------------------------------//
// HELPER FUNCTIONS
double compute_quantile(double quantile, double *array, int sizeArray)
{
  double* array_copy = malloc(sizeof(double) * sizeArray);
  for(int i = 0; i < sizeArray; i++)
     array_copy[i] = array[i];
  double exact_q;
  gsl_sort(array_copy, 1, sizeArray);
  exact_q = gsl_stats_quantile_from_sorted_data(array_copy, 1, sizeArray, quantile);
  free(array_copy);
  return exact_q;
}

double mean(const double *array, int sizeArray)
{
  double sum = 0;
  for(int i=0; i<sizeArray; i++) {
    sum = sum + array[i];
  }
  double mean = sum/sizeArray;
  return mean;
}

double max(const double *array, int sizeArray)
{
  double max = array[0];
  for(int i=0; i<sizeArray; i++) {
    if(array[i] > max) {
        max = array[i];
    }
  }
  return max;
}

double min(const double *array, int sizeArray)
{
  double min = array[0];
  for(int i=0; i<sizeArray; i++) {
    if(array[i] < min) {
        min = array[i];
    }
  }
  return min;
}

void print_vector(double *vector, int size) {
  for (int s=0; s<size; s++)
  {
    printf("%.2f ", vector[s]);
  }
  printf("\n");
}

void print_gsl_vector(gsl_vector *vector, int size) {
  for (int s=0; s<size; s++)
  {
    printf("%.2f ", gsl_vector_get(vector,s));
  }
  printf("\n");
}

void print_matrix(double **matrix, int rows, int columns) {
  for (int r=0; r<rows; r++)
  {
      for(int c=0; c<columns; c++)
      {
         printf("%.2f ", matrix[r][c]);
      }
      printf("\n");
   }
}

void print_gsl_matrix(gsl_matrix *matrix, int rows, int columns) {
  for (int r=0; r<rows; r++)
  {
      for(int c=0; c<columns; c++)
      {
         printf("%.2f ", gsl_matrix_get(matrix,r,c));
      }
      printf("\n");
   }
}

void print_sub_gsl_matrix(gsl_matrix *matrix, int start, int stop) {
  for (int r=start; r<=stop; r++)
  {
      for(int c=start; c<=stop; c++)
      {
         printf("%.2f ", gsl_matrix_get(matrix,r,c));
      }
      printf("\n");
   }
}

gsl_matrix* inverse_matrix(const gsl_matrix *matrix) {
  int s;
  gsl_matrix *return_matrix = gsl_matrix_alloc(matrix->size1,matrix->size2);
  gsl_matrix *temp_matrix = gsl_matrix_alloc(matrix->size1,matrix->size2);
  gsl_matrix_memcpy(temp_matrix, matrix);
  gsl_permutation * p = gsl_permutation_alloc (matrix->size1);

  gsl_linalg_LU_decomp (temp_matrix, p, &s);
  gsl_linalg_LU_invert (temp_matrix, p, return_matrix);
  gsl_matrix_free(temp_matrix);
  gsl_permutation_free (p);
  return return_matrix;
}

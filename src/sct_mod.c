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

void spatial_consistency_test_mod(int *N, int *obs_to_test, double *x, double *y, double *z, double *elevs, double *t, int *nminprof, double *dzmin, double *dhmin, double *dz, double *t2pos, double *t2neg, double *eps2, int *flags, double *sct_out, double *rep_out)
{
    // create box for simplicity
    int n = *N;
    int m = 0;
    for (int i = 0; i < n; i++)
    {
        if (obs_to_test[i] == 1)
            m++;
    }
    int * obs_to_test_indices = malloc(sizeof(int*) * m);
    int counter = 0;
    for(int i = 0; i < n; i++) {
        if(obs_to_test[i] == 1) {
            obs_to_test_indices[counter] = i;
            counter++;
        }
    }

    struct Box currentBox;
    currentBox.x = x;
    currentBox.y = y;
    currentBox.z = elevs;
    currentBox.t = t;
    // no i since no index into global array

    // fill with 0 to keep track of stations that are flagged
    for (int i = 0; i < n; i++)
    {
        flags[i] = 0;
    }

    /*
  Stuff for VP
  */
    double gamma = -0.0065;
    double a = 5.0;
    double meanT = mean(t, n);
    // allocate for output
    double *vp = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++)
    {
        vp[i] = -999;
    }
    double exact_p10 = compute_quantile(0.10, elevs, n);
    double exact_p90 = compute_quantile(0.90, elevs, n);

    // vector (double t0, double gamma, double a, double h0, double h1i)
    // Starting point for optimization
    // printf ("VP input vector set = t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n", meanT, gamma, a, exact_p10, exact_p90);

    // calculate background
    int status = compute_vertical_profile(&currentBox, meanT, gamma, a, exact_p10, exact_p90, nminprof[0], dzmin[0], vp);

    // printf("t0: %.4f gamma: %.4f a: %.4f h0: %.4f h1i: %.4f\n", meanT, gamma, a, exact_p10, exact_p90);

    // printf("status optimizer: %d\n", status);
    // now have temperature profile (vp)
    for (int i = 0; i < n; i++)
    {
        assert(vp[i] != -999);
    }
    // done calculating VP for the first time
    int sizeWhenProfileCalculated = n; // TODO: this number is either redundant or wrong???

    // distance matrices
    double **disth = malloc(sizeof(double *) * n);
    double **distz = malloc(sizeof(double *) * n);
    double *Dh = malloc(sizeof(double) * n);

    // no need to select j since already only have those for a particular box
    // outer product of the matrices
    for (int i = 0; i < n; i++)
    {
        disth[i] = malloc(sizeof(double) * n);
        distz[i] = malloc(sizeof(double) * n);
        double *Dh_vector = malloc(sizeof(double) * (n - 1)); // need to remove one since not considering the diagonal
        for (int j = 0; j < n; j++)
        {
            disth[i][j] = pow((pow((x[i] - x[j]), 2) + pow((y[i] - y[j]), 2) + pow((z[i] - z[j]), 2)), 0.5);
            distz[i][j] = abs(elevs[i] - elevs[j]);
            if (i != j)
            { // do not want to consider the diagonal
                if (i < j)
                {
                    Dh_vector[j - 1] = disth[i][j];
                }
                else if (i > j)
                {
                    Dh_vector[j] = disth[i][j];
                }
            }
        }
        Dh[i] = compute_quantile(0.10, Dh_vector, n - 1);
        free(Dh_vector);
    }

    // set to optimal Dh for the SCT
    // either Dhmin or the average 10-percentile of the set of distances between a
    // station and all the others
    double Dh_mean = mean(Dh, n);
    free(Dh);
    // printf("Dh: %f\n", Dh_mean);
    if (Dh_mean < dhmin[0])
    {
        Dh_mean = dhmin[0];
    }
    // printf("Dh_mean: %f\n", Dh_mean);

    // background error correlation matrix
    float Dz = 200; // good value (use this default)
    gsl_matrix *S, *Sinv;
    S = gsl_matrix_alloc(n, n);
    double **S_global = malloc(sizeof(double *) * n);
    for (int i = 0; i < n; i++)
    {
        S_global[i] = malloc(sizeof(double) * n);
    }
    Sinv = gsl_matrix_alloc(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double value = exp(-.5 * pow((disth[i][j] / Dh_mean), 2) - .5 * pow((distz[i][j] / (dz[0])), 2));
            S_global[i][j] = value;
            if (i == j)
            { // weight the diagonal?? (0.5 default)
                value = value + eps2[i];
            }
            gsl_matrix_set(S, i, j, value);
            //gsl_matrix_set(Sinv,i,j,value); // not yet inverted, but need a copy of S
        }
    }
    // printf("created the S matrix - size1 %lu size2 %lu \n", S->size1, S->size2);
    //print_gsl_matrix(S,n,n);

    // d<-topt-tb
    gsl_vector *d;
    d = gsl_vector_alloc(n);
    double *d_global = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++)
    {
        gsl_vector_set(d, i, (t[i] - vp[i])); // difference between actual temp and temperature from vertical profile
        d_global[i] = (t[i] - vp[i]);
    }

    /* ---------------------------------------------------
  Beginning of real SCT looping
  ------------------------------------------------------*/
    bool first = true;
    int current_n = n;
    int throwOut = 0;
    // loop for SCT
    // note that current_n should be used inside this loop!!!
    while (1)
    {
        if (first)
        {
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
            for (int i = 0; i < current_n; i++)
            {
                double value = gsl_matrix_get(S, i, i) - eps2[i];
                gsl_matrix_set(S, i, i, value);
            }
            //printf("S first\n");
            //print_gsl_matrix(S, current_n, current_n); //(int rows, int columns, gsl_matrix *matrix)
            //printf("Sinv first\n");
            //print_gsl_matrix(Sinv, current_n, current_n);

            // no longer the first iteration
            first = false;
        }
        else
        { // not first time
            if (throwOut > 0)
            { // change to "else if ( >1 ) if implement the shortcut
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
                if (newSize == 0)
                {
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
                for (int i = 0; i < n; i++)
                { // this needs to loop over current_n
                    if (flags[i] == 2)
                    { // stations flagged for removal
                        // printf("Removing column - li: %i, i: %i \n", li, i);
                        flags[i] = 1; // newly removed station now set to 1
                    }
                    else if (flags[i] == 0)
                    { // all rows and columns that we want to keep
                        // update d
                        gsl_vector_set(d_temp, li, d_global[i]);
                        int lj = 0;
                        for (int j = 0; j < n; j++)
                        { // this needs to loop over current_n
                            if (flags[j] == 0)
                            {
                                // update S
                                gsl_matrix_set(s_temp, li, lj, S_global[i][j]);
                                lj++;
                            }
                            else if (flags[j] == 1)
                            {
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
                if (0 && percentageSizeChange > 0.1)
                {
                    printf("size of box has changed significantly (recalculate vp): %f \n", percentageSizeChange);
                    // gsl_vector *vp_input, double *vp
                    int status = compute_vertical_profile(&currentBox, meanT, gamma, a, exact_p10, exact_p90, nminprof[0], dzmin[0], vp);
                    printf("status optimizer: %d\n", status);
                    sizeWhenProfileCalculated = current_n;
                }

                // weight the diagonal again
                li = 0;
                for (int i = 0; i < n; i++)
                {
                    if (flags[i] == 0)
                    {
                        double value = gsl_matrix_get(S, li, li) + eps2[i];
                        gsl_matrix_set(S, li, li, value);
                        li++;
                    }
                }
                gsl_matrix_free(Sinv); // free the old Sinv
                Sinv = gsl_matrix_alloc(current_n, current_n);
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
                for (int i = 0; i < n; i++)
                {
                    if (flags[i] == 0)
                    {
                        double value = gsl_matrix_get(S, li, li) - eps2[i];
                        gsl_matrix_set(S, li, li, value);
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
        for (int i = 0; i < n; i++)
        {
            double acc = 0;
            int lj = 0;
            if (flags[i] == 0)
            {
                for (int j = 0; j < n; j++)
                {
                    if (flags[j] == 0)
                    {
                        acc += gsl_matrix_get(Sinv, li, lj) * gsl_vector_get(d, lj);
                        lj++;
                    }
                }
                gsl_vector_set(Sinv_d, li, acc);
                li++;
            }
        }
        li = 0;
        for (int i = 0; i < n; i++)
        {
            double acc = 0;
            int lj = 0;
            if (flags[i] == 0)
            {
                for (int j = 0; j < n; j++)
                {
                    if (flags[j] == 0)
                    {
                        acc += gsl_matrix_get(S, li, lj) * gsl_vector_get(Sinv_d, lj);
                        lj++;
                    }
                }
                gsl_vector_set(ares_temp, li, acc);
                li++;
            }
        }
        li = 0;
        for (int i = 0; i < n; i++)
        {
            if (flags[i] == 0)
            {
                gsl_vector_set(Zinv, li, (1 / gsl_matrix_get(Sinv, li, li)));                      //Zinv<-1/diag(SRinv)
                gsl_vector_set(ares, li, (gsl_vector_get(ares_temp, li) - gsl_vector_get(d, li))); // ares<-crossprod(S,SRinv.d)-d[sel]
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
        for (int li = 0; li < current_n; li++)
        {
            gsl_vector_set(cvres, li, -1 * gsl_vector_get(Zinv, li));
        }
        gsl_vector_mul(cvres, Sinv_d); // multiplies -Zinv(initial cvres) by Sinv_d (result stored in cvres)
        // print_gsl_vector(cvres, current_n);

        // sig2o<-mean(d[sel]*(-ares))
        gsl_vector *sig2o_temp, *negAres_temp;
        sig2o_temp = gsl_vector_alloc(current_n);
        negAres_temp = gsl_vector_alloc(current_n);
        for (int li = 0; li < current_n; li++)
        {
            gsl_vector_set(negAres_temp, li, -1 * gsl_vector_get(ares, li));
        }
        gsl_vector_memcpy(sig2o_temp, d);         // copies d into sig2o_temp
        gsl_vector_mul(sig2o_temp, negAres_temp); // multiplies d by -ares

        double sig2o = 0;
        for (int li = 0; li < current_n; li++)
        {
            sig2o = sig2o + gsl_vector_get(sig2o_temp, li);
        }
        sig2o = sig2o / current_n;
        // printf("sig2o: %f\n", sig2o);

        gsl_vector_free(sig2o_temp);
        gsl_vector_free(negAres_temp);
        if (sig2o < 0.01)
        {
            sig2o = 0.01;
        }

        // pog[sel]<-(ares*cvres)/sig2o
        gsl_vector *pog, *pog_temp;
        pog = gsl_vector_alloc(current_n);
        pog_temp = gsl_vector_alloc(current_n);
        gsl_vector_memcpy(pog_temp, ares); // copies ares into pog_temp
        gsl_vector_mul(pog_temp, cvres);   // multiplies ares by cvres

        for (int li = 0; li < current_n; li++)
        {
            gsl_vector_set(pog, li, (gsl_vector_get(pog_temp, li) / sig2o));
        }
        gsl_vector_free(pog_temp);

        // figure out if we should flag a station
        throwOut = 0; // reset this number
        li = 0;
        for (int i = 0; i < n; i++)
        {
            if (flags[i] == 0)
            {
                assert(sig2o > 0);
                rep_out[i] = gsl_vector_get(d, li) * gsl_vector_get(ares, li) * -1 / sig2o;
                sct_out[i] = gsl_vector_get(pog, li);
                // does it fail the test
                double curr_cvres = gsl_vector_get(cvres, li);
                if ((curr_cvres < 0 && gsl_vector_get(pog, li) > t2pos[i]) ||
                    (curr_cvres >= 0 && gsl_vector_get(pog, li) > t2neg[i]))
                {
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

        if (throwOut == 0)
        { // no stations to remove, done SCT iterations
            break;
        }
    } // end of while SCT loop
    // FREE ALL THE MEMORY !!!
    free(vp);
    for (int i = 0; i < n; i++)
    {
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

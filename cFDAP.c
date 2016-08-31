/*******************************************/
/* cFDAP.c (c) 2015 Maxim Igaev, Osnabrück */
/*******************************************/

/* Compiling with gsl and blas
   cc/gcc cFDAP.c -o cFDAP -lgsl -lgslcblas -lm */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h> 
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_cdf.h>

/* Global variables */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Default values */
#define DEFAULT_DF 11.0 /* Diffusion constant */
#define DEFAULT_R 3.0 /* The half length of the AR */
#define DEFAULT_T_INI 0.0 /* Time range [T_INI, T_END] */
#define DEFAULT_T_END 112.0
#define DEFAULT_N 113 /* Number of points in the curves */
#define DEFAULT_KON_INIT 0.5 /* Starting value for kon */
#define DEFAULT_KOFF_INIT 0.5 /* Starting value for koff */
#define DEFAULT_XX_INIT 1.0 /* Starting value for xx = kon/koff */
#define DEFAULT_FLAG_WEIGHT 0 /* By default, fitting is unweighted */

/* MACROS */
#define NELEMS_1D(x) (sizeof(x)/sizeof((x)[0]))
#define NROW_2D(x) NELEMS_1D(x)
#define NCOL_2D(x) ((sizeof(x)/sizeof((x)[0][0]))/(sizeof(x)/sizeof((x)[0])))
#define NELEMS_2D(x) (sizeof(x)/sizeof((x)[0][0]))

/* STRUCTURES */
struct data {
    size_t n;
    double Df;
    double R;
    double * time;
    double * y;
    double * sigma;
    char * m;
    size_t p;
    size_t w_flag;
};

/* FUNCTION DECLARATIONS */
/* Functions with the '_x', '_kon' or '_koff' indeces are
 * derivatives of those functions with respect to '_x',
 * '_kon' and '_koff', respectively */
double complex fullModel(double complex s, double kon,
                         double koff, double Df, double R);
double complex fullModel_kon(double complex s, double kon,
                             double koff, double Df, double R);
double complex fullModel_koff(double complex s, double kon,
                              double koff, double Df, double R);
double complex effectiveDiffusion(double complex s, double x,
                                  double Df, double R);
double complex effectiveDiffusion_x(double complex s, double x,
                                    double Df, double R);
double complex reactionDominantPure(double complex s, double kon,
                                    double koff, double Df, double R);
double complex reactionDominantPure_kon(double complex s, double kon,
                                        double koff, double Df, double R);
double complex reactionDominantPure_koff(double complex s, double kon,
                                         double koff, double Df, double R);
double complex hybridModel(double complex s, double kon,
                           double koff, double Df, double R);
double complex hybridModel_kon(double complex s, double kon,
                               double koff, double Df, double R);
double complex hybridModel_koff(double complex s, double kon,
                                double koff, double Df, double R);
double invlap_1(double t, double xx, double Df,
                double R, char *m, int functionOrDerivative);
double invlap_2(double t, double kon, double koff,
                double Df, double R, char *m, int functionOrDerivative);
int model_f(const gsl_vector * x, void *data, gsl_vector * f);
int model_df(const gsl_vector * x, void *data, gsl_matrix * J);
int model_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J);
void print_state (size_t iter, gsl_multifit_fdfsolver * s, size_t p);
void bad_input(void);

/* FUNCTIONS */

/* Laplace images of FDAP(t). Note: These function must take
   a complex argument s and return a complex value */
double complex
fullModel(double complex s, double kon, double koff, double Df, double R) {
    return (1.0/(1.0 + kon/koff)) * (1.0 + kon/(s + koff)) * (1.0/s - 1.0/2.0/s/csqrt(R*R*s*(1.0 + kon/(s + koff))/Df) * (1.0 - cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)))) + (kon/koff/(1.0 + kon/koff))/(s + koff);
}

double complex
fullModel_kon(double complex s, double kon, double koff, double Df, double R) {
    return -koff*(-kon - koff + koff*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + kon*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + 2.0*koff*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + 2.0*kon*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) - 2.0*s + 2.0*s*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)))/4.0/s/(s + koff)/(kon + koff)/(kon + koff)/csqrt(R*R*s*(1.0 + kon/(s + koff))/Df);
}

double complex
fullModel_koff(double complex s, double kon, double koff, double Df, double R) {
    return kon*(-cpow(koff, 2.0) - 2.0*cpow(s, 2.0) - 4.0*s*koff - kon*koff + kon*koff*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + 2.0*koff*koff*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + 2.0*s*s*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + 4.0*s*koff*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + koff*koff*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + 2.0*s*kon*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) + 2.0*kon*koff*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)*cexp(-2.0*csqrt(R*R*s*(1.0 + kon/(s + koff))/Df)) - 2*s*kon)/4.0/s/(s + koff)/(s + koff)/(kon + koff)/(kon + koff)/csqrt(R*R*s*(1.0+kon/(s + koff))/Df);
}

double complex
effectiveDiffusion(double complex s, double xx, double Df, double R) {
    /* xx = kon/koff */
    return 1.0/s - 1.0/2.0/s/csqrt(R*R*s*(1.0 + xx)/Df)*(1.0 - cexp(-2.0*csqrt(R*R*s*(1.0 + xx)/Df)));
}

double complex
effectiveDiffusion_xx(double complex s, double xx, double Df, double R) {
    /* xx = kon/koff */
    return (1.0 - (1.0 + 2.0*csqrt(s*R*R/Df*(1.0 + xx)))*cexp(-2.0*csqrt(s*R*R/Df*(1.0 + xx))))/(4.0*s*(1.0 + xx)*csqrt(s*R*R/Df*(1.0 + xx)));
}

double complex
reactionDominantPure(double complex s, double kon, double koff, double Df, double R) {
    return koff/(kon + koff)*(1.0/s - 1.0/2.0/s/csqrt(R*R*s/Df)*(1.0 - cexp(-2.0*csqrt(R*R*s/Df)))) + kon/(kon + koff)/(s + koff);
}

double complex
reactionDominantPure_kon(double complex s, double kon, double koff, double Df, double R) {
    return (1.0 - kon/(kon + koff))/(kon + koff)/(s + koff) - koff*(1.0/s - 1.0/2.0/s/csqrt(R*R*s/Df)*(1.0 - cexp(-2.0*csqrt(R*R*s/Df))))/cpow(kon + koff, 2.0);
}

double complex
reactionDominantPure_koff(double complex s, double kon, double koff, double Df, double R) {
    return -kon/(kon + koff)/(s + koff)*(1.0/(kon + koff) + 1.0/(s + koff)) + (1.0/s - 1.0/2.0/s/csqrt(R*R*s/Df)*(1.0 - cexp(-2.0*csqrt(R*R*s/Df))))/(kon + koff)*(1.0 - koff/(kon + koff));
}

double complex
hybridModel(double complex s, double kon, double koff, double Df, double R) {
    return (koff/(s + koff))*(1.0/s - 1.0/2.0/s/csqrt(R*R*kon*s/Df/(s + koff))*(1.0 - cexp(-2.0*csqrt(R*R*kon*s/Df/(s + koff))))) + 1.0/(s + koff);
}

double complex
hybridModel_kon(double complex s, double kon, double koff, double Df, double R) {
    return koff/s/(s + koff)*(s*R*R/Df*(1.0 - cexp(-2.0*csqrt(R*R*s/Df*kon/(s + koff))))/4.0/(s + koff)/cpow(R*R*s/Df*kon/(s + koff), 3.0/2.0) - cexp(-2.0*csqrt(R*R*s/Df*kon/(s + koff)))/2.0/kon);
}

double complex
hybridModel_koff(double complex s, double kon, double koff, double Df, double R) {
    return -cexp(-2.0*csqrt(R*R*s/Df*kon/(s + koff)))/(4.0*s*cpow(s + koff, 3)*cpow(R*R*s/Df*kon/(s + koff), 3.0/2.0))*(R*R*s/Df*kon*(koff + 2.0*s)*cexp(2.0*csqrt(R*R*s/Df*kon/(s + koff))) - 2.0*koff*(s + koff)*cpow(R*R*s/Df*kon/(s + koff), 3.0/2.0) - R*R*s/Df*kon*koff - 2.0*kon*R*R*cpow(s, 2.0)/Df);
}

/* Function invlap(t, kon, koff) numerically inverts a Laplace
   image function F(s) into f(t) using the Fast Fourier Transform
   (FFT) algorithm for a specific time moment "t", an upper
   frequency limit "omega", a real parameter "sigma" and the
   number of integration intervals "n_int".
   
   Recommended values: omega > 100, n_int = 50*omega
   Default values:     omega = 200, n_int = 10000
   
   Sigma is a real number which must be a little bit bigger than
   the real part of the rigthmost pole of the function F(s).
   For example, if the rightmost pole is s = 2.0, then sigma could
   be equal to 2.05. This is done to keep all poles at the left of
   the integration area on the complex plane.
  
   Creator:  Fausto Arinos de Almeida Barbuto (Calgary, Canada)
   Date: May 18, 2002
   E-mail: fausto_barbuto@yahoo.ca
  
   Algorithm:
   Huddleston, T. and Byrne, P: "Numerical Inversion of
   Laplace Transforms", University of South Alabama, April
   1999 (found at http://www.eng.usouthal.edu/huddleston/
   SoftwareSupport/Download/Inversion99.doc)
  
   Modified and translated into C code by Maxim Igaev, 2015 */
double
invlap_1(double t, double xx, double Df, double R, char *m, int functionOrDerivative) {
    /* defining constants */
    int i, n_int = 10000;
    double omega = 200.0, sig = 0.05, delta = omega/((double) n_int);
    double sum = 0.0, wi = 0.0, wf, fi, ff;
    double complex witi, wfti; 

    /* loading one of the laplace image functions */
    double complex (*laplace_fun)();
    if (strcmp(m, "effectiveDiffusion") == 0) {
        if (functionOrDerivative == 0) {
            laplace_fun = &effectiveDiffusion;
        }
        else if (functionOrDerivative == 1) {
            laplace_fun = &effectiveDiffusion_xx;
        }
        else {
            fprintf(stderr, "ERROR: in 'invlap_1': functionOrDerivative takes only values 0 or 1 for the model and derivative, respectively.");
        }
    }

    for(i = 0; i < n_int; i++) {
        witi = 0.0 + (wi*t)*I;

        wf = wi + delta;
        wfti = 0.0 + (wf*t)*I;

        fi = creal(cexp(witi)*laplace_fun(sig + wi*I, xx, Df, R));
        ff = creal(cexp(wfti)*laplace_fun(sig + wf*I, xx, Df, R));
        sum += 0.5*(wf - wi)*(fi + ff);
        wi = wf;
    }

    return creal(sum*cexp(sig*t)/M_PI);
}

double
invlap_2(double t, double kon, double koff, double Df, double R, char *m, int functionOrDerivative) {
    /* defining constants */
    int i, n_int = 10000;
    double omega = 200.0, sig = 0.05, delta = omega/((double) n_int);
    double sum = 0.0, wi = 0.0, wf, fi, ff;
    double complex witi, wfti; 

    /* loading one of the laplace image functions */
    double complex (*laplace_fun)();
    if (strcmp(m, "fullModel") == 0) {
        if (functionOrDerivative == 0) {
            laplace_fun = &fullModel;
        }
        else if (functionOrDerivative == 1) {
            laplace_fun = &fullModel_kon;
        }
        else if (functionOrDerivative == 2) {
            laplace_fun = &fullModel_koff;
        }
        else {
            fprintf(stderr, "ERROR: in 'invlap_2': functionOrDerivative takes only values 0, 1 or 2 for the model and derivatives, respectively.");
        }
    }
    else if (strcmp(m, "hybridModel") == 0) {
        if (functionOrDerivative == 0) {
            laplace_fun = &hybridModel;
        }
        else if (functionOrDerivative == 1) {
            laplace_fun = &hybridModel_kon;
        }
        else if (functionOrDerivative == 2) {
            laplace_fun = &hybridModel_koff;
        }
        else {
            fprintf(stderr, "ERROR: in 'invlap_2': functionOrDerivative takes only values 0, 1 or 2 for the model and derivatives, respectively.");
        }
    }
    else if (strcmp(m, "reactionDominantPure") == 0) {
        if (functionOrDerivative == 0) {
            laplace_fun = &reactionDominantPure;
        }
        else if (functionOrDerivative == 1) {
            laplace_fun = &reactionDominantPure_kon;
        }
        else if (functionOrDerivative == 2) {
            laplace_fun = &reactionDominantPure_koff;
        }
        else {
            fprintf(stderr, "ERROR: in 'invlap_2': functionOrDerivative takes only values 0, 1 or 2 for the model and derivatives, respectively.");
        }
    }

    for(i = 0; i < n_int; i++) {
        witi = 0.0 + (wi*t)*I;

        wf = wi + delta;
        wfti = 0.0 + (wf*t)*I;

        fi = creal(cexp(witi)*laplace_fun(sig + wi*I, kon, koff, Df, R));
        ff = creal(cexp(wfti)*laplace_fun(sig + wf*I, kon, koff, Df, R));
        sum += 0.5*(wf - wi)*(fi + ff);
        wi = wf;
    }

    return creal(sum*cexp(sig*t)/M_PI);
}

int
model_f (const gsl_vector * x, void *data, 
        gsl_vector * f) {
    size_t n = ((struct data *)data)->n;
    double Df = ((struct data *)data)->Df;
    double R = ((struct data *)data)->R;
    double *time = ((struct data *)data)->time;
    double *y = ((struct data *)data)->y;
    double *sigma = ((struct data *) data)->sigma;
    char *m = ((struct data *) data)->m;
    size_t p = ((struct data *)data)->p;
    size_t w_flag = ((struct data *)data)->w_flag;

    size_t i;
    double xx, kon, koff;

    if (p == 1) {
        xx = gsl_vector_get (x, 0);

        /* Importing one of the model functions */    
        double (*inverted_fun)();
        inverted_fun = &invlap_1;

        for (i = 0; i < n; i++) {
            /* Model Yi = A * exp(-lambda * i) + b */
            //double Yi = A * exp (-lambda * t) + b;
            
            double Yi = inverted_fun(time[i], xx, Df, R, m, 0);
            if (w_flag == 0) {
                gsl_vector_set (f, i, (Yi - y[i]));
            }
            else if (w_flag == 1) {
                gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
            }
            else {
                fprintf(stderr, "ERROR: in 'model_f': Parameter w_flag is neither 0 nor 1.\n");
                exit(1);
            }
        }
    }
    else if (p == 2) {
        kon = gsl_vector_get (x, 0);
        koff = gsl_vector_get (x, 1);

        /* Importing one of the model functions */    
        double (*inverted_fun)();
        inverted_fun = &invlap_2;

        for (i = 0; i < n; i++) {
            /* Model Yi = A * exp(-lambda * i) + b */
            //double Yi = A * exp (-lambda * t) + b;

            double Yi = inverted_fun(time[i], kon, koff, Df, R, m, 0);
            if (w_flag == 0) {
                gsl_vector_set (f, i, (Yi - y[i]));
            }
            else if (w_flag == 1) {
                gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
            }
            else {
                fprintf(stderr, "ERROR: in 'model_f': Parameter w_flag is neither 0 nor 1.\n");
                exit(1);
            }
        }
    }
    else {
        fprintf(stderr, "ERROR: in 'model_f': Parameter p is neither 1 nor 2.\n");
        exit(1);
    }

    return GSL_SUCCESS;
}

int
model_df(const gsl_vector * x, void *data,
        gsl_matrix * J) {
    size_t n = ((struct data *)data)->n;
    double Df = ((struct data *)data)->Df;
    double R = ((struct data *)data)->R;
    double *time = ((struct data *) data)->time;
    double *sigma = ((struct data *) data)->sigma;
    char *m = ((struct data *) data)->m;
    size_t p = ((struct data *)data)->p;
    size_t w_flag = ((struct data *)data)->w_flag;

    size_t i;
    double xx, kon, koff;

    if (p == 1) {
        xx = gsl_vector_get (x, 0); 

        /* Importing one of the model functions */    
        double (*inverted_fun)();
        inverted_fun = &invlap_1;

        for (i = 0; i < n; i++) {
            /* Jacobian matrix J(i,j) = dfi / dxj, */
            /* where fi = (Yi - yi)/sigma[i],      */
            /*       Yi = A * exp(-lambda * i) + b  */
            /* and the xj are the parameters (A,lambda,b) */

            if (w_flag == 0) {
                gsl_matrix_set (J, i, 0, inverted_fun(time[i], xx, Df, R, m, 1));
            }
            else if (w_flag == 1)
            {
                gsl_matrix_set (J, i, 0, inverted_fun(time[i], xx, Df, R, m, 1)/sigma[i]);
            }
            else {
                fprintf(stderr, "ERROR: in 'model_df': Parameter w_flag is neither 0 nor 1.\n");
                exit(1);
            }
        }
    }
    else if (p == 2) {
        kon = gsl_vector_get (x, 0);
        koff = gsl_vector_get (x, 1);

        /* Importing one of the model functions */    
        double (*inverted_fun)();
        inverted_fun = &invlap_2;

        for (i = 0; i < n; i++) {
            /* Jacobian matrix J(i,j) = dfi / dxj, */
            /* where fi = (Yi - yi)/sigma[i],      */
            /*       Yi = A * exp(-lambda * i) + b  */
            /* and the xj are the parameters (A,lambda,b) */

            if (w_flag == 0) {
                gsl_matrix_set (J, i, 0, inverted_fun(time[i], kon, koff, Df, R, m, 1));
                gsl_matrix_set (J, i, 1, inverted_fun(time[i], kon, koff, Df, R, m, 2));
            }
            else if (w_flag == 1) {
                gsl_matrix_set (J, i, 0, inverted_fun(time[i], kon, koff, Df, R, m, 1)/sigma[i]);
                gsl_matrix_set (J, i, 1, inverted_fun(time[i], kon, koff, Df, R, m, 2)/sigma[i]);
            }
            else {
                fprintf(stderr, "ERROR: in 'model_df': Parameter w_flag is neither 0 nor 1.\n");
                exit(1);
            }
        }
    }
    else {
        fprintf(stderr, "ERROR: in 'model_f': Parameter p is neither 1 nor 2.\n");
        exit(1);
    }

    return GSL_SUCCESS;
}

int
model_fdf(const gsl_vector * x, void *data,
         gsl_vector * f, gsl_matrix * J) {
    model_f (x, data, f);
    model_df (x, data, J);

    return GSL_SUCCESS;
}

void
print_state (size_t iter, gsl_multifit_fdfsolver * s, size_t p) {
    if (p == 1) {
        printf ("iter: %3u xx = % 15.8f "
                "|f(x)| = %g\n",
                iter,
                gsl_vector_get (s->x, 0),
                gsl_blas_dnrm2 (s->f));
    }
    else if (p == 2) {
        printf ("iter: %3u x = % 15.8f % 15.8f "
                "|f(x)| = %g\n",
                iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                gsl_blas_dnrm2 (s->f));
    }
    else {
        fprintf(stderr, "ERROR: in 'print_state': Parameter p is neither 1 nor 2.\n");
        exit(1);
    }
}

void
bad_input(void) {
    fprintf(stderr, "Usage: cFDAP [-m model_type] [-d diffusion_constant]\n");
    fprintf(stderr, "             [-r2 half_activation_area] [-tini initial_time]\n");
    fprintf(stderr, "             [-tend end_time] [-n numsteps]\n");
    fprintf(stderr, "             [-kon0 initial_kon] [-koff0 initial_koff]\n");
    fprintf(stderr, "             [-x0 initial_x] [-w weights] [-i input]\n");
    fprintf(stderr, "             [-sd standard_deviation] [-o output]\n\n");
    fprintf(stderr, "  model_type:             reaction-diffusion model to fit with (mandatory parameter):\n");
    fprintf(stderr, "                          fullModel, hybridModel, reactionDominantPure\n");
    fprintf(stderr, "                          effectiveDiffusion\n");
    fprintf(stderr, "  diffusion_constant:     diffusion constant of unbound proteins (default: 11.0 µm2/s)\n");
    fprintf(stderr, "  half_activation_area:   half length of the activation area (default: 3.0 µm)\n");
    fprintf(stderr, "  initial_time:           initial time in the curve duration range (default: 0.0 s)\n");
    fprintf(stderr, "  end_time:               end time in the curve duration range (default: 112.0 s)\n");
    fprintf(stderr, "  numsteps:               number of steps in the FDAP curve (default: 113)\n");
    fprintf(stderr, "  initial_x:              starting value for x = kon/koff (default: 1.0)\n");
    fprintf(stderr, "                          IMPORTANT: use this parameter only with effectiveDiffusion\n");
    fprintf(stderr, "  initial_kon:            starting value for kon (default: 0.5)\n");
    fprintf(stderr, "  initial_koff:           starting value for koff (default: 0.5)\n");
    fprintf(stderr, "  weights:                whether to use weiths (0 - no, 1 - yes, default: no)\n");
    fprintf(stderr, "  input:                  name of input curve file (mandatory)\n");
    fprintf(stderr, "  standard_error:         name of input SD file (mandatory if weights = yes)\n");
    fprintf(stderr, "  output:                 prefix name of output file (Example: -o tau441wt\n");
    fprintf(stderr, "                          makes cFDAP output 'tau441wt_fit_parameters.dat'\n");
    fprintf(stderr, "                          and 'tau441wt_fit_curve.dat')\n");
    fprintf(stderr, "\n\n");
    exit(1);
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/* MAIN */
int
main(int argc, char *argv[]) {

    int i, status;
    unsigned int iter = 0;
    double chi, chi0;

    /* DEFAULTS */
    char m[80];
    char weights_name[1], curve_name[80], std_name[80], output_prefix[80];
    curve_name[0] = 0; std_name[0] = 0; output_prefix[0] = 0;
    size_t p;
    size_t w_flag = DEFAULT_FLAG_WEIGHT;
    size_t n = DEFAULT_N;
    double Df = DEFAULT_DF, R = DEFAULT_R;
    double t_ini = DEFAULT_T_INI, t_end = DEFAULT_T_END;
    double x_init_1[1] = { DEFAULT_XX_INIT };
    double x_init_2[2] = { DEFAULT_KON_INIT, DEFAULT_KOFF_INIT };

    fprintf(stderr, "\n");
    fprintf(stderr, "  --------------   cFDAP 0.1.0 (C) 2015\n");
    fprintf(stderr, "  |*    cFDAP  |   Author: Maxim Igaev\n");
    fprintf(stderr, "  | *          |   cFDAP is a fitting program for FDAP data\n");
    fprintf(stderr, "  |  ***       |   http://www.neurobiologie.uni-osnabrueck.de/\n");
    fprintf(stderr, "  |     *******|   https://github.com/moozzz\n");
    fprintf(stderr, "  --------------   Email: max_igaev@yahoo.com\n");
    fprintf(stderr, "\n");

    if ((argc < 2) || (argc > 25)) {
        bad_input();
    }

    /* First, a model must be chosen */
    if(strcmp(argv[1], "-m") != 0) {
        fprintf(stderr, "ERROR: First, a model must be chosen.\n\n");
        exit(1);
    }
    else {
        if(argc == 2) {
            fprintf(stderr, "ERROR: Specify the model's name.\n\n");
            exit(1);
        }
        else {
            if(strcmp(argv[2], "fullModel") == 0 ||
               strcmp(argv[2], "hybridModel") == 0 ||
               strcmp(argv[2], "reactionDominantPure") == 0) {
                strcpy(m, argv[2]);
                p = 2;
            }
            else if(strcmp(argv[2], "effectiveDiffusion") == 0) {
                    strcpy(m, argv[2]);
                p = 1;
            }
            else {
                fprintf(stderr, "ERROR: Unknown model '%s'\n\n", argv[2]);
                exit(1);
            }
        }
    }

    /* Parameters */
    for(i = 2; i < argc; i++) {
        if(strcmp(argv[i], "-d") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing diffusion coefficient.\n\n");
                exit(1);
            }
            Df = atof(argv[i + 1]);
            i++;
            if(Df <= 0.0) {
                fprintf(stderr, "ERROR: Would a zero or negative diffusion constant make sense?\n\n");  
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-r2") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing half area size.\n\n");
                exit(1);
            }
            R = atof(argv[i + 1]);
            i++;
            if(R <= 0.0) {
                fprintf(stderr, "ERROR: Would a zero or negative half area size make sense?\n\n");
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-x0") == 0) {
            if(p != 1) {
                fprintf(stderr, "ERROR: The model you chose has two fit parameters.\n\n");
                exit(1);
            }
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing initial value for x0.\n\n");
                exit(1);
            }
            x_init_1[0] = atof(argv[i + 1]);
            i++;
            if(x_init_1[0] < 0.0) {
                fprintf(stderr, "ERROR: Would a negative x = kon/koff make sense?\n\n");
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-kon0") == 0) {
            if(p != 2) {
                fprintf(stderr, "ERROR: The model you chose has one fit parameter.\n\n");
                exit(1);
            }
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing initial value for kon0.\n\n");
                exit(1);
            }
            x_init_2[0] = atof(argv[i + 1]);
            i++;
            if(x_init_2[0] < 0.0) {
                fprintf(stderr, "ERROR: Would a negative kon make sense?\n\n");
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-koff0") == 0) {
            if(p != 2) {
                fprintf(stderr, "ERROR: The model you chose has one fit parameter.\n\n");
                exit(1);
            }
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing initial value for koff0.\n\n");
                exit(1);
            }
            x_init_2[1] = atof(argv[i + 1]);
            i++;
            if(x_init_2[1] < 0.0) {
                fprintf(stderr, "ERROR: Would a negative koff make sense?\n\n");
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-tini") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing initial time.\n\n");
                exit(1);
            }
            t_ini = atof(argv[i + 1]);
            i++;
            if(t_ini < 0.0) {
                fprintf(stderr, "ERROR: Would a negative initial time make sense?\n\n");
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-tend") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing end time.\n\n");
                exit(1);
            }
            t_end = atof(argv[i + 1]);
            i++;
            if(t_end < t_ini) {
                fprintf(stderr, "ERROR: Would t_end < t_ini make sense?\n\n");
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-n") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: Missing number of time points.\n\n");
                exit(1);
            }
            n = atoi(argv[i + 1]);
            i++;
            if(n < 3) {
                fprintf(stderr, "ERROR: Your curve contatins less than 3 points? Are you kidding?\n\n");
                exit(1);
            }
        }
        else if(strcmp(argv[i], "-w") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: So, are you gonna use weights or not?\n\n");
                exit(1);
            }
            w_flag = atof(argv[i + 1]);
            i++;
            if ( !(w_flag == 0 || w_flag == 1) ) {
                fprintf(stderr, "ERROR: -w accepts only 0 (no) or 1 (yes) as arguments.\n\n");
            }
        }
        else if(strcmp(argv[i], "-i") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: No input curve files given.\n\n");
                exit(1);
            }
            strcpy(curve_name, argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-sd") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: No input SD files given.\n\n");
                exit(1);
            }
            strcpy(std_name, argv[i + 1]);
            i++;
        }
        else if(strcmp(argv[i], "-o") == 0) {
            if(i == argc - 1) {
                fprintf(stderr, "ERROR: No output prefix name given.\n\n");
                exit(1);
            }
            strcpy(output_prefix, argv[i + 1]);
            i++;
        }
        else {
            if(strncmp(argv[i], "-", 1) == 0) {
                fprintf(stderr, "\nERROR: Illegal option %s.\n\n", argv[i]);
                bad_input();
            }
        }
    }

    /* Checking whether input and output file names were given */
    char output_prefix_copy[80];
    strcpy(output_prefix_copy, output_prefix);
    if (curve_name[0] == 0 || output_prefix[0] == 0) {
        fprintf(stderr, "ERROR: File input/output is not defined correctly.\n\n");
        exit(1);
    }

    /* Checking whether std name was given if w_flag == 1 */
    if (w_flag == 1 && std_name[0] == 0) {
        fprintf(stderr, "ERROR: You want to use weighted fitting. Specify the std file.\n\n");
        exit(1);
    }

    /* Solver initialization */
    double time[n], y[n], sigma[n], best_fit[n];
    double stepSize = (t_end - t_ini)/(double) (n - 1);
    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    
    gsl_multifit_fdfsolver *s;
    gsl_vector *res_f;
    gsl_matrix *J = gsl_matrix_alloc(n, p); /* Jacobian matrix */
    gsl_matrix *covar = gsl_matrix_alloc (p, p); /* Covariance matrix */

    struct data d = { n, Df, R, time, y, sigma, m, p, w_flag };

    gsl_multifit_function_fdf f;
    gsl_vector_view x;
    if (p == 1) {
        x = gsl_vector_view_array (x_init_1, p);
    }
    else if (p == 2)
    {
        x = gsl_vector_view_array (x_init_2, p);
    }
    else {
        fprintf(stderr, "ERROR: in main: Parameter p is neither 1 nor 2.\n");
        exit(1);
    }

    f.f = &model_f;
    f.df = &model_df;
    f.fdf = &model_fdf;
    f.n = n;
    f.p = p;
    f.params = &d;

    /* Importing the FDAP curve to be fitted */
    double temp;
    FILE *input_curve = fopen(curve_name, "r");
    printf("Opening the curve file...\n");
    if(input_curve == NULL) {
        fprintf(stderr, "ERROR: input_curve file cannot be opened.\n");
        fprintf(stderr, "ERROR: Probably, cFDAP and input_curve must be in the same folder.\n");
        exit(1);
    }
    for(i = 0; i < NELEMS_1D(y); i++) {
        time[i] = t_ini + (double) i*stepSize;
        fscanf(input_curve, "%lf", &temp);
        y[i] = temp;
        if (i < 5) printf("data: %f %g\n", time[i], y[i]);
    }
    if(time[0] == 0.0) time[0] = 0.01;
    printf("Curve file has been successfully read in.\n\n");
    fclose(input_curve);

    /* Importing the errors for the FDAP curve if w_flag == 1 */
    if (w_flag == 1) {
        FILE *error_curve = fopen(std_name, "r");
        printf("Opening the error file...\n");
        if(error_curve == NULL) {
            fprintf(stderr, "ERROR: error_curve file cannot be opened.\n");
            fprintf(stderr, "ERROR: Probably, cFDAP and error_curve must be in the same folder.\n");
            exit(1);
        }
        for(i = 0; i < NELEMS_1D(y); i++) {
            fscanf(error_curve, "%lf", &temp);
            sigma[i] = temp;
            if (i < 5) printf("data: %g\n", sigma[i]);
        }
        printf("Error file has been successfully read in.\n\n");
        fclose(error_curve);
    }

    /* Allocating a new instance for the solver */
    s = gsl_multifit_fdfsolver_alloc (T, n, p);
    if (w_flag == 0) {
        printf("Initializing '%s' solver with NO weights...\n\n", gsl_multifit_fdfsolver_name(s));
    }
    else {
        printf("Initializing '%s' solver with weights...\n\n", gsl_multifit_fdfsolver_name(s));
    }

    /* Initializing a solver with a starting point x */
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

    /* Computing the initial residual norm */
    res_f = gsl_multifit_fdfsolver_residual(s);
    chi0 = gsl_blas_dnrm2(res_f);

    /* Solving the system with a maximum of 500 iterations */
    print_state (iter, s, p);
    do {
        iter++;
        status = gsl_multifit_fdfsolver_iterate (s);

        printf ("current status = %s\n", gsl_strerror (status));

        print_state (iter, s, p);

        if (status)
            break;

        status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
    }
    while (status == GSL_CONTINUE && iter < 500);

    /* Computing the Jacobian and covariace matrix */
    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar (J, 0.0, covar);

    /* Computing the final residual norm */
    chi = gsl_blas_dnrm2(res_f);
    
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    printf("\nSummary from method '%s':\n", gsl_multifit_fdfsolver_name(s));
    printf("Number of iterations done: %zu\n", gsl_multifit_fdfsolver_niter(s));
    printf("Function evaluations: %zu\n", f.nevalf);
    printf("Jacobian evaluations: %zu\n", f.nevaldf);
    printf("Initial |f(x)| = %g\n", chi0);
    printf("Final |f(x)| = %g\n", chi);

    {
        double dof = n - p;
        double c = GSL_MAX_DBL(1, chi / sqrt(dof));

        printf("chisq/dof = %g\n", pow(chi, 2.0)/dof);

        if (p == 1) {
            printf ("x          = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
            printf ("x conf     = +/- %.5f %.5f %.5f\n", c*ERR(0)*gsl_cdf_tdist_Pinv(0.95, dof), c*ERR(0)*gsl_cdf_tdist_Pinv(0.975, dof), c*ERR(0)*gsl_cdf_tdist_Pinv(0.99, dof));
            printf ("bound      = %.5f +/- %.5f\n", 100.0 - 100.0/(1.0 + FIT(0)), 1.0);
        }
        else if (p == 2) {
            printf ("kon        = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
            printf ("kon conf   = +/- %.5f %.5f %.5f\n", c*ERR(0)*gsl_cdf_tdist_Pinv(0.95, dof), c*ERR(0)*gsl_cdf_tdist_Pinv(0.975, dof), c*ERR(0)*gsl_cdf_tdist_Pinv(0.99, dof));
            printf ("koff       = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
            printf ("koff conf  = +/- %.5f %.5f %.5f\n", c*ERR(1)*gsl_cdf_tdist_Pinv(0.95, dof), c*ERR(1)*gsl_cdf_tdist_Pinv(0.975, dof), c*ERR(1)*gsl_cdf_tdist_Pinv(0.99, dof));
            printf ("bound      = %.5f +/- %.5f\n", 100.0 - 100.0/(1.0 + FIT(0)/FIT(1)), 100.0*(c*ERR(0)/FIT(1) - FIT(0)*c*ERR(1)/FIT(1)/FIT(1))/(1.0 + FIT(0)/FIT(1))/(1.0 + FIT(0)/FIT(1)));
        }
        else {
            fprintf(stderr, "ERROR: in main: Parameter p is neither 1 nor 2.\n");
            exit(1);
        }

        /* Writing the fit parameters */
        FILE *fit_params = fopen(strcat(output_prefix, "_fit_params.dat"), "w");
        if (p == 1) {
            fprintf(fit_params, "chisq/dof %g\n", pow(chi, 2.0) / dof);
            fprintf(fit_params, "x_fit %.5f\n", FIT(0));
            fprintf(fit_params, "x_error %.5f\n", c*ERR(0));
            fprintf(fit_params, "bound %.5f\n", 100.0 - 100.0/(1.0 + FIT(0)));
        }
        else if (p == 2) {
            fprintf(fit_params, "chisq/dof %g\n", pow(chi, 2.0) / dof);
            fprintf(fit_params, "kon_fit %.5f\n", FIT(0));
            fprintf(fit_params, "kon_error %.5f\n", c*ERR(0));
            fprintf(fit_params, "koff_fit %.5f\n", FIT(1));
            fprintf(fit_params, "koff_error %.5f\n", c*ERR(1));
            fprintf(fit_params, "bound %.5f\n", 100.0 - 100.0/(1.0 + FIT(0)/FIT(1)));
            fprintf(fit_params, "bound_error %.5f\n", 100.0*(c*ERR(0)/FIT(1) - FIT(0)*c*ERR(1)/FIT(1)/FIT(1))/(1.0 + FIT(0)/FIT(1))/(1.0 + FIT(0)/FIT(1)));
        }
        else {
            fprintf(stderr, "ERROR: in main: Parameter p is neither 1 nor 2.\n");
            exit(1);
        }
        fclose(fit_params);
    }

    printf ("\nSTATUS = %s\n\n", gsl_strerror (status));

    /* Writing the best fit */
    FILE *fit_curve = fopen(strcat(output_prefix_copy, "_best_fit.dat"), "w");
    if (p == 1) {
        for(i = 0; i < NELEMS_1D(best_fit); i++) {
            fprintf(fit_curve, "%f\n", invlap_1(time[i], FIT(0), Df, R, m, 0));
        }
    }
    else if (p == 2) {
        for(i = 0; i < NELEMS_1D(best_fit); i++) {
            fprintf(fit_curve, "%f\n", invlap_2(time[i], FIT(0), FIT(1), Df, R, m, 0));
        }
    }
    else {
        fprintf(stderr, "ERROR: in main: Parameter p is neither 1 nor 2.\n");
        exit(1);
    }
    fclose(fit_curve);

    gsl_multifit_fdfsolver_free (s);
    gsl_matrix_free (covar);
    gsl_matrix_free (J);

    return 0;
}

#pragma once

 

void testTtt();
void tttomo(double* ttpick, double* ttslns, int nrec, double trsp, double rrsp, double Calpr,
	double dtf, double tlod, double slns, double dr, double vav,
	double* tmogrm, double* Rtmogrm, double* ttfit, double* ssfit,
	double* ss0, double* ssd, double* dop);
float TTcstfunc(float* param);
float tvtm(float V);
float vg2v(float V);
float brent(float ax, float bx, float cx, float (*f)(float), float tol, float* xmin);
void amoeba(float** p, float y[], int ndim, float ftol, float (*funk)(float[]), int* nfunk);
float amotry(float** p, float y[], float psum[], int ndim,
    float (*funk)(float[]), int ihi, float fac);
float* vector(long, long);
void free_vector(float*, long, long);
float** matrix(long, long, long, long);
void free_matrix(float**, long, long, long, long);
#define NR_END 1
#define FREE_ARG char*
static float T2R, TR, RR, nR, sdo, gv, Vf, V0, TTrsd[12], Vdat;

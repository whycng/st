#pragma region 备份
//#pragma once
//
// 
//
//void testTtt();
//void tttomo(double* ttpick, double* ttslns, int nrec, double trsp, double rrsp, double Calpr,
//	double dtf, double tlod, double slns, double dr, double vav,
//	double* tmogrm, double* Rtmogrm, double* ttfit, double* ssfit,
//	double* ss0, double* ssd, double* dop);
//float TTcstfunc(float* param);
//float tvtm(float V);
//float vg2v(float V);
//float brent(float ax, float bx, float cx, float (*f)(float), float tol, float* xmin);
//void amoeba(float** p, float y[], int ndim, float ftol, float (*funk)(float[]), int* nfunk);
//float amotry(float** p, float y[], float psum[], int ndim,
//    float (*funk)(float[]), int ihi, float fac);
//float* vector(long, long);
//void free_vector(float*, long, long);
//float** matrix(long, long, long, long);
//void free_matrix(float**, long, long, long, long);
//#define NR_END 1
//#define FREE_ARG char*
//static float T2R, TR, RR, nR, sdo, gv, Vf, V0, TTrsd[12], Vdat;


//#define NR_END 1
//#define FREE_ARG char*
//class Tttomo {
//public:
//    void testTtt();
//    void tttomo(double* ttpick, double* ttslns, int nrec, double trsp, double rrsp, double Calpr,
//        double dtf, double tlod, double slns, double dr, double vav,
//        double* tmogrm, double* Rtmogrm, double* ttfit, double* ssfit,
//        double* ss0, double* ssd, double* dop);
//    float TTcstfunc(float* param); // amoeba使用
//    float tvtm(float V);
//    float vg2v(float V);
//    float brent(float ax, float bx, float cx, float (*f)(float), float tol, float* xmin);
//    void amoeba(float** p, float y[], int ndim, float ftol, 
//        float (*funk)(float[]), int* nfunk);
//    //void amoeba(float** p, float y[], int ndim, float ftol, float (*funk)(float[]), int* nfunk);
//    float amotry(float** p, float y[], float psum[], int ndim,
//        float (*funk)(float[]), int ihi, float fac);
//    float* vector(long, long);
//    void free_vector(float*, long, long);
//    float** matrix(long, long, long, long);
//    void free_matrix(float**, long, long, long, long);
//
//private:
//    float T2R;
//    float TR;
//    float RR;
//    float nR;
//    float sdo;
//    float gv;
//    float Vf;
//    float V0;
//    float TTrsd[12];
//    float Vdat;
//};


#pragma endregion
#pragma once 


void testTtt();
void tttomo(double* ttpick, double* ttslns, int nrec, double trsp, double rrsp, double Calpr,
	double dtf, double tlod, double slns, double dr, double vav,
	double* tmogrm, double* Rtmogrm, double* ttfit, double* ssfit,
	double* ss0, double* ssd, double* dop);
//float TTcstfunc(float* param);
float TTcstfunc(float* param, float* T2R, float TR, float RR, float nR,
	float sdo, float gv, float Vf, float V0, float TTrsd[12], float Vdat);
//float tvtm(float V);
float tvtm(float V, float* T2R, float Vf, float sdo, float V0);
float vg2v(float V, float* T2R, float Vf, float sdo, float gv, float V0);
//float vg2v(float V);
float brent(float ax, float bx, float cx, float (*f)(float, float* T2R, float Vf, float sdo, float gv, float V0),
	float tol, float* xmin, float* T2R, float Vf, float sdo, float gv, float V0);
//float brent(float ax, float bx, float cx, float (*f)(float), float tol, float* xmin);


//void amoeba(float** p, float y[], int ndim, float ftol, float (*funk)(float[]), int* nfunk);
void amoeba(float** p, float y[], int ndim, float ftol,
	float (*funk)(float[], float* T2R, float TR, float RR, float nR,
		float sdo, float gv, float Vf, float V0, float TTrsd[12], float Vdat), int* nfunk,
	float* T2R, float TR, float RR, float nR,
	float sdo, float gv, float Vf, float V0, float TTrsd[12], float Vdat);

//float amotry(float** p, float y[], float psum[], int ndim,
//    float (*funk)(float[]), int ihi, float fac);
float amotry(float** p, float y[], float psum[], int ndim,
	float (*funk)(float[], float* T2R, float TR, float RR, float nR,
		float sdo, float gv, float Vf, float V0, float TTrsd[12], float Vdat), int ihi, float fac,
	float* T2R, float TR, float RR, float nR,
	float sdo, float gv, float Vf, float V0, float TTrsd[12], float Vdat);
float* vector(long, long);
void free_vector(float*, long, long);
float** matrix(long, long, long, long);
void free_matrix(float**, long, long, long, long);

#define NR_END 1
#define FREE_ARG char* 




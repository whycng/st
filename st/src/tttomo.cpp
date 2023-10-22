/*---------------------------------------------------------------------------
 *
 * 名称：tttomo
 *
 * 目的: 反演等梯度速度模型参数并计算层析结果
 *
 *  source       波形数据
 *  dt           时间间隔
 *  wdwdth       找波至设置的窗长
 *  ftm         波至时间
 *-------------------------------------------------------------------------*/

 // #include "ftonset.h" as in the following
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "tttomo.h"
#include <iostream>

void testTtt() {
	;
}

void tttomo(double* ttpick, double* ttslns, int nrec, double trsp, double rrsp, double Calpr,
	double dtf, double tlod, double slns, double dr, double vav,
	double* tmogrm, double* Rtmogrm, double* ttfit, double* ssfit,
	double* ss0, double* ssd, double* dop)
	// ttpick, 提取波至(1~nwf 道)
	// ttslns, 时差积分(1~nwf 道, 不考虑流体部分)
	// Calpr, 井径, m
	// dtf, 流体时差
	// tlod, 仪器直径
	// slns, 当前深度的时差
	// dr, 成像径向步长(成像深度/60)
	// vav, 平均速度
	// tmogrm, 速度剖面(实际)
	// Rtmogrm, 速度剖面(相对)
	// ttfit, 对ttpick拟合的波至
	// ssfit, 对slns拟合的波至
	// ss0, 井壁上的时差
	// ssd, 最大穿透深度上的时差
	// dop, 最大穿透深度
{

	int i, j, j0, iter;
	float  V, Vdeep, testTT, ttf, slnsav, rd, rb;
	float** p, * y, * param;
	float a, b, c, fmin;

	//Allocate matrix and vector for simplex minimization routine //
	p = matrix(1, 3, 1, 2);
	y = vector(1, 3);
	param = vector(1, 2);
	nR = nrec;
	RR = rrsp;
	TR = trsp;
	Vf = 1e6 / dtf;
	rb = Calpr / 2;

	//Traveltime residue, in seconds //
	testTT = 0.0;
	for (i = 0; i < nrec; i++)
	{
		ttf = ttpick[i] - ttslns[i];
		testTT += ttf;
		TTrsd[i] = ttf * 1e-6;
	}
	testTT /= nR;

	slnsav = ttslns[0] / TR;

	sdo = Calpr - tlod;

	ttf = sdo * dtf * sqrt(1.0 - pow(slnsav / dtf, 2));
	Vdat = 1e6 / slns;

	if (testTT - ttf <= 0.0)
	{
		gv = 0.1;
		V0 = Vdat;
		Vdeep = Vdat;
	}
	else
	{

		//Set initial simplex (3 vertices) to start downhill simplex //


		p[1][1] = 100.0; //Initial velocity gradien ft/s /ft Vertix 1//
		p[1][2] = 0.8 * Vdat; //velocity  ft/s //
		p[2][1] = 600.0; //Vertix 2//
		p[2][2] = Vdat;
		p[3][1] = 300.0; //Vertix 3//
		p[3][2] = Vdat * 1.1;

		for (i = 1; i <= 3; i++)
		{
			for (j = 1; j <= 2; j++) param[j] = p[i][j];
			y[i] = TTcstfunc(param);
		}

		amoeba(p, y, 2, 1.0e-3, &TTcstfunc, &iter); //阿米巴巴//
		for (j = 1; j <= 2; j++) param[j] = p[1][j];

		gv = param[1]; //Model velocity gradient//
		V0 = param[2]; //Model borehole wall velocity//

		//Now, find the deepest penetration model velocity //
		T2R = TR + (nR - 1) * RR;
		a = Vdat;
		b = 1.2 * Vdat;
		c = 2 * Vdat;
		fmin = brent(a, b, c, &vg2v, 1e-5, &Vdeep);

	}

	//Compute Tomogram based on velocity model parameters //

	*dop = (Vdeep - V0) / gv;
	if (*dop < 0.01 * rb) *dop = 0.01 * rb;
	*ss0 = 1e6 / V0;
	*ssd = 1e6 / Vdeep;
	j0 = (int)(rb / dr);// Calpr / 2 / dr

	
	// 这一段从 55-65(大概，根据j0定) vav=10，赋值
	for (j = 0; j <= j0; j++)
	{ 
		tmogrm[60 + j] = 0.7 * vav;
		tmogrm[60 - j] = 0.7 * vav;
		Rtmogrm[60 + j] = 0.0;
		Rtmogrm[60 - j] = 0.0; 
	/*	if (60 + j >= 55)
		{ 
			std::cout << "【Message】 tmogrm[60 + " << j << "]" << tmogrm[60 + j] << std::endl;
			std::cout << "【Message】 tmogrm[60 - " << j << "]" << tmogrm[60 - j] << std::endl;
		}*/
	}

	for (j = j0 + 1; j < 61; j++)
	{
		//计算各个径向点的深度rd
		rd = j * dr;
		
		////exponential variation model   这部分是关键，一个衰减指数实现单极纵波的几个短周期的快速拟合// 
		V = (Vdeep - V0) * exp(-2 * (rd - rb) / (*dop));
	/*	std::cout << "【Message142】 V:" << V << std::endl
			<< "Vdeep - V0:" << Vdeep - V0 << std::endl
			<< "exp(-2 * (rd - rb) / (*dop)): " << exp(-2 * (rd - rb) / (*dop)) << std::endl;*/
 
		tmogrm[60 + j] = (Vdeep - V);
		tmogrm[60 - j] = tmogrm[60 + j];
		Rtmogrm[60 + j] = V / Vdeep * 100;
		Rtmogrm[60 - j] = Rtmogrm[60 + j];
	}

	//Compute model traveltime fit to data //

	for (i = 0; i < nrec; i++)
	{
		T2R = TR + i * RR;
		if (gv < 1)
			V = V0;
		else
		{
			a = V0;
			b = 1.2 * V0;
			c = 2 * V0;
			fmin = brent(a, b, c, &vg2v, 1e-6, &V);
		}
		// std::cout << "V:" << V  << std::endl; -- 8000
		ttfit[i] = tvtm(V);
		double t = ttfit[i];
		//std::cout << "【Message172】ttfit[" << i << "]:" << ttfit[i] << std::endl;
	}

	// Compute model velocity fit to data //
	{
		float avt, avr, st, sr, x2r;
		avt = 0.0;
		for (i = 0; i < nrec; i++) {
			avt += ttfit[i];
		}
		avt /= nR;
		avr = TR + 0.5 * RR * (nR - 1);

		st = 0.0;
		sr = 0.0;
		for (i = 0; i < nrec; i++)
		{
			x2r = TR + RR * i;
			st += (x2r - avr) * (ttfit[i] - avt);
			sr += (x2r - avr) * (x2r - avr);

		}

		V = sr / st;
	}

	for (i = 0; i < nrec; i++)
	{
		double test_a = ttslns[i];
		ttfit[i] = ttslns[i] + (ttfit[i] - (TR + i * RR) / V) * 1e6;
		double test_z = (ttfit[i] - (TR + i * RR) / V);
		double test_y = (ttfit[i] - (TR + i * RR) / V) * 1e6;
		double test_x = ttfit[i];
		//std::cout << "【Message】<tttomo>ttfit[" << i << "]:" << ttfit[i] << std::endl;
	}

	*ssfit = 1.0e6 / V;

	free_matrix(p, 1, 3, 1, 2);
	free_vector(y, 1, 3);
	free_vector(param, 1, 2);

	return;

}

/*--------------------------------------------------------------
 * tvtm:
 *
 * 目的: 计算等速梯度模型的传播时间
 *
 * 输入参数:
 *     T2R =  发射器到接收器距离
 *     V   =  在T2R穿透地层最大速度
 *     VF   = 井眼流体速度
 *     sdo  = tool stand off
 * 输出结果:
 *     计算达到时间TT
 *
 *----------------------------------------------------------------*/

float tvtm(float V)
{
	// function TT=tvtm(V)

	// static float T2R, sdo, Vf, V0//
	float sn, cs, gv1, Z, TT;

	sn = Vf / V;
	cs = sqrt(1.0 - sn * sn);
	Z = T2R - sdo * sn / cs;
	gv1 = 2 * V * sqrt(1.0 - V0 * V0 / V / V + 1e-6) / Z;
	if (gv1 < 100.0)
		TT = sdo / Vf / cs + Z / V;
	else
	{
		//	TT = sdo/Vf/cs + 2.*asinh(0.5*gv1*Z/V0)/gv1//
		TT = 0.5 * gv1 * Z / V0;
		TT = log(TT + sqrt(1.0 + TT * TT));
		TT = sdo / Vf / cs + 2. * TT / gv1;
	}

	//std::cout << "【Message】<tvtm> TT:" << TT << std::endl;
	return(TT);
}

/*--------------------------------------------------------------
 * TTcstfunc:
 *
 * 目的: 构建惩罚函数用于测量数据和模拟到时数据Cost function as misfit error between measured and synthetic traveltime and velocity data
 *
 * 输入:
 *     param[1]   = 速度梯度
 *     param[2]  =  井壁地层速度
 *
 * 输出:
 *     TTcstfunc函数值
 *
 *----------------------------------------------------------------*/
float TTcstfunc(float* param)
{
	// function Terr=TTcstfunc(param)//

 //static float T2R, TR, RR, nR, sdo, gv, Vf, V0, TTrsd[12], Vdat//

	float Rcv[12], Tsyn[12];
	float a, b, c, fmin;
	float vsyn, avt, avr, st, sr, Tscl, Terr;
	int i, nrec;

	gv = param[1];
	V0 = param[2];
	nrec = (int)nR;
	for (i = 0; i < nrec; i++) Rcv[i] = TR + i * RR;


	//Compute model traveltime fit to data //

	for (i = 0; i < nrec; i++)
	{
		T2R = Rcv[i];
		if (gv < 1)
			vsyn = V0;
		else
		{
			a = V0;
			b = 1.2 * V0;
			c = 2 * V0;

			fmin = brent(a, b, c, &vg2v, 1e-5, &vsyn);
		}

		Tsyn[i] = tvtm(vsyn);
	}


	//Linear fitting travel time to get synthetic velocity across array //
	avt = 0.0;
	for (i = 0; i < nrec; i++)
	{
		avt += Tsyn[i];
	}

	avt /= nR;
	avr = TR + 0.5 * (Rcv[nrec - 1] - Rcv[0]);

	st = 0.0;
	sr = 0.0;
	for (i = 0; i < nrec; i++)
	{
		st += (Rcv[i] - avr) * (Tsyn[i] - avt);
		sr += (Rcv[i] - avr) * (Rcv[i] - avr);
	}
	vsyn = sr / st;

	//Synthetic traveltime residue data //
	for (i = 0; i < nrec; i++) Tsyn[i] = Tsyn[i] - Rcv[i] / vsyn;

	Tscl = (sdo / Vf);
	Tscl = Tscl * Tscl;
	Terr = 0.0;
	for (i = 0; i < nrec; i++)
	{
		st = Tsyn[i] - TTrsd[i];
		Terr += st * st;
	}

	Terr = Terr / Tscl / nR;
	st = (Vdat - vsyn) / Vf;
	Terr = Terr + st * st;
	Terr = log(0.000001 + Terr);

	return(Terr);

}


/** ------------------------------------------------------------------------------------------------------------------------------
c  Cost function to compute maximum velocity at T2R from velovity gradient
c 输入:
c        T2R = Transmitter to receiver spacing
c        V   =  maximum velocity of penetrated depth at T2R
c        gv =  velocity gradient of the model
c
-------------------------------------------------------------------------------------------------------------------------------****/
float vg2v(float V)
{
	//function TT=vg2v(V)//
// static float T2R, sdo, gv, Vf, V0//
	float sn, cs, TT, Z;
	sn = Vf / V;
	cs = sqrt(1.0 - sn * sn + 0.0000001);
	Z = T2R - sdo * sn / cs;
	TT = gv * Z - 2 * V * sqrt(1.0 - V0 * V0 / V / V + 0.000001);
	TT = TT * TT;
	return(TT);
}

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
/*
ITMAX is the maximum allowed number of iterations;
CGOLD is the golden ratio;
ZEPS is a small number that protects against trying to achieve fractional accuracy for a minimum that
happens to be exactly zero. */
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
float brent(float ax, float bx, float cx, float (*f)(float), float tol, float* xmin)
/* Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
the minimum is returned as xmin, and the minimum function value is returned as brent, the
returned function value.*/
{
	int iter;
	float a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;

	float e = 0.0; //This will be the distance moved on the step before last//
	a = (ax < cx ? ax : cx); //a and b must be in ascending order//
	b = (ax > cx ? ax : cx); //but input abscissas need not be//
	x = w = v = bx; //Initializations...//
	fw = fv = fx = (*f)(x);
	for (iter = 1; iter <= ITMAX; iter++) { //Main program loop//
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) { //Test for done here//
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) { //Construct a trial parabolic fit//
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) p = -p;
			q = fabs(q);
			etemp = e;
			e = d;
			if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			/*The above conditions determine the acceptability of the parabolic fit. Here we
				take the golden section step into the larger of the two segments****/
			else {
				d = p / q; //Take the parabolic step//
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		}
		else {
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = (*f)(u);
		//This is the one function evaluation per iteration//
		if (fu <= fx) { //Now decide what to do with our function evaluation//
			if (u >= x) a = x; else b = x;
			SHFT(v, w, x, u) //Housekeeping follows//
				SHFT(fv, fw, fx, fu)
		}
		else {
			if (u < x) a = u; else b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			}
			else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		} //Done with housekeeping. Back for another iteration//
	}
	printf("Too many iterations in brent");  //Brent方法快速找根！！wzt//
	*xmin = x; //Never get here//
	return fx;
}


#define NMAX 600  //Maximum allowed number of function evaluations//
#define GET_PSUM \
for (j=1;j<=ndim;j++) {\
for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void amoeba(float** p, float y[], int ndim, float ftol,
	float (*funk)(float[]), int* nfunk)
	/****------------------------------------------------------------------------------------------------------------
	Description:
	Multidimensional minimization of the function funk(x) where x[1..ndim] is a vector in ndim
	dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[1..ndim+1]
	[1..ndim] is input. Its ndim+1 rows are ndim-dimensional vectors which are the vertices of
	the starting simplex. Also input is the vector y[1..ndim+1], whose components must be preinitialized
	to the values of funk evaluated at the ndim+1 vertices (rows) of p; and ftol the
	fractional convergence tolerance to be achieved in the function value (n.b.!). On output, p and
	y will have been reset to ndim+1 new points all within ftol of a minimum function value, and
	nfunk gives the number of function evaluations taken.
	-----------------------------------------------------------------------------------------------------------------*****/
{
	float amotry(float** p, float y[], float psum[], int ndim,
		float (*funk)(float[]), int ihi, float fac);
	int i, ihi, ilo, inhi, j, mpts = ndim + 1;
	float rtol, sum, swap, ysave, ytry, * psum;
	psum = vector(1, ndim);
	*nfunk = 0;
	GET_PSUM
		for (;;) {
			ilo = 1;
			/***------------------------------------------------------------------------------------------------
			First we must determine which point is the highest (worst), next-highest, and lowest
			(best), by looping over the points in the simplex
			---------------------------------------------------------------------------------------------------****/
			ihi = y[1] > y[2] ? (inhi = 2, 1) : (inhi = 1, 2);
			for (i = 1; i <= mpts; i++) {
				if (y[i] <= y[ilo]) ilo = i;
				if (y[i] > y[ihi]) {
					inhi = ihi;
					ihi = i;
				}
				else if (y[i] > y[inhi] && i != ihi) inhi = i;
			}
			rtol = 2.0 * fabs(y[ihi] - y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]));

			//Compute the fractional range from highest to lowest and return if satisfactory//

			if (rtol < ftol) { //If returning, put best point and value in slot 1//
				SWAP(y[1], y[ilo])
					for (i = 1; i <= ndim; i++) SWAP(p[1][i], p[ilo][i])
						break;
			}
			if (*nfunk >= NMAX) {
				// printf("NMAX exceeded");
				return;
			}
			*nfunk += 2;
			/***--------------------------------------------------------------------------------------
			Begin a new iteration. First extrapolate by a factor .1 through the face of the simplex
			across from the high point, i.e., reflect the simplex from the high point.
			-----------------------------------------------------------------------------------------****/
			ytry = amotry(p, y, psum, ndim, funk, ihi, -1.0);
			if (ytry <= y[ilo])

				//Gives a result better than the best point, so try an additional extrapolation by a factor 2//

				ytry = amotry(p, y, psum, ndim, funk, ihi, 2.0);
			else if (ytry >= y[inhi]) {

				/***The reflected point is worse than the second-highest, so look for an intermediate
				lower point, i.e., do a one-dimensional contraction****/

				ysave = y[ihi];
				ytry = amotry(p, y, psum, ndim, funk, ihi, 0.5);
				if (ytry >= ysave) {

					//Can't seem to get rid of that high point. Better contract around //
					for (i = 1; i <= mpts; i++) {
						// the lowest (best) point//

						if (i != ilo) {
							for (j = 1; j <= ndim; j++)
								p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
							y[i] = (*funk)(psum);
						}
					}
					*nfunk += ndim;
					GET_PSUM
				}
			}
			else --(*nfunk);
		}
	free_vector(psum, 1, ndim);
}

float amotry(float** p, float y[], float psum[], int ndim,
	float (*funk)(float[]), int ihi, float fac)
	/**----------------------------------------------------------------------
	Extrapolates by a factor fac through the face of the simplex across from the high point, tries
	it, and replaces the high point if the new point is better.
	----------------------------------------------------------------------**/
{
	int j;
	float fac1, fac2, ytry, * ptry;
	ptry = vector(1, ndim);
	fac1 = (1.0 - fac) / ndim;
	fac2 = fac1 - fac;
	for (j = 1; j <= ndim; j++) ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
	ytry = (*funk)(ptry); //Evaluate the function at the trial point//
	if (ytry < y[ihi]) { //If it's better than the highest, then replace the highest//
		y[ihi] = ytry;
		for (j = 1; j <= ndim; j++) {
			psum[j] += ptry[j] - p[ihi][j];
			p[ihi][j] = ptry[j];
		}
	}
	free_vector(ptry, 1, ndim);
	return ytry;
}

float* vector(long nl, long nh)
//allocate a float vector with subscript range v[nl..nh] //
{
	float* v;
	v = (float*)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
	if (!v) printf("allocation failure in vector()");
	return v - nl + NR_END;
}

float** matrix(long nrl, long nrh, long ncl, long nch)
//allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] //
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float** m;
	//allocate pointers to rows //
	m = (float**)malloc((size_t)((nrow + NR_END) * sizeof(float*)));
	if (!m) printf("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	//allocate rows and set pointers to them//
	m[nrl] = (float*)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
	if (!m[nrl]) printf("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i - 1] + ncol;
	//return pointer to array of pointers to rows//
	return m;
}

void free_vector(float* v, long nl, long nh)
//free a float vector allocated with vector() //
{
	free((FREE_ARG)(v + nl - NR_END));
}

void free_matrix(float** m, long nrl, long nrh, long ncl, long nch)
//free a float matrix allocated by matrix() //
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}


//The gateway routine//
//void mexFunction( int nlhs, mxArray *plhs[],
//                  int nrhs, const mxArray *prhs[])
//{
//  double *ttpick,  *ttslns;
//  int nrec; 
//  double trsp,  rrsp,  Calpr, dtf, tlod, slns, dr,  vav;
//  double *tmogrm, *Rtmogrm, *ttfit, *ssfit, *ss0, *ssd, *dop;
//  
////Check for proper number of arguments//
///* NOTE: You do not need an else statement when using
//     mexErrMsgTxt within an if statement. It will never
//     get to the else statement if mexErrMsgTxt is executed.
//     (mexErrMsgTxt breaks you out of the MEX-file.) */
//   
//  if(nrhs != 11) 
//    mexErrMsgTxt("eleven (11) inputs required.");
//   if(nlhs != 7) 
//     mexErrMsgTxt("seven output required.");
//  
////Check and read-in input data //
//
////Get the scalar input //
// 
//  nrec  = mxGetScalar(prhs[2]);
//  trsp  = mxGetScalar(prhs[3]);
//  rrsp  = mxGetScalar(prhs[4]);
//  Calpr = mxGetScalar(prhs[5]);
//  dtf   = mxGetScalar(prhs[6]);
//  tlod  = mxGetScalar(prhs[7]);
//  slns  = mxGetScalar(prhs[8]);
//  dr    = mxGetScalar(prhs[9]);
//  vav   = mxGetScalar(prhs[10]);
//
//
////Create a pointer to the input vector ttpick and ttslsn//
//  //nrec   = mxGetN(prhs[0])//
//  ttpick = mxGetPr(prhs[0]);
//  ttslns = mxGetPr(prhs[1]);
//
//
////Create a C pointer to a copy of the output matrix//
//  plhs[0] = mxCreateDoubleMatrix(1,121, mxREAL);
//  tmogrm = mxGetPr(plhs[0]); 
//  plhs[1] = mxCreateDoubleMatrix(1,121, mxREAL);
//  Rtmogrm = mxGetPr(plhs[1]); 
//  plhs[2] = mxCreateDoubleMatrix(1,nrec, mxREAL);
//  ttfit = mxGetPr(plhs[2]); 
//  plhs[3] = mxCreateDoubleMatrix(1,1, mxREAL);
//  ssfit   = mxGetPr(plhs[3]); 
//  plhs[4] = mxCreateDoubleMatrix(1,1, mxREAL);
//  ss0     = mxGetPr(plhs[4]); 
//  plhs[5] = mxCreateDoubleMatrix(1,1, mxREAL);
//  ssd     = mxGetPr(plhs[5]); 
//  plhs[6] = mxCreateDoubleMatrix(1,1, mxREAL);
//  dop     = mxGetPr(plhs[6]); 
//
//tttomo(ttpick, ttslns, nrec, trsp, rrsp, Calpr, dtf, tlod, slns, dr, vav, 
//			 tmogrm, Rtmogrm, ttfit, ssfit, ss0, ssd, dop);
//  return;
//}
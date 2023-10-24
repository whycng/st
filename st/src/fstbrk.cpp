/*---------------------------------------------------------------------------
 * fstbrk
 *
 *  程序目的: 得到首波波至，Find the first break.
 *  程序原理: 通过衰减半正弦实现首波的快速模拟和计算
 *        fit the wave's first half cycle with a function.
 *
 *  Arguments     Type      in-out    Description
  -------------------  -------  ---------------------------------------
 *     source            float      input     waveform data
 *     dt                                              time increment.
 *     wdwdth                                    window width of source to find first break.
 *     *ftm                                          output first break time
 *-------------------------------------------------------------------------*/

 // #include "ftonset.h" as in the following
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "fstbrk.h"


//void test1(float slowness, float tt, int nsamp, int nR, float RR, float period,
//	float wdwdth, float thrs, float wvlen, float pkmethod, float* rarvtm,
//	float signal[8 * 672], float tstart[8]) {
//	int x = 0;
//}


// #include "mex.h"     //Matlab调用C的函数
void fstbrkX(float slowness, float tt, int nsamp, int nR, float RR, float period,
	float wdwdth, float thrs, float wvlen, float pkmethod, float* rarvtm,
	float* signal, float* tstart )
{
	int i, j, lwn;
	int istrt;
	float* source, * wvstk, * sptr;
	float** sigary;
	float brktm, x1, x2, wvdc;
	lwn = (int)(wvlen / period);
	if (lwn > nsamp / 2) lwn = nsamp / 2;
	source = (float*)calloc(lwn + 3, sizeof(float));
	wvstk = (float*)calloc(lwn, sizeof(float));
	sigary = (float**)malloc(nR * sizeof(float*));
	for (i = 0; i < nR; i++)
		sigary[i] = (float*)calloc(lwn, sizeof(float));

	//Purpose: Form new array data across array along slowness moveout, remove DC ahead of first arrival//

	for (i = 0; i < nR; i++) {
		x1 = (tt + slowness * i * RR - tstart[i]) / period;
		istrt = (int)x1;
		if (istrt < 0) istrt = 0;
		if (istrt > lwn) istrt = lwn;
		x1 = x1 - istrt;
		wvdc = 0.0;
		for (j = 0; j <= istrt; j++) wvdc += signal[nsamp * i + j];
		wvdc = wvdc / (istrt + 1);
		for (j = 0; j < lwn; j++)
		{
			sptr = &signal[nsamp * i + istrt + j];
			sigary[i][j] = *sptr * (1.0 - x1) + *(sptr + 1) * x1 - wvdc;
			int temp = sigary[i][j];
		}
	}
	if (pkmethod > 1)
	{
		// Purpose: find first arrival for every receiver//
		for (i = 0; i < nR; i++)
		{
			for (j = 0; j < 3; j++) source[j] = 0;
			for (j = 0; j < lwn; j++)
			{
				source[j + 3] = sigary[i][j];
				int temp = source[j];
			}
			brktm = pkfstbrk(source, period, wdwdth, thrs);
			rarvtm[i] = brktm - 3 * period + tt + slowness * i * RR;
		}
	}
	else
	{

		/*Pick first arrival only on the stacked trace;
		use cross-correlation to find relative time difference between receivers */

		//Form stacked trace//   
		for (j = 0; j < lwn; j++)
		{
			source[j + 3] = 0;
			for (i = 0; i < nR; i++) source[j + 3] += sigary[i][j];
			source[j + 3] /= (float)nR;
		}

		//Weight data toward front//
		for (j = 0; j < lwn; j++)
		{
			x1 = (float)(lwn - j + 1) / (lwn + 1);
			x1 = x1 * x1;
			wvstk[j] = source[j + 3] * x1;
			for (i = 0; i < nR; i++) sigary[i][j] *= x1;
		}

		//Crosscorrelate to find relative delay among receivers//
		for (i = 0; i < nR; i++)
		{
			sptr = sigary[i];
			xcorr(lwn, wvstk, sptr, 5, &x1, &x2);
			rarvtm[i] = x2 * period;
		}
		// Calibrate to first arrival time and add slowness moveout//
		for (i = 0; i < nR; i++) rarvtm[i] = rarvtm[i] + tt + slowness * i * RR;

		brktm = pkfstbrk(source, period, wdwdth, thrs);

		for (i = 0; i < nR; i++) rarvtm[i] += brktm - 3 * period;
	}
	free(source);
	free(wvstk);
	for (i = 0; i < nR; i++) {
		free(sigary[i]);
		sigary[i] = NULL;
	}
	return;
}

float pkfstbrk(float* source, float period, float wdwdth, float thrs)
{

	int i, j, k, l, nk, np, nn_k, lwn;
	int istrt, iend;
	float tm, period_n, * e1;
	float v, v1, wv0, wv1, pt0, pt;
	float tpk, tzx, tstrt_n, dlt, dtzx;
	float dx, x1, x2, xmid, f, fmid, xroot;
	float* newdata;
	float wvdc;

	//period = controlP->data.trwav.period; /* Sample interval //
	lwn = (int)(wdwdth / period); //Window length (in samples)//
	//Find the first peak/trough time in mumber of samples//
	e1 = (float*)calloc(lwn, sizeof(float));
	np = lwn;
	tm = fstpktime(source, np, thrs); //峰值或谷值时间//
	np = (int)tm;
	if (np < 2) np = 2;

	//first, find the start and end time of the data around the picked peak.  Also condiser the possibility that ttp is the zero-crossing time //
	v = fabs(source[np]);
	if (source[np] * source[np + 1] > 0.0 || source[np] * source[np - 1] > 0.0)
	{      //ttp is a peak time //
		if (v < fabs(source[np - 1])) v = fabs(source[np - 1]);
		if (v < fabs(source[np - 2])) v = fabs(source[np - 2]);
		l = np;
		k = 0;
	}
	else { //ttp is a zero-crossing time, no picking //
		return (tm * period);
	}
	wv0 = source[np];
	i = 1;
	istrt = 0;
	while (i <= l) { //find rough first break left of peak //
		wv1 = *(source + l - i);
		v1 = fabs(wv1);
		if (v1 <= 0.1 * v || wv1 * wv0 <= 0.0 || v1 > fabs(wv0)) {
			istrt = i;
			break;
		}
		wv0 = wv1;
		i++;
	}

	if (istrt == 0)
	{
		return(tm * period);
	}
	else istrt++;

	//compute wave DC ahead of onset //
	wvdc = 0.0;

	/**	 k=0;  DC has been removed in data stacking stage
		 for(i=0; i<=l-istrt; i++)
		 {
			 wvdc += *(source+i);
			 k++;
		 }
		  if (k > 0) wvdc = wvdc/k;  测试之用，可删除**/

	wv0 = source[l] - wvdc;
	v = 0.0;
	i = 0;
	iend = 2;
	while (i < lwn - np) { //find zero-crossing right of peak//
		wv1 = *(source + np + i) - wvdc;
		if (wv1 * wv0 <= 0.0) {
			iend = i;
			break;
		}
		v1 = fabs(wv1);
		if (v <= v1) v = v1;
		wv0 = wv1;
		i++;
	}

	//ill in the data buffer//
	for (i = -istrt; i <= iend; i++) e1[i + istrt] = *(source + np + i) - wvdc;

	//Then upsample the data by a factor of 5//
	k = iend + istrt + 1;
	nn_k = 5 * k;
	period_n = period / 5.0;
	newdata = (float*)calloc(nn_k + 2, sizeof(float));
	resample(e1, k, newdata, nn_k);

	//find the peak time, zero-crossing time of the new data //
	v = 0.0;
	l = 1;
	wv0 = newdata[5 * istrt + 1];
	for (i = 0; i < nn_k; i++) {
		wv1 = newdata[i];
		v1 = fabs(wv1);
		if (v < v1 && wv1 * wv0 > 0.) {
			l = i;
			v = v1;
		}
	}

	//quadratic fitting of the three points around peak//
	wv0 = newdata[l];
	wv1 = newdata[l - 1];
	v1 = newdata[l + 1];
	v = (v1 - wv1) / (2. * wv0 - wv1 - v1);
	tpk = (l + 0.5 * v) * period_n;
	tstrt_n = period * (np - istrt);
	v = wv0 + (v1 - wv1) * v / 8.0;

	//normalize the wave data //
	for (i = 0; i < nn_k + 2; i++) newdata[i] /= v;

	//zero-crossing right of peak of new data //
	wv0 = 1.0;
	i = l;
	iend = l + 1;
	while (i < nn_k) {
		wv1 = newdata[i];
		if (wv0 * wv1 <= 0.0) {
			iend = i;
			break;
		}
		wv0 = wv1;
		i++;
	}

	wv0 = fabs(newdata[iend - 1]);
	wv1 = fabs(newdata[iend]);
	dtzx = wv0 / (wv0 + wv1) * period_n;
	tzx = (iend - 1) * period_n + dtzx;
	iend = iend - 1;

	//find a rough first break of the new data //
	wv0 = newdata[l];
	i = l;
	istrt = 0;
	while (i >= 0) {
		wv1 = newdata[i];
		v1 = fabs(wv1);
		if (wv1 * wv0 <= 0.0 || v1 <= 0.1 || v1 > fabs(wv0)) {
			istrt = i;
			break;
		}
		wv0 = wv1;
		i--;
	}

	dlt = tzx - tpk;
	pt0 = 2. * (tzx - istrt * period_n);

	//find the best fitting fuction parameter pt using bisection method //
	x1 = 0.9;
	x2 = 3.5;
	pt = x1 * pt0;
	f = ft_onset(pt, period_n, dlt, dtzx, newdata, istrt, iend);
	while (x2 > x1) {
		pt = x2 * pt0;
		fmid = ft_onset(pt, period_n, dlt, dtzx, newdata, istrt, iend);
		if (f * fmid < 0.0) break;
		x2 -= 0.2;
	}

	if (f * fmid >= 0.0) {
		xroot = 1.0;
		goto L400; // no root found, go to use the rough 1st break //
	}

	if (f < 0.0) {
		xroot = x1;
		dx = x2 - x1;
	}
	else {
		xroot = x2;
		dx = x1 - x2;
	}

	i = 1;
	while (i < 40) {
		dx = 0.5 * dx;
		xmid = xroot + dx;
		pt = xmid * pt0;
		fmid = ft_onset(pt, period_n, dlt, dtzx, newdata, istrt, iend);
		if (fmid <= 0.0) xroot = xmid;
		if (fabs(dx) < 0.001 || fmid == 0.0) break;
		i++;
	}     // bisection done //

	if (i > 38) xroot = 1.0;  // no root found //

L400:    pt = xroot * pt0;


	//picking on the fitted curve //
	i = 0;
	while (i < nn_k) {
		tm = i * period_n;
		xroot = func_onset(pt, tm, dlt);
		if (xroot >= 0.1) break;
		i++;
	}
	tm -= period_n;
	tm = tstrt_n + (tzx - (0.5 * pt - tm));
	free(e1);
	free(newdata);
	return(tm);
}

/*--------------------------------------------------------------
 * ft_onset:  快速拟合处理
 *
 * 目的: Fit the wave oneset with the function  using the least-squares method
 *
 * 输入参数:
 *     wvonset = array containing wave onset
 *     dt   = sampling period of wvonset
 *     pt   = fitted wave onset period T
 *     dlt  = time difference bewteen wave peak and zero-crossing right of peak
 *     dtzx = time difference tzx - (int)(tzx)
 *     istrt= start time index in wvonset
 *     iend = end time index in wvonset
 *
 * 输出参数:
 *     wave residue between wvonset and the theoretical function
 *
 * 需要调用函数:  func_onset，如下
 *----------------------------------------------------------------*/
float ft_onset(float pt, float period, float dlt, float dtzx,
	float* wvonset, int istrt, int iend)
{
	int i;
	float pt1, wvrsd, y0, y1, dy, tm, dpt;

	wvrsd = 0.0;
	for (i = iend, i >= istrt; i--;)
	{
		tm = 0.5 * pt - (dtzx + (iend - i) * period);
		y0 = func_onset(pt, tm, dlt);
		if (y0 == 0.0) continue;
		dpt = 0.001 * pt;
		pt1 = pt + dpt;
		y1 = func_onset(pt1, tm, dlt);
		dy = (y1 - y0) / dpt;
		wvrsd += dy * (y0 - wvonset[i]);
	}
	return(wvrsd);
}

/*-----------------------------------------------------------------------------------
 * func_onset:
 *
 * 目的: computes the wave oneset fitting function
 *
 * 输入参数:
 *     pt   = fitted wave onset period T
 *     dlt  = time difference bewteen wave peak and zero-crossing right of peak
 *     t    = time
 * 输出参数:
 *     function value
 *
 * 需要调用函数: 无
 *-------------------------------------------------------------------------------------*/
float func_onset(float pt, float t, float dlt)
{
	int i;
	float pi = 4. * atan(1.0), pi2;
	float t0, x, y, fval;
	double tmp, alfa;
	pi2 = pi + pi;
	t0 = 0.5 * pt - dlt;
	x = pi2 * t0 / pt;
	if (fabs(2. * x / pi - 1.) < 0.001)
	{
		fval = sin(pi2 * t / pt);
	}
	else {
		y = pi2 * t0;
		alfa = -y / tan(x) / pt;
		x = alfa * pt / y;
		tmp = t / t0;
		fval = pow(tmp, alfa) * sqrt(1. + x * x) * sin(pi2 * t / pt);
	}
	return(fval);
}


/*-----------------------------------------------------------------------------------
c  目的: routine to register positive and negative peaks in a waveform temp(nt)
c
c  输入参数:
c        *wvin = waveform of length ntpts
c        thrs = threshhold to pick, it's the fraction of the peak amplitude within npts
c 输出参数:
c        rgst = registred pick location in # of samples
*-------------------------------------------------------------------------------------*/

float fstpktime(float* wvin, int npts, float thrs)
{
	int nrg;
	int rgst[10] = {0};
	int j, j0, j1, jmx;
	float  amx, amin, alclmx, alclmin, cnst, ap0;

	amx = 0.;
	cnst = 1e20;
	amin = cnst;
	for (j = 0; j < npts; j++)
	{
		if (amx <= *(wvin + j)) amx = *(wvin + j);
		if (amin >= *(wvin + j)) amin = *(wvin + j);
	}

	alclmx = 0.0;
	alclmin = cnst;
	nrg = 0;

	//register peaks //
	for (j = 0; j < npts; j++)
	{
		ap0 = wvin[j];
		if (ap0 > thrs * amx && alclmx <= ap0)
		{
			alclmx = ap0;
			j0 = j;

			if (alclmin != cnst)
			{
				rgst[nrg] = j1;
				nrg = nrg + 1;

				if (nrg > 9) goto L100;
				alclmin = cnst;
			}
		}

		if (ap0 < thrs * amin && alclmin >= ap0)
		{
			alclmin = ap0;
			j1 = j;
			if (alclmx != 0.0) {
				rgst[nrg] = j0;
				nrg = nrg + 1;

				if (nrg > 9) goto L100;
				alclmx = 0.0;
			}
		}
	}

L100:    return(rgst[0]);
}

/*-----------------------------------------------------------------------------------
$START_DOC
void resample ( float *data1, int ndata1, float *data2,
				int ndata2, int *icode)
目的:  差值处理， Interpolate data using cubic spline.

参数说明:
 Name                 In/Out   Description
 -------------------  -------  ---------------------------------------
 data1                Input    Original data of size ndata1.
 ndata1               Input    Number of samples in data1.
 data2                Output   New data of size ndata2.
 ndata2               Output   Size of new data.

RETURN VALUE: 不返回值

METHOD: 重新抽样波形

NOTES:

  1) Include resample.h
	 #include "resample.h"
$END_DOC
 *************************************************************/

void   resample(float* data1, int ndata1, float* data2, int ndata2)  //define the resample function//
{
	float* x, * y, * b, * c, * d;
	int i, j;
	float fnw1, fnw2;
	float u;

	if (ndata1 == ndata2)
	{
		for (i = 0; i < ndata1; i++)
		{
			data2[i] = data1[i];
		}
		return;
	}

	if (ndata1 == 0)
	{
		for (i = 0; i < ndata2; i++)
		{
			data2[i] = 0;
		}
		return;
	}
	if (ndata1 == 1)
	{
		for (i = 0; i < ndata2; i++)
		{
			data2[i] = data1[0];
		}
		return;
	}

	x = (float*)calloc(sizeof(float), ndata1);
	y = (float*)calloc(sizeof(float), ndata1);
	b = (float*)calloc(sizeof(float), ndata1);
	c = (float*)calloc(sizeof(float), ndata1);
	d = (float*)calloc(sizeof(float), ndata1);

	fnw1 = (float)(ndata1 - 1);
	for (i = 1; i <= ndata1; i++)
	{
		x[i - 1] = (float)(i - 1) / fnw1;
		y[i - 1] = data1[i - 1];
	}

	//construct spline of the input data//
	SPline(x, y, ndata1, b, c, d);  //调用函数如下//
	//resample// 
	fnw2 = (float)(ndata2 - 1);
	for (i = 1; i <= ndata2; i++)
	{
		u = (float)(i - 1) / fnw2;
		data2[i - 1] = SPeval(x, y, ndata1, b, c, d, u);
	}
	free(x);
	free(y);
	free(b);
	free(c);
	free(d);
	return;
}

/******************************************************************
$START_DOC
void SPline(float *x,float *y, int n,float *b,float *c, float *d)

目的:    This function compute the natural spline polynomial coefficients given the n data points.

参数说明:

 Name                 In/Out   Description
 -------------------  -------  ---------------------------------------

 x                    Input    Data points for polynomial on x-axis.
 y                    Input    Data points for polynomial on y-axis.
 n                    Input    Total number of data point.
 b                    Input    Data points for polynomial.
 c                    Input    Data points for polynomial.
 d                    Input    Data points for polynomial.

RETURN VALUE: NO
NOTES:       #include "resample.h"
$END_DOC
 ******************************************************************/

void   SPline(float* x, float* y, int n, float* b, float* c, float* d)
{
	int i, j, nm1;
	float p, fi;
	i = 0;
	b[i] = 0.0;
	d[i] = 0.0;
	nm1 = n - 1;

	//set matrix //

	for (i = 1; i <= nm1; i++)
	{
		j = i + 1;
		c[i - 1] = 1.0 / (x[j - 1] - x[i - 1]);
		b[j - 1] = 3.0 * (y[j - 1] - y[i - 1]) * c[i - 1] * c[i - 1];
		d[j - 1] = 2.0 * c[i - 1];
		b[i - 1] = b[i - 1] + b[j - 1];
		d[i - 1] = d[i - 1] + d[j - 1];
	}
	//forward substitution//
	for (i = 2; i <= n; i++)
	{
		j = i - 1;
		p = c[j - 1] / d[j - 1];
		d[i - 1] = d[i - 1] - c[j - 1] * p;
		b[i - 1] = b[i - 1] - b[j - 1] * p;
	}
	//back substitution //
	b[n - 1] = b[n - 1] / d[n - 1];
	for (i = 1; i <= nm1; i++)
	{
		j = n - i;
		b[j - 1] = (b[j - 1] - c[j - 1] * b[j]) / d[j - 1];
	}
	//b now contains the constants b(i)//
	for (i = 1; i <= nm1; i++)
	{
		j = i + 1;
		fi = (y[j - 1] - y[i - 1]) * c[i - 1];
		d[i - 1] = (-2.0 * fi + b[j - 1] + b[i - 1]) * c[i - 1] * c[i - 1];
		c[i - 1] = (3.0 * fi - b[j - 1] - 2.0 * b[i - 1]) * c[i - 1];
	}
	return;
}

/******************************************************************
$START_DOC

float SPeval( float *x, float *y, int n, float *b,
			  float *c, float *d, float u)

目的:    This function is to evaluate the spline function  value at u (v=y(u)).

参数说明:

 Name                 In/Out   Description
 -------------------  -------  ---------------------------------------

 x                    Input    Data points for polynomial on x-axis.
 y                    Input    Data points for polynomial on y-axis.
 n                    Input    Total number of data point.
 b                    Input    Data points for polynomial.
 c                    Input    Data points for polynomial.
 d                    Input    Data points for polynomial.
 u                    Input    Data points for polynomial.

RETURN VALUE:  Returns a float value.
NOTES:      #include <MyHeaderFile.h>
$END_DOC
 ******************************************************************/
float  SPeval(float* x, float* y, int n, float* b, float* c, float* d,
	float u)
{
	int i, j, k;
	float z, w, v;
	i = 1;
	j = n;
L10:	if (j <= (i + 1)) goto L50;
	k = (i + j) / 2;

	if ((u - x[k - 1]) < 0) goto  L20;
	if ((u - x[k - 1]) == 0) goto  L40;
	if ((u - x[k - 1]) > 0) goto  L30;
L20:    j = k;
	goto L10;
L30:    i = k;
	goto L10;
L40:    i = k;
L50:    z = u - x[i - 1];

	w = z * d[i - 1];
	v = c[i - 1] + w;
	w = v + w;
	v = b[i - 1] + z * v;
	w = v + z * w;
	v = y[i - 1] + z * v;
	return (v);
}

//-------------------------------------------- function cross -------------------------------------------- */

void xcorr(int np, float* x, float* y, int nxmax,
	float* xcormax, float* xlag)
	/***
		   correlation *pxcormax, lag *plag, and reflection coefficient 	pscale[].

	*    y    -  time series 2 for correlation
   *    np   - number of samples in both x and y
   *    xcor - correlaiton function between x and y
   *    plag - lag of x relative to y
   *    pxcormax - Maximum crosscorrelation coefficient
   *    nxmax  - number of correlation steps in positive and negative shifts***/

{
	int    i, i0, j, npts, nxcor, n;
	int    xindex, yindex, length;
	float  xsum, ysum, tot, xnorm;
	float xx[6], newdata[32], * xcor;
	double prod;

	npts = np;
	nxcor = 2 * nxmax;
	xcor = (float*)calloc(nxcor + 1, sizeof(float));
	xsum = 0.0;
	ysum = 0.0;
	tot = 0.0;
	for (i = 0; i < npts; i++) {
		xsum += (x[i] * x[i]);
		ysum += (y[i] * y[i]);
	}

	prod = xsum * ysum;
	xnorm = sqrt(prod);
	for (i = 0; i < nxcor; i++) xcor[i] = 0.0;
	for (i = 0; i < nxmax; i++) {
		xindex = 0;
		yindex = nxmax - i - 1;
		length = npts - nxmax + i + 1;
		tot = 0.0;
		for (j = 0; j < length; j++)
		{
			tot = tot + x[xindex] * y[yindex];
			xindex++;
			yindex++;
		}
		xcor[i] = tot / xnorm;
	}

	for (i = nxmax; i < nxcor; i++) {
		yindex = 0;
		xindex = i - nxmax + 1;
		length = npts - (i - nxmax + 1);
		tot = 0.0;
		for (j = 0; j < length; j++)
		{
			tot = tot + x[xindex] * y[yindex];
			xindex++;
			yindex++;
		}
		xcor[i] = tot / xnorm;
	}

	*xcormax = 0.0;
	i0 = 2;
	for (i = 2; i < nxcor; i++)
	{
		if (*xcormax <= xcor[i])
		{
			*xcormax = xcor[i];
			i0 = i;
		}
	}

	npts = 6;
	for (i = 0; i < npts; i++) xx[i] = xcor[i0 - 2 + i];
	/** Then upsample the data by a factor of 5 **/

	n = 5 * npts;

	// newdata = (float *) calloc ( nn+2, sizeof ( float) );
	resample(xx, npts, newdata, n);

	*xcormax = 0.0;
	for (i = 0; i < n; i++)
	{
		if (*xcormax <= newdata[i])
		{
			*xcormax = newdata[i];
			*xlag = i / 5.0;
		}
	}

	*xlag = i0 - 2 - nxmax + 1 + *xlag;

	free(xcor);

	return;
}

//
////The gateway routine //
//void mexFunction( int nlhs, mxArray *plhs[],
//                  int nrhs, const mxArray *prhs[])
//{
//  int    i, nt,nr;
//  double *sig, *tst, *arvtime;
//  float  dt, tt, tlength, rr, slns, thrs, wvlen, pkmethod, *fstm;
//  float  *sigin;
//  
////Check for proper number of arguments//
//   
//  if(nrhs != 12) 
//    mexErrMsgTxt("yi-shi-san (12) inputs required.");
//   if(nlhs != 1) 
//     mexErrMsgTxt("one output required.");
//  
////Check and read-in input data/ /
//
////Get the scalar input  //
// 
//  slns     = mxGetScalar(prhs[2]);
//  tt       = mxGetScalar(prhs[3]);
//  nt       = mxGetScalar(prhs[4]);
//  nr       = mxGetScalar(prhs[5]);
//  rr       = mxGetScalar(prhs[6]);
//  dt       = mxGetScalar(prhs[7]);
//  tlength  = mxGetScalar(prhs[8]);
//  thrs     = mxGetScalar(prhs[9]);
//  wvlen    = mxGetScalar(prhs[10]);
//  pkmethod = mxGetScalar(prhs[11]);
//
//  
////Create a pointer to the input matrix sig //
//   
////Create a pointer to the input start time array//
//  tst = mxGetPr(prhs[1]);
////Create a pointer to the input matrix sig //
//  nt = mxGetM(prhs[0]);
//  nr = mxGetN(prhs[0]);
//  sig = mxGetPr(prhs[0]);
// //printf("see  %d  %f %f\n", nt, dt, tlength); 
// {
//   int i;
//    sigin = (float *) calloc(nt*nr, sizeof(float));
//    for (i=0; i<nr*nt; i++)  *(sigin + i)=(float)sig[i];
//  }
//
////Create a C pointer to a copy of the output matrix. //
//  plhs[0] = mxCreateDoubleMatrix(1,nr, mxREAL);
//  arvtime = mxGetPr(plhs[0]);  
//
//// Call the C subroutine//
//  fstm  = (float *) calloc (nr, sizeof (float));
//
//  fstbrk( sigin, tst, slns, tt,  nt,  nr, rr, dt, tlength,  thrs, wvlen, pkmethod, fstm );
//
// //printf("returned  %f %f\n", fstm[0], fstm[nr-1]); 
//
//  for(i=0; i<nr; i++) arvtime[i] = fstm[i];
//  free(sigin);
//  free(fstm);
//  return;
//}
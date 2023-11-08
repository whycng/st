#pragma once
 
//void test1(float slowness, float tt, int nsamp, int nR, float RR, float period,
//    float wdwdth, float thrs, float wvlen, float pkmethod, float* rarvtm,
//    float signal[8 * 672], float tstart[8]);--tmp
void fstbrkX(float slowness, float tt, int nsamp, int nR, float RR, float period,
    float wdwdth, float thrs, float wvlen, float pkmethod, float* rarvtm,
    float* signal, float* tstart);


float pkfstbrk(float* source, float period, float wdwdth, float thrs);
float ft_onset(float pt, float period, float dlt, float dtzx,
    float* wvonset, int istrt, int iend);
float func_onset(float pt, float t, float dlt);
float fstpktime(float* wvin, int npts, float thrs);

//extern void qip(int *n, float *x, float *y, float *u, float *f);
float  SPeval(float* x, float* y, int n, float* b, float* c, float* d,
    float u);
void   SPline(float* x, float* y, int n, float* b, float* c, float* d);
void   resample(float* data1, int ndata1, float* data2, int ndata2);
void xcorr(int np, float* x, float* y, int nxmax, float* xcormax, float* plag);
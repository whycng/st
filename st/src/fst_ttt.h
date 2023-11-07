#pragma once
#include <fstream>
#include <malloc.h>
#include <iostream>
#include <sstream>

// TFWV01.dat头部信息
struct tfwv_firstLine {
	double depth_begin;
	double depth_end;
	int nR;//8个接收器
	int nsamp;//672个数据
	double period;//us -- 每个接收点之间的时间间隔是12us 0.000012
	double RR; 
	/*RR = input('输入接收器间距, R-to-R spacing (e.g., 0.5ft & 0.1524m)');*/
	double RR2;// 深度间距
	double TR;// 发射器T距离第一个接收器R的距离 3.81 (tfw_a改为TR)
};

class Fst_ttt
{
public:
	double* m_tmogrm;// 输出
	double* m_Rtmogrm; 
	double* m_ttfit;
	double m_ssfit;
	double m_ss0;
	double m_ssd;
	double m_dop;
	// 内部数据
private:
	bool b_Caliper = false;
	tfwv_firstLine tfwv_f;// TFWI 头部信息

	float** Aals;//测试
	float* DTC;// slowness.txt DTC->slowness
	float* m_dept;// 深度值
	float* TTAV_tt; 
	float** m_allData;// 每个深度点的所有数据
	float** m_rarvtm_all;//所有的TTT
	float* m_rarvtm;// fstbrk的输出，用于tttomo的输入
	float* m_TT1;
	float* m_tstart;

	double* m_ttpick;// tttomo的输入，实际就是m_rarvtm
	double* m_ttslns;// 0
	double* Caliper;//井径

	double m_dtf;
	double m_dr;/* dr = Rdis/60; */
	double m_Vav;
	float m_pkmethod;
	float m_wdwdth;/*wdlen = input('For first arrival picking, enter window length (us)');
	%2-5个波长周期就行，200-500，是调试*/
	float m_thrs;/* 百分数？0.01 0.02 --? 这个需要测试 */
	float m_wvlen;
	float m_tlod_diaTOOL;/*两者是一回事，有可能单位不统一，现在改为英寸*/
	int m_num_depth;// 深度个数 
	int m_size_tmp;
	int m_MFilter;//中值滤波参数 3?


#pragma region 原始测试数据
	//原始测试数据
//float ret[100];
//float scard[8] = { 0 };//
//float depth;
//float data[8 * 672] = { 0 };// 此数据应由read_TFWV_dat读取到每一个点的数据
//// char cc[50]; 
//float t;
//float tstart[8];
//float slowness;// = 120.0f;
//float tt;//= 0.0f;
//int nsamp;// = 672;
//int nR;// = 8;
//float TR;/* 外部数据, 3.81 
//TR之间的距离一般也就是3.8-4.2米之间，不同公司的仪器这个距离不同
//TR = input('输入源距, T-to-R spacing (e.g., 12ft & 3.6576m)');
//*/
//float RR;// = 0.1524f;
//float period;//= 12.0f;
//float wdwdth;// = 300.0f;
//float thrs;// = 0.5f;// 150.0f;  --?
//float wvlen;// = 1200.0f;
//float pkmethod;//= 2.0f;
//float* rarvtm;// = (float*)alloca(nR * sizeof(float)); // 长度为 nR
//// tttomo
//FILE* fp1;

//// tttomo返回的数值
//int size_tmp;//= 120;//tmogrm[60 + j] = 0.7 * vav;   120--对称的
//double* tmogrm;//= (double*)alloca(size_tmp * sizeof(double));// 60是因为tmogrm[60 + j] = 0.7 * vav;
//double* Rtmogrm;//= (double*)alloca(size_tmp * sizeof(double));
//double* ttfit;// = (double*)alloca(nR * sizeof(double));//   ttfit[i] = tvtm(V); 
//double ssfit;  //	*ssfit = 1.0e6 / V;
//double ss0; //*ss0 = 1e6 / V0;
//double ssd; //*ssd = 1e6 / Vdeep;;
//double dop; //*dop = (Vdeep - V0) / gv;;
//// tttomo函数需要的值
//double* ttpick;// = (double*)alloca(nR * sizeof(double)); // nR 为数组长度
//double* ttslns;// = (double*)alloca(nR * sizeof(double));
//double Calpr;//  
//double dtf;//起步200;
//double tlod;//0.08;
//double dr;// 
//double vav;/*vav这个就不好说了，要看实际地层，可以从slowness里面找一个数就行，
//那个数的分之一再乘以1000000就是速度的单位了
//1/155.9999 * 1000000    */
//int index_t;//= 0;
//int index_tk;//= 0;

#pragma endregion

public:
	// 方法
	void ft_main();
	void read_DTC(const char* filename , std::string mode);
	void read_TFWV_dat(const char* filename );
	//void handle_fst(const char* outFielName);
	void handle_fst(std::string& path,
		std::string file_tmogrm = "/tmogrm.txt",
		std::string file_Rtmogrm = "/Rtmogrm.txt",
		std::string file_others = "/others.txt");
	void init_par();
	void compute_par();
	int GetSamp();
	void test_read_real_fstbrk(const char* filename);
};

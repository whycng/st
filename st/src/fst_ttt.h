#pragma once
#include <fstream>
#include <malloc.h>
#include <iostream>
#include <sstream>

// TFWV01.datͷ����Ϣ
struct tfwv_firstLine {
	double depth_begin;
	double depth_end;
	int nR;//8��������
	int nsamp;//672������
	double period;//us -- ÿ�����յ�֮���ʱ������12us 0.000012
	double RR; 
	/*RR = input('������������, R-to-R spacing (e.g., 0.5ft & 0.1524m)');*/
	double RR2;// ��ȼ��
	double TR;// ������T�����һ��������R�ľ��� 3.81 (tfw_a��ΪTR)
};

class Fst_ttt
{
public:
	double* m_tmogrm;// ���
	double* m_Rtmogrm; 
	double* m_ttfit;
	double m_ssfit;
	double m_ss0;
	double m_ssd;
	double m_dop;
	// �ڲ�����
private:
	bool b_Caliper = false;
	tfwv_firstLine tfwv_f;// TFWI ͷ����Ϣ

	float** Aals;//����
	float* DTC;// slowness.txt DTC->slowness
	float* m_dept;// ���ֵ
	float* TTAV_tt; 
	float** m_allData;// ÿ����ȵ����������
	float** m_rarvtm_all;//���е�TTT
	float* m_rarvtm;// fstbrk�����������tttomo������
	float* m_TT1;
	float* m_tstart;

	double* m_ttpick;// tttomo�����룬ʵ�ʾ���m_rarvtm
	double* m_ttslns;// 0
	double* Caliper;//����

	double m_dtf;
	double m_dr;/* dr = Rdis/60; */
	double m_Vav;
	float m_pkmethod;
	float m_wdwdth;/*wdlen = input('For first arrival picking, enter window length (us)');
	%2-5���������ھ��У�200-500���ǵ���*/
	float m_thrs;/* �ٷ�����0.01 0.02 --? �����Ҫ���� */
	float m_wvlen;
	float m_tlod_diaTOOL;/*������һ���£��п��ܵ�λ��ͳһ�����ڸ�ΪӢ��*/
	int m_num_depth;// ��ȸ��� 
	int m_size_tmp;
	int m_MFilter;//��ֵ�˲����� 3?


#pragma region ԭʼ��������
	//ԭʼ��������
//float ret[100];
//float scard[8] = { 0 };//
//float depth;
//float data[8 * 672] = { 0 };// ������Ӧ��read_TFWV_dat��ȡ��ÿһ���������
//// char cc[50]; 
//float t;
//float tstart[8];
//float slowness;// = 120.0f;
//float tt;//= 0.0f;
//int nsamp;// = 672;
//int nR;// = 8;
//float TR;/* �ⲿ����, 3.81 
//TR֮��ľ���һ��Ҳ����3.8-4.2��֮�䣬��ͬ��˾������������벻ͬ
//TR = input('����Դ��, T-to-R spacing (e.g., 12ft & 3.6576m)');
//*/
//float RR;// = 0.1524f;
//float period;//= 12.0f;
//float wdwdth;// = 300.0f;
//float thrs;// = 0.5f;// 150.0f;  --?
//float wvlen;// = 1200.0f;
//float pkmethod;//= 2.0f;
//float* rarvtm;// = (float*)alloca(nR * sizeof(float)); // ����Ϊ nR
//// tttomo
//FILE* fp1;

//// tttomo���ص���ֵ
//int size_tmp;//= 120;//tmogrm[60 + j] = 0.7 * vav;   120--�ԳƵ�
//double* tmogrm;//= (double*)alloca(size_tmp * sizeof(double));// 60����Ϊtmogrm[60 + j] = 0.7 * vav;
//double* Rtmogrm;//= (double*)alloca(size_tmp * sizeof(double));
//double* ttfit;// = (double*)alloca(nR * sizeof(double));//   ttfit[i] = tvtm(V); 
//double ssfit;  //	*ssfit = 1.0e6 / V;
//double ss0; //*ss0 = 1e6 / V0;
//double ssd; //*ssd = 1e6 / Vdeep;;
//double dop; //*dop = (Vdeep - V0) / gv;;
//// tttomo������Ҫ��ֵ
//double* ttpick;// = (double*)alloca(nR * sizeof(double)); // nR Ϊ���鳤��
//double* ttslns;// = (double*)alloca(nR * sizeof(double));
//double Calpr;//  
//double dtf;//��200;
//double tlod;//0.08;
//double dr;// 
//double vav;/*vav����Ͳ���˵�ˣ�Ҫ��ʵ�ʵز㣬���Դ�slowness������һ�������У�
//�Ǹ����ķ�֮һ�ٳ���1000000�����ٶȵĵ�λ��
//1/155.9999 * 1000000    */
//int index_t;//= 0;
//int index_tk;//= 0;

#pragma endregion

public:
	// ����
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

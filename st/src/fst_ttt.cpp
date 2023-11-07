 
#include <iomanip>
#include "fst_ttt.h"
#include "fstbrk.h"
#include "tttomo.h" 
#include "base.h"

#pragma region ԭʼ���Ժ���
//void Fst_ttt::ft_main()
//{
//    //
//    //slowness = 120.0f;
//    // tt = 0.0f;
//    // nsamp = 672;// 672������
//    // nR = 8; // 8��������
//    // RR = 0.1524f; // ���������
//    // period = 12.0f;// ��12us ?
//    // wdwdth = 300.0f;
//    // thrs = 0.5f;// 150.0f;  --?
//    // wvlen = 1200.0f;
//    // pkmethod = 2.0f;
//    // rarvtm = (float*)alloca(nR * sizeof(float)); // ����Ϊ nR
//    //// tttomo���ص���ֵ
//    //  size_tmp = 120;//tmogrm[60 + j] = 0.7 * vav;  ��������µ�120
//    // tmogrm = (double*)malloc(size_tmp * sizeof(double));// 60����Ϊtmogrm[60 + j] = 0.7 * vav;
//    // Rtmogrm = (double*)malloc(size_tmp * sizeof(double));
//    //  ttfit = (double*)malloc(nR * sizeof(double));//   ttfit[i] = tvtm(V); 
// 
//    //// tttomo������Ҫ��ֵ
//    //  ttpick = (double*)malloc(nR * sizeof(double)); // nR Ϊ���鳤��
//    //  ttslns = (double*)malloc(nR * sizeof(double));
//    // Calpr = 0.2;
//    // TR = 3.81; // 3.8 -- 4.2
//    // dtf = 200;  // 5����
//    // tlod = 0.08;// --> 0.09
//    // dr = 0.015;// ??? ��İ� 5/60=0.0833333
//    // vav = 1 / 155.9999 * 1000000;
//    //  index_t = 0;
//    //  index_tk = 0;
//
//    
//     
//
//    //"E:\Proj\vsProj\st_FileSave\wtest.dat"
//    //fopen_s(&fp1, "E:\\Proj\\vsProj\\st\\wtest.dat", "rb");
//    fopen_s(&fp1, "E:\\Proj\\vsProj\\st_FileSave\\wtest.dat", "rb");
//
//    for (int i = 0; i < 8; i++)
//    {
//        tstart[i] = 0;
//    }
//    if (fp1 == NULL) {
//        puts("���ļ�ʧ��");
//        return ;
//    }
//
//    fseek(fp1, 0, SEEK_SET);
//    while (fscanf_s(fp1, "%f", &t) != EOF) {
//        //process data
//            if (index_t < 8)
//            {
//                scard[index_t] = t;
//            }
//            else if (index_t == 8)
//            {
//                depth = t;
//            }
//            else
//            {
//                data[index_tk] = t;
//                index_tk++;
//            }
//        index_t++;
//    }
//    fclose(fp1);
//
//  /*  float m = 120.0f;
//    float ma[3] = { 2.3f, 3.4f, 4.5f };*/
//
//    //����ֵrarvtm
//    /*
//    tshft=(1:1:8);   1��2��3 ... 8
//    tshft = fstbrk(...)
//    */
//    fstbrkX(slowness, tt, nsamp, nR, RR, period, wdwdth, thrs, wvlen, pkmethod, rarvtm, data, tstart);
//
//    for (int i = 0; i < nR; i++)
//    {
//        //ttpick[i] = 0;
//        ttpick[i] = rarvtm[i]; //ttpick = rarvtm; 
//        ttslns[i] = 0;// ��0��Ӱ��
//    }
//
//    /*
//    ttpick = rarvtm�� ttslns=0�� nR=8��TR=3.81 ,  RR=0.1524, Calpr=0.015, dtf=200�𲽣�
//    ���ܸ�5�������˼����ˮ���ٶȵĵ�������Լ��200us/ft, �����us/m�ĵ�λ�Ǿ���
//    200us/0.3048��=656us/m����slowness��λһ���ˡ�
//    tlod = 0.08m�� slowness = 120�� dr=0.02�� vav=1/155.9999 * 1000000 
//    */
//    tttomo(ttpick, ttslns, nR, TR, RR, Calpr, dtf,
//        tlod, slowness, dr, vav,
//        tmogrm, Rtmogrm, ttfit, &ssfit,
//        &ss0, &ssd, &dop);
//      
//    for (int i = 0; i < size_tmp; i++)
//    {
//        std::cout << "  tmogrm[" << i << "]:" << tmogrm[i]
//            << "  Rtmogrm[" << i << "]:" << Rtmogrm[i] << std::endl;
//    }
//
//    for (int i = 0; i < nR; i++)
//    {
//        std::cout << "  ttfit[" << i << "]:" << ttfit[i] << std::endl;
//    }
//
//    std::cout << " ssfit:" << ssfit << " ,"
//        << " ss0:" << ss0 << " ,"
//        << " ssd:" << ssd << " ,"
//        << " dop:" << dop << std::endl;
//
//
//    printf("\n---------���rarvtm--------\n" );
//    for (int i = 0; i < nR; i++)
//    {
//        printf("rarvtm[%d]: %f\n", i, rarvtm[i]);
//    }
//}

#pragma endregion


// ��ȡSlowness.txt����ȡ����DTC���� slowness�� mode=Caliper ʱЯ����ֱ������
void Fst_ttt::read_DTC(const char* filename, std::string mode = "")
{
    std::ifstream ifs(filename);
    // ��� ifs����
    if (!ifs)
    {
        std::cout << "[ERROR]<Fst_ttt::read_DTC> ifs NOT EXIST" << std::endl;
        return;
    }
    std::string line;
    int line_count = 0;
    int index = 0;  
    //float* dept;
    //float* DTC;
    std::ifstream file(filename);
    while (getline(file, line)) {
        line_count++;
    }
    file.close();
    file.open("filename");
    line_count--;
    //*num_depth = line_count;
    m_num_depth = line_count;
    std::cout << "[Message]line_count:" << line_count << std::endl;
    m_dept = (float*)malloc(line_count * (sizeof(float)));
    DTC = (float*)malloc(line_count * (sizeof(float)));
    Caliper = (double*)malloc(line_count * (sizeof(double)));
    if (!mode.compare("Caliper"))//���Я��������
    {
        std::cout << "[Message] Caliper" << std::endl;
        b_Caliper = true;
        
    }

    std::getline(ifs, line);
    std::cout << std::fixed << std::setprecision(4);
    while (std::getline(ifs, line)) {
        std::stringstream word(line);  /*�и������׵Ĵ���������ɶ*/
        word >> m_dept[index];
        word >> DTC[index]; 
        if (!mode.compare("Caliper"))
        {
            word >> Caliper[index];
            //Caliper[index] *= 0.0254;/*Ϊ�˷��㣬��Ӣ��in תΪ��m*/
        }
       /* std::cout << "dept[index]:" << dept[index] << std::endl
            << DTC[index] << std::endl;*/
      
        index++;
    } 
}
 
// ��ȡ
void Fst_ttt::read_TFWV_dat(const char* filename )
{
    std::ifstream ifs(filename);
    // ��� ifs����
    if (!ifs)
    {
        std::cout << "[ERROR]<Fst_ttt::read_DTC> ifs NOT EXIST" << std::endl;
        return;
    }
    std::string line; 
    int index = 0;  
    
    m_allData = (float**)malloc(m_num_depth * sizeof(float*));
    //std::cout << "[Message]" << "num_depth:" << num_depth << std::endl;
    // ��һ��
    std::getline(ifs, line); 
    std::stringstream word(line);
    word >> tfwv_f.depth_begin >> tfwv_f.depth_end >> tfwv_f.nR >>
        tfwv_f.nsamp >> tfwv_f.period >> tfwv_f.RR >> tfwv_f.RR2 >>
        tfwv_f.TR;

    tfwv_f.period *= 1000000; /*   6λ ��������*/

    std::cout << "[Message]" << "tfwv_f.depth_begin:" << tfwv_f.depth_begin
        << "     tfwv_f.depth_end:" << tfwv_f.depth_end
        << "      tfwv_f.nR:" << tfwv_f.nR 
        << "      tfwv_f.nsamp:" << tfwv_f.nsamp
        << "      tfwv_f.period * 100000�Ľ��:" << tfwv_f.period
        << "      tfwv_f.RR:" << tfwv_f.RR
        << "      tfwv_f.RR2:" << tfwv_f.RR2
        << "      tfwv_f.TR:" << tfwv_f.TR 
        << std::endl;

    while (!ifs.eof()) {
       
       /* if (index == 1)
            return;*/

        ifs >> line;// ���,���� 
        //m_allData[index] = (float*)malloc(8 * tfwv_f.nsamp * sizeof(float));
        m_allData[index] = (float*)malloc(tfwv_f.nR * tfwv_f.nsamp * sizeof(float));
        int index_in = 0;
        while (index_in < tfwv_f.nsamp * tfwv_f.nR)
        {
            // 8��������������� 8 * 672�����ݣ�һ����ȵ�
            ifs >> m_allData[index][index_in];
            /*std::cout << "[Message]" << "m_allData[" << index << "][" << index_in << "]:" << m_allData[index][index_in]
                << std::endl;*/
            index_in++;
        }
        if( index % 20 == 0)
            std::cout << "[Message]idnex:" << index <<  "       depth-���:" << line << std::endl;
          //  << "    index_in:" << index_in <<  std::endl;

 

        index++;
    }
}

// ������Ҫ����Ĳ���
void Fst_ttt::compute_par()
{
    double s2r;
    double SAV;
    int iavg; 
    int jstrt;
    int jend;
    float* TT0 = (float*)malloc(m_num_depth * sizeof(float));
    /* ��һά��ʾ��ά:   nR1���֣�nR2���֣�...*/
    m_TT1 = (float*)malloc(m_num_depth * tfwv_f.nR * sizeof(float));
    /*TTAV_tt = (float*)malloc(m_num_depth * sizeof(float));*/
    for (int i = 0; i < tfwv_f.nR; i++)
    {
        s2r = tfwv_f.TR + i * tfwv_f.RR;/*Դ�㵽�������ľ���*/
         
        iavg = s2r / abs(tfwv_f.RR2);/* fix��βȡ��*/
        /*ddep = XTF.Curve(nWave).C_SPACING;  Ӧ���ǲ������*/
        /*ddep = ddep/0.3048;  %ת��ΪӢ��,Ҳ���Բ�ת��������ʵ���������*/

        /*SSS=zeros(1,EndRecord - StartRecord + 1);*/
        /*SSS(i)=slowness;   %��ȡ�õ���ʱ��ݲ��������P��
        SSS ���� DTC*/
        for (int j = 0; j < m_num_depth; j++)
        {
            if (j == 213 && i == 7)
                int tp = 3;
            jstrt = j ;/*jstrt = max(1,j+ishft); ishft=0*/
            jend =  m_num_depth  < (j + iavg) ? m_num_depth  : (j + iavg); /*jend = min(length(SSS),j+ishft+iavg);*/
            jend--;/* �±�ע�� */
            /*SAV = mean(SSS(jstrt:jend));*/
         /*   if (j == 327)
                int tp = 3;*/
            SAV = abs(mean_d(DTC, jstrt, jend));
            /*TT0(j) = (Calpr(j)-diaTOOL)*dtF.*sqrt(1-SAV.^2/dtF^2); */
            TT0[j] = ( Caliper[j] - m_tlod_diaTOOL) * m_dtf * sqrt(1.0 - pow(SAV / m_dtf, 2));
            //double a = (Caliper[j] - m_tlod_diaTOOL);
            //double b = sqrt(1.0 - pow(SAV / m_dtf, 2));/*m_dtf ֻ��200����SAV�����ܴ���200*/
            //double c = a * b * m_dtf;
            /*m_TT1(j,i)= s2r*SAV; %Purely based on slowness*/
            
            m_TT1[i * m_num_depth + j] = s2r * SAV;
        } 
    }

    for (int i = 0; i < m_num_depth; i++)
    {
        /*TTAV = m_TT1(i,1) + TT0(i);*/
        if (i >= 213)
            int tp = 3;
        double x1 = m_TT1[i];
        double x2 = TT0[i];
        TTAV_tt[i] = m_TT1[i] + TT0[i];
    }
}

// �ڶ�ȡ������֮�󣬲��ֲ�����Ҫ��ʼ��
void  Fst_ttt::init_par()
{
    m_MFilter = 11;
    /*MFilter = 11;      
XFactor = 3; 
YFactor = 3;   
light_filter = 1;  
NFilter = 7;  */

    m_dtf = 650;// us/m
    m_wdwdth = 300.0f; //300
    m_thrs = 0.01f;//  �ٷ�����0.01 0.02 --? 0.5
    m_wvlen = 1200.0f;// 800 �� 1000 ������Χ us -- matlab�ｨ��1200
    m_pkmethod = 1.0f;// 1 ÿ��������������matlab
    m_rarvtm = (float*)malloc(tfwv_f.nR * sizeof(float)); // ����Ϊ nR
    m_rarvtm_all = (float**)malloc(m_num_depth * sizeof(float*));
    // tttomo���ص���ֵ
    m_size_tmp = 121;//  ��121
    m_tmogrm = (double*)malloc(m_size_tmp * sizeof(double));// 60����Ϊtmogrm[60 + j] = 0.7 * vav;
    m_Rtmogrm = (double*)malloc(m_size_tmp * sizeof(double));
    m_ttfit = (double*)malloc(tfwv_f.nR * sizeof(double));//   ttfit[i] = tvtm(V); 
    // tttomo������Ҫ��ֵ
    m_ttpick = (double*)malloc(tfwv_f.nR * sizeof(double)); // nR Ϊ���鳤��
    m_ttslns = (double*)malloc(tfwv_f.nR * sizeof(double));
    // diaTOOL ���� tlod
    m_tlod_diaTOOL = 0.09;//
    /* x39.37תΪӢ�磬/12תΪӢ�� 
     diaTOOL = input('������Դ/������ֱ��(inches), enter tool source/receiver diameter (inches)');  %����ʵ������ֱ��ת��ΪӢ�磬����8cm=3.1496in
    diaTOOL = diaTOOL/12.0;   %��ת��ΪӢ�� */
    TTAV_tt = (float*)malloc(m_num_depth * sizeof(float));

    m_dr = 5.0f / 60.0f; 

    m_tstart = (float*)malloc(tfwv_f.nR * sizeof(float));
    memset(m_tstart, 0, tfwv_f.nR * sizeof(float)); 

    /*��λת��*/
    // תΪӢ�� slowness
    for (int i = 0; i < m_num_depth; i++)
    {
        // Ӧ��û�仯�ģ������ԭ��û��˵����
         DTC[i] *= 0.3048; // �Ƚ�ȷ����ת��--- ԭ�����double b = sqrt(1.0 - pow(SAV / m_dtf, 2));
         /* SAV ����С�� m_dtf�� SAV = abs(mean_d(DTC, jstrt, jend)); */
         
        // ע�͵�����
        //Caliper[i] /= 0.3048; // ��tttomo����ת
        Caliper[i] *= 0.08333; // Ӣ��תӢ��
    }
    tfwv_f.RR /= 0.3048; //Ӣ��
     
    tfwv_f.RR2 /= 0.3048;  
    tfwv_f.TR /= 0.3048;
    m_dtf *= 0.3048;  // us/m -- us/ft  //  
    m_tlod_diaTOOL /= 0.3048;

    // ��Ҫת��ߴ�����
    /*Vav = 1e6/mean(SSS(find(SSS))); %find()���ҷ���Ԫ�ص�������ֵ*/
    double sum = 0;
    int count = 0;
    for (int i = 0; i < m_num_depth; i++) {
        if (DTC[i] != 0) {
            sum += DTC[i];
            count++;
        }
    }
    double mean = sum / count;
    m_Vav = 1e6 / mean;

    if (b_Caliper)//���������
    {
        compute_par();
    }
    else
    { 
        memset(TTAV_tt, 0, tfwv_f.nR * sizeof(float)); 
        for (int i = 0; i < m_num_depth; i++)
            Caliper[i] = 0.2f;
    }


#pragma region ԭʼ��������
    //// ԭʼ���Բ���--���Ե����ʱ��ֵ
    //slowness = 120.0f;
    //tt = 0.0f;// xx
    //nsamp = 672;// 672������
    //nR = 8; // 8��������
    //RR = 0.1524f; // ���������
    //period = 12.0f;// ��12us ? ������֮��ʱ���
    //wdwdth = 300.0f;
    //thrs = 0.5f;// 150.0f;  --?
    //wvlen = 1200.0f;
    //pkmethod = 2.0f;
    //rarvtm = (float*)alloca(nR * sizeof(float)); // ����Ϊ nR
    //// tttomo���ص���ֵ
    //size_tmp = 120;//tmogrm[60 + j] = 0.7 * vav;  ��������µ�120
    //tmogrm = (double*)malloc(size_tmp * sizeof(double));// 60����Ϊtmogrm[60 + j] = 0.7 * vav;
    //Rtmogrm = (double*)malloc(size_tmp * sizeof(double));
    //ttfit = (double*)malloc(nR * sizeof(double));//   ttfit[i] = tvtm(V); 

    //// tttomo������Ҫ��ֵ
    //ttpick = (double*)malloc(nR * sizeof(double)); // nR Ϊ���鳤��
    //ttslns = (double*)malloc(nR * sizeof(double));
    //Calpr = 0.2; // ��λӢ�磿
    //TR = 3.81; // 3.8 -- 4.2
    //dtf = 650;  // 5����   200us/ft  650us/m
    //tlod = 0.09;// m  ������ֱ��
    //dr = 0.015; // 1/60 
    //vav = 1 / 155.9999 * 1000000;
    //i = 0;// ��Ҫ�޸�
    //k = 0;
#pragma endregion

    
}

// ��ȡÿ����ȵ����
int Fst_ttt::GetSamp()
{ 
    return tfwv_f.nR * tfwv_f.nsamp;
}



// ���м��㣬��Ҫ����, path:·������������ļ���Ĭ��������
void  Fst_ttt::handle_fst(std::string& path, 
    std::string file_tmogrm  ,  
    std::string file_Rtmogrm ,  
    std::string file_others  )
{  
    std::string outFilePath_tmogrm = path + file_tmogrm;
    std::string outFilePath_Rtmogrm = path + file_Rtmogrm;
    std::string outFilePath_others = path + file_others;
    std::string outFilePath_tmpfile = path + "/fst_out.txt";
    std::string outFilePath_tmpfile2 = path + "/fst_same.txt";
    std::string outFilePath_tmpfile_ttlsns = path + "/fst_ttlsns.txt";
    // ��txt�ļ�д��
    FILE* output_tmogrm = NULL;
    FILE* output_Rtmogrm = NULL;
    FILE* output_others = NULL;
    FILE* tmp_out = NULL;// ����fstbrk���
    FILE* tmp_save = NULL;// ����fstbrk���
    FILE* tmp_ttlsns = NULL;// ����fstbrk���
    
    fopen_s(&tmp_out, outFilePath_tmpfile.c_str(), "w");
    fopen_s(&tmp_save, outFilePath_tmpfile2.c_str(), "w");
    fopen_s(&tmp_ttlsns, outFilePath_tmpfile_ttlsns.c_str(), "w");
    //fopen_s(&output_tmogrm, "outFielName.txt", "w");
    //fopen_s(&output_tmogrm, outFielName, "w");
    fopen_s(&output_tmogrm, outFilePath_tmogrm.c_str(), "w");
    fopen_s(&output_Rtmogrm, outFilePath_Rtmogrm.c_str(), "w");
    fopen_s(&output_others, outFilePath_others.c_str(), "w");
    if (output_tmogrm == NULL)
    {
        /* �����ͬ���ļ��Ḳ�ǣ�����--�� Ӧ�ø���ϵͳʱ�������ļ��У�
        ��дһ��������������Ҫ���ļ���*/
        std::cout << "[ERROR]<handle_fst> outFielName ·�����󣡣�" << std::endl;
        return;
    }
    
    /*SSS = dtfilter(SSS); %ʱ��ƽ������һ�㶼��Ҫ��ʱ���ƽ������
    ����Ҳ�ǿ��Եģ��Ƽ���ƽ��*/
    //dtfilter(DTC, m_num_depth, 5);

    // �����������
    for (int i = 0; i < m_num_depth; i++)
    {
        // ���ڲ���
        /*if (i > 2)
            break;*/

        if (i >= 213)
            int tp = 3;

        // tt->TTAV_tt ,pkmethod��  m_thrs��Ҫ����
        /* RR (e.g., 0.5ft & 0.1524m)
        * wdwdth, 2-5���������ھ���,  200-500,  �ǵ���,us  --300
        *   m_tstart = 0 0 0 
        */
     /*   fstbrkX(DTC[i], TTAV_tt[i], tfwv_f.nsamp, tfwv_f.nR, tfwv_f.RR,
            tfwv_f.period, m_wdwdth, m_thrs, m_wvlen, m_pkmethod, m_rarvtm, m_allData[i], m_tstart);
       */
        
        if (i >= 5)
            int debug = 1;

        // ʹ����������� --tmp
        m_rarvtm_all[i] = (float*)malloc(tfwv_f.nR * sizeof(float));
        fstbrkX(DTC[i], TTAV_tt[i], tfwv_f.nsamp, tfwv_f.nR, tfwv_f.RR,
            tfwv_f.period, m_wdwdth, m_thrs, m_wvlen, m_pkmethod, m_rarvtm_all[i], m_allData[i], m_tstart);

        //period 1.2e-5
        if (i >= 5)
            int debug = 1;

        //for(int index_i = 0; index_i < 7; index_i++)
        //{ 
        //    if (m_rarvtm[index_i] > 20000 || m_rarvtm[index_i] < -20000) //Ϊ�˻�ͼ����
        //        m_rarvtm[index_i] = 0;
        //    fprintf(tmp_out, "%.5lf ", m_rarvtm[index_i]);
        //}
        //if (m_rarvtm[7] > 20000 || m_rarvtm[7] < -20000) //Ϊ�˻�ͼ����
        //    m_rarvtm[7] = 0;
        //fprintf(tmp_out, "%.5lf\n", m_rarvtm[7]);

       
#pragma region ע��tttomo


 //       for (int j = 0; j < tfwv_f.nR; j++)// float ת double
 //       {
 //           if (b_Caliper)
 //               m_ttslns[j] = m_TT1[j * m_num_depth + i];/*  nR1, nR2, ... */
 //           else
 //               m_ttslns[j] = 0;
 //           m_ttpick[j] = m_rarvtm[j]; //ttpick = rarvtm; 
 //           //std::cout << "[Message] m_rarvtm :" << m_rarvtm[j] << std::endl;

 //           // ����
 //           //m_ttpick[j] = Aals[i][j];


 //           // std::cout << "Aal:" << Aals[i][j];
 //           //std::cout << "m_ttslns:" << Aals[i][j];
 //          
 //         /*  if( i >= 200 && i <= 250)
 //               std::cout << "Aal:" << Aals[i][j];
 //           std::cout << "\n";*/

 //           /*if(j < 7)
 //               fprintf(tmp_save, "%.5lf ", Aals[i][j]);
 //           else
 //               fprintf(tmp_save, "%.5lf\n", Aals[i][j]);*/
 //       }
 //       //std::cout << std::endl;
 //       if (i >= 53)
 //           int tp = 3;

 //       if (  i >= 10)
 //           int tp = 3;
 //       // NMAX exceeded
 //       /* ��ss0���յó���ֵ��֪������ֻ�ܵ�slnsӰ�죬��DTC[i]����� DTC��λΪ��
 //          ǰ��Ƚ�ȷ��DTCת�˵�λft ???*/
 //       
 ///*       std::cout << "tfwv_f.nR:" << tfwv_f.nR
 //           << "  tfwv_f.TR:" << tfwv_f.TR
 //           << "  Caliper[i]:" << Caliper[i]
 //           << "  m_dtf:" << m_dtf
 //           << "  m_tlod_diaTOOL:" << m_tlod_diaTOOL
 //           << "  DTC[i]:" << DTC[i]
 //           << "  m_dr:" << m_dr
 //           << "  m_Vav:" << m_Vav << std::endl;*/


 //       tttomo(m_ttpick, m_ttslns, tfwv_f.nR, tfwv_f.TR, tfwv_f.RR, Caliper[i], m_dtf ,
 //           m_tlod_diaTOOL, DTC[i], m_dr, m_Vav,
 //           m_tmogrm, m_Rtmogrm, m_ttfit, &m_ssfit,
 //           &m_ss0, &m_ssd, &m_dop);

 //       
 ///*       tttomo(m_ttpick, m_ttslns, tfwv_f.nR, tfwv_f.TR * 0.3048, tfwv_f.RR * 0.3048, Caliper[i] * 0.3048, m_dtf / 0.3048,
 //           m_tlod_diaTOOL * 0.3048, DTC[i], m_dr, m_Vav / 0.3048,
 //           m_tmogrm, m_Rtmogrm, m_ttfit, &m_ssfit,
 //           &m_ss0, &m_ssd, &m_dop);*/

 //       if (i % 70 == 0)
 //       {
 //           std::cout << "[Message] ����д��ڣ�i:" << i << "������" << std::endl;
 //       }
 //       
 //       // ���
 //       fprintf(output_tmogrm, "%.5lf ", m_dept[i]);
 //       fprintf(output_Rtmogrm, "%.5lf ", m_dept[i]);
 //       fprintf(output_others, "%.5lf ", m_dept[i]);
 //       for (int index_i = 0; index_i < m_size_tmp - 1; index_i++)
 //       { 
 //           //std::cout << "��Message��m_tmogrm[" << index_i << "]:" << m_tmogrm[index_i] << std::endl;
 //           // �����ֵĴ�
 //           //fprintf(output_tmogrm, "m_tmogrm[%d]:%.5lf ", index_i, m_tmogrm[index_i]);
 //           // ֻ������ 
 //           fprintf(output_tmogrm, "%.5lf ", m_tmogrm[index_i]);
 //           fprintf(output_Rtmogrm, "%.5lf ", m_Rtmogrm[index_i]);
 //       }
 //       fprintf(output_tmogrm, "%.5lf\n", m_tmogrm[m_size_tmp - 1]);
 //       fprintf(output_Rtmogrm, "%.5lf\n", m_Rtmogrm[m_size_tmp - 1]);


 //       for (int index_i = 0; index_i < tfwv_f.nR; index_i++)
 //       {
 //           //std::cout << "��Message��m_ttfit[" << index_i << "]:" << m_ttfit[index_i] << std::endl;
 //           //fprintf(output_tmogrm, "m_ttfit[%d]:%.5lf ", index_i, m_ttfit[index_i]);
 //           fprintf(output_others, "%.5lf ", m_ttfit[index_i]);
 //       }
 //       fprintf(output_others, "%.5lf ", m_ssfit);
 //       fprintf(output_others, "%.5lf ", m_ss0);
 //       fprintf(output_others, "%.5lf ", m_ssd);
 //       //fprintf(output_others, "m_dop:%.5lf ", m_dop);
 //       fprintf(output_others, "\n");
#pragma endregion

    }

    // ��ֵ�˲� 
    //medianFilter(m_rarvtm, tfwv_f.nR, m_MFilter);
    float* ttt_tp = (float*)alloca(m_num_depth * sizeof(float));
  
    for (int j = 0; j < tfwv_f.nR; j++)
    {
        //float sum = 0;
        for (int i = 0; i < m_num_depth; i++)
        { 
          /*  if (std::isnan(m_rarvtm_all[i][j]))
            {
                m_rarvtm_all[i][j] = 0;
            }*/
            ttt_tp[i] = m_rarvtm_all[i][j];
           /* sum += ttt_tp[i];*/
        }

       /* for (int i = 0; i < m_num_depth; i++)
        {
            if (ttt_tp[i] == 0)
                ttt_tp[i] = sum / m_num_depth;
        }*/


        medianFilter(ttt_tp, m_num_depth, m_MFilter); // ��ֵ�˲�

        dtfilter(ttt_tp, m_num_depth, 200); // ƽ���˲�


        for (int i = 0; i < m_num_depth; i++)
        {
            m_rarvtm_all[i][j] = ttt_tp[i];
        } 
    }
   
    //��ӡTTT
    for (int i = 0; i < m_num_depth; i++)
    {
        for (int j = 0; j < tfwv_f.nR - 1; j++)
        { 
            fprintf(tmp_out, "%.5lf ", m_rarvtm_all[i][j]);
        }
        fprintf(tmp_out, "%.5lf\n", m_rarvtm_all[i][tfwv_f.nR - 1]);
    }

    // �ó���TTT����ʼ����
    for (int i = 0; i < m_num_depth; i++)
    {
        for (int j = 0; j < tfwv_f.nR; j++)// float ת double
        {
            if (b_Caliper)
            {
                m_ttslns[j] = m_TT1[j * m_num_depth + i];/*  nR1, nR2, ... */
                
            } 
            else
                m_ttslns[j] = 0;

            if( j < tfwv_f.nR - 1)
                fprintf(tmp_ttlsns, "%.5lf ", m_ttslns[j]);
            // ԭʼ����
            m_ttpick[j] = m_rarvtm_all[i][j]; //ttpick = rarvtm; 
            // ����,��ȷ���� --tmp
            //m_ttpick[j] = Aals[i][j] + 25;
            if(j < 7)
               fprintf(tmp_save, "%.5lf ", Aals[i][j]);
            else
               fprintf(tmp_save, "%.5lf\n", Aals[i][j]);
 
            // NMAX exceeded 
        }
        fprintf(tmp_ttlsns, "%.5lf\n", m_ttslns[tfwv_f.nR - 1]);

        tttomo(m_ttpick, m_ttslns, tfwv_f.nR, tfwv_f.TR, tfwv_f.RR, Caliper[i], m_dtf,
            m_tlod_diaTOOL, DTC[i], m_dr, m_Vav,
            m_tmogrm, m_Rtmogrm, m_ttfit, &m_ssfit,
            &m_ss0, &m_ssd, &m_dop);
         

        if (i % 70 == 0)
        {
            std::cout << "[Message] ����д��ڣ�i:" << i << "������" << std::endl;
        }

        // ���
        fprintf(output_tmogrm, "%.5lf ", m_dept[i]);
        fprintf(output_Rtmogrm, "%.5lf ", m_dept[i]);
        fprintf(output_others, "%.5lf ", m_dept[i]);
        for (int index_i = 0; index_i < m_size_tmp - 1; index_i++)
        { 
            fprintf(output_tmogrm, "%.5lf ", m_tmogrm[index_i]);
            fprintf(output_Rtmogrm, "%.5lf ", m_Rtmogrm[index_i]);
        }
        fprintf(output_tmogrm, "%.5lf\n", m_tmogrm[m_size_tmp - 1]);
        fprintf(output_Rtmogrm, "%.5lf\n", m_Rtmogrm[m_size_tmp - 1]);


        for (int index_i = 0; index_i < tfwv_f.nR; index_i++)
        {
             fprintf(output_others, "%.5lf ", m_ttfit[index_i]);
        }
        fprintf(output_others, "%.5lf ", m_ssfit);
        fprintf(output_others, "%.5lf ", m_ss0);
        fprintf(output_others, "%.5lf ", m_ssd); 
        fprintf(output_others, "\n");
        
    }
    fclose(output_tmogrm);
    fclose(output_Rtmogrm);
    fclose(output_others);
    fclose(tmp_ttlsns);
}



//���� ֱ�Ӷ�ȡfstbrk�������
// ��ȡSlowness.txt����ȡ����DTC���� slowness�� mode=Caliper ʱЯ����ֱ������
void Fst_ttt::test_read_real_fstbrk(const char* filename )
{
    std::ifstream ifs(filename);
    // ��� ifs����
    if (!ifs)
    {
        std::cout << "[ERROR]<Fst_ttt::read_DTC> ifs NOT EXIST" << std::endl;
        return;
    }
    std::string line;
    int line_count = 0;
    int index = 0; 
    std::ifstream file(filename);
    while (getline(file, line)) {
        line_count++;
    }
    file.close();
    file.open("filename");
    line_count--;  
    std::cout << "[Message]line_count:" << line_count << std::endl;
    Aals = (float**)malloc(m_num_depth * (sizeof(float*)));  

    std::getline(ifs, line); 
    while (std::getline(ifs, line)) {
        double dept;
        Aals[index] = (float*)malloc(8 * (sizeof(float)));
        std::stringstream word(line);   
        word >> dept;
        for (int i = 0; i < 8; i++)
        {
            word >> Aals[index][i];
        }
        index++;  
    }
    Aals[m_num_depth - 1 ] = (float*)malloc(8 * (sizeof(float)));
    for (int i = 0; i < 8; i++)
    {
        Aals[m_num_depth - 1][i]=0;
    }
}



#include "test.h"
#include "base.h"

#include <fstream>
#include <iostream>

/*

����cpp����Ҫ������ʱ���ԣ���Ӱ������ʹ��


*/

 
// ��ȡ����--��������� --tmp
void read_test_b(const char* fileName, const char* outFielName, int mode )
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    FILE* output = NULL;
    fopen_s(&output, outFielName, "w");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    long remaining = fsize - ftell(fp);
    const int size = remaining;
    char* buff = (char*)alloca(remaining * sizeof(char));

    int len = fread(buff, 1, remaining, fp);

    for (int i = 0; i < len; i++) {
        printf("%02x ", buff[i]);
    }
    // ֱ�Ӵ�ӡ�ֽ�
    char t;
    while (fscanf_s(fp, "%c", &t) != EOF) {
        std::cout << "." << t;
    }
}
// ��ȡ����--ȫ����float4����double8 --tmp
void read_test_fd(const char* fileName, const char* outFielName, int mode )
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
    // ȫ����float4��ȡ
    int count = 0;
    while (!feof(fp)) {
        //float val = readfloat4(fp);
        double val = readdouble8(fp);


        std::cout << count << ":" << val << " ";
        //fprintf(output, "%f ", val);

        if (int(val) == 7425)
        {
            std::cout << "\n\n\n����:" << count << std::endl;
        }

        count++;
        if (count >= 2000)
            return;
    }

}



#pragma region 23.10.25
//���ԣ����ڼ��㵥����ȵ����
//fst_ttt.ft_main(); //���Է���û�д��� NMAX

/*  ��λ���⣺--> �ڲ���λȫΪӢ�� -
������Ӣ��  ��תӢ�� double -
ʱ�λ���� us/m Ҳ������ us/ft(Ӣ��)  -
    dr:3ft/60 ������ 1m/60, -
    ���0.1524m��ΪӢ�� -
    TR,RR תӢ�� -
    W1 W2 kHz

�����ٶȿ��ǣ�10e6 -

*/

/* �������̣�
�˲�  37 -
���� ttt   81   249 SAV=mean(...) �����ʱ��������� -

m_ttpick ����ȥ֮ǰ�����˲�-һ����ֵ�˲� medfilt1 -
m_ttslns ��Ҫ���� -
vav ��Ҫ���� -


�����ļ����棬����tmg001,tmg002 ...
Rtmg001 Rtmg002
4��ֵ + 8��ttfit  ����dop
���λתΪm
*/
#pragma endregion


#pragma region 23.10.22 
/*
ImgRes_Pre_XX.dat��513 (double����)
���ࣺ501
*/

/*����ļ�����Դ�ļ���һ�����ɣ�
��reflection_image.dat--> reflection_image.txt������������ʶ�� */

// ������float4-VeloSection-501
// "E:\Proj\vsProj\st_FileSave\�½��ļ��� (2)\ImgRes_Pre_XX.dat"
// ����VeloSection dat "E:\Proj\vsProj\st_FileSave\FY3-H3.dat"
//read_VeloSection_dat("E:\\Proj\\vsProj\\Desktop\\VeloSection.dat", "output_VeloSection.txt");
//read_VeloSection_dat("E:\\Proj\\vsProj\\st_FileSave\\FY3-H3.dat", "output_FY3-H3.txt");

//���ݶ���double8-reflection-513
//read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\ImgRes_Pre_XX.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\output_ImgRes_Pre_XX.txt",513);//  1026/2
//read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\reflection_image.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\output_reflection_image.txt", 513);
// �ٲ���
//read_VeloSection_dat("E:\\Proj\\vsProj\\st_FileSave\\datas_all\\��������\\����705\\velocity_profile.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\datas_all\\��������\\����705\\output_velocity_profile.txt", 501);//



//�����ƶ�ȡ���� read_as_WIS_d
//read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\reflection_image.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\output_reflection_image_b2.dat", 513);//

//read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\reflection_image.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\output_reflection_image_WIS.dat", 513);//

// ��ȡ�ı��ļ�תΪ WIS
/*   read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\reflection_image.dat",
       "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\SWI.txt",513);*/


       // ������������˲���mode������������mode=1 SWI ��mode=2 VEP 
       /*read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\reflection_image.txt",
           "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\SWI.txt",513,1);*/

           /*  read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\velocity_profile.txt",
                 "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\VEP.txt", 501,2);*/

                 // �����޸�  
                 //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\����-ԭʼ���� - ����)", ".txt");
                 //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\datas_all\��������)", ".txt");
                 //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\datas_all\�½��ļ���)", ".dat");

                 // ��������û��
                 /*batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\����23-10-20)", ".dat");*/

                 //read_as_WIS_d(R"(E:\Proj\vsProj\st_FileSave\����23-10-20\SWI_LR.dat)",
                 //    R"(E:\Proj\vsProj\st_FileSave\����23-10-20\SWI_LR.txt)", 513, 1);
                 //read_as_WIS_d(R"(E:\Proj\vsProj\st_FileSave\����23-10-20\SWI_UD.dat)",
                 //    R"(E:\Proj\vsProj\st_FileSave\����23-10-20\SWI_UD.txt)", 513, 1);

                /* read_VeloSection_dat(R"(E:\Proj\vsProj\st_FileSave\����23-10-20\velocity_profile.dat)",
                     R"(E:\Proj\vsProj\st_FileSave\����23-10-20\output_velocity_profile.txt)", 501);*/

                     // float
//read_as_WIS(R"(E:\Proj\vsProj\st_FileSave\����23-10-20\velocity_profile.dat)",
//    R"(E:\Proj\vsProj\st_FileSave\����23-10-20\VEP.txt)", 501, 2);//VEP.txt
#pragma endregion


#pragma region ��ʼ�汾

//// ����ָ���ļ����������ļ���
//void get_sub_folders(std::string path, std::vector<std::string>& folders) {
//
//    intptr_t handle = 0;
//    struct _finddata_t file_info;
//
//    std::string search_path = path + "/*";
//
//    if ((handle = _findfirst(search_path.c_str(), &file_info)) != -1) {
//
//        do {
//
//            if (file_info.attrib & _A_SUBDIR && strcmp(file_info.name, ".") && strcmp(file_info.name, "..")) {
//
//                std::string folder_path = path + "/" + file_info.name;
//                folders.push_back(folder_path);
//
//                // �ݹ�������ļ���
//                get_sub_folders(folder_path, folders);
//
//            }
//
//        } while (_findnext(handle, &file_info) == 0);
//
//        _findclose(handle);
//
//    }
//
//}
//
//// ��ȡĳһ�ļ����£�ָ���ļ���׺���ļ�
//void get_need_file(std::string path, std::vector<std::string>& file, std::string ext)
//{
//    intptr_t file_handle = 0;
//    struct _finddata_t file_info;
//    std::string temp;
//    if ((file_handle = _findfirst(temp.assign(path).append("/*" + ext).c_str(), &file_info)) != -1)
//    {
//        do
//        {
//            file.push_back(temp.assign(path).append("/").append(file_info.name));
//        } while (_findnext(file_handle, &file_info) == 0);
//        _findclose(file_handle);
//    }
//}
// 


// ���Ժ���
/*
// ��ȡ����--��������� --tmp
void read_test_b(const char* fileName, const char* outFielName, int mode = 0)
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    FILE* output = NULL;
    fopen_s(&output, outFielName, "w");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    long remaining = fsize - ftell(fp);
    const int size = remaining;
    char* buff = (char*)alloca(remaining * sizeof(char));

    int len = fread(buff, 1, remaining, fp);

    for (int i = 0; i < len; i++) {
        printf("%02x ", buff[i]);
    }
    // ֱ�Ӵ�ӡ�ֽ�
    char t;
    while (fscanf_s(fp, "%c", &t) != EOF) {
        std::cout << "." << t;
    }
}
// ��ȡ����--ȫ����float4����double8 --tmp
void read_test_fd(const char* fileName, const char* outFielName, int mode=0)
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
    // ȫ����float4��ȡ
    int count = 0;
    while (!feof(fp)) {
        //float val = readfloat4(fp);
        double val = readdouble8(fp);


        std::cout << count << ":" << val << " ";
        //fprintf(output, "%f ", val);

        if (int(val) == 7425)
        {
            std::cout << "\n\n\n����:" << count << std::endl;
        }

        count++;
        if (count >= 2000)
            return;
    }

}
*/


// ����base.cpp
 /*
double readdouble8(FILE* fpin)
{
    double data;
    char* str;
    str = (char*)&data;
    for (int i = 7; i >= 0; i--)
        fread(str + i, 1, 1, fpin);

    return (data);
}

short int readint2(FILE* fpin)
{
    short int data;
    char* str;
    str = (char*)&data;
    for (int i = 1; i >= 0; i--)
        fread(str + i, 1, 1, fpin);

    return (data);
}

float readfloat4(FILE* fpin)
{
    float data;
    char* str;
    str = (char*)&data;
    for (int i = 3; i >= 0; i--)
        fread(str + i, 1, 1, fpin);

    return (data);
}

int readint4(FILE* fpin)
{
    int data;
    char* str;
    str = (char*)&data;
    for (int i = 3; i >= 0; i--)
        fread(str + i, 1, 1, fpin);

    return (data);
}
*/



// �����ļ��У���txt
/*
    std::vector<std::string> my_file;
    std::string need_extension = ".txt";


    // ���ԣ������ļ����������ļ���
    std::vector<std::string> folders;
    std::string path = R"(E:\Proj\vsProj\st_FileSave\����-ԭʼ���� - ����)";

    get_sub_folders(path, folders);

    // ��ӡ���
    for (std::string folder : folders) {
        std::cout << folder << std::endl;
        my_file.clear();
        get_need_file(folder, my_file, need_extension);
        for (int i = 0; i < my_file.size(); i++)
        {
            std::cout << "File " << i + 1 << " is:" << std::endl;
            std::cout << my_file[i] << std::endl;
        }
        if (my_file.size() == 0)
        {
            std::cout << "No file can be found!" << std::endl;
        }
        else
        {
            std::cout << std::endl << "Find " << my_file.size() << " file(s)." << std::endl;
        }
    }

*/

 

//���ԣ�����ĳ�ļ���������txt�ļ� 
/*
std::string file_path = R"(E:\Proj\vsProj\st_FileSave\����-ԭʼ���� - ����\1-��Դ3-H3)";
std::vector<std::string> my_file;
std::string need_extension = ".txt";
   get_need_file(file_path, my_file, need_extension);
    for (int i = 0; i < my_file.size(); i++)
    {
        std::cout << "File " << i + 1 << " is:" << std::endl;
        std::cout << my_file[i] << std::endl;
    }
    if (my_file.size() == 0)
    {
        std::cout << "No file can be found!" << std::endl;
    }
    else
    {
        std::cout << std::endl << "Find " << my_file.size() << " file(s)." << std::endl;
    } */
 

// ����main
    //float ret[100];
    //float scard[8] = {0};//
    //float depth;
    //float data[8 * 672] = {0};
    //// char cc[50]; 
    //float t;
    //float tstart[8];
    //float slowness = 120.0f;
    //float tt = 0.0f;
    //int nsamp = 672;
    //int nR = 8;
    //float RR = 0.1524f;
    //float period = 12.0f;
    //float wdwdth = 300.0f;
    //float thrs = 0.5f;// 150.0f;  --?
    //float wvlen = 1200.0f;
    //float pkmethod = 2.0f;
    //float* rarvtm = (float*)alloca(nR * sizeof(float)); // ����Ϊ nR
    //// tttomo
    //FILE* fp1;
    //fopen_s(&fp1, "E:\\Proj\\vsProj\\st\\wtest.dat", "rb");
    //// tttomo���ص���ֵ
    //int size_tmp = 121;//tmogrm[60 + j] = 0.7 * vav;  ��������µ�120
    //double* tmogrm = (double*)alloca(size_tmp * sizeof(double));// 60����Ϊtmogrm[60 + j] = 0.7 * vav;
    //double* Rtmogrm = (double*)alloca(size_tmp * sizeof(double));
    //double* ttfit = (double*)alloca(nR * sizeof(double));//   ttfit[i] = tvtm(V); 
    //double ssfit;  //	*ssfit = 1.0e6 / V;
    //double ss0; //*ss0 = 1e6 / V0;
    //double ssd; //*ssd = 1e6 / Vdeep;;
    //double dop; //*dop = (Vdeep - V0) / gv;;
    //// tttomo������Ҫ��ֵ
    //double* ttpick = (double*)alloca(nR * sizeof(double)); // nR Ϊ���鳤��
    //double* ttslns = (double*)alloca(nR * sizeof(double));
    //double Calpr = 10;
    //double dtf = 5;
    //double tlod = 10;
    //double dr = 10;
    //double vav = 10;
    //int i = 0;
    //int k = 0;

    /*for (int i = 0; i < 8; i++)
    {
        tstart[i] = 0;
    }
    if (fp1 == NULL) {
        puts("���ļ�ʧ��");
        return 0;
    }

    fseek(fp1, 0, SEEK_SET);
    while (fscanf_s(fp1, "%f", &t) != EOF) {
         process data
        if (i < 8)
        {
            scard[i] = t;
        }
        else if (i == 8)
        {
            depth = t;
        }
        else
        {
            data[k] = t;
            k++;
        }
        i++;
    }
    fclose(fp1);

    float m = 120.0f;
    float ma[3] = { 2.3f, 3.4f, 4.5f };

     ����ֵrarvtm
    fstbrkX(slowness,tt, nsamp,nR,RR,period,wdwdth,thrs,wvlen,pkmethod,rarvtm, data, tstart);

    for (int i = 0; i < nR; i++)
    {
        ttpick[i] = 0;
        ttslns[i] = 0;
    }



    tttomo(ttpick,ttslns, nR, TR, RR, Calpr, dtf,
        tlod, slowness, dr, vav,
        tmogrm, Rtmogrm, ttfit, &ssfit,
        &ss0, &ssd, &dop);

    for (int i = 0; i < size_tmp; i++)
    {
        std::cout << "  tmogrm[" << i << "]:" << tmogrm[i]
            << "  Rtmogrm[" << i << "]:" << Rtmogrm[i]  << std::endl;
    }

    for (int i = 0; i < nR; i++)
    {
        std::cout << "  ttfit[" << i << "]:" << ttfit[i] << std::endl;
    }

    std::cout << " ssfit:" << ssfit << " ,"
        << " ss0:" << ss0 << " ,"
        << " ssd:" << ssd << " ,"
        << " dop:" << dop << std::endl;*/
        //void tttomo(double* ttpick, double* ttslns, int nrec, double trsp, double rrsp, double Calpr,
        //    double dtf, double tlod, double slns, double dr, double vav,
        //    double* tmogrm, double* Rtmogrm, double* ttfit, double* ssfit,
        //    double* ss0, double* ssd, double* dop);

            /*

            [tmogrm, Rtmogrm, tfit, slsfit, sls0, slsd, dpop] =  tttomo(ttpick,...   %�ٶȲ����������tttomo()����
      ttslns,...
      nR,...nrec
      TR,...trsp
      RR,... rrsp
      bhod,... Calpr
      dtF,... dtf
      diaTOOL,... tlod
      slowness,... slns
      dr,...dr
      Vav); vav
         // ttpick, ��ȡ����(1~nwf ��)
                 // ttslns, ʱ�����(1~nwf ��, ���������岿��)
                 // Calpr, ����, m
                 // dtf, ����ʱ��
                 // tlod, ����ֱ��
                 // slns, ��ǰ��ȵ�ʱ��
                 // dr, �����򲽳�(�������/60)
                 // vav, ƽ���ٶ�
                 // tmogrm, �ٶ�����(ʵ��)
                 // Rtmogrm, �ٶ�����(���)
                 // ttfit, ��ttpick��ϵĲ���
                 // ssfit, ��slns��ϵĲ���
                 // ss0, �����ϵ�ʱ��
                 // ssd, ���͸����ϵ�ʱ��
                 // dop, ���͸���*/
                 //printf("\n---------���rarvtm--------\n" );
                 //for (int i = 0; i < nR; i++)
                 //{
                 //    printf("rarvtm[%d]: %f\n", i, rarvtm[i]);
                 //}

#pragma endregion
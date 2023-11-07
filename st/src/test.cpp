
#include "test.h"
#include "base.h"

#include <fstream>
#include <iostream>

/*

测试cpp，主要用于临时测试，不影响正常使用


*/

 
// 读取测试--二进制输出 --tmp
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
    // 直接打印字节
    char t;
    while (fscanf_s(fp, "%c", &t) != EOF) {
        std::cout << "." << t;
    }
}
// 读取测试--全部以float4或者double8 --tmp
void read_test_fd(const char* fileName, const char* outFielName, int mode )
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
    // 全部以float4读取
    int count = 0;
    while (!feof(fp)) {
        //float val = readfloat4(fp);
        double val = readdouble8(fp);


        std::cout << count << ":" << val << " ";
        //fprintf(output, "%f ", val);

        if (int(val) == 7425)
        {
            std::cout << "\n\n\n长度:" << count << std::endl;
        }

        count++;
        if (count >= 2000)
            return;
    }

}



#pragma region 23.10.25
//测试，用于计算单个深度点测试
//fst_ttt.ft_main(); //测试发现没有触发 NMAX

/*  单位问题：--> 内部单位全为英尺 -
井径是英寸  ，转英尺 double -
时差单位可能 us/m 也可能是 us/ft(英尺)  -
    dr:3ft/60 而不是 1m/60, -
    间隔0.1524m改为英尺 -
    TR,RR 转英尺 -
    W1 W2 kHz

计算速度考虑：10e6 -

*/

/* 处理流程：
滤波  37 -
计算 ttt   81   249 SAV=mean(...) 计算的时候处理非零项 -

m_ttpick 传进去之前进行滤波-一阶中值滤波 medfilt1 -
m_ttslns 需要计算 -
vav 需要计算 -


三个文件保存，列名tmg001,tmg002 ...
Rtmg001 Rtmg002
4个值 + 8个ttfit  不存dop
最后单位转为m
*/
#pragma endregion


#pragma region 23.10.22 
/*
ImgRes_Pre_XX.dat：513 (double类型)
其余：501
*/

/*输出文件名和源文件名一样即可，
如reflection_image.dat--> reflection_image.txt，便于批处理识别 */

// 数据是float4-VeloSection-501
// "E:\Proj\vsProj\st_FileSave\新建文件夹 (2)\ImgRes_Pre_XX.dat"
// 测试VeloSection dat "E:\Proj\vsProj\st_FileSave\FY3-H3.dat"
//read_VeloSection_dat("E:\\Proj\\vsProj\\Desktop\\VeloSection.dat", "output_VeloSection.txt");
//read_VeloSection_dat("E:\\Proj\\vsProj\\st_FileSave\\FY3-H3.dat", "output_FY3-H3.txt");

//数据都是double8-reflection-513
//read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\ImgRes_Pre_XX.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\output_ImgRes_Pre_XX.txt",513);//  1026/2
//read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\reflection_image.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\output_reflection_image.txt", 513);
// 再测试
//read_VeloSection_dat("E:\\Proj\\vsProj\\st_FileSave\\datas_all\\更新数据\\满深705\\velocity_profile.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\datas_all\\更新数据\\满深705\\output_velocity_profile.txt", 501);//



//二进制读取测试 read_as_WIS_d
//read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\reflection_image.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\output_reflection_image_b2.dat", 513);//

//read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\reflection_image.dat",
//    "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\output_reflection_image_WIS.dat", 513);//

// 读取文本文件转为 WIS
/*   read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\reflection_image.dat",
       "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\SWI.txt",513);*/


       // 下面两个添加了参数mode：控制列名，mode=1 SWI ，mode=2 VEP 
       /*read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\reflection_image.txt",
           "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\SWI.txt",513,1);*/

           /*  read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\velocity_profile.txt",
                 "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\VEP.txt", 501,2);*/

                 // 遍历修改  
                 //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\数据-原始备份 - 副本)", ".txt");
                 //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\datas_all\更新数据)", ".txt");
                 //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\datas_all\新建文件夹)", ".dat");

                 // 遍历单层没有
                 /*batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\数据23-10-20)", ".dat");*/

                 //read_as_WIS_d(R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\SWI_LR.dat)",
                 //    R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\SWI_LR.txt)", 513, 1);
                 //read_as_WIS_d(R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\SWI_UD.dat)",
                 //    R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\SWI_UD.txt)", 513, 1);

                /* read_VeloSection_dat(R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\velocity_profile.dat)",
                     R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\output_velocity_profile.txt)", 501);*/

                     // float
//read_as_WIS(R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\velocity_profile.dat)",
//    R"(E:\Proj\vsProj\st_FileSave\数据23-10-20\VEP.txt)", 501, 2);//VEP.txt
#pragma endregion


#pragma region 初始版本

//// 遍历指定文件夹下所有文件夹
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
//                // 递归遍历子文件夹
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
//// 获取某一文件夹下，指定文件后缀的文件
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


// 测试函数
/*
// 读取测试--二进制输出 --tmp
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
    // 直接打印字节
    char t;
    while (fscanf_s(fp, "%c", &t) != EOF) {
        std::cout << "." << t;
    }
}
// 读取测试--全部以float4或者double8 --tmp
void read_test_fd(const char* fileName, const char* outFielName, int mode=0)
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
    // 全部以float4读取
    int count = 0;
    while (!feof(fp)) {
        //float val = readfloat4(fp);
        double val = readdouble8(fp);


        std::cout << count << ":" << val << " ";
        //fprintf(output, "%f ", val);

        if (int(val) == 7425)
        {
            std::cout << "\n\n\n长度:" << count << std::endl;
        }

        count++;
        if (count >= 2000)
            return;
    }

}
*/


// 引入base.cpp
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



// 遍历文件夹，及txt
/*
    std::vector<std::string> my_file;
    std::string need_extension = ".txt";


    // 测试，遍历文件夹下所有文件夹
    std::vector<std::string> folders;
    std::string path = R"(E:\Proj\vsProj\st_FileSave\数据-原始备份 - 副本)";

    get_sub_folders(path, folders);

    // 打印结果
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

 

//测试，遍历某文件夹下所有txt文件 
/*
std::string file_path = R"(E:\Proj\vsProj\st_FileSave\数据-原始备份 - 副本\1-富源3-H3)";
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
 

// 早期main
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
    //float* rarvtm = (float*)alloca(nR * sizeof(float)); // 长度为 nR
    //// tttomo
    //FILE* fp1;
    //fopen_s(&fp1, "E:\\Proj\\vsProj\\st\\wtest.dat", "rb");
    //// tttomo返回的数值
    //int size_tmp = 121;//tmogrm[60 + j] = 0.7 * vav;  根据这个猜的120
    //double* tmogrm = (double*)alloca(size_tmp * sizeof(double));// 60是因为tmogrm[60 + j] = 0.7 * vav;
    //double* Rtmogrm = (double*)alloca(size_tmp * sizeof(double));
    //double* ttfit = (double*)alloca(nR * sizeof(double));//   ttfit[i] = tvtm(V); 
    //double ssfit;  //	*ssfit = 1.0e6 / V;
    //double ss0; //*ss0 = 1e6 / V0;
    //double ssd; //*ssd = 1e6 / Vdeep;;
    //double dop; //*dop = (Vdeep - V0) / gv;;
    //// tttomo函数需要的值
    //double* ttpick = (double*)alloca(nR * sizeof(double)); // nR 为数组长度
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
        puts("打开文件失败");
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

     返回值rarvtm
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

            [tmogrm, Rtmogrm, tfit, slsfit, sls0, slsd, dpop] =  tttomo(ttpick,...   %速度层析成像调用tttomo()函数
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
                 // dop, 最大穿透深度*/
                 //printf("\n---------输出rarvtm--------\n" );
                 //for (int i = 0; i < nR; i++)
                 //{
                 //    printf("rarvtm[%d]: %f\n", i, rarvtm[i]);
                 //}

#pragma endregion
// #include<stdio.h> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <string> 
#include <vector>
#include <io.h>

#include <malloc.h>

#include "fstbrk.h"
#include "tttomo.h" 

 
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

//WIS 读取文件头信息 --tmp
void Read_FileInfor( FILE* fp)
{
    float m_S_Depth = readfloat4(fp);
    float m_E_Depth = readfloat4(fp);
    short int m_nwf = readint2(fp);
    short int m_nwf_L = readint2(fp);
    float m_dt = (float)(readfloat4(fp) * 1000000.0);
    float m_ds = readfloat4(fp);
    float m_dr = readfloat4(fp);
    float m_tr = readfloat4(fp);

    std::cout << "[M]" << m_S_Depth << "  "
        << m_E_Depth << "  "
        << m_nwf << "  "
        << m_nwf_L << "  "
        << m_dt << "  "
        << m_ds << "  "
        << m_dr << "  "
        << m_tr << "  "
        << std::endl;
    // 写入txt
    
}

//WIS 读取二进制文件主程序 --tmp
void read_dat(const char* fileName)
{   

    FILE* fp = NULL; 
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
    char t; 
    int index = 0;
    fseek(fp, 0, SEEK_SET);

    std::cout << std::endl;
    // 直接打印字节
    //while (fscanf_s(fp, "%c", &t) != EOF) {
    //    if (index > 10000)
    //        break;
    //    std::cout << " " << t;
    //    //std::cout << t;
    //    index++;
    //} 
    const int size = 1024 * 10;
    char buff[size];

    int len = fread(buff, 1, size, fp);

    for (int i = 0; i < len; i++) {
        printf("%02x ", buff[i]);
    }

    fclose(fp);
   
    // 先读取文件头,并写入txt
    // Read_FileInfor(fp);

    std::cout <<std::endl;
    std::cout << "[M] read_dat over" << std::endl;
}

// Velo Section 

struct VelHeader
{
    float Arrival_1st;
    int Arrival_Index;
    float W_Len ;
    float TimeBond;
    float Step_t;
    float FrontWidth;
    float BackWidth;
    int FirFilterTag;
    float S_Freq;
    float E_Freq;
    int UseEneTag;
    float EneCoef;
    float WLenCoef;
    float S_DT_1D;
    float E_DT_1D;
    float Step_DT_1D;
    int PointProcess; 
};

// 读取头信息函数 --头文件格式不明确
void read_VelHeadInfo(FILE* fp, FILE* output, VelHeader& velHeader)
{
    int PREFIX_LEN = 266; 
    char* prefix = (char*)alloca((PREFIX_LEN + 1) * sizeof(float));

    fread(prefix, 1, PREFIX_LEN, fp);
    prefix[PREFIX_LEN] = '\0'; //添加字符串结束符
    for (int i = 0; i < PREFIX_LEN + 1; i++)
    {
        std::cout << "prefix:" << prefix[i] << std::endl;
    }
    std::cout << "over" << std::endl;

}

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
// velo Section 读取 数据是float4
void read_VeloSection_dat(const char* fileName, const char* outFielName, int nwf_L_set = 501)
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
  
    //int nwf_L = 1026, Depth_n;// 501
    float S_depth, E_depth, ds, S_Value, Step_Value;

    // read_VelHeadInfo(fp, nwf_L, Depth_n, S_depth, E_depth, ds, S_Value, Step_Value);

    // 打开txt文件写入
    FILE* output = NULL;
    //fopen_s(&output, "output.txt", "w");
    fopen_s(&output, outFielName, "w");
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    VelHeader velHeaderInfo;
    int index = 0;

    double dep = readdouble8(fp);
    fprintf(output, "%.5lf ", dep);
    while (!feof(fp))
    {   
        if (index != 0)// 去除首行空行
        {
            double dep = readdouble8(fp);
            fprintf(output, "\n%.5lf ", dep);
        }
        
        for (int j = 0; j < nwf_L_set; j++)//nwf_L改为外部输入
        {
            long remaining = fsize - ftell(fp) ;
            // std::cout << remaining << std::endl;
            // 
            // V_Section_HeadInfor.Description = "Velo_Profile"; 这里可以看出来
            // 文件末尾确实打印了  Velo_Profile
            if (remaining <= 266) // 文件末尾有信息，但是实现被封装了，无法查看,具体多长不知道
            {// 文件末尾
               //read_VelHeadInfo(fp, output, velHeaderInfo);
                 
                const int size = remaining;
                char* buff = (char*)alloca(remaining * sizeof(char)) ;

                int len = fread(buff, 1, remaining, fp);

                for (int i = 0; i < len; i++) {
                    printf("%02x ", buff[i]);
                }
               // 直接打印字节
                char t;
                while (fscanf_s(fp, "%c", &t) != EOF) { 
                    std::cout << "." << t; 
                } 
                 
                return;
            }
            else
            {
                float velo = readfloat4(fp);
                if (j == nwf_L_set - 1)// 行尾，去空格
                {
                    fprintf(output, "%.5f", velo);
                }
                else
                {
                    fprintf(output, "%.5f ", velo);
                }
                
            }
        }
        index++;
    }
    fclose(fp);
    fclose(output);


}
// velo Section 全部double型
void read_VeloSection_dat_d(const char* fileName, const char* outFielName, int nwf_L_set = 501)
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    } 
    // 打开txt文件写入
    FILE* output = NULL;
    //fopen_s(&output, "output.txt", "w");
    fopen_s(&output, outFielName, "w");
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    VelHeader velHeaderInfo;
    int index = 0;

    double dep = readdouble8(fp);
    fprintf(output, "%.5lf ", dep);
    while (!feof(fp))
    {
        if (index != 0)// 去除首行空行
        {
            double dep = readdouble8(fp);
            fprintf(output, "\n%.5lf ", dep);
        }

        for (int j = 0; j < nwf_L_set; j++)//nwf_L改为外部输入
        {
            long remaining = fsize - ftell(fp);
            // std::cout << remaining << std::endl;
            // 
            // V_Section_HeadInfor.Description = "Velo_Profile"; 这里可以看出来
            // 文件末尾确实打印了  Velo_Profile
            if (remaining <= 266) // 文件末尾有信息，但是实现被封装了，无法查看,具体多长不知道
            {// 文件末尾
               //read_VelHeadInfo(fp, output, velHeaderInfo);

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

                return;
            }
            else
            {
                //float velo = readfloat4(fp);
                double velo = readdouble8(fp);
                if (j == nwf_L_set - 1)// 行尾，去空格
                {
                    fprintf(output, "%.5lf", velo);
                }
                else
                {
                    fprintf(output, "%.5lf ", velo);
                }

            }
        }
        index++;
    }
    fclose(fp);
    fclose(output);


}

// 读取文本数据，保存为WIS格式
void read_txt_as_WIS_d(const char* fileName, const char* outFielName, int nwf_L_set = 501, int mode = 1)
{ 
    std::ifstream ifs(fileName);
    // 检查 ifs存在
    if (!ifs)
    {
        std::cout << "[ERROR] ifs NOT EXIST" << std::endl;
        return;
    }

    // 打开txt文件写入
    FILE* output = NULL;
    //fopen_s(&output, "output.txt", "w");
    //fopen_s(&output, outFielName, "wb");  
    if (fopen_s(&output, outFielName, "wb") != 0) {
        // 打开文件失败,报错
        std::cout << "[ERROR] open file ERROR ： " << outFielName << std::endl;
        return;
    }
    int index = 0;
    double tp_d = 0.0;//临时赋值
    double STDEP, ENDEP, RLEV, num;
    //char line[1024] = {0};
    std::string line;
    double dep;
 
    long currPos[10] = { 0 }; long curr_index = 0;

    fprintf_s(output, "FORWARD_TEXT_FORMAT_1.0\n");
    currPos[curr_index] = ftell(output);
    curr_index++;
    char buffer[50];
    int len = sprintf_s(buffer, "STDEP =     %5.4lf\n", tp_d);
    fwrite(buffer, len, 1, output); 
    currPos[curr_index] = ftell(output);
    curr_index++;
    len = sprintf_s(buffer, "ENDEP =     %5.4lf\n", tp_d);
    fwrite(buffer, len, 1, output); 
    currPos[curr_index] = ftell(output);
    curr_index++;
    len = sprintf_s(buffer, "RLEV  =     %5.4lf\n", tp_d);
    fwrite(buffer, len, 1, output); 

    fprintf_s(output, "CURVENAME ="); // 

    char name[11];
    if (mode == 1)
    {
        for (int i = 0; i < nwf_L_set - 1; i++)
        {
            sprintf_s(name, " SWIS%04d,", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " SWIS%04d", nwf_L_set - 1); //使用sprintf函数来格式化字符串 
        fprintf_s(output, " %s", name);
    }
    else if(mode == 2)
    {
        for (int i = 0; i < nwf_L_set - 1; i++)
        {
            sprintf_s(name, " VEP%04d,", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " VEP%04d", nwf_L_set - 1); //使用sprintf函数来格式化字符串 
        fprintf_s(output, " %s", name);
    }
    else
    {
        std::cout << "[ERROR] mode" << std::endl;
        return;
    }


    fprintf_s(output, "\nEND\n");
    fprintf_s(output, "#DEPTH    ");
    if (mode == 1)
    {
        for (int i = 0; i < nwf_L_set; i++)
        {
            char name[10];
            sprintf_s(name, " SWIS%04d", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
    }
    else
    {
        for (int i = 0; i < nwf_L_set; i++)
        {
            char name[10];
            sprintf_s(name, " VEP%04d", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
    }
 

  
    index = 0;
    while (std::getline(ifs, line)){
        std::stringstream word(line);  /*有个很离谱的错误，忘了是啥*/
        word >> dep;
        if (index == 0)
        {
            STDEP = dep;
        }
        //std::cout << "dep:" << dep << "   ";
        fprintf(output, "\n%.5lf ", dep);
        ENDEP = dep;
        while (word >> num)
        {
            if (word.eof())// 行尾，去空格
            {
                fprintf(output, "%.5lf", num);
            }
            else
            {
                fprintf(output, "%.5lf ", num);
            }
        }
        index++; 
    }
     
   
    RLEV = (ENDEP - STDEP) / index;
    // 定位回变量定义行前


    fseek(output, currPos[0], SEEK_SET); 
    fprintf_s(output, "STDEP = %5.4lf", STDEP);

    fseek(output, currPos[1], SEEK_SET);
    fprintf_s(output, "ENDEP = %5.4lf", ENDEP);

    fseek(output, currPos[2], SEEK_SET);
    fprintf_s(output, "RLEV  = %5.4lf", RLEV);
 

 
    fclose(output);


}

// 读取dat数据，保存为WIS格式 数据是double8类型
void read_as_WIS_d(const char* fileName, const char* outFielName, int nwf_L_set = 501 , int mode = 1)
{
    FILE* fp = NULL;
  
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
    // 打开txt文件写入
    FILE* output = NULL;
    //fopen_s(&output, "output.txt", "w");
    fopen_s(&output, outFielName, "wb");
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    VelHeader velHeaderInfo;
    int index = 0;
    double tp_d = 0.0;//临时赋值
    double STDEP, ENDEP, RLEV;


  /*  char str[50]; 
    sprintf(str, "%f", value);
    int decimal = strcspn(str, ".");
    int integer_digits = decimal;*/


    //char** curve_name = (char**)alloca(nwf_L_set * sizeof(char*));
    long currPos[10] = { 0 }; long curr_index = 0;

    fprintf_s(output, "FORWARD_TEXT_FORMAT_1.0\n");
    currPos[curr_index] = ftell(output);
    curr_index++;
    char buffer[50];
    int len = sprintf_s(buffer, "STDEP =     %5.4lf\n", tp_d);
    fwrite(buffer, len, 1, output);
    //fprintf_s(output, "STDEP = %5.4lf\n", tp_d); // 7305.0698
  /*  char format[20] = "STDEP = %5.4lf\n";
    fprintf_s(output, format, tp_d);*/
    //fprintf_s(output, "STDEP = %5.4lf\n",tp_d);  
    currPos[curr_index] = ftell(output);
    curr_index++;
    len = sprintf_s(buffer, "ENDEP =     %5.4lf\n", tp_d);
    fwrite(buffer, len, 1, output);
    //fprintf_s(output, "ENDEP = %5.4lf\n", tp_d); // 7886.4800
 /*   char format2[20] = "ENDEP = %5.4lf\n";
    fprintf_s(output, format2, tp_d);*/
    //fprintf_s(output, "ENDEP = %5.4lf\n", tp_d);
    currPos[curr_index] = ftell(output);
    curr_index++;
    len = sprintf_s(buffer, "RLEV  =     %5.4lf\n", tp_d);
    fwrite(buffer, len, 1, output);
    //fprintf_s(output, "RLEV  = %5.4lf\n", tp_d); // RLEV  =     0.0762
    //fprintf_s(output, "RLEV  = %5.4lf\n", tp_d);
    
    fprintf_s(output, "CURVENAME ="); // 
    //for (int i = 0; i < nwf_L_set; i++)
    //{
    //    char name[11];
    //    sprintf_s(name, " SWIS%04d,", i); //使用sprintf函数来格式化字符串 
    //    fprintf_s(output, " %s", name);
    //}
    if (mode == 1)
    {
        char name[11];
        for (int i = 0; i < nwf_L_set - 1; i++)
        {
            sprintf_s(name, " SWIS%04d,", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " SWIS%04d", nwf_L_set - 1); //使用sprintf函数来格式化字符串 
        fprintf_s(output, " %s", name);
    }
    else if (mode == 2)
    {
        char name[11];
        for (int i = 0; i < nwf_L_set - 1; i++)
        {
            sprintf_s(name, " VEP%04d,", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " VEP%04d", nwf_L_set - 1); //使用sprintf函数来格式化字符串 
        fprintf_s(output, " %s", name);
    }
    else
    {
        std::cout << "[ERROR] mode" << std::endl;
        return;
    }



    fprintf_s(output, "\nEND\n");
    fprintf_s(output, "#DEPTH    ");
    if (mode == 1)
    {
        for (int i = 0; i < nwf_L_set; i++)
        {
            char name[10];
            sprintf_s(name, " SWIS%04d", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
    }
    else if (mode == 2)
    {
        for (int i = 0; i < nwf_L_set; i++)
        {
            char name[10];
            sprintf_s(name, " VEP%04d", i); //使用sprintf函数来格式化字符串 
            fprintf_s(output, " %s", name);
        }
    }


    //第一个深度
    double dep = readdouble8(fp);
    STDEP = dep;
    fprintf(output, "\n%.5lf ", dep);
    while (!feof(fp))
    {
        std::cout << "dep: " << dep << std::endl;
        if (index != 0)// 去除首行空行
        {
            //double dep = readdouble8(fp);
            dep = readdouble8(fp);
            ENDEP = dep;
            fprintf(output, "\n%.5lf ", dep);
        }

        for (int j = 0; j < nwf_L_set; j++)//nwf_L改为外部输入
        {
            long remaining = fsize - ftell(fp);
            // std::cout << remaining << std::endl;
            // 
            // V_Section_HeadInfor.Description = "Velo_Profile"; 这里可以看出来
            // 文件末尾确实打印了  Velo_Profile
            if (remaining <= 266) // 文件末尾有信息，但是实现被封装了，无法查看,具体多长不知道
            {// 文件末尾
               //read_VelHeadInfo(fp, output, velHeaderInfo);

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

                //return;  没有关闭文件
                break;
            }
            else
            {
                //float velo = readfloat4(fp);
                double velo = readdouble8(fp);
                if (j == nwf_L_set - 1)// 行尾，去空格
                {
                    fprintf(output, "%.5lf", velo);
                }
                else
                {
                    fprintf(output, "%.5lf ", velo);
                }

            }
        }
        index++;
    }

    RLEV = (ENDEP - STDEP) / index;
    // 定位回变量定义行前

    
    fseek(output, currPos[0], SEEK_SET); 
    //fflush(output);
    fprintf_s(output, "STDEP = %5.4lf", STDEP);  
  
    fseek(output, currPos[1], SEEK_SET);
    fprintf_s(output, "ENDEP = %5.4lf", ENDEP); 

    fseek(output, currPos[2], SEEK_SET);
    fprintf_s(output, "RLEV  = %5.4lf", RLEV);  
    
    

    //fseek(output, currPos[0], SEEK_SET);

    //len = sprintf_s(buffer, "STDEP = %20.4lf\n", STDEP);
    //std::cout << "\n\n\n\n LEN:::::" << len << "\n\n\n\n\n\n";
    //fwrite(buffer, len, 1, output);

    //fseek(output, currPos[1], SEEK_SET);

    //len = sprintf_s(buffer, "ENDEP = %5.4lf\n", ENDEP);
    //fwrite(buffer, len, 1, output);

    //fseek(output, currPos[2], SEEK_SET);

    //len = sprintf_s(buffer, "RLEV = %5.4lf\n", RLEV);
    //fwrite(buffer, len, 1, output);

    


    fclose(fp);
    fclose(output);


}



 
// 遍历指定文件夹下所有文件夹
void get_sub_folders(std::string path, std::vector<std::string>& folders) {

    intptr_t handle = 0;
    struct _finddata_t file_info;

    std::string search_path = path + "/*";

    if ((handle = _findfirst(search_path.c_str(), &file_info)) != -1) {

        do {

            if (file_info.attrib & _A_SUBDIR && strcmp(file_info.name, ".") && strcmp(file_info.name, "..")) {

                std::string folder_path = path + "/" + file_info.name;
                folders.push_back(folder_path);

                // 递归遍历子文件夹
                get_sub_folders(folder_path, folders);

            }

        } while (_findnext(handle, &file_info) == 0);

        _findclose(handle);

    }

}

// 获取某一文件夹下，指定文件后缀的文件
void get_need_file(std::string path, std::vector<std::string>& file, std::string ext)
{
    intptr_t file_handle = 0;
    struct _finddata_t file_info;
    std::string temp;
    if ((file_handle = _findfirst(temp.assign(path).append("/*" + ext).c_str(), &file_info)) != -1)
    {
        do
        {
            file.push_back(temp.assign(path).append("/").append(file_info.name));
        } while (_findnext(file_handle, &file_info) == 0);
        _findclose(file_handle);
    }
}

 

// 遍历文件 批量处理 read_as_WIS_d；read_txt_as_WIS_d
void batch_txt2WIS(std::string path, std::string need_extension, int nwf_L_set = 501)
{
    int index = 0;
    std::vector<std::string> my_file;
    //std::string need_extension = ".txt"; 
    std::vector<std::string> folders;
    //std::string path = R"(E:\Proj\vsProj\st_FileSave\数据-原始备份 - 副本)";

    get_sub_folders(path, folders);

    // 打印结果
    for (std::string folder : folders) {
        std::cout << "index:" << index << "      " << folder << std::endl;
        my_file.clear();
        get_need_file(folder, my_file, need_extension);
        for (int i = 0; i < my_file.size(); i++)
        {
            std::cout << "File " << i + 1 << " is:" << std::endl;
            std::cout << my_file[i] << std::endl;

            std::string::size_type pos = my_file[i].find_last_of('/'); 
            std::string filename = my_file[i].substr(pos + 1); 
            std::string foldername = my_file[i].substr(0, pos);
            std::string outfile_SWI = foldername + "/SWI.txt";  //SWI VEP
            std::string outfile_VEP = foldername + "/VEP.txt";  //SWI VEP
            std::cout << filename << std::endl << outfile_SWI 
                << std::endl << outfile_VEP << std::endl;

            if (!filename.compare("reflection_image.txt"))
            {
               //std::cout << "[Test] 相等" << std::endl;
                read_txt_as_WIS_d(my_file[i].c_str(), outfile_SWI.c_str(),513,1);
                // 测试
              /*  if (index == 1)
                {
                    read_txt_as_WIS_d(my_file[i].c_str(), outfile_SWI.c_str());
                }*/
            }
            else if(!filename.compare("velocity_profile.txt"))
            {
                //std::cout << "[Test] 不等" << std::endl;
                read_txt_as_WIS_d(my_file[i].c_str(), outfile_VEP.c_str(), 501,2);
            }
            else
            {
                std::cout << "[Warning] 既不是velocity_profile.txt 也不是reflection_image.txt" << std::endl;
            }
        }
        if (my_file.size() == 0)
        {
            std::cout << "No file can be found!" << std::endl;
        }
        else
        {
            std::cout << "Find " << my_file.size() << " file(s)." << std::endl << std::endl << std::endl;
        }
        index++;
    }

}
int main()
{ 
    /*  
    ImgRes_Pre_XX.dat：513 (double类型)
    其余：501 
    */

    // 数据是float4
    // "E:\Proj\vsProj\st_FileSave\新建文件夹 (2)\ImgRes_Pre_XX.dat"
    // 测试VeloSection dat "E:\Proj\vsProj\st_FileSave\FY3-H3.dat"
    //read_VeloSection_dat("E:\\Proj\\vsProj\\Desktop\\VeloSection.dat", "output_VeloSection.txt");
    //read_VeloSection_dat("E:\\Proj\\vsProj\\st_FileSave\\FY3-H3.dat", "output_FY3-H3.txt");
    
    //数据都是double8
    //read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\ImgRes_Pre_XX.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\output_ImgRes_Pre_XX.txt",513);//  1026/2
    //read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\reflection_image.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\output_reflection_image.txt", 513);//
  
    //二进制读取测试 read_as_WIS_d
    //read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\reflection_image.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\数据\\9-哈得302-H2\\output_reflection_image_b2.dat", 513);//

    //read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\reflection_image.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\output_reflection_image_WIS.dat", 513);//

    // 读取文本文件转为 WIS
 /*   read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\reflection_image.dat",
        "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\SWI.txt",513);*/


    // 下面两个添加了参数mode
    /*read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\reflection_image.txt",
        "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\SWI.txt",513,1);*/

    read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\velocity_profile.txt",
        "E:\\Proj\\vsProj\\st_FileSave\\数据-原始备份 - 副本\\1-富源3-H3\\VEP.txt", 501,2);

    // 遍历修改 
    //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\数据-原始备份 - 副本)", ".txt");

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



   /*
    for (int i = 0; i < 8; i++)
    {
        tstart[i] = 0;
    }
    if (fp1 == NULL) {
        puts("打开文件失败");
        return 0;
    }
    
    fseek(fp1, 0, SEEK_SET); 
    while (fscanf_s(fp1, "%f", &t) != EOF) {
        // process data
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
   
    // 返回值rarvtm
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
        << " dop:" << dop << std::endl;
    */
    // --



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
   
 
      
	return 0;
}
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

//WIS ��ȡ�ļ�ͷ��Ϣ --tmp
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
    // д��txt
    
}

//WIS ��ȡ�������ļ������� --tmp
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
    // ֱ�Ӵ�ӡ�ֽ�
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
   
    // �ȶ�ȡ�ļ�ͷ,��д��txt
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

// ��ȡͷ��Ϣ���� --ͷ�ļ���ʽ����ȷ
void read_VelHeadInfo(FILE* fp, FILE* output, VelHeader& velHeader)
{
    int PREFIX_LEN = 266; 
    char* prefix = (char*)alloca((PREFIX_LEN + 1) * sizeof(float));

    fread(prefix, 1, PREFIX_LEN, fp);
    prefix[PREFIX_LEN] = '\0'; //����ַ���������
    for (int i = 0; i < PREFIX_LEN + 1; i++)
    {
        std::cout << "prefix:" << prefix[i] << std::endl;
    }
    std::cout << "over" << std::endl;

}

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
// velo Section ��ȡ ������float4
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

    // ��txt�ļ�д��
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
        if (index != 0)// ȥ�����п���
        {
            double dep = readdouble8(fp);
            fprintf(output, "\n%.5lf ", dep);
        }
        
        for (int j = 0; j < nwf_L_set; j++)//nwf_L��Ϊ�ⲿ����
        {
            long remaining = fsize - ftell(fp) ;
            // std::cout << remaining << std::endl;
            // 
            // V_Section_HeadInfor.Description = "Velo_Profile"; ������Կ�����
            // �ļ�ĩβȷʵ��ӡ��  Velo_Profile
            if (remaining <= 266) // �ļ�ĩβ����Ϣ������ʵ�ֱ���װ�ˣ��޷��鿴,����೤��֪��
            {// �ļ�ĩβ
               //read_VelHeadInfo(fp, output, velHeaderInfo);
                 
                const int size = remaining;
                char* buff = (char*)alloca(remaining * sizeof(char)) ;

                int len = fread(buff, 1, remaining, fp);

                for (int i = 0; i < len; i++) {
                    printf("%02x ", buff[i]);
                }
               // ֱ�Ӵ�ӡ�ֽ�
                char t;
                while (fscanf_s(fp, "%c", &t) != EOF) { 
                    std::cout << "." << t; 
                } 
                 
                return;
            }
            else
            {
                float velo = readfloat4(fp);
                if (j == nwf_L_set - 1)// ��β��ȥ�ո�
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
// velo Section ȫ��double��
void read_VeloSection_dat_d(const char* fileName, const char* outFielName, int nwf_L_set = 501)
{
    FILE* fp = NULL;
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    } 
    // ��txt�ļ�д��
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
        if (index != 0)// ȥ�����п���
        {
            double dep = readdouble8(fp);
            fprintf(output, "\n%.5lf ", dep);
        }

        for (int j = 0; j < nwf_L_set; j++)//nwf_L��Ϊ�ⲿ����
        {
            long remaining = fsize - ftell(fp);
            // std::cout << remaining << std::endl;
            // 
            // V_Section_HeadInfor.Description = "Velo_Profile"; ������Կ�����
            // �ļ�ĩβȷʵ��ӡ��  Velo_Profile
            if (remaining <= 266) // �ļ�ĩβ����Ϣ������ʵ�ֱ���װ�ˣ��޷��鿴,����೤��֪��
            {// �ļ�ĩβ
               //read_VelHeadInfo(fp, output, velHeaderInfo);

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

                return;
            }
            else
            {
                //float velo = readfloat4(fp);
                double velo = readdouble8(fp);
                if (j == nwf_L_set - 1)// ��β��ȥ�ո�
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

// ��ȡ�ı����ݣ�����ΪWIS��ʽ
void read_txt_as_WIS_d(const char* fileName, const char* outFielName, int nwf_L_set = 501, int mode = 1)
{ 
    std::ifstream ifs(fileName);
    // ��� ifs����
    if (!ifs)
    {
        std::cout << "[ERROR] ifs NOT EXIST" << std::endl;
        return;
    }

    // ��txt�ļ�д��
    FILE* output = NULL;
    //fopen_s(&output, "output.txt", "w");
    //fopen_s(&output, outFielName, "wb");  
    if (fopen_s(&output, outFielName, "wb") != 0) {
        // ���ļ�ʧ��,����
        std::cout << "[ERROR] open file ERROR �� " << outFielName << std::endl;
        return;
    }
    int index = 0;
    double tp_d = 0.0;//��ʱ��ֵ
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
            sprintf_s(name, " SWIS%04d,", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " SWIS%04d", nwf_L_set - 1); //ʹ��sprintf��������ʽ���ַ��� 
        fprintf_s(output, " %s", name);
    }
    else if(mode == 2)
    {
        for (int i = 0; i < nwf_L_set - 1; i++)
        {
            sprintf_s(name, " VEP%04d,", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " VEP%04d", nwf_L_set - 1); //ʹ��sprintf��������ʽ���ַ��� 
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
            sprintf_s(name, " SWIS%04d", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
    }
    else
    {
        for (int i = 0; i < nwf_L_set; i++)
        {
            char name[10];
            sprintf_s(name, " VEP%04d", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
    }
 

  
    index = 0;
    while (std::getline(ifs, line)){
        std::stringstream word(line);  /*�и������׵Ĵ���������ɶ*/
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
            if (word.eof())// ��β��ȥ�ո�
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
    // ��λ�ر���������ǰ


    fseek(output, currPos[0], SEEK_SET); 
    fprintf_s(output, "STDEP = %5.4lf", STDEP);

    fseek(output, currPos[1], SEEK_SET);
    fprintf_s(output, "ENDEP = %5.4lf", ENDEP);

    fseek(output, currPos[2], SEEK_SET);
    fprintf_s(output, "RLEV  = %5.4lf", RLEV);
 

 
    fclose(output);


}

// ��ȡdat���ݣ�����ΪWIS��ʽ ������double8����
void read_as_WIS_d(const char* fileName, const char* outFielName, int nwf_L_set = 501 , int mode = 1)
{
    FILE* fp = NULL;
  
    fopen_s(&fp, fileName, "rb");
    if (fp == NULL)
    {
        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
        return;
    }
    // ��txt�ļ�д��
    FILE* output = NULL;
    //fopen_s(&output, "output.txt", "w");
    fopen_s(&output, outFielName, "wb");
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    VelHeader velHeaderInfo;
    int index = 0;
    double tp_d = 0.0;//��ʱ��ֵ
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
    //    sprintf_s(name, " SWIS%04d,", i); //ʹ��sprintf��������ʽ���ַ��� 
    //    fprintf_s(output, " %s", name);
    //}
    if (mode == 1)
    {
        char name[11];
        for (int i = 0; i < nwf_L_set - 1; i++)
        {
            sprintf_s(name, " SWIS%04d,", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " SWIS%04d", nwf_L_set - 1); //ʹ��sprintf��������ʽ���ַ��� 
        fprintf_s(output, " %s", name);
    }
    else if (mode == 2)
    {
        char name[11];
        for (int i = 0; i < nwf_L_set - 1; i++)
        {
            sprintf_s(name, " VEP%04d,", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
        sprintf_s(name, " VEP%04d", nwf_L_set - 1); //ʹ��sprintf��������ʽ���ַ��� 
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
            sprintf_s(name, " SWIS%04d", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
    }
    else if (mode == 2)
    {
        for (int i = 0; i < nwf_L_set; i++)
        {
            char name[10];
            sprintf_s(name, " VEP%04d", i); //ʹ��sprintf��������ʽ���ַ��� 
            fprintf_s(output, " %s", name);
        }
    }


    //��һ�����
    double dep = readdouble8(fp);
    STDEP = dep;
    fprintf(output, "\n%.5lf ", dep);
    while (!feof(fp))
    {
        std::cout << "dep: " << dep << std::endl;
        if (index != 0)// ȥ�����п���
        {
            //double dep = readdouble8(fp);
            dep = readdouble8(fp);
            ENDEP = dep;
            fprintf(output, "\n%.5lf ", dep);
        }

        for (int j = 0; j < nwf_L_set; j++)//nwf_L��Ϊ�ⲿ����
        {
            long remaining = fsize - ftell(fp);
            // std::cout << remaining << std::endl;
            // 
            // V_Section_HeadInfor.Description = "Velo_Profile"; ������Կ�����
            // �ļ�ĩβȷʵ��ӡ��  Velo_Profile
            if (remaining <= 266) // �ļ�ĩβ����Ϣ������ʵ�ֱ���װ�ˣ��޷��鿴,����೤��֪��
            {// �ļ�ĩβ
               //read_VelHeadInfo(fp, output, velHeaderInfo);

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

                //return;  û�йر��ļ�
                break;
            }
            else
            {
                //float velo = readfloat4(fp);
                double velo = readdouble8(fp);
                if (j == nwf_L_set - 1)// ��β��ȥ�ո�
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
    // ��λ�ر���������ǰ

    
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



 
// ����ָ���ļ����������ļ���
void get_sub_folders(std::string path, std::vector<std::string>& folders) {

    intptr_t handle = 0;
    struct _finddata_t file_info;

    std::string search_path = path + "/*";

    if ((handle = _findfirst(search_path.c_str(), &file_info)) != -1) {

        do {

            if (file_info.attrib & _A_SUBDIR && strcmp(file_info.name, ".") && strcmp(file_info.name, "..")) {

                std::string folder_path = path + "/" + file_info.name;
                folders.push_back(folder_path);

                // �ݹ�������ļ���
                get_sub_folders(folder_path, folders);

            }

        } while (_findnext(handle, &file_info) == 0);

        _findclose(handle);

    }

}

// ��ȡĳһ�ļ����£�ָ���ļ���׺���ļ�
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

 

// �����ļ� �������� read_as_WIS_d��read_txt_as_WIS_d
void batch_txt2WIS(std::string path, std::string need_extension, int nwf_L_set = 501)
{
    int index = 0;
    std::vector<std::string> my_file;
    //std::string need_extension = ".txt"; 
    std::vector<std::string> folders;
    //std::string path = R"(E:\Proj\vsProj\st_FileSave\����-ԭʼ���� - ����)";

    get_sub_folders(path, folders);

    // ��ӡ���
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
               //std::cout << "[Test] ���" << std::endl;
                read_txt_as_WIS_d(my_file[i].c_str(), outfile_SWI.c_str(),513,1);
                // ����
              /*  if (index == 1)
                {
                    read_txt_as_WIS_d(my_file[i].c_str(), outfile_SWI.c_str());
                }*/
            }
            else if(!filename.compare("velocity_profile.txt"))
            {
                //std::cout << "[Test] ����" << std::endl;
                read_txt_as_WIS_d(my_file[i].c_str(), outfile_VEP.c_str(), 501,2);
            }
            else
            {
                std::cout << "[Warning] �Ȳ���velocity_profile.txt Ҳ����reflection_image.txt" << std::endl;
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
    ImgRes_Pre_XX.dat��513 (double����)
    ���ࣺ501 
    */

    // ������float4
    // "E:\Proj\vsProj\st_FileSave\�½��ļ��� (2)\ImgRes_Pre_XX.dat"
    // ����VeloSection dat "E:\Proj\vsProj\st_FileSave\FY3-H3.dat"
    //read_VeloSection_dat("E:\\Proj\\vsProj\\Desktop\\VeloSection.dat", "output_VeloSection.txt");
    //read_VeloSection_dat("E:\\Proj\\vsProj\\st_FileSave\\FY3-H3.dat", "output_FY3-H3.txt");
    
    //���ݶ���double8
    //read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\ImgRes_Pre_XX.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\dat_Files\\output_ImgRes_Pre_XX.txt",513);//  1026/2
    //read_VeloSection_dat_d("E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\reflection_image.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\output_reflection_image.txt", 513);//
  
    //�����ƶ�ȡ���� read_as_WIS_d
    //read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\reflection_image.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\����\\9-����302-H2\\output_reflection_image_b2.dat", 513);//

    //read_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\reflection_image.dat",
    //    "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\output_reflection_image_WIS.dat", 513);//

    // ��ȡ�ı��ļ�תΪ WIS
 /*   read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\reflection_image.dat",
        "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\SWI.txt",513);*/


    // ������������˲���mode
    /*read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\reflection_image.txt",
        "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\SWI.txt",513,1);*/

    read_txt_as_WIS_d("E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\velocity_profile.txt",
        "E:\\Proj\\vsProj\\st_FileSave\\����-ԭʼ���� - ����\\1-��Դ3-H3\\VEP.txt", 501,2);

    // �����޸� 
    //batch_txt2WIS(R"(E:\Proj\vsProj\st_FileSave\����-ԭʼ���� - ����)", ".txt");

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



   /*
    for (int i = 0; i < 8; i++)
    {
        tstart[i] = 0;
    }
    if (fp1 == NULL) {
        puts("���ļ�ʧ��");
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
   
    // ����ֵrarvtm
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
   
 
      
	return 0;
}
// #include<stdio.h> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <string> 
#include <vector>
#include <io.h>

#include <malloc.h>

#include "base.h" // ��������
#include "fst_ttt.h" // ���� fstbrk tttomo
#include "test.h"// ���ڲ���

/*��Ҫ����ʵ�ּ� fst_ttt.cpp �� fst_ttt.h 
�ļ��ṹ��fst_ttt������� fstbrk.cpp �� tttomo.cpp
          test.cpp ���ڲ���
          base.cpp ʵ�ֳ��õĻ�������*/

// ֱ���۵��ⲿ�ִ��뼴�ɣ����������޹�
#pragma region Velo,SWI,WIS ���ݶ�ȡ�ȣ����������޹�
 
//WIS ��ȡ�ļ�ͷ��Ϣ --��Ϣ����ȷ������
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

//WIS ��ȡ�������ļ�������  
//void read_dat(const char* fileName)
//{   
//
//    FILE* fp = NULL; 
//    fopen_s(&fp, fileName, "rb");
//    if (fp == NULL)
//    {
//        std::cout << "[ERROR] fopen_s ERROR" << std::endl;
//        return;
//    }
//    char t; 
//    int index = 0;
//    fseek(fp, 0, SEEK_SET);
//
//    std::cout << std::endl;
//    // ֱ�Ӵ�ӡ�ֽ�
//    //while (fscanf_s(fp, "%c", &t) != EOF) {
//    //    if (index > 10000)
//    //        break;
//    //    std::cout << " " << t;
//    //    //std::cout << t;
//    //    index++;
//    //} 
//    const int size = 1024 * 10;
//    char buff[size];
//
//    int len = fread(buff, 1, size, fp);
//
//    for (int i = 0; i < len; i++) {
//        printf("%02x ", buff[i]);
//    }
//
//    fclose(fp);
//   
//    // �ȶ�ȡ�ļ�ͷ,��д��txt
//    // Read_FileInfor(fp);
//
//    std::cout <<std::endl;
//    std::cout << "[M] read_dat over" << std::endl;
//}

// Velo Section ͷ�ļ���Ϣ
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

// ��ȡͷ��Ϣ���� --ͷ�ļ�д���ʽ����ȷ������
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
 
 
/*
velo Section ���ݶ�ȡ��������float4����
fileName��velo Section�ļ���
outFielName������ļ���
nwf_L_set��ÿ����ȵ����ݸ���
*/
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
/*
��������β _d ����������Ϊ double8��
*/
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

 
/*
* ��ȡ�ı����ݣ�����ΪWIS��ʽ
mode=1 ����SWI
mode=2 ����VEP
*/
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
        // std::cout << "dep: " << dep << std::endl;
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

// ��ȡdat���ݣ�����ΪWIS��ʽ  
void read_as_WIS(const char* fileName, const char* outFielName, int nwf_L_set = 501, int mode = 1)
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
        // std::cout << "dep: " << dep << std::endl;
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
                float velo = readfloat4(fp);
                //double velo = readdouble8(fp);
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
 



    fclose(fp);
    fclose(output);


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
            std::string outfile_SWI_LR = foldername + "/SWI_LR.txt";  //SWI VEP
            std::string outfile_SWI_UD = foldername + "/SWI_UD.txt";  //SWI VEP
            std::string outfile_VEP = foldername + "/VEP.txt";  //SWI VEP
            std::cout << filename << std::endl << outfile_SWI 
                << std::endl << outfile_VEP << std::endl;

   /*         std::string del1 = foldername + "/SWI.txt";
            std::string del2 = foldername + "/SWI.txt_LR";
            std::string del3 = foldername + "/SWI.txt_UD";

            std::cout << "del1:" << del1 << std::endl;

            std::remove(del1.c_str());
            std::remove(del2.c_str());
            std::remove(del3.c_str());*/

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
            else if (!filename.compare("SWI_LR.dat"))
            {  
                read_as_WIS_d(my_file[i].c_str(), outfile_SWI_LR.c_str(), 513, 1);
            }
            else if (!filename.compare("SWI_UD.dat"))
            { 
                read_as_WIS_d(my_file[i].c_str(), outfile_SWI_UD.c_str(), 513, 1);
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
#pragma endregion

int main()
{ 
    // 1.
    // ���� Velo, SWI, WIS���� �����ƻ����ı����ͣ����磺
    // read_as_WIS_d("E:\\����\\reflection_image.dat",
//    "E:����\\output_reflection_image_b2.dat", 513);//
    // ʾ���� test.cpp 

    // 2.��Ҫ����
    //------���� fstbrk tttomo----------------------------------------------------------------------
    /*ע��㣺NMAX*/
    Fst_ttt fst_ttt; 
    
    int test_fa = 2;/*�����������ݣ���һ�鲻����ֱ���ģ��ڶ������ֱ��*/
    switch (test_fa) 
    {
    case 1:
    {
        std::string outFilePath = R"(E:\Proj\vsProj\st_FileSave\��������)";
        // ��ȡ��DTC-slowness���ݣ����� --Caliper
        fst_ttt.read_DTC(R"(E:\Proj\vsProj\st_FileSave\��������\Slowness.txt)", "");
        // ��ȡ���ļ�ͷ��Ϣ���ṹ�� tfw_firstLine
        fst_ttt.read_TFWV_dat(R"(E:\Proj\vsProj\st_FileSave\��������\TFWV01.dat)");
        // ������ʼ��
        fst_ttt.init_par();
        // ��������-�洢ʵ���������ݵ�txt,��������δ����
        fst_ttt.handle_fst(outFilePath);
        break;
    }
    case 2: 
    {
        std::string outFilePath = R"(E:\Proj\vsProj\st_FileSave\��������\file2Capi)";
        fst_ttt.read_DTC(R"(E:\Proj\vsProj\st_FileSave\��������\file2Capi\SLW-Caliper.dat)", "Caliper"); 
        fst_ttt.read_TFWV_dat(R"(E:\Proj\vsProj\st_FileSave\��������\file2Capi\TFWV01.dat)");
        fst_ttt.init_par();
        fst_ttt.handle_fst(outFilePath);
        break;
    }
    default:
        break;
    } 

     
	return 0;
}
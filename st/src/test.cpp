
#include "test.h"
#include "base.h"

#include <fstream>
#include <iostream>

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
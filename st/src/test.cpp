
#include "test.h"
#include "base.h"

#include <fstream>
#include <iostream>

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
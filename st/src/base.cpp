#include "base.h"
#include <vector>
#include <io.h>

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

// 求平均值
double mean_d(const float* array, int& start, int& end)
{
    double sum = 0;
    for (int i = start; i <= end; i++) {
        sum += array[i];
    }
    return sum / (end - start + 1);
}

// 中值滤波
void medianFilter(float* Data , int length, int MFilter) {

    int halfFilterSize = MFilter / 2;
    float* tempData = (float*)malloc(length * sizeof(float));

    float sum = 0;
    for (int i = 0; i < length; i++) 
    {
        if (std::isnan(Data[i]))
        {
            Data[i] = 0;
        } 
        sum += Data[i];
    }
    for (int i = 0; i < length; i++)
    {
        if (Data[i]==0)
        {
            Data[i] = sum/ length;
        }
    }
  

    /* 不改就越界了
    for (int i = 0 + halfFilterSize + 1; i < length - halfFilterSize - 1; i++)
    */
    for (int i = 0; i < length ; i++) { 
        // 复制待滤波数据到临时数组
        for (int j = 0; j < length; j++) {
         /*   if (std::isnan(Data[j]))
                tempData[j] = 0;
            else*/
                tempData[j] = Data[j];
        }

        // 对滑动窗口内的数据 
        /*for (int j = i - halfFilterSize; j <= i + halfFilterSize; j++) {
            if (j < 0) {
                tempData[j + length] = Data[0];
            }
            else if (j >= length) {
                tempData[j - length] = Data[length - 1];
            }
            else {
                tempData[j] = Data[j];
            }
        }*/

        for (int j = 0; j <= MFilter; j++) {
            if (i - halfFilterSize + j < 0) {
                tempData[j] = Data[0];
            }
            else if (i - halfFilterSize + j >= length) {
                tempData[j] = Data[length - 1];
            }
            else {
                tempData[j] = Data[i - halfFilterSize + j];
            }

        }
        // 对滑动窗口内的数据进行冒泡排序
        for (int j = i - halfFilterSize; j <= i + halfFilterSize; j++) {
            for (int k = i - halfFilterSize; k <= i + halfFilterSize; k++) {
       
                if (tempData[(k+length)%length] > tempData[(k + 1 + length)%length]) {
                    float temp = tempData[(k + length) % length];
                    tempData[(k + length) % length] = tempData[(k + 1 + length) % length];
                    tempData[(k + 1 + length) % length] = temp;
                }
            }
        }

        if (i >= 325)
            int tp = 3;
        // 获取中值
        Data[i] = tempData[i];
    }

    free(tempData);
}


//// 中值滤波
//void medianFilter(float* Data, int length, int MFilter) {
//
//    int halfFilterSize = MFilter / 2;
//    float* tempData = (float*)malloc(length * sizeof(float));
//
//    /* 不改就越界了
//    for (int i = 0 + halfFilterSize + 1; i < length - halfFilterSize - 1; i++)
//    */
//    for (int i = 0 + halfFilterSize + 1; i < length - halfFilterSize - 1; i++) {
//        if (i == length - 1)
//            int tp = 2;
//        // 复制待滤波数据到临时数组
//        for (int j = 0; j < length; j++) {
//            tempData[j] = Data[j];
//        }
//
//        // 对滑动窗口内的数据进行排序
//        for (int j = i - halfFilterSize; j <= i + halfFilterSize; j++) {
//            if (j < 0) {
//                tempData[j + length] = Data[0];
//            }
//            else if (j >= length) {
//                tempData[j - length] = Data[length - 1];
//            }
//            else {
//                tempData[j] = Data[j];
//            }
//        }
//
//        // 对滑动窗口内的数据进行冒泡排序
//        for (int j = i - halfFilterSize; j <= i + halfFilterSize; j++) {
//            for (int k = i - halfFilterSize; k <= i + halfFilterSize; k++) {
//                /*        int tk;
//                        if (k >= length)
//                        {
//                            tk = k - length;
//                        }
//                        else if (k < 0)
//                        {
//                            tk = k + length;
//                        }
//                        if (tempData[k] > tempData[k + 1]) {
//                            float temp = tempData[k];
//                            tempData[k] = tempData[k + 1];
//                            tempData[k + 1] = temp;
//                        }*/
//
//
//                if (tempData[k] > tempData[k + 1]) {
//                    float temp = tempData[k];
//                    tempData[k] = tempData[k + 1];
//                    tempData[k + 1] = temp;
//                }
//            }
//        }
//
//        // 获取中值
//        Data[i] = tempData[i];
//    }
//
//    free(tempData);
//}

// 平滑滤波 
// 任意点数加权平均,对各种曲线可以进行平滑处理
void dtfilter(float* Data , int length, int numPoints) {
    float* fltr = (float*)malloc(numPoints * sizeof(float));
    float wt = 0.0;
    for (int i = 0; i < numPoints; i++) {
        fltr[i] = 1.0;  // 每个点的权重都设置为1
        wt += fltr[i];
    }
    float* dtout = (float*)malloc(length * sizeof(float));

    // 初始化滑动窗口
    float sumwt = 0.0;
    for (int i = 0; i < numPoints; i++) {
        sumwt += fltr[i] * Data[i];
    }
    dtout[numPoints / 2] = sumwt / wt;

    // 滑动窗口循环移位更新加权平均结果
    for (int i = numPoints / 2 + 1; i < length - numPoints / 2; i++) {
        sumwt -= fltr[0] * Data[i - numPoints / 2 - 1];
        sumwt += fltr[numPoints - 1] * Data[i + numPoints / 2];
        dtout[i] = sumwt / wt;
    }

    free(fltr);
    free(dtout);
}

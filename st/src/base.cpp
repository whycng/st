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

// ��ƽ��ֵ
double mean_d(const float* array, int& start, int& end)
{
    double sum = 0;
    for (int i = start; i <= end; i++) {
        sum += array[i];
    }
    return sum / (end - start + 1);
}

// ��ֵ�˲�
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
  

    /* ���ľ�Խ����
    for (int i = 0 + halfFilterSize + 1; i < length - halfFilterSize - 1; i++)
    */
    for (int i = 0; i < length ; i++) { 
        // ���ƴ��˲����ݵ���ʱ����
        for (int j = 0; j < length; j++) {
         /*   if (std::isnan(Data[j]))
                tempData[j] = 0;
            else*/
                tempData[j] = Data[j];
        }

        // �Ի��������ڵ����� 
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
        // �Ի��������ڵ����ݽ���ð������
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
        // ��ȡ��ֵ
        Data[i] = tempData[i];
    }

    free(tempData);
}


//// ��ֵ�˲�
//void medianFilter(float* Data, int length, int MFilter) {
//
//    int halfFilterSize = MFilter / 2;
//    float* tempData = (float*)malloc(length * sizeof(float));
//
//    /* ���ľ�Խ����
//    for (int i = 0 + halfFilterSize + 1; i < length - halfFilterSize - 1; i++)
//    */
//    for (int i = 0 + halfFilterSize + 1; i < length - halfFilterSize - 1; i++) {
//        if (i == length - 1)
//            int tp = 2;
//        // ���ƴ��˲����ݵ���ʱ����
//        for (int j = 0; j < length; j++) {
//            tempData[j] = Data[j];
//        }
//
//        // �Ի��������ڵ����ݽ�������
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
//        // �Ի��������ڵ����ݽ���ð������
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
//        // ��ȡ��ֵ
//        Data[i] = tempData[i];
//    }
//
//    free(tempData);
//}

// ƽ���˲� 
// ���������Ȩƽ��,�Ը������߿��Խ���ƽ������
void dtfilter(float* Data , int length, int numPoints) {
    float* fltr = (float*)malloc(numPoints * sizeof(float));
    float wt = 0.0;
    for (int i = 0; i < numPoints; i++) {
        fltr[i] = 1.0;  // ÿ�����Ȩ�ض�����Ϊ1
        wt += fltr[i];
    }
    float* dtout = (float*)malloc(length * sizeof(float));

    // ��ʼ����������
    float sumwt = 0.0;
    for (int i = 0; i < numPoints; i++) {
        sumwt += fltr[i] * Data[i];
    }
    dtout[numPoints / 2] = sumwt / wt;

    // ��������ѭ����λ���¼�Ȩƽ�����
    for (int i = numPoints / 2 + 1; i < length - numPoints / 2; i++) {
        sumwt -= fltr[0] * Data[i - numPoints / 2 - 1];
        sumwt += fltr[numPoints - 1] * Data[i + numPoints / 2];
        dtout[i] = sumwt / wt;
    }

    free(fltr);
    free(dtout);
}

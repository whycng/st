#pragma once
#include <fstream>
#include <vector>
double readdouble8(FILE* fpin);
short int readint2(FILE* fpin);
float readfloat4(FILE* fpin);
int readint4(FILE* fpin);

void get_sub_folders(std::string path, std::vector<std::string>& folders);

void get_need_file(std::string path, std::vector<std::string>& file, std::string ext);

double mean_d(const float* array, int& start, int& end);

void medianFilter(float* Data , int length, int MFilter);

void dtfilter(float Data[], int length, int numPoints);
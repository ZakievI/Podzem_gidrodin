#pragma once
#include <vector>
#include<fstream>

using namespace std;
void MKV(const double R, const double r, const double p_0, const double p_1, const int N);
void out_file(vector<double>& r, vector<double>& p);
void print(vector<double>& r, vector<double>& p);
void print(vector<double>& tetta);
#pragma once
#include<iostream>
#include<iomanip>
#include<string>
#include<fstream>
class Grid
{
public:
	Grid(int row, int col);
	~Grid();
	int index(int row, int col);
	void print();
	void setIMAX(int);
	void setJMAX(int);
	void transfiniteInterpolation();
	void metrics();
	void record(std::string);
	void gridInit();
	int getRow();
	int getCol();
protected:
	int m_row;
	int m_col;

public:
	double* x;
	double* y;
	double* dxdeta;
	double* dydeta;
	double* dxdksi;
	double* dydksi;
	double* J;
	double deta;
	double dksi;
};

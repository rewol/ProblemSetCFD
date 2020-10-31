/****************************************************/
/* Generates a flow square flow field.				*/
/* u, v, p, T, rho, e, Et as primitive properties	*/
/* mu and k as thermophysical properties			*/		
/****************************************************/

#pragma once
#include<string>
#include<fstream>
#include"Altitude.h"
#include<iostream>
#include<iomanip>
class Field2D
{
	friend class Solver;
public:
	Field2D(int, int, float, float, float, float); // Ctor list: nx, ny, L, M, Tw/Tinf, altitude
	~Field2D();
	int index(int, int) const;					// Access matrix elements
	void plotPrimitive() const;				
	void plotThermophysical() const;
	void printPrimitive() const;

protected:
	void fieldInit();					// Initialize the flow field, all zero
	void fieldSet();					// Sets BC and IC for primitives
	void thermophysical();				// Sets viscosity and conductivity

protected:
	double height;
	int m_row;
	int m_col;
	float m_LHORI;
	float m_LVERT;
	float m_dx;
	float m_dy;
	const int size = 5;

public:
	Altitude A;									// Altitude object to set ref values	
	float p_inf;
	float rho_inf;
	float u_inf;
	float T_inf;
	float mu_inf;
	float e_inf;
	float Re_inf;
	float sigma;

	float M;
	float tr;									// Twall/Tinf ratio to calculate wall temperature
	const float m_cv = 287.0f / (1.4f - 1.0f);
	const float m_cp = 1.4f * m_cv;

	float Pr		= 0.71f;
	float* m_u		= nullptr;
	float* m_v		= nullptr;
	float* m_V		= nullptr;
	float* m_T		= nullptr;
	float* m_p		= nullptr;
	float* m_rho	= nullptr;
	float* m_e		= nullptr;
	float* m_Et		= nullptr;
	float* m_mu		= nullptr;
	float* m_k		= nullptr;
	float* x		= nullptr;
	float* y		= nullptr;
};


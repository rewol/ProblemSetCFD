#pragma once
#include "Field2D.h"
#include<iomanip>
#include<iostream>
class Solver
{
public:
	Solver(Field2D* F, float, float);			// Field object, residual, mdot residual
	~Solver();									// Dcon
	int index(int, int) const;	
	void MAC();									// McCormack Scheme

public:											// Printing Functions.
	void printU() const;
	void printE() const;
	void printF() const;
	void plotViscous() const;					// writes down shear stresses and heat conductions
	void printIteration() const;				// Prints iteration data	

private:
	void solverInit();							// Initialize solver parameters
	void setTimestep();							// Sets timestep dt	
	void setTimeStepModified();					// Modified time step !
	void BC();									// Sets boundary conditions
	void viscousSetPredictE();					// Sets TauXX, TauXY, qX for PREDICTION
	void viscousSetPredictF();					// Sets TauXY, TauYY, qY for PREDICTION
	void viscousSetCorrectE();					// Sets TauXX, TauXY, qX for CORRECTION
	void viscousSetCorrectF();					// Sets TauXY, TauYY, qY for CORRECTION
	void subViscous();							// Calculates viscous terms in boundaries except outlet.
	void decodePredicted();						// Decodes primitive props from predicted U
	void decodeCorrected();						// Decodes primitive props from corrected U
	void vectorSetE();							// Sets vector E using U
	void vectorSetEfromPred();					// Sets vector E using predicted U
	void vectorSetFfromPred();					// Sets vector F using predicted U
	void vectorSetF();							// Sets vector F using U
	void vectorSetU();							// Sets vector U using primitives
	void cellReynolds();						// Sets cell Reynolds numbers
	float maxFinder(float*);					// Calculates max value of provided array
	float minFinder(float*);					// Calculates min value of provided array
	void convergence();							// Finds the maximum rho difference between iterations

private:
	Field2D* f = nullptr;						// Pointer to primitive field
	const int size = 5;							// Print distance size
	int m = m_row * m_col;						// Array size
	float m_res;								// Requested Residual NOT IMPLEMENTED !
	float m_mdot_res = 1;						// Mass Flow Check NOT IMPLEMENTED !
	float m_rho_res = 1;						// RHO residual
	int m_col;									// Number of columns
	int m_row;									// Number of raws
	int iteration = 0;
	float m_error;
	float m_mdot_error;	
	float dx, dy;
	float time = 0;
	float dt = 0;
	float K = 0.6;								// Courant Number
	float vprime = 0;	
	
	

private:
	float* m_U1		= nullptr;
	float* m_U2		= nullptr;
	float* m_U3		= nullptr;
	float* m_U5		= nullptr;
	float* m_E1		= nullptr;
	float* m_E2		= nullptr;
	float* m_E3		= nullptr;
	float* m_E5		= nullptr;
	float* m_F1		= nullptr;
	float* m_F2		= nullptr;
	float* m_F3		= nullptr;
	float* m_F5		= nullptr;
	float* m_txx	= nullptr;
	float* m_txy	= nullptr;
	float* m_tyy	= nullptr;
	float* m_qx		= nullptr;
	float* m_qy		= nullptr;
	float* m_Rex	= nullptr;
	float* m_Rey	= nullptr;
	float* a		= nullptr;
	float* vp		= nullptr;
	float* vp_mod	= nullptr;
	float* dt_CFL	= nullptr;
	float* m_U1P	= nullptr;
	float* m_U2P	= nullptr;
	float* m_U3P	= nullptr;
	float* m_U5P	= nullptr;
	float* dU1dt	= nullptr;
	float* dU2dt	= nullptr;
	float* dU3dt	= nullptr;
	float* dU5dt	= nullptr;
	float* dU1Pdt	= nullptr;
	float* dU2Pdt	= nullptr;
	float* dU3Pdt	= nullptr;
	float* dU5Pdt	= nullptr;
	float* dU1avdt	= nullptr;
	float* dU2avdt	= nullptr;
	float* dU3avdt	= nullptr;
	float* dU5avdt	= nullptr;
	float* m_U1n	= nullptr;
	float* m_U2n	= nullptr;
	float* m_U3n	= nullptr;
	float* m_U5n	= nullptr;
	float* m_rho_temp = nullptr;

};


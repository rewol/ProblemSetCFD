#include "Solver.h"
#include<algorithm>
#include<math.h>

Solver::Solver(Field2D* F, float e1, float e2) : m_res(e1), m_mdot_res(e2), m_row(F->m_row), m_col(F->m_col), f(F), dx(F->m_dx), dy(F->m_dy)
{
	iteration		= 0;						// Set Iteration number
	m_error			= 1;
	m_mdot_error		= 1;	
	solverInit();								// Initialize solution field
}

int Solver::index(int row, int col) const
{
	return row + col * m_col;
}

void Solver::solverInit()
{
	m_U1		= new float[m_row * m_col];
	m_U2		= new float[m_row * m_col];
	m_U3		= new float[m_row * m_col];
	m_U5		= new float[m_row * m_col];
	m_E1		= new float[m_row * m_col];
	m_E2		= new float[m_row * m_col];
	m_E3		= new float[m_row * m_col];
	m_E5		= new float[m_row * m_col];
	m_F1		= new float[m_row * m_col];
	m_F2		= new float[m_row * m_col];
	m_F3		= new float[m_row * m_col];
	m_F5		= new float[m_row * m_col];
	m_txx		= new float[m_row * m_col];
	m_tyy		= new float[m_row * m_col];
	m_txy		= new float[m_row * m_col];
	m_qx		= new float[m_row * m_col];
	m_qy		= new float[m_row * m_col];
	m_Rex		= new float[m_row * m_col];
	m_Rey		= new float[m_row * m_col];
	a		= new float[m_row * m_col];
	vp		= new float[m_row * m_col];
	vp_mod		= new float[m_row * m_col];
	dt_CFL		= new float[m_row * m_col];
	m_U1P		= new float[m_row * m_col];
	m_U2P		= new float[m_row * m_col];
	m_U3P		= new float[m_row * m_col];
	m_U5P		= new float[m_row * m_col];
	dU1dt		= new float[m_row * m_col];
	dU2dt		= new float[m_row * m_col];
	dU3dt		= new float[m_row * m_col];
	dU5dt		= new float[m_row * m_col];
	dU1Pdt		= new float[m_row * m_col];
	dU2Pdt		= new float[m_row * m_col];
	dU3Pdt		= new float[m_row * m_col];
	dU5Pdt		= new float[m_row * m_col];
	dU1avdt		= new float[m_row * m_col];
	dU2avdt		= new float[m_row * m_col];
	dU3avdt		= new float[m_row * m_col];
	dU5avdt		= new float[m_row * m_col];
	m_U1n		= new float[m_row * m_col];
	m_U2n		= new float[m_row * m_col];
	m_U3n		= new float[m_row * m_col];
	m_U5n		= new float[m_row * m_col];
	m_rho_temp	= new float[m_row * m_col];
	
	std::fill_n(m_U1, m_row * m_col, 0);
	std::fill_n(m_U2, m_row * m_col, 0);
	std::fill_n(m_U3, m_row * m_col, 0);
	std::fill_n(m_U5, m_row * m_col, 0);
	std::fill_n(m_E1, m_row * m_col, 0);
	std::fill_n(m_E2, m_row * m_col, 0);
	std::fill_n(m_E3, m_row * m_col, 0);
	std::fill_n(m_E5, m_row * m_col, 0);
	std::fill_n(m_F1, m_row * m_col, 0);
	std::fill_n(m_F2, m_row * m_col, 0);
	std::fill_n(m_F3, m_row * m_col, 0);
	std::fill_n(m_F5, m_row * m_col, 0);
	std::fill_n(m_Rex, m_row * m_col, 0);
	std::fill_n(m_Rey, m_row * m_col, 0);
	std::fill_n(m_txy, m_row * m_col, 0);
	std::fill_n(m_txx, m_row * m_col, 0);
	std::fill_n(m_tyy, m_row * m_col, 0);
	std::fill_n(m_qx, m_row * m_col, 0);
	std::fill_n(m_qy, m_row * m_col, 0);
	std::fill_n(a, m_row * m_col, 0);
	std::fill_n(vp, m_row * m_col, 0);
	std::fill_n(vp_mod, m_row * m_col, 0);
	std::fill_n(dt_CFL, m_row * m_col, 1);
	std::fill_n(m_U1P, m_row * m_col, 0);
	std::fill_n(m_U2P, m_row * m_col, 0);
	std::fill_n(m_U3P, m_row * m_col, 0);
	std::fill_n(m_U5P, m_row * m_col, 0);
	std::fill_n(dU1dt, m_row * m_col, 0);
	std::fill_n(dU2dt, m_row * m_col, 0);
	std::fill_n(dU3dt, m_row * m_col, 0);
	std::fill_n(dU5dt, m_row * m_col, 0);
	std::fill_n(dU1Pdt, m_row * m_col, 0);
	std::fill_n(dU2Pdt, m_row * m_col, 0);
	std::fill_n(dU3Pdt, m_row * m_col, 0);
	std::fill_n(dU5Pdt, m_row * m_col, 0);
	std::fill_n(dU1avdt, m_row * m_col, 0);
	std::fill_n(dU2avdt, m_row * m_col, 0);
	std::fill_n(dU3avdt, m_row * m_col, 0);
	std::fill_n(dU5avdt, m_row * m_col, 0);
	std::fill_n(m_U1n, m_row * m_col, 0);
	std::fill_n(m_U2n, m_row * m_col, 0);
	std::fill_n(m_U3n, m_row * m_col, 0);
	std::fill_n(m_U5n, m_row * m_col, 0);
	std::fill_n(m_rho_temp, m_row * m_col, 0);
	
}

Solver::~Solver()
{
	delete[] m_U1;
	delete[] m_U2;
	delete[] m_U3;
	delete[] m_U5;
	delete[] m_E1;
	delete[] m_E2;
	delete[] m_E3;
	delete[] m_E5;
	delete[] m_F1;
	delete[] m_F2;
	delete[] m_F3;
	delete[] m_F5;
	delete[] m_txx;
	delete[] m_tyy;
	delete[] m_txy;
	delete[] m_qx;
	delete[] m_qy;
	delete[] m_Rex;
	delete[] m_Rey;
	delete[] a;
	delete[] vp;
	delete[] vp_mod;
	delete[] dt_CFL;
	delete[] m_U1P;
	delete[] m_U2P;
	delete[] m_U3P;
	delete[] m_U5P;
	delete[] dU1dt;
	delete[] dU2dt;
	delete[] dU3dt;
	delete[] dU5dt;
	delete[] dU1Pdt;
	delete[] dU2Pdt;
	delete[] dU3Pdt;
	delete[] dU5Pdt;
	delete[] dU1avdt;
	delete[] dU2avdt;
	delete[] dU3avdt;
	delete[] dU5avdt;
	delete[] m_U1n;
	delete[] m_U2n;
	delete[] m_U3n;
	delete[] m_U5n;

}

void Solver::MAC()
{
	int maxiter = 10000;
	f->thermophysical();		// Set mu and k
	vectorSetU();				// Set U for all field
	//f->plotPrimitive();

	while (iteration < maxiter && m_rho_res > 1e-6)
	{
		// Set timestep 
		setTimestep();
		//std::cout << dt << std::endl;
		// Copy density field for convergence check
		for (int i = 0; i < m_col; i++)
		{
			for (int j = 0; j < m_row; j++)
			{
				m_rho_temp[index(i, j)] = f->m_rho[index(i, j)];
			}
		}
		viscousSetPredictE();
		vectorSetE();
		viscousSetPredictF();
		vectorSetF();		

		//1. Calculate dUdt for interior points. Exterior dUdt values will be zero and not be utilized in the solution.
		for (int i = 1; i < m_col - 1; i++)
		{
			for (int j = 1; j < m_row - 1; j++)
			{
				dU1dt[index(i, j)] = -1 * (m_E1[index(i + 1, j)] - m_E1[index(i, j)]) / dx - (m_F1[index(i, j + 1)] - m_F1[index(i, j)]) / dy;
				dU2dt[index(i, j)] = -1 * (m_E2[index(i + 1, j)] - m_E2[index(i, j)]) / dx - (m_F2[index(i, j + 1)] - m_F2[index(i, j)]) / dy;
				dU3dt[index(i, j)] = -1 * (m_E3[index(i + 1, j)] - m_E3[index(i, j)]) / dx - (m_F3[index(i, j + 1)] - m_F3[index(i, j)]) / dy;
				dU5dt[index(i, j)] = -1 * (m_E5[index(i + 1, j)] - m_E5[index(i, j)]) / dx - (m_F5[index(i, j + 1)] - m_F5[index(i, j)]) / dy;
			}
		}

		// Calculate inner predicted U vectors (ALL INTERIOR PREDICTED)
		for (int i = 1; i < m_col - 1; i++)
		{
			for (int j = 1; j < m_row - 1; j++)
			{
				m_U1P[index(i, j)] = m_U1[index(i, j)] + dU1dt[index(i, j)] * dt;
				m_U2P[index(i, j)] = m_U2[index(i, j)] + dU2dt[index(i, j)] * dt;
				m_U3P[index(i, j)] = m_U3[index(i, j)] + dU3dt[index(i, j)] * dt;
				m_U5P[index(i, j)] = m_U5[index(i, j)] + dU5dt[index(i, j)] * dt;
			}
		}

		// Decode for predicted primitives from predicted vectors in internal field.
		decodePredicted();
		// Set Boundary conditions, ALL PRIMITIVE FIELD IS PREDICTED. NOT VECTORS.
		BC();		
		// Set viscosity and thermal conductivity
		f->thermophysical();		
		
		// Calculate viscous terms with predicted primitives
		viscousSetCorrectE();
		vectorSetEfromPred();
		viscousSetCorrectF();
		vectorSetFfromPred();
				
		// 2. Calculate predicted dUdt for interior points.
		for (int i = 1; i < m_col - 1; i++)
		{
			for (int j = 1; j < m_row - 1; j++)
			{
				dU1Pdt[index(i, j)] = -1 * (m_E1[index(i, j)] - m_E1[index(i - 1, j)]) / dx - (m_F1[index(i, j)] - m_F1[index(i, j - 1)]) / dy;
				dU2Pdt[index(i, j)] = -1 * (m_E2[index(i, j)] - m_E2[index(i - 1, j)]) / dx - (m_F2[index(i, j)] - m_F2[index(i, j - 1)]) / dy;
				dU3Pdt[index(i, j)] = -1 * (m_E3[index(i, j)] - m_E3[index(i - 1, j)]) / dx - (m_F3[index(i, j)] - m_F3[index(i, j - 1)]) / dy;
				dU5Pdt[index(i, j)] = -1 * (m_E5[index(i, j)] - m_E5[index(i - 1, j)]) / dx - (m_F5[index(i, j)] - m_F5[index(i, j - 1)]) / dy;
			}
		}

		// 3. Calculate averaged dUdt to calculate U at t + dt
		for (int i = 1; i < m_col - 1; i++)
		{
			for (int j = 1; j < m_row - 1; j++)
			{
				dU1avdt[index(i, j)] = 0.5f * (dU1dt[index(i, j)] + dU1Pdt[index(i, j)]);
				dU2avdt[index(i, j)] = 0.5f * (dU2dt[index(i, j)] + dU2Pdt[index(i, j)]);
				dU3avdt[index(i, j)] = 0.5f * (dU3dt[index(i, j)] + dU3Pdt[index(i, j)]);
				dU5avdt[index(i, j)] = 0.5f * (dU5dt[index(i, j)] + dU5Pdt[index(i, j)]);
			}
		}

		// 4. Calculate corrected U vector at time d + dt
		// INTERIOR U IS CORRECTED.
		for (int i = 1; i < m_col - 1; i++)
		{
			for (int j = 1; j < m_row - 1; j++)
			{
				m_U1n[index(i, j)] = m_U1[index(i, j)] + dU1avdt[index(i, j)] * dt;
				m_U2n[index(i, j)] = m_U2[index(i, j)] + dU2avdt[index(i, j)] * dt;
				m_U3n[index(i, j)] = m_U3[index(i, j)] + dU3avdt[index(i, j)] * dt;
				m_U5n[index(i, j)] = m_U5[index(i, j)] + dU5avdt[index(i, j)] * dt;
			}
		}

		decodeCorrected();
		BC();		
		convergence();
		f->thermophysical();

		//Copy new U vector values into old values.
		for (int i = 1; i < m_col - 1; i++)
		{
			for (int j = 1; j < m_row - 1; j++)
			{
				m_U1[index(i, j)] = m_U1n[index(i, j)];
				m_U2[index(i, j)] = m_U2n[index(i, j)];
				m_U3[index(i, j)] = m_U3n[index(i, j)];
				m_U5[index(i, j)] = m_U5n[index(i, j)];
			}
		}
		
		// Set CORRECTED U for BOUNDARIES
		// LEFT
		for (int j = 0; j < m_row; j++)
		{
			m_U1[index(0, j)] = f->m_rho[index(0, j)];
			m_U2[index(0, j)] = f->m_rho[index(0, j)] * f->m_u[index(0, j)];
			m_U3[index(0, j)] = f->m_rho[index(0, j)] * f->m_v[index(0, j)];
			m_U5[index(0, j)] = f->m_Et[index(0, j)];
		}
		// RIGHT
		for (int j = 0; j < m_row; j++)
		{
			m_U1[index(m_col - 1, j)] = f->m_rho[index(m_col - 1, j)];
			m_U2[index(m_col - 1, j)] = f->m_rho[index(m_col - 1, j)] * f->m_u[index(m_col - 1, j)];
			m_U3[index(m_col - 1, j)] = f->m_rho[index(m_col - 1, j)] * f->m_v[index(m_col - 1, j)];
			m_U5[index(m_col - 1, j)] = f->m_Et[index(m_col - 1, j)];
		}
		// DOWN
		for (int i = 1; i < m_col - 1; i++)
		{
			m_U1[index(i, 0)] = f->m_rho[index(i, 0)];
			m_U2[index(i, 0)] = f->m_rho[index(i, 0)] * f->m_u[index(i, 0)];
			m_U3[index(i, 0)] = f->m_rho[index(i, 0)] * f->m_v[index(i, 0)];
			m_U5[index(i, 0)] = f->m_Et[index(i, 0)];
		}
		// UP
		for (int i = 1; i < m_col - 1; i++)
		{
			m_U1[index(i, m_row - 1)] = f->m_rho[index(i, m_row - 1)];
			m_U2[index(i, m_row - 1)] = f->m_rho[index(i, m_row - 1)] * f->m_u[index(i, m_row - 1)];
			m_U3[index(i, m_row - 1)] = f->m_rho[index(i, m_row - 1)] * f->m_v[index(i, m_row - 1)];
			m_U5[index(i, m_row - 1)] = f->m_Et[index(i, m_row - 1)];
		}
		
		
		// Update iteration
		iteration += 1;
		time += dt;
		// Call print functions
		printIteration();		
	}
	
}

float Solver::maxFinder(float* residual)
{
	float maxVal = 0;
	for (int j = 1; j < m_row - 1; j++)
	{
		for (int i = 1; i < m_col - 1; i++)
		{
			if(residual[index(i, j)] > maxVal)
				{
					maxVal = residual[index(i, j)];
				}
		}
	}

	return maxVal;
}

float Solver::minFinder(float* timestep)
{
	float minVal = 1;
	for (int j = 1; j < m_row - 1; j++)
	{
		for (int i = 1; i < m_col - 1; i++)
		{
			if (timestep[index(i, j)] < minVal)
			{
				minVal = timestep[index(i, j)];
			}
		}
	}

	return minVal;
}

void Solver::decodePredicted()
{
	//decode for internal fields 
	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			f->m_rho[index(i, j)] = m_U1P[index(i, j)];
			f->m_u[index(i, j)] = m_U2P[index(i, j)] / m_U1P[index(i, j)];
			f->m_v[index(i, j)] = m_U3P[index(i, j)] / m_U1P[index(i, j)];
			f->m_V[index(i, j)] = pow(f->m_u[index(i, j)] * f->m_u[index(i, j)] + f->m_v[index(i, j)] * f->m_v[index(i, j)], 0.5);
			f->m_e[index(i, j)] = m_U5P[index(i, j)] / m_U1P[index(i, j)] - pow(f->m_V[index(i, j)], 2) / 2;
			f->m_Et[index(i, j)] = m_U5P[index(i, j)];
			f->m_T[index(i, j)] = f->m_e[index(i, j)] / f->m_cv;
			f->m_p[index(i, j)] = f->m_rho[index(i, j)] * 287 * f->m_T[index(i, j)];
		}
	}
}

void Solver::decodeCorrected()
{
	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			f->m_rho[index(i, j)] = m_U1n[index(i, j)];
			f->m_u[index(i, j)] = m_U2n[index(i, j)] / m_U1n[index(i, j)];
			f->m_v[index(i, j)] = m_U3n[index(i, j)] / m_U1n[index(i, j)];
			f->m_V[index(i, j)] = pow(f->m_u[index(i, j)] * f->m_u[index(i, j)] + f->m_v[index(i, j)] * f->m_v[index(i, j)], 0.5);
			f->m_e[index(i, j)] = m_U5n[index(i, j)] / m_U1n[index(i, j)] - pow(f->m_V[index(i, j)], 2) / 2;
			f->m_Et[index(i, j)] = m_U5n[index(i, j)];
			f->m_T[index(i, j)] = f->m_e[index(i, j)] / f->m_cv;
			f->m_p[index(i, j)] = f->m_rho[index(i, j)] * 287 * f->m_T[index(i, j)];
		}
	}
}

void Solver::BC()
{
	// First set the primitives on outlet (except wall and upper boundary)
	for (int j = 1; j < m_row - 1; j++)
	{
		f->m_u[index(m_col - 1, j)] = 2 * f->m_u[index(m_col - 2, j)] - f->m_u[index(m_col - 3, j)];
		f->m_v[index(m_col - 1, j)] = 2 * f->m_v[index(m_col - 2, j)] - f->m_v[index(m_col - 3, j)];
		f->m_p[index(m_col - 1, j)] = 2 * f->m_p[index(m_col - 2, j)] - f->m_p[index(m_col - 3, j)];
		f->m_T[index(m_col - 1, j)] = 2 * f->m_T[index(m_col - 2, j)] - f->m_T[index(m_col - 3, j)];
		f->m_rho[index(m_col - 1, j)] = f->m_p[index(m_col - 1, j)] / (287.0f * f->m_T[index(m_col - 1, j)]);
		f->m_e[index(m_col - 1, j)] = f->m_T[index(m_col - 1, j)] * f->m_cv;
		f->m_V[index(m_col - 1, j)] = pow(f->m_u[index(m_col - 1, j)] * f->m_u[index(m_col - 1, j)] + f->m_v[index(m_col - 1, j)] * f->m_v[index(m_col - 1, j)], 0.5);
		f->m_Et[index(m_col - 1, j)] = f->m_rho[index(m_col - 1, j)] * (f->m_e[index(m_col - 1, j)] + pow(f->m_V[index(m_col - 1, j)], 2) / 2);
	}
	// Now set the primitives on the wall except Leading Edge. Set only p and and calculate rho with ideal gas relation.
	for (int i = 1; i < m_col; i++)
	{
		f->m_p[index(i, 0)] = 2 * f->m_p[index(i, 2)] - f->m_p[index(i, 1)];
		f->m_rho[index(i, 0)] = f->m_p[index(i, 0)] / (287.0f * f->m_T[index(i, 0)]);
		f->m_Et[index(i, 0)] = f->m_rho[index(i, 0)] * f->m_e[index(i, 0)];
	}
}

void Solver::setTimestep()
{
	// First local speed of sound in the flow field is calculated.
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			a[index(i, j)] = pow(1.4 * 287 * f->m_T[index(i, j)], 0.5);
		}
	}

	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			vp[index(i, j)] = 1.33 * f->m_mu[index(i, j)] * (1.4 * f->m_mu[index(i, j)] / f->Pr) / f->m_rho[index(i, j)];
		}
	}

	//vprime = *std::max_element(vp, vp + n);
	vprime = maxFinder(vp);
	//std::cout << "vprime :" << vprime << std::endl;

	float term1 = pow((1 / (dx * dx)) + (1 / (dy * dy)), 0.5);
	//std::cout << "term1 :" << term1 << std::endl;
	float term2 = (1 / (dx * dx)) + (1 / (dy * dy));
	//std::cout << "term2 :" << term2 << std::endl;

	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			dt_CFL[index(i, j)] = K * pow(fabs(f->m_u[index(i, j)] / dx) + fabs(f->m_v[index(i, j)] / dy) + a[index(i, j)] * term1 + 2 * vprime * term2, -1);
			//std::cout << dt_CFL[index(i, j)] << std::endl;
		}
	}

	//dt = *std::min_element(dt_CFL, dt_CFL + n);
	dt = minFinder(dt_CFL);
	//std::cout << "Timestep dt is: " << dt << std::endl;

		
}

void Solver::setTimeStepModified() // DO NOT USE THIS ONE
{
	
	float max;
	// First local speed of sound in the flow field is calculated.
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			a[index(i, j)] = pow(1.4 * 287 * f->m_T[index(i, j)], 0.5);
		}
	}
}

void Solver::viscousSetPredictE()
{
	subViscous();
	// Internal Field, where x derivatives will be BWD and y derivates will be CNT differenced.
	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			m_txy[index(i, j)] = f->m_mu[index(i, j)] * ((f->m_u[index(i, j + 1)] - f->m_u[index(i, j - 1)]) / (2 * dy)
				+ (f->m_v[index(i, j)] - f->m_v[index(i - 1, j)]) / dx);

			m_txx[index(i, j)] = (4 / 3) * f->m_mu[index(i, j)] * (f->m_u[index(i, j)] - f->m_u[index(i - 1, j)]) / dx
				- (2 / 3) * f->m_mu[index(i, j)] * (f->m_v[index(i, j + 1)] - f->m_v[index(i, j - 1)]) / (2 * dy);

			m_qx[index(i, j)] = -1.0f * f->m_k[index(i, j)] * (f->m_T[index(i, j)] - f->m_T[index(i - 1, j)]) / dx;
		}
	}

	// OUTLET
	for (int j = 1; j < m_row - 1; j++)
	{
		m_txy[index(m_col - 1, j)] = f->m_mu[index(m_col - 1, j)] * ((f->m_v[index(m_col - 1, j)] - f->m_v[index(m_col - 2, j)]) / dx
			+ (f->m_u[index(m_col - 1, j + 1)] - f->m_u[index(m_col - 1, j - 1)]) / (2 * dy));
		m_txx[index(m_col - 1, j)] = (4 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_u[index(m_col - 1, j)] - f->m_u[index(m_col-2, j)]) / dx
			- (2 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_v[index(m_col - 1, j + 1)] - f->m_v[index(m_col - 1, j - 1)]) / (2 * dy);
		m_qx[index(m_col - 1, j)] = -1.0f * f->m_k[index(m_col - 1, j)] * (f->m_T[index(m_col - 1, j)] - f->m_T[index(m_col - 2, j)]) / dx;
	}

}

void Solver::viscousSetPredictF()
{
	//subViscous();
	// Internal Field, where x in central, y in BWD
	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			m_txy[index(i, j)] = f->m_mu[index(i, j)] * ((f->m_u[index(i, j)] - f->m_u[index(i, j - 1)]) / dy
				+ (f->m_v[index(i + 1, j)] - f->m_v[index(i - 1, j)]) / (2 * dx));

			m_tyy[index(i, j)] = (4 / 3) * f->m_mu[index(i, j)] * (f->m_v[index(i, j)] - f->m_v[index(i, j - 1)]) / dy
				- (2 / 3) * f->m_mu[index(i, j)] * (f->m_u[index(i + 1, j)] - f->m_u[index(i - 1, j)]) / (2 * dx);

			m_qy[index(i, j)] = -1.0f * f->m_k[index(i, j)] * (f->m_T[index(i, j)] - f->m_T[index(i, j - 1)]) / dy;
		}
	}

	// OUTLET
	for (int j = 1; j < m_row - 1; j++)
	{
		m_txy[index(m_col - 1, j)] = f->m_mu[index(m_col - 1, j)] * ((f->m_v[index(m_col - 1, j)] - f->m_v[index(m_col - 2, j)]) / dx
			+ (f->m_u[index(m_col - 1, j)] - f->m_u[index(m_col - 1, j - 1)]) / dy);
		m_tyy[index(m_col - 1, j)] = (4 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_v[index(m_col - 1, j)] - f->m_v[index(m_col - 1, j - 1)]) / dy
			- (2 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_u[index(m_col - 1, j)] - f->m_u[index(m_col - 2, j)]) / dx;
		m_qy[index(m_col - 1, j)] = -1.0f * f->m_k[index(m_col - 1, j)] * (f->m_T[index(m_col - 1, j)] - f->m_T[index(m_col - 1, j - 1)]) / dy;
	}

}

void Solver::viscousSetCorrectE()
{
	subViscous();
	// Internal Field, where x derivatives will be FWD and y derivates will be CNT differenced.
	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			m_txy[index(i, j)] = f->m_mu[index(i, j)] * ((f->m_u[index(i, j + 1)] - f->m_u[index(i, j - 1)]) / (2 * dy)
				+ (f->m_v[index(i + 1, j)] - f->m_v[index(i, j)]) / dx);

			m_txx[index(i, j)] = (4 / 3) * f->m_mu[index(i, j)] * (f->m_u[index(i + 1, j)] - f->m_u[index(i, j)]) / dx
				- (2 / 3) * f->m_mu[index(i, j)] * (f->m_v[index(i, j + 1)] - f->m_v[index(i, j - 1)]) / (2 * dy);

			m_qx[index(i, j)] = -1.0f * f->m_k[index(i, j)] * (f->m_T[index(i + 1, j)] - f->m_T[index(i, j)]) / dx;
		}
	}

	// OUTLET
	for (int j = 1; j < m_row - 1; j++)
	{
		m_txy[index(m_col - 1, j)] = f->m_mu[index(m_col - 1, j)] * ((f->m_v[index(m_col - 1, j)] - f->m_v[index(m_col - 2, j)]) / dx
			+ (f->m_u[index(m_col - 1, j + 1)] - f->m_u[index(m_col - 1, j - 1)]) / (2 * dy));
		m_txx[index(m_col - 1, j)] = (4 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_u[index(m_col - 1, j)] - f->m_u[index(m_col - 2, j)]) / dx
			- (2 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_v[index(m_col - 1, j + 1)] - f->m_v[index(m_col - 1, j - 1)]) / (2 * dy);
		m_qx[index(m_col - 1, j)] = -1.0f * f->m_k[index(m_col - 1, j)] * (f->m_T[index(m_col - 1, j)] - f->m_T[index(m_col - 2, j)]) / dx;
	}
	
}

void Solver::viscousSetCorrectF() 
{
	//subViscous();
	// Internal Field, where x in central, y in FWD
	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			m_txy[index(i, j)] = f->m_mu[index(i, j)] * ((f->m_u[index(i, j + 1)] - f->m_u[index(i, j)]) / dy
				+ (f->m_v[index(i + 1, j)] - f->m_v[index(i - 1, j)]) / (2 * dx));

			m_tyy[index(i, j)] = (4 / 3) * f->m_mu[index(i, j)] * (f->m_v[index(i, j + 1)] - f->m_v[index(i, j)]) / dy
				- (2 / 3) * f->m_mu[index(i, j)] * (f->m_u[index(i + 1, j)] - f->m_u[index(i - 1, j)]) / (2 * dx);

			m_qy[index(i, j)] = -1.0f * f->m_k[index(i, j)] * (f->m_T[index(i, j + 1)] - f->m_T[index(i, j)]) / dy;
		}
	}

	// OUTLET
	for (int j = 1; j < m_row - 1; j++)
	{
		m_txy[index(m_col - 1, j)] = f->m_mu[index(m_col - 1, j)] * ((f->m_v[index(m_col - 1, j)] - f->m_v[index(m_col - 2, j)]) / dx
			+ (f->m_u[index(m_col - 1, j + 1)] - f->m_u[index(m_col - 1, j)]) / dy);
		m_tyy[index(m_col - 1, j)] = (4 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_v[index(m_col - 1, j + 1)] - f->m_v[index(m_col - 1, j)]) / dy
			- (2 / 3) * f->m_mu[index(m_col - 1, j)] * (f->m_u[index(m_col - 1, j)] - f->m_u[index(m_col - 2, j)]) / dx;
		m_qy[index(m_col - 1, j)] = -1.0f * f->m_k[index(m_col - 1, j)] * (f->m_T[index(m_col - 1, j + 1)] - f->m_T[index(m_col - 1, j)]) / dy;
	}
}

void Solver::subViscous()
{
	// LOWER LEFT CORNER: LEADING EDGE
	// dv/dx = 0
	m_txy[index(0, 0)] = f->m_mu[index(0, 0)] * (f->m_u[index(0, 1)] - f->m_u[index(0, 0)]) / dy;
	// du/dx = 0, dv/dy = 0
	m_txx[index(0, 0)] = 0;
	m_tyy[index(0, 0)] = 0;
	m_qx[index(0, 0)] = 0;
	m_qy[index(0, 0)] = 0;
	// INLET
	for (int j = 1; j < m_row - 1; j++)
	{
		// du/dy = 0
		m_txy[index(0, j)] = f->m_mu[index(0, j)] * (f->m_v[index(1, j)] - f->m_v[index(0, j)]) / dx;
		// dv/dy = 0
		m_txx[index(0, j)] = (4/3) * f->m_mu[index(0, j)] * (f->m_u[index(1, j)] - f->m_u[index(0, j)]) / dx;
		m_tyy[index(0, j)] = (-2/3) * f->m_mu[index(0, j)] * (f->m_u[index(1, j)] - f->m_u[index(0, j)]) / dx;
		m_qx[index(0, j)] = -f->m_k[index(0, j)] * (f->m_T[index(1, j)] - f->m_T[index(0, j)]) / dx;
		m_qy[index(0, j)] = 0;
	}
	// UPPER LEFT CORNER
	// du/dy = 0, dv/dx = 0, du/dx = 0 dv/dy = 0
	m_txy[index(0, m_row - 1)] = 0;
	m_txx[index(0, m_row - 1)] = 0;
	m_tyy[index(0, m_row - 1)] = 0;
	m_qy[index(0, m_row - 1)] = 0;
	m_qy[index(0, m_row - 1)] = 0;
	// WALL
	for (int i = 1; i < m_col - 1; i++)
	{
		// dv/dx = 0
		m_txy[index(i, 0)] = f->m_mu[index(i, 0)] * (f->m_u[index(i, 1)] - f->m_u[index(i, 0)]) / dy;
		// du/dx = 0
		m_txx[index(i, 0)] = (-2 / 3) * f->m_mu[index(i, 0)] * (f->m_v[index(i, 1)] - f->m_v[index(i, 0)]) / dy;
		m_tyy[index(i, 0)] = (4 / 3) * f->m_mu[index(i, 0)] * (f->m_v[index(i, 1)] - f->m_v[index(i, 0)]) / dy;
		m_qx[index(i, 0)] = 0;
		m_qy[index(i, 0)] = -f->m_k[index(i, 0)] * (f->m_T[index(i, 1)] - f->m_T[index(i, 0)]) / dy;
	}
	// UPPER SIDE
	for (int i = 1; i < m_col - 1; i++)
	{
		// dv/dx = 0
		m_txy[index(i, m_row - 1)] = f->m_mu[index(i, m_row - 1)]  * (f->m_u[index(i, m_row - 1)] - f->m_u[index(i, m_row - 2)]) / dy;
		// du/dx = 0
		m_txx[index(i, m_row - 1)] = (-2 / 3) * f->m_mu[index(i, m_row - 1)] * (f->m_v[index(i, m_row - 1)] - f->m_v[index(i, m_row - 2)]) / dy;
		m_tyy[index(i, m_row - 1)] = (4 / 3) * f->m_mu[index(i, m_row - 1)] * (f->m_v[index(i, m_row - 1)] - f->m_v[index(i, m_row - 2)]) / dy;
		m_qx[index(i, m_row - 1)] = 0;
		m_qy[index(i, m_row - 1)] = -f->m_k[index(i, m_row - 1)] * (f->m_T[index(i, m_row - 1)] - f->m_T[index(i, m_row - 2)]) / dy;
	}
	// UPPER RIGHT CORNER
	m_txy[index(m_col - 1, m_row - 1)] = f->m_mu[index(m_col - 1, m_row - 1)] * (f->m_u[index(m_col - 1, m_row - 1)] - f->m_u[index(m_col - 1, m_row - 2)]) / dy;
	m_txx[index(m_col - 1, m_row - 1)] = (-2 / 3) * f->m_mu[index(m_col - 1, m_row - 1)] * (f->m_v[index(m_col - 1, m_row - 1)] - f->m_v[index(m_col - 1, m_row - 2)]) / dy;
	m_tyy[index(m_col - 1, m_row - 1)] = (4 / 3) * f->m_mu[index(m_col - 1, m_row - 1)] * (f->m_v[index(m_col - 1, m_row - 1)] - f->m_v[index(m_col - 1, m_row - 2)]) / dy;
	m_qx[index(m_col - 1, m_row - 1)] = 0;
	m_qy[index(m_col - 1, m_row - 1)] = -f->m_k[index(m_col - 1, m_row - 1)] * (f->m_T[index(m_col - 1, m_row - 1)] - f->m_T[index(m_col - 1, m_row - 2)]) / dy;
	// LOWER RIGHT CORNER
	m_txy[index(m_col - 1, 0)] = f->m_mu[index(m_col - 1, 0)] * (f->m_u[index(m_col - 1, 1)] - f->m_u[index(m_col - 1, 0)]) / dy;
	m_txx[index(m_col - 1, 0)] = (-2 / 3) * f->m_mu[index(m_col - 1, 0)] * (f->m_v[index(m_col - 1, 1)] - f->m_v[index(m_col - 1, 0)]) / dy;
	m_tyy[index(m_col - 1, 0)] = (4 / 3) * f->m_mu[index(m_col - 1, 0)] * (f->m_v[index(m_col - 1, 1)] - f->m_v[index(m_col - 1, 0)]) / dy;
	m_qx[index(m_col - 1, 0)] = 0;
	m_qy[index(m_col - 1, 0)] = -f->m_k[index(m_col - 1, 0)] * (f->m_T[index(m_col - 1, 1)] - f->m_T[index(m_col - 1, 0)]) / dy;

}

void Solver::vectorSetU()
{
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_U1[index(i, j)] = f->m_rho[index(i, j)];
			m_U2[index(i, j)] = f->m_rho[index(i, j)] * f->m_u[index(i, j)];
			m_U3[index(i, j)] = f->m_rho[index(i, j)] * f->m_v[index(i, j)];
			m_U5[index(i, j)] = f->m_Et[index(i, j)];
		}
	}
}

void Solver::vectorSetE()
{
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_E1[index(i, j)] = m_U2[index(i, j)];
			m_E2[index(i, j)] = m_U2[index(i, j)] * (m_U2[index(i, j)] / m_U1[index(i, j)])
				+ ((m_U5[index(i, j)] / m_U1[index(i, j)]) - (pow((m_U2[index(i, j)] / m_U1[index(i, j)]), 2) + pow((m_U3[index(i, j)] / m_U1[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1[index(i, j)] * 287
				- m_txx[index(i, j)];
			m_E3[index(i, j)] = m_U2[index(i, j)] * (m_U3[index(i, j)] / m_U1[index(i, j)]) - m_txy[index(i, j)];
			m_E5[index(i, j)] = (m_U5[index(i, j)] + ((m_U5[index(i, j)] / m_U1[index(i, j)]) - (pow((m_U2[index(i, j)] / m_U1[index(i, j)]), 2) + pow((m_U3[index(i, j)] / m_U1[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1[index(i, j)] * 287)
				* (m_U2[index(i, j)] / m_U1[index(i, j)]) - (m_U2[index(i, j)] / m_U1[index(i, j)]) * m_txx[index(i, j)]
				- (m_U3[index(i, j)] / m_U1[index(i, j)]) * m_txy[index(i, j)] + m_qx[index(i, j)];
		}
	}
}

void Solver::vectorSetEfromPred()
{
	// Set PREDICTED U for BOUNDARIES
	// LEFT
	for (int j = 0; j < m_row; j++)
	{
		m_U1P[index(0, j)] = f->m_rho[index(0, j)];
		m_U2P[index(0, j)] = f->m_rho[index(0, j)] * f->m_u[index(0, j)];
		m_U3P[index(0, j)] = f->m_rho[index(0, j)] * f->m_v[index(0, j)];
		m_U5P[index(0, j)] = f->m_Et[index(0, j)];
	}
	// RIGHT
	for (int j = 0; j < m_row; j++)
	{
		m_U1P[index(m_col - 1, j)] = f->m_rho[index(m_col - 1, j)];
		m_U2P[index(m_col - 1, j)] = f->m_rho[index(m_col - 1, j)] * f->m_u[index(m_col - 1, j)];
		m_U3P[index(m_col - 1, j)] = f->m_rho[index(m_col - 1, j)] * f->m_v[index(m_col - 1, j)];
		m_U5P[index(m_col - 1, j)] = f->m_Et[index(m_col - 1, j)];
	}
	// DOWN
	for (int i = 1; i < m_col - 1; i++)
	{
		m_U1P[index(i, 0)] = f->m_rho[index(i, 0)];
		m_U2P[index(i, 0)] = f->m_rho[index(i, 0)] * f->m_u[index(i, 0)];
		m_U3P[index(i, 0)] = f->m_rho[index(i, 0)] * f->m_v[index(i, 0)];
		m_U5P[index(i, 0)] = f->m_Et[index(i, 0)];
	}
	// UP
	for (int i = 1; i < m_col - 1; i++)
	{
		m_U1P[index(i, m_row - 1)] = f->m_rho[index(i, m_row - 1)];
		m_U2P[index(i, m_row - 1)] = f->m_rho[index(i, m_row - 1)] * f->m_u[index(i, m_row - 1)];
		m_U3P[index(i, m_row - 1)] = f->m_rho[index(i, m_row - 1)] * f->m_v[index(i, m_row - 1)];
		m_U5P[index(i, m_row - 1)] = f->m_Et[index(i, m_row - 1)];
	}
	
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_E1[index(i, j)] = m_U2P[index(i, j)];
			m_E2[index(i, j)] = m_U2P[index(i, j)] * (m_U2P[index(i, j)] / m_U1P[index(i, j)])
				+ ((m_U5[index(i, j)] / m_U1P[index(i, j)]) - (pow((m_U2P[index(i, j)] / m_U1P[index(i, j)]), 2) + pow((m_U3P[index(i, j)] / m_U1P[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1P[index(i, j)] * 287
				- m_txx[index(i, j)];
			m_E3[index(i, j)] = m_U2P[index(i, j)] * (m_U3P[index(i, j)] / m_U1P[index(i, j)]) - m_txy[index(i, j)];
			m_E5[index(i, j)] = (m_U5P[index(i, j)] + ((m_U5P[index(i, j)] / m_U1P[index(i, j)]) - (pow((m_U2P[index(i, j)] / m_U1P[index(i, j)]), 2) + pow((m_U3P[index(i, j)] / m_U1P[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1P[index(i, j)] * 287)
				* (m_U2P[index(i, j)] / m_U1P[index(i, j)]) - (m_U2P[index(i, j)] / m_U1P[index(i, j)]) * m_txx[index(i, j)]
				- (m_U3P[index(i, j)] / m_U1P[index(i, j)]) * m_txy[index(i, j)] + m_qx[index(i, j)];
		}
	}
}

void Solver::vectorSetF()
{
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_F1[index(i, j)] = m_U3[index(i, j)];
			m_F2[index(i, j)] = m_U3[index(i, j)] * (m_U2[index(i, j)] / m_U1[index(i, j)]) - m_txy[index(i, j)];
			m_F3[index(i, j)] = m_U3[index(i, j)] * (m_U3[index(i, j)] / m_U1[index(i, j)])
				+ ((m_U5[index(i, j)] / m_U1[index(i, j)]) - (pow((m_U2[index(i, j)] / m_U1[index(i, j)]), 2) + pow((m_U3[index(i, j)] / m_U1[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1[index(i, j)] * 287
				- m_tyy[index(i, j)];
			m_F5[index(i, j)] = (m_U5[index(i, j)] + ((m_U5[index(i, j)] / m_U1[index(i, j)]) - (pow((m_U2[index(i, j)] / m_U1[index(i, j)]), 2) + pow((m_U3[index(i, j)] / m_U1[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1[index(i, j)] * 287)
				* (m_U3[index(i, j)] / m_U1[index(i, j)]) - (m_U2[index(i, j)] / m_U1[index(i, j)]) * m_txy[index(i, j)]
				- (m_U3[index(i, j)] / m_U1[index(i, j)]) * m_tyy[index(i, j)] + m_qy[index(i, j)];
		}
	}
}

void Solver::vectorSetFfromPred()
{
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_F1[index(i, j)] = m_U3P[index(i, j)];
			m_F2[index(i, j)] = m_U3P[index(i, j)] * (m_U2P[index(i, j)] / m_U1P[index(i, j)]) - m_txy[index(i, j)];
			m_F3[index(i, j)] = m_U3P[index(i, j)] * (m_U3P[index(i, j)] / m_U1P[index(i, j)])
				+ ((m_U5P[index(i, j)] / m_U1P[index(i, j)]) - (pow((m_U2P[index(i, j)] / m_U1P[index(i, j)]), 2) + pow((m_U3P[index(i, j)] / m_U1P[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1P[index(i, j)] * 287
				- m_tyy[index(i, j)];
			m_F5[index(i, j)] = (m_U5P[index(i, j)] + ((m_U5P[index(i, j)] / m_U1P[index(i, j)]) - (pow((m_U2P[index(i, j)] / m_U1P[index(i, j)]), 2) + pow((m_U3P[index(i, j)] / m_U1P[index(i, j)]), 2)) / 2) / (f->m_cv) * m_U1P[index(i, j)] * 287)
				* (m_U3P[index(i, j)] / m_U1P[index(i, j)]) - (m_U2P[index(i, j)] / m_U1P[index(i, j)]) * m_txy[index(i, j)]
				- (m_U3P[index(i, j)] / m_U1P[index(i, j)]) * m_tyy[index(i, j)] + m_qy[index(i, j)];
		}
	}
}

void Solver::plotViscous() const
{
	std::cout << "Tau_xy" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(5) << std::setw(12) << m_txy[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Tau_xx" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(5) << std::setw(12) << m_txx[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Tau_yy" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(5) << std::setw(12) << m_tyy[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "q_x" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(5) << std::setw(12) << m_qx[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "q_y" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(5) << std::setw(12) << m_qy[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;
}

void Solver::printU() const
{
	// Prints U vector
	
	std::cout << std::endl;
	std::cout << "Vector U1" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_U1[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector U2" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_U2[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector U3" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_U3[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector U5" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_U5[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

void Solver::printE() const
{
	// Prints E vector
	std::cout << std::endl;
	std::cout << "Vector E1" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_E1[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector E2" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_E2[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector E3" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_E3[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector E5" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_E5[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

void Solver::printF() const
{
	// Prints F vector
	std::cout << std::endl;
	std::cout << "Vector F1" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_F1[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector F2" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_F2[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector F3" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_F3[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
	std::cout << "Vector F5" << std::endl;
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_F5[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}

void Solver::printIteration() const
{
	std::cout << "*********************************************************" << std::endl;
	std::cout << "Iteration No: " << iteration << std::endl;
	std::cout << "*********************************************************" << std::endl;
	std::cout << "rho Residual: " << m_rho_res << std::endl;
	std::cout << "*********************************************************" << std::endl;
	std::cout << "Time step : " << dt << std::endl;
	std::cout << "Time elapsed : " << time << std::endl;
	std::cout << "*********************************************************" << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;	
}

void Solver::cellReynolds()
{
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_Rex[index(i, j)] = f->m_rho[index(i, j)] * f->m_u[index(i, j)] * dx / f->m_mu[index(i, j)];
			m_Rey[index(i, j)] = f->m_rho[index(i, j)] * f->m_v[index(i, j)] * dy / f->m_mu[index(i, j)];
		}
	}
}

void Solver::convergence()
{
	float temp = 0;
	m_rho_res = 0;
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			temp = fabs(f->m_rho[index(i, j)] - m_rho_temp[index(i, j)]) / f->m_rho[index(i, j)];
			if (temp > m_rho_res)
			{
				m_rho_res = temp;
			}
				
		}
	}
}


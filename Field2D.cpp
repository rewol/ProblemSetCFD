#include "Field2D.h"
Field2D::Field2D(int row, int col, float L, float mach, float ratio, float alt) : m_row(row), m_col(col), m_LHORI(L), M(mach), tr(ratio), A(alt)
{
	p_inf = A.getPres();
	rho_inf = A.getrho();
	u_inf = M * A.getSOS();
	T_inf = A.getTemp();
	mu_inf = A.getVis();
	e_inf = m_cv * T_inf;
	Re_inf = (rho_inf * u_inf * m_LHORI) / mu_inf;
	sigma = (5.0f * m_LHORI) / pow(Re_inf, 0.5f);
	m_LVERT = 5.0f * sigma;
	m_dx = m_LHORI / (m_col - 1);
	m_dy = m_LVERT / (m_row - 1);
	fieldInit();
	fieldSet();	
}

Field2D::~Field2D()
{
	delete[] m_u;
	delete[] m_v;
	delete[] m_V;
	delete[] m_rho;
	delete[] m_e;
	delete[] m_Et;
	delete[] m_p;
	delete[] m_mu;
	delete[] m_T;
	delete[] m_k;
	delete[] x;
	delete[] y;
}

void Field2D::fieldInit()
{
	m_u		= new float[m_row * m_col];
	m_v		= new float[m_row * m_col];
	m_V		= new float[m_row * m_col];
	m_T		= new float[m_row * m_col];
	m_p		= new float[m_row * m_col];
	m_e		= new float[m_row * m_col];
	m_Et	= new float[m_row * m_col];
	m_rho	= new float[m_row * m_col];
	m_mu	= new float[m_row * m_col];
	m_k		= new float[m_row * m_col];
	x		= new float[m_row * m_col];
	y		= new float[m_row * m_col];

	std::fill_n(m_u, m_row * m_col, 0);
	std::fill_n(m_v, m_row * m_col, 0);
	std::fill_n(m_rho, m_row * m_col, 0);
	std::fill_n(m_e, m_row * m_col, 0);
	std::fill_n(m_p, m_row * m_col, 0);
	std::fill_n(m_V, m_row * m_col, 0);
	std::fill_n(m_mu, m_row * m_col, 0);
	std::fill_n(m_k, m_row * m_col, 0);
	std::fill_n(m_Et, m_row * m_col, 0);
	std::fill_n(x, m_row * m_col, 0);
	std::fill_n(y, m_row * m_col, 0);
}

int Field2D::index(int row, int col) const
{
	return row + col * m_col;
}

void Field2D::fieldSet()
{
	// Set x and y coordinates !
	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			x[index(i, j)] =  i * m_dx;
		}
	}
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			y[index(i, j)] = j * m_dy;
		}
	}	
	
	// LEFT
	for (int j = 1; j < m_row - 1; j++)
	{
		m_u[index(0, j)] = u_inf;
		m_v[index(0, j)] = 0;
		m_p[index(0, j)] = p_inf;
		m_T[index(0, j)] = T_inf;
	}

	// LEADING EDGE
	m_u[index(0, 0)] = 0;
	m_v[index(0, 0)] = 0;
	m_p[index(0, 0)] = p_inf;
	m_T[index(0, 0)] = T_inf;
	
	// INTERNAL FIELD
	for (int i = 1; i < m_col - 1; i++)
	{
		for (int j = 1; j < m_row - 1; j++)
		{
			m_u[index(i, j)] = u_inf;
			m_v[index(i, j)] = 0;
			m_p[index(i, j)] = p_inf;
			m_T[index(i, j)] = T_inf;
		}
	}

	// UPPER
	for (int i = 0; i < m_col; i++)
	{
		m_u[index(i, m_row - 1)] = u_inf;
		m_v[index(i, m_row - 1)] = 0;
		m_p[index(i, m_row - 1)] = p_inf;
		m_T[index(i, m_row - 1)] = T_inf;
	}

	// RIGHT
	for (int j = 1; j < m_row - 1; j++)
	{
		m_u[index(m_col - 1, j)] = 
			2 * m_u[index(m_col - 2, j)] - m_u[index(m_col - 3, j)];
		m_v[index(m_col - 1, j)] =
			2 * m_v[index(m_col - 2, j)] - m_v[index(m_col - 3, j)];
		m_p[index(m_col - 1, j)] =
			2 * m_p[index(m_col - 2, j)] - m_p[index(m_col - 3, j)];
		m_T[index(m_col - 1, j)] =
			2 * m_T[index(m_col - 2, j)] - m_T[index(m_col - 3, j)];
	}

	// WALL
	for (int i = 1; i < m_col; i++)
	{
		m_u[index(i, 0)] = 0;
		m_v[index(i, 0)] = 0;
		m_T[index(i, 0)] = tr * T_inf;
		m_p[index(i, 0)] = 2 * m_p[index(i, 1)] - m_p[index(i, 2)];
	}

	//
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_rho[index(i, j)] = m_p[index(i, j)] / (287 * m_T[index(i, j)]);
			m_e[index(i, j)] = m_cv * m_T[index(i, j)];
			m_Et[index(i, j)] = m_rho[index(i, j)] *
				(m_e[index(i, j)] + pow(m_u[index(i, j)], 2) / 2);
			m_V[index(i, j)] = m_u[index(i, j)];
		}
	}

	thermophysical();
}

void Field2D::thermophysical()
{
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			m_mu[index(i, j)] = 1.7894e-5f * ((288.16f + 110.0f) / (m_T[index(i, j)] + 110.0f)) * pow(m_T[index(i, j)] / 288.16f, 1.5f);
			m_k[index(i, j)] = m_mu[index(i, j)] * m_cp / Pr;
		}
	}
}

void Field2D::plotPrimitive() const
{
	
	std::cout << "X-Velocity [m / s]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_u[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Y-Velocity [m / s]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_v[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Pressure [Pa]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_p[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Temperature [K]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_T[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Density [kg/m3]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_rho[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Energy [J]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_e[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Total Energy [J]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_Et[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
}

void Field2D::plotThermophysical() const
{
	std::cout << "Dynamic Viscosity [kg / s]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_mu[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;

	std::cout << "Thermal cond. [W / m.K]" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(size) << std::setw(12) << m_k[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}

	std::cout << std::endl;
}

void Field2D::printPrimitive() const
{
	int font = 20;
	// Print X Velocity
	std::fstream file;
	std::string filename = "x_vel.txt";
	file.open(filename, std::fstream::out);

	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			file << std::left << std::setw(font) << x[index(i, j)] << std::setw(font) << y[index(i, j)] << std::setw(font) << m_u[index(i, j)] << std::endl;
			if (j == m_row - 1)
				file << std::endl;
		}
	}
	file.close();

	// Print Y Velocity
	std::fstream file1;
	std::string filename1 = "y_vel.txt";
	file1.open(filename1, std::fstream::out);

	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			file1 << std::left << std::setw(font) << x[index(i, j)] << std::setw(font) << y[index(i, j)] << std::setw(font) << m_v[index(i, j)] << std::endl;
			if (j == m_row - 1)
				file1 << std::endl;
		}
	}
	file1.close();

	// Print Density
	std::fstream file2;
	std::string filename2 = "rho.txt";
	file2.open(filename2, std::fstream::out);

	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			file2 << std::left << std::setw(font) << x[index(i, j)] << std::setw(font) << y[index(i, j)] << std::setw(font) << m_rho[index(i, j)] << std::endl;
			if (j == m_row - 1)
				file2 << std::endl;
		}
	}
	file2.close();

	// Print Pressure
	std::fstream file3;
	std::string filename3 = "p.txt";
	file3.open(filename3, std::fstream::out);

	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			file3 << std::left << std::setw(font) << x[index(i, j)] << std::setw(font) << y[index(i, j)] << std::setw(font) << m_p[index(i, j)] << std::endl;
			if (j == m_row - 1)
				file3 << std::endl;
		}
	}
	file3.close();

	// Print Temperature
	// Print Density
	std::fstream file4;
	std::string filename4 = "T.txt";
	file4.open(filename4, std::fstream::out);

	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			file4 << std::left << std::setw(font) << x[index(i, j)] << std::setw(font) << y[index(i, j)] << std::setw(font) << m_T[index(i, j)] << std::endl;
			if (j == m_row - 1)
				file4 << std::endl;
		}
	}
	file4.close();
}




#include "Grid.h"

Grid::Grid(int row, int col) : m_row(row), m_col(col), deta(1 / ((double)m_row - 1)), dksi(1 / ((double)m_col - 1))
{
	gridInit();
}

Grid::~Grid()
{
	delete[] x;
	delete[] y;
	delete[] dxdeta;
	delete[] dydeta;
	delete[] dxdksi;
	delete[] dydksi;
	delete[] J;
	std::cout << "Grid dtor is called !" << std::endl;
}

int Grid::index(int row, int col)
{
	return row + col * m_col;
}

void Grid::gridInit()
{
	x = new double[m_row * m_col];
	y = new double[m_row * m_col];
	dxdeta = new double[m_row * m_col];
	dydeta = new double[m_row * m_col];
	dxdksi = new double[m_row * m_col];
	dydksi = new double[m_row * m_col];
	J = new double[m_row * m_col];

	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			x[index(i, j)] = 0;
			y[index(i, j)] = 0;
			dxdeta[index(i, j)] = 0;
			dydeta[index(i, j)] = 0;
			dxdksi[index(i, j)] = 0;
			dydksi[index(i, j)] = 0;
			J[index(i, j)] = 0;
		}
	}

}

void Grid::setIMAX(int row)
{
	m_row = row;
	gridInit();
}

void Grid::setJMAX(int col)
{
	m_col = col;
	gridInit();
}

void Grid::print()
{
	std::cout << std::left << "x-coordinates" << std::endl;
	std::cout << "-------------" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(3) << std::setw(8) << x[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
	std::cout << std::endl;
	std::cout << std::left << "y-coordinates" << std::endl;
	std::cout << "-------------" << std::endl;
	for (int j = 0; j < m_row; j++)
		for (int i = 0; i < m_col; i++)
		{
			std::cout << std::setprecision(3) << std::setw(8) << y[index(i, j)];
			if (i == m_col - 1)
				std::cout << std::endl;
		}
}

void Grid::transfiniteInterpolation()
{
	int M = m_col - 1;
	int N = m_row - 1;
	for (int i = 1; i < M; i++)
	{
		for (int j = 1; j < N; j++)
		{
			x[index(i, j)] = (i / static_cast<double>(M)) * x[index(M, j)] + ((static_cast<double>(M) - i) / static_cast<double>(M)) * x[index(0, j)] + (j / static_cast<double>(N)) * x[index(i, N)]
				+ ((static_cast<double>(N) - j) / static_cast<double>(N)) * x[index(i, 0)] - (i / static_cast<double>(M)) * (j / static_cast<double>(N)) * x[index(M, N)]
				- (i / static_cast<double>(M)) * ((static_cast<double>(N) - j) / static_cast<double>(N)) * x[index(M, 0)] - ((static_cast<double>(M) - i) / static_cast<double>(M)) * (j / static_cast<double>(N)) * x[index(0, N)]
				- ((static_cast<double>(M) - i) / static_cast<double>(M)) * ((static_cast<double>(N) - j) / static_cast<double>(N)) * x[index(0, 0)];

			y[index(i, j)] = (i / static_cast<double>(M)) * y[index(M, j)] + ((static_cast<double>(M) - i) / static_cast<double>(M)) * y[index(0, j)] + (j / static_cast<double>(N)) * y[index(i, N)]
				+ ((static_cast<double>(N) - j) / static_cast<double>(N)) * y[index(i, 0)] - (i / static_cast<double>(M)) * (j / static_cast<double>(N)) * y[index(M, N)]
				- (i / static_cast<double>(M)) * ((static_cast<double>(N) - j) / static_cast<double>(N)) * y[index(M, 0)] - ((static_cast<double>(M) - i) / static_cast<double>(M)) * (j / static_cast<double>(N)) * y[index(0, N)]
				- ((static_cast<double>(M) - i) / static_cast<double>(M)) * ((static_cast<double>(N) - j) / static_cast<double>(N)) * y[index(0, 0)];
		}
	}
}

void Grid::record(std::string filename)
{
	filename = "C:/Users/Refik Alper/Desktop/plot/" + filename;
	std::fstream file;
	file.open(filename, std::fstream::out);
	for (int i = 0; i < m_col; i++)
		for (int j = 0; j < m_row; j++)
		{
			file << std::left << std::setw(20) << x[index(i, j)] << std::setw(20) << y[index(i, j)] << std::setw(20) << "0" << std::endl;
			if (j == m_row - 1)
				file << std::endl;
		}
	file.close();
}

void Grid::metrics()
{
	for (int i = 0; i < m_col; i++)
	{
		for (int j = 0; j < m_row; j++)
		{
			if (j == 0)
			{
				dxdeta[index(i, j)] = (-3 * x[index(i, j)] + 4 * x[index(i, j + 1)] - x[index(i, j + 2)]) / (2 * deta);
				dydeta[index(i, j)] = (-3 * y[index(i, j)] + 4 * y[index(i, j + 1)] - y[index(i, j + 2)]) / (2 * deta);
			}

			else if (j == m_row - 1)
			{
				dxdeta[index(i, j)] = (-3 * x[index(i, j)] + 4 * x[index(i, j - 1)] - x[index(i, j - 2)]) / (2 * deta);
				dydeta[index(i, j)] = (-3 * y[index(i, j)] + 4 * y[index(i, j - 1)] - y[index(i, j - 2)]) / (2 * deta);
			}

			else
			{
				dxdeta[index(i, j)] = (x[index(i, j + 1)] - x[index(i, j - 1)]) / (2 * deta);
				dydeta[index(i, j)] = (y[index(i, j + 1)] - y[index(i, j - 1)]) / (2 * deta);
			}

		}
	}

	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			if (i == 0)
			{
				dxdksi[index(i, j)] = (-3 * x[index(i, j)] + 4 * x[index(i + 1, j)] - x[index(i + 2, j)]) / (2 * dksi);
				dydksi[index(i, j)] = (-3 * y[index(i, j)] + 4 * y[index(i + 1, j)] - y[index(i + 2, j)]) / (2 * dksi);
			}
			else if (i == m_col - 1)
			{
				dxdksi[index(i, j)] = (-3 * x[index(i, j)] + 4 * x[index(i - 1, j)] - x[index(i - 2, j)]) / (2 * dksi);
				dydksi[index(i, j)] = (-3 * y[index(i, j)] + 4 * y[index(i - 1, j)] - y[index(i - 2, j)]) / (2 * dksi);
			}
			else
			{
				dxdksi[index(i, j)] = (x[index(i + 1, j)] - x[index(i - 1, j)]) / (2 * dksi);
				dydksi[index(i, j)] = (y[index(i + 1, j)] - y[index(i - 1, j)]) / (2 * dksi);
			}
		}
	}

	for (int j = 0; j < m_row; j++)
	{
		for (int i = 0; i < m_col; i++)
		{
			J[index(i, j)] = (dxdksi[index(i, j)] * dydeta[index(i, j)] - dydksi[index(i, j)] * dxdeta[index(i, j)]);
		}
	}
}

int Grid::getRow()
{
	return m_row;
}

int Grid::getCol()
{
	return m_col;
}

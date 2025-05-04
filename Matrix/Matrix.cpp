#include <iostream>
#include "Matrix.h"
#include <cmath>
#include <vector>
#include <chrono>
#include <thread>

using namespace std;

matrix::matrix(int n, int m)
{
	this->n = n;
	this->m = m;
	try
	{
		cells = new double* [n];
		for (int i = 0; i < n; i++)
		{
			cells[i] = new double[m];
			for (int j = 0; j < m; j++)
			{
				cells[i][j] = 0.0;
			}
		}
	}
	catch (const bad_alloc& e)
	{
		cerr << "Allocation failed: " << e.what() << '\n';
	}
}
matrix::matrix()
{
	this->n = 0;
	this->m = 0;
	try
	{
		cells = new double* [n];
		for (int i = 0; i < n; i++)
		{
			cells[i] = new double[m];
			for (int j = 0; j < m; j++)
			{
				cells[i][j] = 0.0;
			}
		}
	}
	catch (const bad_alloc& e)
	{
		cerr << "Allocation failed: " << e.what() << '\n';
	}
}
matrix::matrix(const matrix& other)
{
	n = other.n;
	m = other.m;
	cells = new double* [n];
	for (int i = 0; i < n; i++)
	{
		cells[i] = new double[m];
		for (int j = 0; j < m; j++)
		{
			cells[i][j] = other.cells[i][j];
		}
	}
}
matrix::~matrix()
{
	for (int i = 0; i < n; i++)
	{
		delete[] cells[i];
	}
	delete[] cells;
}

void matrix::setDiag(int diag, double value)
{
	int i = 0, j = 0;
	if (diag < 0)
		i -= diag;
	else
		j += diag;
	while (i < n && j < m)
	{
		cells[i++][j++] = value;
	}
}
void matrix::setValue(double value)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cells[i][j] = value;
		}
	}
}
void matrix::setValue(double (*f)(int))
{
	int element_number = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cells[i][j] = f(element_number++);
		}
	}
}

matrix matrix::residual(const matrix& A, const matrix& b, const matrix& x)
{
	return A * x - b;
}
double matrix::norm(const matrix& r)
{
	double sum = 0.0;
	for (int i = 0; i < r.n; i++)
	{
		for (int j = 0; j < r.m; j++)
		{
			sum += r.cells[i][j] * r.cells[i][j];
		}
	}
	return sqrt(sum);
}

solution matrix::jacobi(const matrix& A, const matrix& b, double epsilon)
{
	matrix x(A.m, 1);
	solution sol;
	double n = norm(residual(A, b, x));
	int iterations = 0;

	x.setValue(1.0);

	sol.norms.push_back(n);

	auto start = chrono::high_resolution_clock::now();

	while(n > epsilon && iterations <= MAX_ITERATIONS)
	{
		matrix x_new(A.m, 1);

		for (int i = 0; i < A.m; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < A.m; j++)
			{
				if (i != j)
				{
					sum += A.cells[i][j] * x.cells[j][0];
				}
			}
			x_new.cells[i][0] = (b.cells[i][0] - sum) / A.cells[i][i];
		}
		
		x = x_new;
		n = norm(residual(A, b, x));
		sol.norms.push_back(n);
		iterations++;
	}

	auto end = chrono::high_resolution_clock::now();

	sol.time = chrono::duration<double>(end - start).count();
	sol.x = x;
	sol.iterations = iterations;
	return sol;
}
solution matrix::gaussSeidel(const matrix& A, const matrix& b, double epsilon)
{
	matrix x(A.m, 1);
	solution sol;
	double n = norm(residual(A, b, x));
	int iterations = 0;

	x.setValue(1.0);

	sol.norms.push_back(n);

	auto start = chrono::high_resolution_clock::now();

	while (n > epsilon && iterations <= MAX_ITERATIONS)
	{
		for (int i = 0; i < A.m; i++)
		{
			double sum = 0.0;
			for (int j = 0; j < A.m; j++)
			{
				if (i != j)
				{
					sum += A.cells[i][j] * x.cells[j][0];
				}
			}
			x.cells[i][0] = (b.cells[i][0] - sum) / A.cells[i][i];
		}
		n = norm(residual(A, b, x));
		sol.norms.push_back(n);
		iterations++;
	}

	auto end = chrono::high_resolution_clock::now();

	sol.time = chrono::duration<double>(end - start).count();
	sol.x = x;
	sol.iterations = iterations;
	return sol;
}
solution matrix::lu(const matrix& A, const matrix& b)
{
	matrix x(A.m, 1);
	solution sol;

	auto start = chrono::high_resolution_clock::now();

	matrix U = A;
	matrix L(A.n, A.m);
	L.setDiag(0, 1.0);
	for (int i = 2; i <= A.m; i++)
	{
		for (int j = 1; j < i; j++)
		{
			L.cells[i-1][j-1] = U.cells[i-1][j-1] / U.cells[j-1][j-1];
			for (int k = j; k < A.m; k++)
			{
				U.cells[i-1][k] -= L.cells[i-1][j-1] * U.cells[j-1][k];
			}
		}
	}

	//Ly = b
	matrix y(A.m, 1);
	for (int i = 0; i < A.m; i++)
	{
		double sum = 0.0;
		for (int j = 0; j < i; j++)
		{
			sum += L.cells[i][j] * y.cells[j][0];
		}
		y.cells[i][0] = (b.cells[i][0] - sum) / L.cells[i][i];
	}
	//Ux = y
	for (int i = A.m - 1; i >= 0; i--)
	{
		double sum = 0.0;
		for (int j = A.m - 1; j > i; j--)
		{
			sum += U.cells[i][j] * x.cells[j][0];
		}
		x.cells[i][0] = (y.cells[i][0] - sum) / U.cells[i][i];
	}

	auto end = chrono::high_resolution_clock::now();

	sol.time = chrono::duration<double>(end - start).count();
	sol.norms.push_back(norm(residual(A, b, x)));
	sol.x = x;
	return sol;
}

matrix matrix::operator+(const matrix& other) const
{
	if (n == other.n && m == other.m)
	{
		matrix result(n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				result.cells[i][j] = cells[i][j] + other.cells[i][j];
			}
		}
		return result;
	}
	else
	{
		cerr << "Matrix dimensions do not match\n";
		return matrix(0, 0);
	}
}
matrix matrix::operator-(const matrix& other) const
{
	if (n == other.n && m == other.m)
	{
		matrix result(n, m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				result.cells[i][j] = cells[i][j] - other.cells[i][j];
			}
		}
		return result;
	}
	else
	{
		cerr << "Matrix dimensions do not match\n";
		return matrix(0, 0);
	}
}
matrix matrix::operator-() const
{
	matrix result(n, m);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			result.cells[i][j] = -cells[i][j];
		}
	}
	return result;
}
matrix matrix::operator*(const matrix& other) const
{
	if (m == other.n)
	{
		matrix result(n, other.m);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < other.m; j++)
			{
				for (int k = 0; k < m; k++)
				{
					result.cells[i][j] += cells[i][k] * other.cells[k][j];
				}
			}
		}
		return result;
	}
	else
	{
		cerr << "Matrix dimensions do not match\n";
		return matrix(0, 0);
	}
}
matrix& matrix::operator=(const matrix& other)
{
	if (n != other.n || m != other.m)
	{
		for (int i = 0; i < n; i++)
		{
			delete[] cells[i];
		}
		delete[] cells;
		n = other.n;
		m = other.m;
		cells = new double* [other.n];
		for (int i = 0; i < n; i++)
		{
			cells[i] = new double[other.m];
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cells[i][j] = other.cells[i][j];
		}
	}
	return *this;
}
bool matrix::operator!=(const matrix& other) const
{
	if (n != other.n || m != other.m)
	{
		return true;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (cells[i][j] != other.cells[i][j])
			{
				return true;
			}
		}
	}
	return false;
}
ostream& operator<<(ostream& os, const matrix& mat)
{
	for (int i = 0; i < mat.n; i++)
	{
		for (int j = 0; j < mat.m; j++)
		{
			os << mat.cells[i][j] << "\t";
		}
		os << "\n";
	}
	return os;
}

solution::solution()
{
	x = matrix(0, 0);
	time = 0.0;
	iterations = 0;
	norms = vector<double>();
}
solution::~solution()
{

}
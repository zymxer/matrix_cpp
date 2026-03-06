#include <iostream>
#include "Matrix.h"
#include <cmath>
#include <vector>
#include <chrono>
#include <thread>
#include <algorithm>
#include <execution> // для параллельного execution policy

using namespace std;

matrix::matrix(size_t rows, size_t cols, double value)
	: n(rows), m(cols), cells(rows* cols, value) {
}

void matrix::setDiag(size_t diag, double value)
{
	size_t i = 0, j = 0;
	if (diag < 0)
		i -= diag;
	else
		j += diag;
	while (i < n && j < m)
	{
		(*this)(i++, j++) = value;
	}
}

void matrix::setValue(double value)
{
	for (auto& cell : cells)
	{
		cell = value;
	}
}

void matrix::setValue(double (*f)(int))
{
	int element_number = 0;
	for (auto& cell : cells)
	{
		cell = f(element_number++);
	}
}

matrix matrix::residual(const matrix& A, const matrix& b, const matrix& x)
{
	return A * x - b;
}

double matrix::norm(const matrix& r)
{
	double sum = 0.0;
	for (auto& cell : r.cells)
	{
		sum += cell * cell;
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
	matrix x_new(A.m, 1);
	while (n > epsilon && iterations <= MAX_ITERATIONS)
	{
		for (size_t i = 0; i < A.m; i++)
		{
			double sum = 0.0;
			for (size_t j = 0; j < A.m; j++)
			{
				if (i != j)
				{
					sum += A(i, j) * x(j, 0);
				}
			}
			x_new(i, 0) = (b(i, 0) - sum) / A(i, i);
		}

		x = x_new;
		n = norm(residual(A, b, x));
		sol.norms.push_back(n);
		iterations++;
	}

	auto end = chrono::high_resolution_clock::now();

	sol.time = chrono::duration<double>(end - start).count();
	sol.x = std::move(x);
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
		for (size_t i = 0; i < A.m; i++)
		{
			double sum = 0.0;
			for (size_t j = 0; j < A.m; j++)
			{
				if (i != j)
				{
					sum += A(i, j) * x(j, 0);
				}
			}
			x(i, 0) = (b(i, 0) - sum) / A(i, i);
		}
		n = norm(residual(A, b, x));
		sol.norms.push_back(n);
		iterations++;
	}

	auto end = chrono::high_resolution_clock::now();

	sol.time = chrono::duration<double>(end - start).count();
	sol.x = std::move(x);
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
			L(i - 1, j - 1) = U(i - 1, j - 1) / U(j - 1, j - 1);
			for (int k = j; k < A.m; k++)
			{
				U(i - 1, k) -= L(i - 1, j - 1) * U(j - 1, k);
			}
		}
	}

	//Ly = b
	matrix y(A.m, 1);
	for (size_t i = 0; i < A.m; i++)
	{
		double sum = 0.0;
		for (size_t j = 0; j < i; j++)
		{
			sum += L(i, j) * y(j, 0);
		}
		y(i, 0) = (b(i, 0) - sum) / L(i, i);
	}
	//Ux = y
	for (int i = A.m - 1; i >= 0; i--)
	{
		double sum = 0.0;
		for (size_t j = A.m - 1; j > i; j--)
		{
			sum += U(i, j) * x(j, 0);
		}
		x(i, 0) = (y(i, 0) - sum) / U(i, i);
	}

	auto end = chrono::high_resolution_clock::now();

	sol.time = chrono::duration<double>(end - start).count();
	sol.norms.push_back(norm(residual(A, b, x)));
	sol.x = std::move(x);
	return sol;
}

matrix matrix::operator+(const matrix& other) const
{
	if (n != other.n || m != other.m)
	{
		cerr << "Matrix dimensions do not match\n";
		return matrix(0, 0);
	}
	matrix result(n, m);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < m; j++)
		{
			result(i, j) = (*this)(i, j) + other(i, j);
		}
	}
	return result;
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
				result(i, j) = (*this)(i, j) - other(i, j);
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
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < m; j++)
		{
			result(i, j) = -(*this)(i, j);
		}
	}
	return result;
}
matrix matrix::operator*(const matrix& other) const
{
	if (m == other.n)
	{
		matrix result(n, other.m);
		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < other.m; j++)
			{
				for (size_t k = 0; k < m; k++)
				{
					result(i, j) += (*this)(i, k) * other(k, j);
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

bool matrix::operator!=(const matrix& other) const
{
	if (n != other.n || m != other.m)
	{
		return true;
	}
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < m; j++)
		{
			if ((*this)(i, j) != other(i, j))
			{
				return true;
			}
		}
	}
	return false;
}
ostream& operator<<(ostream& os, const matrix& mat)
{
	for (size_t i = 0; i < mat.n; i++)
	{
		for (size_t j = 0; j < mat.m; j++)
		{
			os << mat(i, j) << "\t";
		}
		os << "\n";
	}
	return os;
}

solution::solution()
	: x(std::move(matrix(0, 0))), time(0.0), iterations(0), norms(vector<double>())
{

}

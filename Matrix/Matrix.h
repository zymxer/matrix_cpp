#pragma once
#include <cmath>
#include <vector>

#define MAX_ITERATIONS 500


class solution;

class matrix
{
public:
	matrix(size_t rows = 0, size_t cols = 0, double value = 0.0);

	void setDiag(size_t diag, double value);

	size_t rows() inline const { return n; }
	size_t cols() inline const { return m; }

	double& operator()(size_t i, size_t j) {
		return cells[i * m + j];
	}

	double operator()(size_t i, size_t j) const {
		return cells[i * m + j];
	}

	matrix(matrix&& other) noexcept
		: n(other.n), m(other.m), cells(std::move(other.cells))
	{
		other.n = 0;
		other.m = 0;
	}

	matrix& operator=(matrix&& other)
	{
		if (this != &other) {
			n = other.n; m = other.m;
			cells = std::move(other.cells);
			other.n = other.m = 0;
		}
		return *this;
	}

	matrix(const matrix& other)
		: n(other.n), m(other.m), cells(other.cells)
	{
	}

	matrix& operator=(const matrix& other) {
		if (this != &other) {
			n = other.n;
			m = other.m;
			cells = other.cells;
		}
		return *this;
	}

	void setValue(double value);
	void setValue(double (*f)(int));

	static matrix residual(const matrix& A, const matrix& b, const matrix& x);
	static double norm(const matrix& r);

	static solution jacobi(const matrix& A, const matrix& b, double epsilon);
	static solution gaussSeidel(const matrix& A, const matrix& b, double epsilon);
	static solution lu(const matrix& A, const matrix& b);


	matrix operator+(const matrix& other) const;
	matrix operator-(const matrix& other) const;
	matrix operator-() const;
	matrix operator*(const matrix& other) const;
	bool operator!=(const matrix& other) const;
	friend std::ostream& operator<<(std::ostream& os, const matrix& m);


private:
	size_t n;
	size_t m;
	std::vector<double> cells;

};

class solution
{
public:
	solution();
	matrix x;
	double time;
	int iterations;
	std::vector<double> norms;
};
#pragma once
#include <cmath>
#include <vector>

#define MAX_ITERATIONS 500

using namespace std;

class solution;

class matrix
{
public:
	matrix(int n, int m);
	matrix();
	matrix(const matrix& other);
	~matrix();


	void setDiag(int diag, double value);
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
	matrix& operator=(const matrix& other);
	bool operator!=(const matrix& other) const;
	friend std::ostream& operator<<(std::ostream& os, const matrix& m);


private:
	int n;
	int m;
	double** cells;

};

class solution
{
public:
	solution();
	~solution();
	matrix x;
	double time;
	int iterations;
	vector<double> norms;
};
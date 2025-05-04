#include <iostream>
#include <cmath>
#include "Matrix.h"
#include <fstream>

#define E 1e-9

using namespace std;

double matrix_value(int n)
{
	return sin(n * (6 + 1));
}

void normsToCSV(const string& filename, const vector<double>& jacobiNorms,
	const vector<double>& gaussSeidelNorms) {
	
	ofstream file(filename);

	
	file << "Jacobi Norms,Gauss-Seidel Norms\n";
 
	int maxSize = max(jacobiNorms.size(), gaussSeidelNorms.size());

	for (size_t i = 0; i < maxSize; ++i) {

		if (i < jacobiNorms.size()) {
			file << jacobiNorms[i] << ",";
		}
		else {
			file << ",";
		}

		if (i < gaussSeidelNorms.size()) {
			file << gaussSeidelNorms[i] << "\n";
		}
		else {
			file << "\n";
		}
	}

	file.close();
}

void timesToCSV(const string& filename, const vector<double>& jacobiTimes,
	const vector<double>& gaussSeidelTimes, const vector<double>& LUTimes) {

	ofstream file(filename);


	file << "Jacobi Times,Gauss-Seidel Times,LU Times\n";

	int size = jacobiTimes.size();

	for (size_t i = 0; i < size; ++i) {

		file << jacobiTimes[i] << "," << gaussSeidelTimes[i] << "," << LUTimes[i] << "\n";
	}

	file.close();
}

int main() {

	int N = 923;

	int a1 = 12, a2 = -1, a3 = -1;
	matrix A(N, N);

	A.setDiag(0, a1);
	A.setDiag(1, a2);
	A.setDiag(2, a3);
	A.setDiag(-1, a2);
	A.setDiag(-2, a3);

	matrix b(N, 1);
	b.setValue(matrix_value);

	solution s = matrix::jacobi(A, b, E);
	solution s2 = matrix::gaussSeidel(A, b, E);
	cout << "Jacobi:" << endl;
	cout << "iterations: " << s.iterations << endl;
	cout << "time: " << s.time << endl;
	cout << "Gauss-Seidel: " << endl;
	cout << "iterations: " << s2.iterations << endl;
	cout << "time: " << s2.time << endl;
	//normsToCSV("norms.csv", s.norms, s2.norms);

	a1 = 3.0;
	a2 = -1;
	a3 = -1;
	A = matrix(N, N);

	A.setDiag(0, a1);
	A.setDiag(1, a2);
	A.setDiag(2, a3);
	A.setDiag(-1, a2);
	A.setDiag(-2, a3);

	solution s3 = matrix::jacobi(A, b, E);
	solution s4 = matrix::gaussSeidel(A, b, E);
	//normsToCSV("norms2.csv", s3.norms, s4.norms);


	solution s5 = matrix::lu(A, b);
	//cout << s5.norms[0] << endl;


	vector<int> sizes = { 100, 500, 1000, 2000, 3000, 4000, 5000};
	vector<double> timesJacobi;
	vector<double> timesGaussSeidel;
	vector<double> timesLU;
	a1 = 12;
	a2 = -1;
	a3 = -1;

	cout << "\nTime to solve:" << endl;

	for (int i = 0; i < sizes.size(); i++)
	{
		matrix B(sizes[i], sizes[i]);
		B.setDiag(0, a1);
		B.setDiag(1, a2);
		B.setDiag(2, a3);
		B.setDiag(-1, a2);
		B.setDiag(-2, a3);

		matrix c(sizes[i], 1);
		c.setValue(matrix_value);

		solution j = matrix::jacobi(B, c, E);
		solution gs = matrix::gaussSeidel(B, c, E);
		solution lu = matrix::lu(B, c);
		timesJacobi.push_back(j.time);
		timesGaussSeidel.push_back(gs.time);
		timesLU.push_back(lu.time);

		cout << "---------------------------------" << endl;
		cout << "Matrix of size " << sizes[i] << ":" << endl;
		cout << "Jacobi: " << j.time << endl;
		cout << "Gauss-Seidel: " << gs.time << endl;
		cout << "LU: " << lu.time << endl;
	}
	//timesToCSV("times.csv", timesJacobi, timesGaussSeidel, timesLU);	

	return 0;
}
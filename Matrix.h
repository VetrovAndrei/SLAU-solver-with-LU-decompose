#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;


typedef float real;
typedef double dubl;

class Matrix
{
private:
	int n;
	vector< vector<real> > matrix;
public:
	Matrix(int x);
	Matrix(void);
	~Matrix(void);
	void setMatrix(vector< vector<real> > A, int x);
	void ToProf(MatrixProf &A);
	void Gauss(vector<real> *X, vector<real> B)
};

class MatrixProf
{
private:
	int n;
	int col;
	vector<real> di;
	vector<int> ia;
	vector<real> al;
	vector<real> au;
	vector<real> F;
	
public:
	MatrixProf(void);
	MatrixProf(int x, int c);
	~MatrixProf(void);
	void load (ifstream &size, ifstream &matrix, ifstream &vect);
	void save (ofstream &output);
	void LUDec();
	vector<dubl> Direct();
	void Reverse(vector<real> F, vector<real> *X);
	void ToTight(Matrix *A);
};



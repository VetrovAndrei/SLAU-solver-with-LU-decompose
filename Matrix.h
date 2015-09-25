#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <locale>
#include <iostream>
using namespace std;

class MatrixProf;

typedef float real;
typedef float dubl;

class Matrix
{
private:
	int n;
	int col;
	vector< vector<real> > matrix;
	vector<real> F;
public:
	Matrix(int x, ifstream &vect);
	Matrix(int x);
	Matrix(void);
	~Matrix(void);
	void getCol();
	void setMatrix(vector< vector<real> > A, int x, vector<real> B);
	void ToProf(MatrixProf *A);
	vector<dubl> Gauss();
	vector<real> GoodGauss();
	void Gilbert();
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
	void load (ifstream &matrix, ifstream &vect);
	void save (ofstream &output);
	void LUDec();
	vector<dubl> Direct();
	vector<dubl> Reverse(vector<dubl> F);
	void ToTight(Matrix *A);
	vector<dubl> SLAU();
	void setProf(int bn, int bcol, vector<real> bdi, vector<int> bia, vector<real> bal, vector<real> bau, vector<real> bF);
};



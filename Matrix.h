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


typedef float real;
typedef double dubl;

class Matrix
{
private:
	int n;
	int col;
	vector< vector<real> > matrix;
	vector<real> F;
public:
	Matrix(int x);
	Matrix(void);
	~Matrix(void);
	void getCol();
	void setMatrix(vector< vector<real> > A, int x, vector<real> B);
	void ToProf(MatrixProf *A);
	vector<real> Gauss(vector<real> B);
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
};



#include "Matrix.h"


Matrix::Matrix(int x)
{
	this->n = x;
	for (int i = 0; i < n; i++)
	{
		vector<real> buf;
		for (int j = 0; j < n; j++)
			buf.push_back(0);
		matrix.push_back(buf);
	}
}

void Matrix::setMatrix(vector< vector<real> > A, int x)
{
	n = x;
	for (int i = 0; i < n; i++)
	{
		vector<real> buf;
		real b;
		for (int j = 0; j < n; j++)
		{
			b = A[i][j];
			buf.push_back(b);
		}
		matrix.push_back(buf);
	}
}

Matrix::Matrix(void)
{
}


Matrix::~Matrix(void)
{
}

void Matrix::Gauss(vector<real> *X, vector<real> B)
{
	vector<dubl> bufer;
	for(int i = 0; i < n; i++)
		bufer.push_back(0);
	for (int i = 1; i<n; i++)
		for (int j=i; j<n; j++)
		{
			real m = matrix[j][i-1]/matrix[i-1][i-1];
			for (int k = 0; k < n; k++)
				matrix[j][k]=matrix[j][k]-m*matrix[i-1][k];
			B[j] = B[j] - m * B[i-1]; 
		}
	for (int i=n-1; i>=0; i--)
		{
			for (int j = i+1; j < n; j++)
			{
				bufer[i] -= matrix[i][j]*bufer[j];
			}
			bufer[i] = bufer[i]/matrix[i][i];
		}
	for(int i = 0; i < n; i++)
		X->push_back(bufer[i]);
}

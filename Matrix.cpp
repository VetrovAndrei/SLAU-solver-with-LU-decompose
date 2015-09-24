#include "Matrix.h"


Matrix::Matrix(int x)
{
	this->n = x;
	col = 0;
	matrix.resize(n);
	F.resize(n);
}

void Matrix::setMatrix(vector< vector<real> > A, int x, vector<real> B)
{
	n = x;
	matrix = A;
	F = B;
}

Matrix::Matrix(void)
{
}


Matrix::~Matrix(void)
{
}

vector<real> Matrix::Gauss( vector<real> B)
{
	
	for (int i = 1; i<n; i++)
		for (int j=i; j<n; j++)
		{
			real m = matrix[j][i-1]/matrix[i-1][i-1];
			for (int k = 0; k < n; k++)
				matrix[j][k]=matrix[j][k]-m*matrix[i-1][k];
			B[j] = B[j] - m * B[i-1]; 
		}
	for (int i = n-1; i >= 0; i--)
		{
			dubl buf = 0;
			for (int j = i+1; j < n; j++)
			{
				buf += matrix[i][j]*B[j];
			}
			B[i] = B[i] - buf;
			B[i] = B[i]/matrix[i][i];
		}
	return B;
}

void Matrix::Gilbert()
{
	for (int i = 0; i < n; i++)
		for (int j = 0; i < n; j++)
			matrix[i][j] = 1.0/(i+j+1);
}

void Matrix::ToProf(MatrixProf *A)
{
	getCol();
	vector <real> bdi(n);
	vector <int> bia(n+1);
	vector <real> bau(col);
	vector <real> bal(col);
	int s = 0;
	int flag;
	for (int i = 0; i < n; i++)
	{
		bdi[i] = matrix[i][i];
		bia[i] = s;
		flag = 0;
		for (int j = 0; j < i; j++)
		{
			if (matrix[i][j] != 0 || matrix[j][i] != 0)
			{
				flag = 1;
			}
			if (flag == 1)
			{
				bau[s] = matrix[i][j];
				bal[s] = matrix[j][i];
				s++;
			}
		}
	}
	bia[n] = s;

}


void Matrix::getCol()
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
		{
			if (matrix[i][j] != 0 || matrix[j][i] != 0)
				col++;
		}
}


#include "Matrix.h"


MatrixProf::MatrixProf(void)
{
	n = 0;
	col = 0;
	vector<real> di(n);
	vector<int> ia(n+1);
	vector<real> al(col);
	vector<real> au(col);
	vector<real> F(n);
}

MatrixProf::MatrixProf(int x, int c)
{
	n = x;
	col = c;
	vector<real> di(n);
	vector<int> ia(n+1);
	vector<real> al(col);
	vector<real> au(col);
	vector<real> F(n);
}


MatrixProf::~MatrixProf(void)
{
}

void MatrixProf::load(ifstream &size, ifstream &matrix, ifstream &vect)
{
	size >> n >> col;	
	for( int i = 0; i < n+1; i++)
	{
		matrix >> ia[i];
	}
	for( int i = 0; i < n; i++)
	{
		matrix >> di[i];
	}
	for( int i = 0; i < col; i++)
	{
		matrix >> al[i];
	}
	for( int i = 0; i < col; i++)
	{
		matrix >> au[i];
	}
	for( int i = 0; i < n; i++)
	{
		vect >> F[i];
	}
}

void MatrixProf::LUDec()
{
	for(int i = 0; i < n; i++){
		int i0 = ia[i]; // индекс первого элемента строки 
		int i1 = ia[i+1]; // индекс первого элемента следующей строки
		int j = i - (ia[i+1]-ia[i]); // номер первого ненулевого элемента в реальной матрице или количество нулевых элементов в строке
		for(int k = ia[i]; k < ia[i+1]; k++,j++)
		{ 
			int ki = ia[i]; 
			int kj = ia[j];
			int dif = k - ia[i] - ia[j+1] + ia[j]; 
			if(dif < 0)
				kj += abs(dif);
			else
				ki += dif;
			for(ki; ki<k; ki++,kj++){
				al[k] -= al[ki]*au[kj];
				au[k] -= au[ki]*al[kj];
			}
			au[k] = au[k] / di[j];
			di[i] -= al[k]*au[k];
		}
		if (di[i] == 0)
			throw 1;
	}
}

void MatrixProf::ToTight(Matrix *A)
{
	vector< vector< real > > B;
	for (int i = 0; i < n; i++)
	{
		vector<real> buf;
		for (int j = 0; j < n; j++)
			buf.push_back(0);
		B.push_back(buf);
	}
	for(int i = 0; i < n; i++)
	{
		int j = i - (ia[i+1] - ia[i]);
		for (int k = ia[i]; k < ia[i+1]; k++, j++)
		{
			B[i][j] = al[k];
			B[j][i] = au[k];
		}
		B[i][i] = di[i];
	}
	A->setMatrix(B,n);

}

vector<dubl> MatrixProf::Direct()
{
	vector<dubl> bufer(n);
	for (int i = 0; i < n; i++)
	{
		int j = i - (ia[i+1]-ia[i]);
		for (int k = ia[i]; k < ia[i+1]; k++, j++)
		{
			bufer[i] += bufer[j] * al[k];
		}
		bufer[i] = (F[i] - bufer[i]) / di[i];
	}
	return bufer;
}

void MatrixProf::Reverse(vector<real> F, vector<real> *X)
{
	vector<real> bufer(n);
	bufer = F;
	for(int i = n - 1; i > 0; i--)
	{
		int j = i - (ia[i+1]-ia[i]);
		for (int k = ia[i]; ia[i] < ia[i+1]; k++, j++)
		{
			bufer[j] -= bufer[i]*au[k];
		}
	}
	*X = bufer;
}


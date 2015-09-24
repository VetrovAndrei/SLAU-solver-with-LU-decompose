#include "Matrix.h"

void main()
{
	setlocale (LC_CTYPE, "Russian");
	ifstream size ("size.txt");
	ifstream matr("matrix.txt");
	ifstream vect("vector.txt");
	int n, col;
	size >> n >> col;
	MatrixProf myMat(n, col);
	vector<dubl> X(n);
	Matrix tightMat(n);
	myMat.load(matr,vect);
	try
	{
		myMat.ToTight(&tightMat);
		X = myMat.SLAU();
	}
	catch(int error)
	{
		switch (error)
		{
		case 1:
			{
				cout << "что-то пошло не так";
				system("pause");
				break;
			}
		}
	}
}
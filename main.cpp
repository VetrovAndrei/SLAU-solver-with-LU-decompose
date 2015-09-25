#include "Matrix.h"

void main()
{
	setlocale (LC_CTYPE, "Russian");
	ifstream size ("size.txt");
	ifstream matr("matrix.txt");
	ifstream vect("vector.txt");
	ofstream proffile("prof.txt");
	ofstream tightfile("tight.txt");
	int n, col;
	size >> n >> col;
	MatrixProf myMat(n, col);
	vector<dubl> X(n);
	vector<real> Y(n);
	Matrix tightMat(n, vect);
	int flag = 0;
	try
	{
		cout << "1 - LU из файла. 2 - гильберт в плотной и в LU. 3" << endl;
		cin >> flag;
		switch(flag)
		{
		case 1:
			{
				myMat.load(matr,vect);
				myMat.ToTight(&tightMat);
				X = myMat.SLAU();
				Y = tightMat.Gauss();
				tightfile.precision(7);
				proffile.precision(7);
				for (int i = 0; i < n; i++)
				{
					tightfile << Y[i] << endl;
					proffile << X[i] << endl;
				}

				break;
			}
		case 2:
			{
				tightMat.Gilbert();
				tightMat.ToProf(&myMat);
				X = myMat.SLAU();
				Y = tightMat.Gauss();
				tightfile.precision(7);
				proffile.precision(7);
				for (int i = 0; i < n; i++)
				{
					tightfile << Y[i] << endl;
					proffile << X[i] << endl;
				}
				break;
			}
		default:
			{

			}



		}

		
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
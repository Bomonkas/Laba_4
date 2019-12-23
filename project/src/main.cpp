#include "eigenvalue.h"

int main()
{
	TYPE **A;
	int n = file_input(&A, "tests/D1.txt");	

	if (n)
	{
		cout << "Source matrix:" << endl
			 << endl;
		print_m(A, n);

		TYPE **H = new TYPE *[n];
		TYPE **A1 = new TYPE *[n];
		for (int i = 0; i < n; i++)
		{
			H[i] = new TYPE[n];
			A1[i] = new TYPE[n];
		}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A1[i][j] = A[i][j];

		cout << endl;
		//приводим матрицу а к форме ’ессенберга
		H = hessenberg(A1, n);
		cout << "Hessenberg matrix:" << endl
			 << endl;
		print_m(H, n);

		TYPE *lambda; //массив собственных значений
		lambda = new TYPE[n];

		TYPE **x1; //матрица, столбцами которой €вл€ютс€ собственные векторы
		x1 = new TYPE *[n];
		for (int i = 0; i < n; i++)
			x1[i] = new TYPE[n];

		//поиск собвтвенных значений
		lambda = EigenValues(H, n);

		//поиск собственных векторов методом обратных итераций
		x1 = EigenVectors(A, lambda, n);

		H = hessenberg(A1, n);

		//проверка
		chek(A, H, lambda, x1, n);

		TYPE *lambdaR; //массив собственных значений (поиск отношением –еле€)
		lambdaR = new TYPE[n];

		//поиск собственных значений и собственных векторов с помощью отношени€ –эле€
		RayleighRatio(A, x1, lambdaR, n);
		//проверка
		chek(A, H, lambda, x1, n);

		delete_m(A, n);
		delete_m(A1, n);
		delete_m(H, n);
		delete_m(x1, n);
		delete[] lambda;

	} //end if (flag)
	return (0);
}
#include "eigenvalue.h"

TYPE *gauss(TYPE **a1, TYPE *b1, int p) //реализация функции double *gauss (метода Гаусса)
{
	TYPE *x, buf, t;
	int flag;
	x = new TYPE[p]; //массив неизвестных

	//cout<<"Метод Гаусса:"<<endl;

	//выбор ведущего элемента
	for (int k = 0; k < p; k++)
	{
		//поиск ведущего элемента по столбцу
		int imax = k;
		for (int i = k; i < p; ++i)
			if (fabs(a1[i][k]) > fabs(a1[imax][k]))
				imax = i;

		flag = 1;
		if (imax != k)
		{
			for (int j = k; j < p + 1; j++)
			{
				buf = a1[k][j];
				a1[k][j] = a1[imax][j];
				a1[imax][j] = buf;
			} //end for(int j = k; j < n + 1; j++)
			t = b1[k];
			b1[k] = b1[imax];
			b1[imax] = t;
		} //end if (imax != k)

		//прямой ход метода Гаусса
		if (flag)
		{
			for (int i = k + 1; i < p; i++)
			{
				buf = a1[i][k] / a1[k][k];
				b1[i] -= b1[k] * buf;
				for (int j = k; j < p; j++)
					a1[i][j] -= a1[k][j] * buf;
			} //end for (int i = k + 1; i < n; i++)
			x = obrhod(a1, b1, x, p);

		} //end if (flag)
	}	 //end for (int k = 0; k < n; k++)

	return x;
} //end double * gauss(double **a, int n, int m)

int qrmeth(TYPE **a, TYPE **Q, TYPE **R, TYPE **T, int m) //метод QR-разложения
{
	TYPE
		**E, //единичная матрица
		s,   //переменная для матрицы поворота
		c,   //переменная для матрицы поворота
		v, w, znam;
	E = new TYPE *[m];

	for (int i = 0; i < m; i++)
		E[i] = new TYPE[m];

	for (int i = 0; i < m; i++) //задание единичной матрицы
		for (int j = 0; j < m; j++)
			if (i == j)
				E[i][j] = 1.0;
			else
				E[i][j] = 0.0;

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
		{
			R[i][j] = a[i][j];
			T[i][j] = E[i][j];
			Q[i][j] = E[i][j];
		} //end for (int j=0; j<m; j++)

	//прямой ход QR-разложения
	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			T[i][j] = E[i][j];

	for (int i = 0; i < m - 1; i++)
	{
		for (int j = i + 1; j < m; j++)
		{
			znam = sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
			if (znam > 0.1e-010)
			{
				c = R[i][i] / znam;
				s = R[j][i] / znam;
				R[i][i] = znam; //текущий элемент на диагонали равен знаменателю коэффициентов c и s
				R[j][i] = 0.0;  //элемент, стойщий в i-ом столбце, j-ой строке, под текущим диагональным, обнуляется

				//перемножаем матрицу T и матрицу R
				for (int k = i + 1; k < m; ++k)
				{
					v = R[i][k];
					w = R[j][k];
					R[i][k] = c * v + s * w;
					R[j][k] = -s * v + c * w;
				} //end for (int k=i+1; k<m; ++k)

				for (int k = 0; k < m; ++k)
				{
					v = T[i][k];
					w = T[j][k];
					T[i][k] = c * v + s * w;
					T[j][k] = -s * v + c * w;
				} //end for (int k=0; k<m; ++k)

			} //end (znam>0.1e-010)
		}	 //end for for (int j=i+1; j<m; j++)
	}		  //end for (int i = 0; i< m-1; ++i)

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			Q[i][j] = T[j][i];

	for (int i = 0; i < m; i++)
		delete[] E[i];
	delete[] E;

	return 0;
} //end double *qrmeth(double **a, double *b, int n)

TYPE *obrhod(TYPE **a1, TYPE *y1, TYPE *x2, int n) //обратный ход метода Гаусса
{
	TYPE t;

	for (int j = n - 1; j >= 0; j--)
	{
		t = y1[j];
		for (int i = j + 1; i < n; i++)
			t = t - a1[j][i] * x2[i];
		x2[j] = t / a1[j][j];
		if (fabs(x2[j]) < eps)
			x2[j] = 0.0;
	} //for

	return x2;
} //end TYPE *obrhod (TYPE **a1, TYPE *y1, TYPE *x2, int n)

// TYPE **inverse_matrix(TYPE **a1, TYPE **Q, TYPE **R, TYPE **T, TYPE *b1, int m) //вычисляет обратную матрицу
// {
// 	TYPE **E,   //единичная матрица
// 		**a1_obr, //обратная матрица
// 		*x,		  //неизвестные при решении СЛАУ R*x=p
// 		*e,		  //столбец единичной матрицы
// 		*p;		  //вектор равный T*e

// 	E = new TYPE *[m];
// 	a1_obr = new TYPE *[m];
// 	x = new TYPE[m];
// 	p = new TYPE[m];
// 	e = new TYPE[m];

// 	for (int i = 0; i < m; i++)
// 	{
// 		E[i] = new TYPE[m];
// 		a1_obr[i] = new TYPE[m];
// 	} //end for (int i=0; i<m; i++)

// 	qrmeth(a1, Q, R, T, m);

// 	for (int i = 0; i < m; i++) //задаем единичную матрицу
// 		for (int j = 0; j < m; j++)
// 		{
// 			if (i == j)
// 				E[i][j] = 1.0;
// 			else
// 				E[i][j] = 0.0;
// 		} //end for (int j=0; j<m; j++)

// 	for (int i = 0; i < m; i++) //ищем обратную матрицу, решая систему Rx=p
// 	{
// 		for (int j = 0; j < m; j++)
// 			e[j] = E[j][i];

// 		for (int k = 0; k < m; k++) //умножем столбец единичной матрицы слева на T
// 		{
// 			TYPE t = 0.0;
// 			t = T[k][i] * e[i];
// 			p[k] = t;

// 			/*for (int j=0;j<m; j++)
// 				    t+=T[j][k]*e[k];
// 			    p[k]=t;*/
// 		} //end for (int k=0; k<m; k++)

// 		x = obrhod(R, p, x, m);
// 		for (int j = 0; j < m; j++)
// 			a1_obr[j][i] = x[j];
// 	} //end for (int i=0; i<m; i++)

// 	for (int i = 0; i < m; i++)
// 		delete[] E[i];
// 	delete[] E;
// 	delete[] e;
// 	delete[] x;
// 	delete[] p;

// 	return a1_obr;
// } //end TYPE **obr_matrix(TYPE **a1, TYPE **Q, TYPE **R, TYPE **T, int m)

TYPE **hessenberg(TYPE **a1, int m) //приведение матрицы к форме Хессенберга
{
	TYPE
		**a2, //матрица в форме Хессенберга
		s,	//переменная для матрицы поворота
		c,	//переменная для матрицы поворота
		v, w, znam;
	//int r;
	//q = new TYPE[m];
	//x = new TYPE[m];
	a2 = new TYPE *[m];

	for (int i = 0; i < m; i++)
		a2[i] = new TYPE[m];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			a2[i][j] = a1[i][j];

	for (int i = 1; i < m - 1; i++)
	{
		for (int j = i + 1; j < m; j++)
		{
			znam = sqrt(a2[i][i - 1] * a2[i][i - 1] + a2[j][i - 1] * a2[j][i - 1]);

			//прямой ход QR-разложения
			if (znam > 0.1e-010)
			{
				c = a2[i][i - 1] / znam;
				s = a2[j][i - 1] / znam;
				a2[i][i - 1] = znam; //текущий элемент на диагонали равен знаменателю коэффициентов c и s
				a2[j][i - 1] = 0.0;  //элемент, стойщий в i-ом столбце, j-ой строке, под текущим диагональным, обнуляется

				for (int k = i; k < m; ++k)
				{
					v = a2[i][k];
					w = a2[j][k];
					a2[i][k] = v * c + w * s;
					a2[j][k] = -v * s + w * c;
				} //end for (int k=0; k<m; ++k)

				a2[i - 1][i] = znam;
				a2[i - 1][j] = 0.0;

				for (int k = i; k < m; ++k)
				{
					v = a2[k][i];
					w = a2[k][j];
					a2[k][i] = v * c + w * s;
					a2[k][j] = -v * s + w * c;
				} //end for (int k=0; k<m; ++k)

			} //end (znam>0.1e-010)
		}	 //end for for (int j=i+1; j<m; j++)
	}		  //end for (int i = 0; i< m-1; ++i)

	return a2;

} //end TYPE **hessenberg(TYPE **a, int m)

TYPE *EigenValues(TYPE **H1, int n) //поиск собственных значений
{
	int *iter; //массив итераций для поиска собственных значений
	TYPE **a1_t, **Q, **R, **T, d;

	iter = new int[n];
	iter[0] = 0;

	for (int i = n - 1; i > 0; i--)
	{
		TYPE lambda_t = 0.0;
		a1_t = new TYPE *[i + 1];
		for (int j = 0; j < i + 1; j++)
			a1_t[j] = new TYPE[n];

		for (int j = 0; j < i + 1; j++)
			for (int p = 0; p < i + 1; p++)
				a1_t[j][p] = H1[j][p];

		iter[i] = 0;
		do
		{
			if (iter[i] > 0)
				lambda_t = a1_t[i][i];
			iter[i]++;
			d = /*0.0;*/ a1_t[i][i];
			Q = new TYPE *[i + 1];
			R = new TYPE *[i + 1];
			T = new TYPE *[i + 1];
			for (int j = 0; j < i + 1; j++)
			{
				Q[j] = new TYPE[i + 1];
				R[j] = new TYPE[i + 1];
				T[j] = new TYPE[i + 1];
			} //end for (int j=0; j<i+1; j++)
			qrmeth(shift_minus(a1_t, d, i + 1), Q, R, T, i + 1);
			a1_t = shift_plus(mult_m_m(R, Q, i + 1), d, i + 1);
		} //end do
		while (fabs(a1_t[i][i] - lambda_t) >= eps);
		for (int j = 0; j < i + 1; j++)
			for (int k = 0; k < i + 1; k++)
				H1[j][k] = a1_t[j][k];
	} //end for (int i=n-1; i>0; i--)
	cout << "Матрица с собственными значениями на диагонали:" << endl
		 << endl;
	print_m(H1, n);
	cout << endl;

	TYPE *lambda1;
	lambda1 = new TYPE[n];

	for (int i = 0; i < n; i++)
		lambda1[i] = H1[i][i];
	cout << endl
		 << endl;

	cout << "Приближенные собственные зачения:" << endl
		 << endl;
	for (int i = 0; i < n; i++)
		cout << i + 1 << ". " << lambda1[i] << ";   "
			 << "число итераций: " << iter[i] << endl;

	delete[] iter;

	return lambda1;
} //end TYPE *EigenValues(TYPE **H1, int n)

TYPE **EigenVectors(TYPE **a, TYPE *lambda, int n) //поиск собственных векторов
{
	int *iter; //массив итераций для поиска собственных значений
	iter = new int[n];

	//ищем собственные векторы методом обратных итераций
	TYPE *z, *y, *x, **a2, **x1;
	z = new TYPE[n];
	y = new TYPE[n];
	x = new TYPE[n];
	a2 = new TYPE *[n];
	x1 = new TYPE *[n];

	for (int i = 0; i < n; i++)
	{
		a2[i] = new TYPE[n];
		x1[i] = new TYPE[n];
	} //end for (int i=0; i<n; i++)

	for (int i = 0; i < n; i++)
	{
		iter[i] = 0;
		for (int j = 0; j < n; j++)
			if (i == j)
				y[j] = 1.0;
			else
				y[j] = 0.0;
		/*y[0]=-0.0119071;
			y[1]=0.709759;
			y[2]=-0.607555;
			y[3]=0.356337;*/

		a2 = shift_minus(a, lambda[i], n);

		do
		{
			if (iter[i] > 0)
				for (int j = 0; j < n; j++)
					y[j] = x[j];
			iter[i] += 1;
			z = gauss(a2, y, n);

			for (int j = 0; j < n; j++)
				x[j] = z[j] / sqrt(EuclidMult(z, z, n));
		} //end do
		while ((1 - fabs(EuclidMult(x, y, n) / (sqrt(EuclidMult(x, x, n)) * sqrt(EuclidMult(y, y, n))))) >= eps);

		for (int j = 0; j < n; j++)
			x1[i][j] = x[j];
	} //for (int i=0; i<n; i++)

	cout << endl;

	cout << "Приближенные собственные векторы: " << endl
		 << endl;
	for (int i = 0; i < n; i++)
	{
		cout << i + 1 << ". Число итераций: " << iter[i] << endl;
		print_v(x1[i], n);
	} //end for (int i=0; i<n; i++)

	delete[] x;
	delete[] y;
	delete[] z;
	delete[] iter;
	for (int i = 0; i < n; i++)
		delete[] a2[i];
	delete[] a2;

	return x1;
} //end TYPE *EigenVectors(TYPE **a, TYPE *lambda, int n)

TYPE **shift_minus(TYPE **a1, TYPE d1, int m) //сдвиг диагонального элемента матрицы на -d
{
	TYPE **tmp;
	tmp = new TYPE *[m];
	for (int i = 0; i < m; i++)
		tmp[i] = new TYPE[m];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
		{
			if (i == j)
				tmp[i][j] = a1[i][j] - d1;
			else
				tmp[i][j] = a1[i][j];
		} //end for (int j=0; j<m; j++)

	return tmp;
} //end TYPE **shift_minus(TYPE **a1, TYPE d1, int m)

TYPE **shift_plus(TYPE **a1, TYPE d1, int m) //сдвиг диагонального элемента матрицы на +d
{
	TYPE **tmp;
	tmp = new TYPE *[m];
	for (int i = 0; i < m; i++)
		tmp[i] = new TYPE[m];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
		{
			if (i == j)
				tmp[i][j] = a1[i][j] + d1;
			else
				tmp[i][j] = a1[i][j];
		} //end for (int j=0; j<m; j++)

	return tmp;
} //end TYPE **shift_plus(TYPE **a1, TYPE d1, int m)

TYPE EuclidMult(TYPE *b1, TYPE *b2, int m) //Евклидово скалярное произведение
{
	TYPE sc = 0.0;

	for (int i = 0; i < m; i++)
		sc += b1[i] * b2[i];

	return sc;
} //end TYPE EuclidMult(TYPE *b1, TYPE *b2, m)

TYPE **det_gauss(TYPE **a1, int m) //вычисление определителя методом Гаусса
{
	TYPE **temp, buf;
	temp = new TYPE *[m];
	for (int i = 0; i < m; i++)
		temp[i] = new TYPE[m];

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			temp[i][j] = a1[i][j];

	for (int i = 0; i < m - 1; i++)
		for (int j = i + 1; j < m; j++)
		{
			buf = temp[j][i] / temp[i][i];
			for (int k = 0; k < m; k++)
				temp[j][k] = temp[i][k] * buf - temp[j][k];
		} //end for (int j=i+1, j<n; j++)

	return temp;
} //end TYPE det_gauss(TYPE **a1, int m)

TYPE trace(TYPE **a1, int m) //след матрицы
{
	TYPE tr = 1.0;
	for (int i = 0; i < m; i++)
		tr *= a1[i][i];

	return tr;
} //end TYPE trace(TYPE **a1, int m)

void chek(TYPE **a1, TYPE **H1, TYPE *lambd, TYPE **x1, int m) //проверка
{
	TYPE **x, *temp, **a2;
	x = new TYPE *[m];
	a2 = new TYPE *[m];
	temp = new TYPE[m];

	for (int i = 0; i < m; i++)
	{
		x[i] = new TYPE[m];
		a2[i] = new TYPE[m];
	} //end for (int i=0; i<m; i++)

	for (int i = 0; i < m; i++)
		for (int j = 0; j < m; j++)
			a2[i][j] = a1[i][j];

	cout << "Проверка:" << endl
		 << endl;
	for (int i = 0; i < m; i++)
	{
		cout << "det(A-lambda_" << i + 1 << " E)=" << fabs(trace(det_gauss(shift_minus(H1, lambd[i], m), m), m)) << endl
			 << endl;
		temp = mult_v_n(x1[i], lambd[i], m);
		x[i] = gauss(a2, temp, m);
		for (int j = 0; j < m; j++)
			for (int k = 0; k < m; k++)
				a2[j][k] = a1[j][k];
	} //end for (int i=0; i<m; i++)

	for (int i = 0; i < m; i++)
	{
		cout << "Вектор ошибки " << i + 1 << " собственного вектора: " << endl;
		print_v(diff_v(x1[i], x[i], m), m);
	} //end for (int i=0; i<m; i++)
	cout << endl
		 << endl;
} //end  TYPE *chek(TYPE **a1, TYPE **H1, TYPE *lambd, TYPE *x1, int m)

TYPE RayleighRatio(TYPE **a, TYPE **x1, TYPE *lambdaR, int m)
{
	TYPE *z, *y, *x, **a2;
	z = new TYPE[m];
	y = new TYPE[m];
	x = new TYPE[m];
	a2 = new TYPE *[m];

	int *iter; //массив итераций для поиска собственных значений
	iter = new int[m];

	for (int i = 0; i < m; i++)
		a2[i] = new TYPE[m];

	cout << "Поиск собственных значений и собственных векторов с помощью отношения Рэлея" << endl
		 << endl;
	//TYPE *lambdaR;
	//lambdaR=new TYPE [m];
	TYPE norm_R;

	for (int i = 0; i < m; i++)
	{
		iter[i] = 0;
		TYPE lambda_t = 0.0;
		//приближенное собственное значение
		/*for (int j=0; j<m; j++)
				y[j]=rand()%3 - 1;*/
		if (fil == 1)
		{
			if (i == 0)
			{
				y[0] = 1.0;
				y[1] = 0.0;
				y[2] = 0.0;
				y[3] = 0.5;
				/* y[0]=-0.0119071;
			       y[1]=0.709759;
			       y[2]=-0.607555;
			       y[3]=0.356337;*/
			} //end if (i == 0)
			if (i == 1)
			{
				y[0] = 0.5;
				y[1] = 1.0;
				y[2] = -1.0;
				y[3] = 0.0;
			} //end if (i == 1)
			if (i == 2)
			{
				y[0] = 0.0;
				y[1] = 1.0;
				y[2] = 1.0;
				y[3] = 1.0;
			} //end if (i == 2)
			if (i == 3)
			{
				y[0] = 0.0;
				y[1] = -1.0;
				y[2] = 0.0;
				y[3] = 1.0;
			} //end if (i == 3)
		}	 //end if (fil==1)
		else if (fil == 2)
		{
			if (i == 0)
			{
				y[0] = -1.0;
				y[1] = 0.0;
				y[2] = 0.0;
				y[3] = 0.0;
			} //end if (i == 0)
			if (i == 1)
			{
				y[0] = 0.0;
				y[1] = -1.0;
				y[2] = -0.1;
				y[3] = 0.3;
			} //end if (i == 1)
			if (i == 2)
			{
				y[0] = 0.0;
				y[1] = 0.0;
				y[2] = -1.0;
				y[3] = 0.0;
			} //end if (i == 2)
			if (i == 3)
			{
				y[0] = 0.0;
				y[1] = 0.0;
				y[2] = 0.5;
				y[3] = 1.0;
			} //end if (i == 3)
		}	 //end else if (fil==2)
		else
		{
			if (i == 3)
			{
				y[0] = 0.0;
				y[1] = 1.0;
				y[2] = 0.0;
				y[3] = -0.2;
			} //end if (i == 0)
			if (i == 2)
			{
				y[0] = 0.0;
				y[1] = 0.2;
				y[2] = 0.0;
				y[3] = 1.0;
			} //end if (i == 1)
			if (i == 1)
			{
				y[0] = 0.0;
				y[1] = 0.0;
				y[2] = -1.0;
				y[3] = 0.0;
			} //end if (i == 2)
			if (i == 0)
			{
				y[0] = 1.0;
				y[1] = 0.0;
				y[2] = 0.0;
				y[3] = 0.0;
			} //end if (i == 3)
		}	 //end else if (fil==2)

		norm_R = sqrt(EuclidMult(y, y, m));
		do
		{
			if (iter[i] > 0)
			{
				for (int j = 0; j < m; j++)
					y[j] = x[j];
				lambda_t = lambdaR[i];
				norm_R = sqrt(EuclidMult(z, z, m));
			} //end if (iter[i]>0)
			iter[i] += 1;
			lambdaR[i] = EuclidMult(mult_m_v(a, y, m), y, m);

			a2 = shift_minus(a, lambdaR[i], m);
			z = gauss(a2, y, m);

			for (int j = 0; j < m; j++)
				x[j] = z[j] / sqrt(EuclidMult(z, z, m));
		} //end do
		while (((1 - fabs(EuclidMult(x, y, m) / (sqrt(EuclidMult(x, x, m)) * sqrt(EuclidMult(y, y, m))))) >= eps) && (fabs(lambda_t - lambdaR[i]) >= eps));

		for (int j = 0; j < m; j++)
			x1[i][j] = x[j];
	} //end for (int i=0; i<m; i++)

	//cout<<"Приближенные собственные зачения и собственные вектора:"<<endl<<endl;
	for (int i = 0; i < m; i++)
	{
		cout << i + 1 << ". "
			 << "Собственное значение: " << lambdaR[i] << ";   "
			 << "число итераций: " << iter[i] << endl
			 << endl;
		cout << "Собственный вектор: " << endl;
		print_v(x1[i], m);
	}

	cout << "norm_R = " << norm_R << endl;
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] iter;
	for (int i = 0; i < m; i++)
		delete[] a2[i];
	delete[] a2;

	return 0;
} //end TYPE RayleighRatio(TYPE **a, TYPE **x1, TYPE *lambdaR, int m)
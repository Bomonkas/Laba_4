#include "matrix.h"

int     file_input(TYPE*** A, string name_file)
{
	ifstream file(name_file);
	if (!file)
	{
		cerr << "Error with opening file\n";
		return (-1);
	}
	string h;
	int size;
	file >> size;
	*A = new TYPE *[size];
	for (int i = 0; i < size; i++)
		(*A)[i] = new TYPE[size];
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			file >> (*A)[i][j];
	file.close();
	return (size);
}

TYPE    **copy_m(TYPE **matr, const int size)
{
	TYPE **copy_matr = new TYPE*[size];
	for (int i = 0; i < size; i++) {
		copy_matr[i] = new TYPE[size];
		for (int j = 0; j < size; j++) {
			copy_matr[i][j] = matr[i][j];
		}
	}
	return copy_matr;
}

TYPE    *copy_v(TYPE *vec, const int size)
{
	TYPE *copy_vec = new TYPE[size];
	for (int i = 0; i < size; i++) {
		copy_vec[i] = vec[i];
	}
	return copy_vec;
}

void    delete_m(TYPE** matr, const int size)
{
	if (!matr)
	{
		cout << "Matrix is not exist\n";
		return;
	}
	for (int i = 0; i < size; i++)
		delete[] matr[i];
	delete[] matr;
}

TYPE    norm_inf_v(TYPE *vec, const int size)
{
	TYPE max = fabs(vec[0]);

	for (int i = 1; i < size; i++)
		if (max < fabs(vec[i]))
			max = fabs(vec[i]);
	return (max);
}

TYPE    norm_1_v(TYPE *vec, const int size)
{
	TYPE sum = 0.0;
	for (int i = 0; i < size; i++)
		sum += fabs(vec[i]);
	return (sum);
}

TYPE    norm_inf_m(TYPE **A, const int size)
{
	if (!A)
	{
		cout << "Matrix is not exist\n";
		return 0;
	}
	double norm = 0, sum = 0;

	for (int i = 0; i < size; i++)
	{
		sum = 0;
		for (int j = 0; j < size; j++)
			sum += fabs(A[i][j]);
		if (sum > norm)
			norm = sum;
	}
	return (norm);
}

TYPE    norm_1_m(TYPE **A, const int size)
{
	if (!A)
	{
		cout << "Matrix is not exist\n";
		return (0);
	}
	double norm = 0, sum = 0;
	for (int i = 0; i < size; i++)
	{
		sum = 0;
		for (int j = 0; j < size; j++)
			sum += fabs(A[j][i]);
		if (sum > norm)
			norm = sum;
	}
	return (norm);
}

int     is_zero(TYPE elem)
{
	if (fabs(elem) <= eps)
		return (1);
	return (0);
}

TYPE    **transp_m(TYPE** matr, int size)
{
	TYPE tmp;
	for (int i = 0; i < size; i++)
	{
		for (int j = i; j < size; j++)
		{
			tmp = matr[i][j];
			matr[i][j] = matr[j][i];
			matr[j][i] = tmp;
		}
	}
	return (matr);
}

TYPE	*mult_v_n(TYPE *vec, TYPE n, int size)
{
	TYPE *new_vec = new TYPE[size];
	for (int i = 0; i < size; i++)
		new_vec[i] = vec[i] * n;
	return (new_vec);
}

TYPE	**mult_m_n(TYPE **matr, TYPE n, const int size)
{
	TYPE **new_matr = new TYPE *[size];
	for (int i = 0; i < size; i++)
	{
		new_matr[i] = new TYPE[size];
		for (int j = 0; j < size; j++)
			new_matr[i][j] = matr[i][j] * n;
	}
	return (new_matr);
}

TYPE 	*mult_m_v(TYPE **matr, TYPE *vec, const int size)
{
	TYPE *res_vec = NULL;
	res_vec = new TYPE[size];

	for (int i = 0; i < size; i++)
	{
		TYPE sum = 0;
		for (int j = 0; j < size; j++)
			sum += matr[i][j] * vec[j];
		if (is_zero(sum))
			res_vec[i] = 0;
		else
			res_vec[i] = sum;
	}
	return (res_vec);
}

TYPE 	**mult_m_m(TYPE **l_matr, TYPE **r_matr, const int size)
{
	TYPE **res_matr = NULL;
	res_matr = new TYPE *[size];
	for (int i = 0; i < size; i++)
		res_matr[i] = new TYPE[size];

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			TYPE sum = 0;
			for (int k = 0; k < size; k++)
				sum += l_matr[i][k] * r_matr[k][j];
			if (is_zero(sum))
				res_matr[i][j] = 0;
			else
				res_matr[i][j] = sum;
		}
	}
	return (res_matr);
}

TYPE    *diff_v(TYPE *x, TYPE *y, const int size)
{
	TYPE *diff = new TYPE[size];

	for (int i = 0; i < size; i++)
		diff[i] = x[i] - y[i];
	return (diff);
}

void    print_v(TYPE* vec, const int size)
{
	if (!vec)
	{
		cout << "Matrix is degenerated" << endl;
		return;
	}
	for (int i = 0; i < size; i++)
		cout << vec[i] << " ";
	cout << endl;
}

void    print_m(TYPE** matr, const int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			cout << fixed << setprecision(4) << setw(13) << matr[i][j] << ' ';
		cout << endl;
	}
	cout << endl;
}
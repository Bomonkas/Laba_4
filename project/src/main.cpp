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
		//�������� ������� � � ����� �����������
		H = hessenberg(A1, n);
		cout << "Hessenberg matrix:" << endl
			 << endl;
		print_m(H, n);

		TYPE *lambda; //������ ����������� ��������
		lambda = new TYPE[n];

		TYPE **x1; //�������, ��������� ������� �������� ����������� �������
		x1 = new TYPE *[n];
		for (int i = 0; i < n; i++)
			x1[i] = new TYPE[n];

		//����� ����������� ��������
		lambda = EigenValues(H, n);

		//����� ����������� �������� ������� �������� ��������
		x1 = EigenVectors(A, lambda, n);

		H = hessenberg(A1, n);

		//��������
		chek(A, H, lambda, x1, n);

		TYPE *lambdaR; //������ ����������� �������� (����� ���������� �����)
		lambdaR = new TYPE[n];

		//����� ����������� �������� � ����������� �������� � ������� ��������� �����
		RayleighRatio(A, x1, lambdaR, n);
		//��������
		chek(A, H, lambda, x1, n);

		delete_m(A, n);
		delete_m(A1, n);
		delete_m(H, n);
		delete_m(x1, n);
		delete[] lambda;

	} //end if (flag)
	return (0);
}
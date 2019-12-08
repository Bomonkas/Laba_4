#include "eigenvalue.h"

void	Hessenberg(TYPE **A, const int size)
{
	for (int i = 1; i < size - 1; i++)
		for(int j = i + 1; j < size; j++)	
		{
			double sq = sqrt(A[i][i-1] * A[i][i-1] + A[j][i-1] * A[j][i-1]);
			if (sq > eps)
			{
				double c = A[i][i-1] / sq;
				double s = A[j][i-1] / sq;
				double *a = new double [size];
				for (int k = 0; k < size; k++)
				{
					a[k] = A[i][k] * c + A[j][k] * s;
					A[j][k] = c * A[j][k] - s * A[i][k];
				}
				delete[] A[i];
				A[i] = a;
				double exchange;
				for (int k = 0; k < size; k++)
				{
					exchange = A[k][i] * c + A[k][j] * s;
					A[k][j]=c * A[k][j] - s * A[k][i];
					A[k][i]=exchange;
				}
			}
		}
}
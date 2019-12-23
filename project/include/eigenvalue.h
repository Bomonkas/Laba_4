#pragma once

#include <iostream>
#include <cmath>
// #include <iomanip>
#include <fstream>
#include <cstring>
#include "matrix.h"

int fil = 1;

TYPE 	*obrhod(TYPE **a1, TYPE *y1, TYPE *x, int n);     							//�������� ��� ������ ������  
TYPE 	*gauss(TYPE **a, TYPE *b, int n);											//������ ��� ������ ������
int 	qrmeth(TYPE **a, TYPE **Q, TYPE **R, TYPE **T, int n);						//QR-����������
// TYPE 	**inverse_matrix(TYPE **a1, TYPE **Q, TYPE **R, TYPE **T, TYPE *b1, int m); //���������� �������� �������
TYPE 	**hessenberg(TYPE **a1, int m);
TYPE 	*EigenValues(TYPE **H1, int n);												//����� ����������� ��������
TYPE 	**EigenVectors(TYPE **a, TYPE *lambda, int n);								//����� ����������� ��������
TYPE 	**shift_minus(TYPE **a1, TYPE d1, int m);									//����� ������������� �������� ������� �� -d
TYPE 	**shift_plus(TYPE **a1, TYPE d1, int m);									//����� ������������� �������� ������� �� +d
TYPE 	EuclidMult(TYPE *b1, TYPE *b2, int m);										//��������� ��������� ������������
TYPE 	**det_gauss(TYPE **a1, int m);												//���������� ������������ ������� ������
TYPE 	trace(TYPE **a1, int m);													//���� �������
void 	chek(TYPE **a1, TYPE **H1, TYPE *lambd, TYPE **x1, int m);					//�������� ����������� �������� � ����������� �����
TYPE 	RayleighRatio(TYPE **a, TYPE **x1, TYPE *lambdaR, int m);


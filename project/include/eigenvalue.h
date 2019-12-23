#pragma once

#include <iostream>
#include <cmath>
// #include <iomanip>
#include <fstream>
#include <cstring>
#include "matrix.h"

int fil = 1;

TYPE 	*obrhod(TYPE **a1, TYPE *y1, TYPE *x, int n);     							//обратный ход метода Гаусса  
TYPE 	*gauss(TYPE **a, TYPE *b, int n);											//прямой ход метода Гаусса
int 	qrmeth(TYPE **a, TYPE **Q, TYPE **R, TYPE **T, int n);						//QR-разложение
// TYPE 	**inverse_matrix(TYPE **a1, TYPE **Q, TYPE **R, TYPE **T, TYPE *b1, int m); //вычисление обратной матрицы
TYPE 	**hessenberg(TYPE **a1, int m);
TYPE 	*EigenValues(TYPE **H1, int n);												//поиск собственных значений
TYPE 	**EigenVectors(TYPE **a, TYPE *lambda, int n);								//поиск собвтвенных векторов
TYPE 	**shift_minus(TYPE **a1, TYPE d1, int m);									//сдвиг диагонального элемента матрицы на -d
TYPE 	**shift_plus(TYPE **a1, TYPE d1, int m);									//сдвиг диагонального элемента матрицы на +d
TYPE 	EuclidMult(TYPE *b1, TYPE *b2, int m);										//Евклидово скалярное произведение
TYPE 	**det_gauss(TYPE **a1, int m);												//вычисление определителя методом Гаусса
TYPE 	trace(TYPE **a1, int m);													//след матрицы
void 	chek(TYPE **a1, TYPE **H1, TYPE *lambd, TYPE **x1, int m);					//проверка собственных значений и собственных чисел
TYPE 	RayleighRatio(TYPE **a, TYPE **x1, TYPE *lambdaR, int m);


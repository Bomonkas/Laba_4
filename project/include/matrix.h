#pragma once
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>

#define TYPE double

#define eps 1e-7

using namespace std;

int     file_input(TYPE*** A, string name_file);
TYPE    **copy_m(TYPE **matr, const int size);
TYPE    *copy_v(TYPE *vec, const int size);
void    delete_m(TYPE** matr, const int size);

TYPE    norm_inf_v(TYPE *vec, const int size);
TYPE    norm_1_v(TYPE *vec, const int size);
TYPE    norm_inf_m(TYPE **A, const int size);
TYPE    norm_1_m(TYPE **A, const int size);

int     is_zero(TYPE elem);

TYPE    **transp_m(TYPE** matr, int size);
void	mult_v_n(TYPE *vec, TYPE n, int size);
void	mult_m_n(TYPE **matr, TYPE n, const int size);
TYPE    *mult_m_v(TYPE** matr, TYPE* vec, const int size);
TYPE    **mult_m_m(TYPE** l_matr, TYPE** r_matr, const int size);

TYPE    *diff_v(TYPE *x, TYPE *y, const int size);

void    print_v(TYPE* vec, const int size);
void    print_m(TYPE** matr, const int size);

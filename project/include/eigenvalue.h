#pragma once

#include <iostream>
#include <cmath>
// #include <iomanip>
// #include <fstream>
// #include <cstring>
#include "matrix.h"

void	QR(TYPE **A, const int size);
void	Shift(TYPE **A, const int size);
void	Hessenberg(TYPE **A, const int size);

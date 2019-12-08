#include "eigenvalue.h"

int main()
{
	TYPE **A;
	int n = file_input(&A, "tests/D1.txt");
	cout << "Source matrix A:\n";
	print_m(A, n);


	//Reduction of the matrix to the Hessenberg form
	Hessenberg(A, n);
	cout << "Matrix A in Hessenberg form:\n";
	print_m(A, n);

	// //Shift
	// Shift(A, n);
	// cout << "matrix A with shift:\n";
	// print_m(A, n);

	// //QR-decomposition
	// QR(A, n);

	delete_m(A, n);
	return (0);
}
#include <exception>
#include <iostream>
#include <sstream>

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <armadillo>

#include "Jacobi.hpp"

using namespace std;
using namespace arma;

int main(int argc, char const *argv[]) {
	unsigned int N = 10;
	double h =  1. / N;
	double d =  2. / (h * h);
	double a = -1. / (h * h);
	clock_t tick, tock;

	mat A;
	vec eigval;
	mat eigvec;
	mat eigvecs;
	double t1, t2;

	N = 8;
	for (unsigned int i = 0; i < 6; i++) {
		Jacobi J(a, d, N);
		A = J.input();
		eigvecs = mat(N - 1, 2);

		// Køyr Jacobi og ta tida
		tick = clock();
		J.run();
		tock = clock();
		t1 = (double)(tock - tick) / CLOCKS_PER_SEC;

		// Henta ut eigenvektoren høyrande til det minste eigenverdet
		eigvecs.col(0) = J.vectors(true).col(0);

		// Køyr utreknaren til Armadillo
		tick = clock();
		eig_sym(eigval, eigvec, A);
		tock = clock();
		t2 = (double)(tock - tick) / CLOCKS_PER_SEC;

		// Samla løysingane i ei matrise og skriv til fil
		ostringstream os;
		os << "evec_" << N << ".csv";
		eigvecs.col(1) = eigvec.col(0);
		eigvecs.save(os.str(), csv_ascii);

		// Meld ifrå om tida
		cout << "N=" << N << ": " << t1 << " sek. (Jacobi; " << J.iterations() << " iterasjonar), " << t2 << " sek. (Armadillo)" << endl;

		N *= 2;
	}

	return EXIT_SUCCESS;
}

#include <exception>
#include <iostream>
#include <sstream>

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <armadillo>

#include "Tridiagonal.hpp"

using namespace std;
using namespace arma;

double f(double x);
double u(double x);

int main(int argc, char const *argv[]) {
	if (argc < 5) {
		cout << "Programmet løyser eit likningssett på tridiagonal matriseform" << endl << endl
			 << "Bruk:" << endl
			 << "  " << argv[0] << " e a b c" << endl << endl
			 << "der" << endl
			 << "  n=10^e er storleiken på matrisa (n×n) og lengda på løysingsvektoren" << endl
			 << "  a er talet i det nedre bandet i matrisa" << endl
			 << "  b er talet i mellombandet i matrisa" << endl
			 << "  c er talet i det øvre bandet i matrisa" << endl << endl
			 << "Resultata vert skrivne til filene vgen_n.csv (generell algoritme)," << endl
			 << "vspec_n.csv (spesialisert algoritme) og u_n.csv (eksakt løysing)." << endl;
		return argc == 1 ? EXIT_SUCCESS : EXIT_FAILURE;
	}

	// Les kommandolinja
	unsigned int n;
	vec a, b, c;
	try {
		int i = stoi(argv[1]);
		if (i < 1) {
			throw invalid_argument("n må vera 1 eller høgre");
		}
		n = pow(10, i);

		a = vec(n - 1);
		a.fill(stoi(argv[2]));

		b = vec(n);
		b.fill(stoi(argv[3]));

		c = vec(n - 1);
		c.fill(stoi(argv[4]));
	} catch (exception &e) {
		cout << "Greidde ikkje å lesa argument: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	// Sett opp vektorar til diskretiseringa av x og f(x)
	double h = 1 / (n + 1.f);
	vec xi = vec(n);
	vec y = vec(n);
	for (unsigned int i = 0; i < n; i++) {
		xi(i) = (i + 1) * h;
		y(i)  = h*h * f(xi(i));
	}

	// Vektorar til løysingane
	vec vgen, vspec, ui;

	// Sett opp Tridiagonal-klassa og løys likningane
	try {
		Tridiagonal tridiag;
		tridiag.init(n, a, b, c);

		vgen = tridiag.solve(y);

		vspec = tridiag.solve_spec(y);
	} catch (exception &e) {
		cout << "Eitkvart gjekk ikkje rett føre seg: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	// Finn eksakt løysing
	ui = vec(n);
	for (unsigned int i = 0; i < n; i++) {
		ui(i) = u(xi(i));
	}

	#ifdef DEBUG
	ui.print("u=");
	#endif

	// Skriv filer
	ostringstream os;
	os << "vgen_" << n << ".csv";
	vgen.save(os.str(), csv_ascii);

	os.str("");
	os << "vspec_" << n << ".csv";
	vspec.save(os.str(), csv_ascii);

	os.str("");
	os << "u_" << n << ".csv";
	ui.save(os.str(), csv_ascii);

	return EXIT_SUCCESS;

// LU i Armadillo
{
	arma_rng::set_seed_random();
	mat A = randu<mat>(5, 5);
	vec b = randu<vec>(5);

	A.print("A=");
	b.print("b=");

	vec x = solve(A, b);
	x.print("x=");

	mat L, U, P;
	lu(L, U, P, A);
	L.print("L=");
	U.print("U=");
	P.print("P=");
	(P*A-L*U).print("Test of LU");

	U.save("U_mat.txt", csv_ascii);
	b.save("b_vec.txt", csv_ascii);
}
// End LU

	return EXIT_SUCCESS;
}

double f(double x) {
	return 100*exp(-10*x);
}

double u(double x) {
	return 1 - (1 - exp(-10))*x - exp(-10*x);
}

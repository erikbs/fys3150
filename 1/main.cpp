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
		((argc == 1) ? cout : cerr)
			 << "Programmet løyser eit likningssett på tridiagonal matriseform på tri måtar." << endl << endl
			 << "Bruk:" << endl
			 << "  " << argv[0] << " e a b c [LU=1]" << endl << endl
			 << "der" << endl
			 << "  n=10^e er storleiken på matrisa (n×n) og lengda på løysingsvektoren" << endl
			 << "  a er talet i det nedre bandet i matrisa" << endl
			 << "  b er talet i mellombandet i matrisa" << endl
			 << "  c er talet i det øvre bandet i matrisa" << endl
			 << "  LU styrer om det minnekrevjande LU-steget vert køyrt (0: nei, 1: ja)" << endl << endl
			 << "Resultata vert skrivne til filene v-gn_n.csv (generell algoritme)," << endl
			 << "v-sp_n.csv (spesialisert algoritme), v-lu_n.csv (LU) og u_n.csv (eksakt løysing)." << endl;
		return argc == 1 ? EXIT_SUCCESS : EXIT_FAILURE;
	}

	// Slå av vitskapleg talformat på console out
	cout << fixed;

	// Les kommandolinja
	unsigned int n, lu = 1;
	vec a, b, c;
	try {
		int i = stoi(argv[1]);
		if (i < 1) {
			throw invalid_argument("n må vera 1 eller høgre");
		}
		n = (unsigned int)pow(10, i);

		a = vec(n - 1);
		a.fill(stoi(argv[2]));

		b = vec(n);
		b.fill(stoi(argv[3]));

		c = vec(n - 1);
		c.fill(stoi(argv[4]));

		if (argc > 5) {
			lu = stoi(argv[5]);
		}
	} catch (exception &e) {
		cerr << "Greidde ikkje å lesa argument: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	// Set opp vektorar til diskretiseringa av x og f(x)
	double h = 1 / (n + 1.f);
	vec xi = vec(n);
	vec y = vec(n);

	for (unsigned int i = 0; i < n; i++) {
		xi(i) = (i + 1) * h;
		y(i)  = h*h * f(xi(i));
	}

	// Vektorar til løysingane
	vec v_gn, v_sp, v_lu, ui;

	// Set opp Tridiagonal-klassa og løys likningane
	try {
		Tridiagonal tridiag;
		tridiag.init(n, a, b, c);

		v_gn = tridiag.solve_gn(y);
		v_sp = tridiag.solve_sp(y);
		if (lu) {
			v_lu = tridiag.solve_lu(y);
		}
	} catch (exception &e) {
		cerr << "Eitkvart gjekk ikkje rett føre seg: " << e.what() << endl;
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

	// Jamfør løysingane
	double eps_gn, eps_sp, eps_lu;
	eps_gn = abs((v_gn - ui) / ui).max();
	eps_sp = abs((v_sp - ui) / ui).max();
	if (lu) {
		eps_lu = abs((v_lu - ui) / ui).max();
	}

	cout << "Fråvik (relativ feil):" << endl;
	cout << "  Generell algoritme:     " << log10(eps_gn) << endl;
	cout << "  Spesialisert algoritme: " << log10(eps_sp) << endl;
	if (lu) {
		cout << "  Ved LU-faktorisering:   " << log10(eps_lu) << endl;
	}

	// Skriv filer
	cout << "Skriv løysingane til fil ..." << endl;
	ostringstream os;
	os << "v-gn_" << n << ".csv";
	v_gn.save(os.str(), csv_ascii);

	os.str("");
	os << "v-sp_" << n << ".csv";
	v_sp.save(os.str(), csv_ascii);

	if (lu) {
		os.str("");
		os << "v-lu_" << n << ".csv";
		v_lu.save(os.str(), csv_ascii);
	}

	os.str("");
	os << "u_" << n << ".csv";
	ui.save(os.str(), csv_ascii);

	return EXIT_SUCCESS;
}

inline double f(double x) {
	return 100*exp(-10*x);
}

inline double u(double x) {
	return 1 - (1 - exp(-10))*x - exp(-10*x);
}

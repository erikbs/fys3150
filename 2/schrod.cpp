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
	unsigned int N;
	double rho_min = 0, rho_max = 1;

	if (argc < 2) {
		((argc == 1) ? cout : cerr)
			<< "Bruk: " << endl
			<< "  " << argv[0] << " N [ρ]" << endl << endl
			<< "der" << endl
			<< " N er talet på integrasjonspunkt" << endl
			<< " ρ er det høgste verdet til potensialet ρ (førevalt verd: 1)" << endl << endl
			<< "Resultata vert skrivne til fila lambda_N_ρ.csv" << endl;
		return argc == 1 ? EXIT_SUCCESS : EXIT_FAILURE;
	}

	try {
		N = stoi(argv[1]);
		if (N < 3) {
			throw invalid_argument("N lyt vera minst 3");
		}

		if (argc > 2) {
			rho_max = stof(argv[2]);
			if (rho_max <= 0) {
				throw invalid_argument("ρ lyt vera større enn 0");
			}
		}
	} catch (exception &e) {
		cerr << "Greidde ikkje å lesa argument: " << e.what() << endl;
		return EXIT_FAILURE;
	}

	double h =  (rho_max - rho_min) / N;
	double d =  2. / (h * h);
	double a = -1. / (h * h);

	vec v_a(N - 2);
	vec v_d(N - 1);

	v_a.fill(a);
	v_d.fill(d);

	// Legg til potensialet på diagonalen
	for (unsigned int i = 1; i < N; i++) {
		double rho_i = rho_min + i * h;
		v_d(i - 1) += rho_i * rho_i;
	}

	Jacobi J(v_a, v_d);

	clock_t tick, tock;
	tick = clock();
	J.run();
	tock = clock();
	vec eigval = J.values(true);

	// Skriv (inntil) dei fire fyrste tala
	cout << "λ = [";
	for (int i = 0; i < 4 && i < eigval.size(); i++) {
		cout << eigval(i);
		if (i < eigval.size() - 1) {
			cout << ", ";
		}
	}
	cout << ((eigval.size() > 4) ? "…]" : "]") << endl << endl;

	cout << "Utrekningstid: " << (double)(tock - tick) / CLOCKS_PER_SEC << " sekund" << endl;

	// Skriv til fil
	vec l = eigval.head(eigval.size() > 4 ? 4 : eigval.size());
	ostringstream os;
	os << "lambda_" << N << "_" << fixed << rho_max << ".csv";
	l.save(os.str(), csv_ascii);

	return EXIT_SUCCESS;
}

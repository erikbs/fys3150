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

	vec omega {0.01, 0.5, 1., 5.};
	mat eigmat = mat(N - 1, omega.size());
	for (unsigned int i = 0; i < omega.size(); i++) {
		double o = omega(i);

		// Legg til potensialet på diagonalen
		v_d.fill(d);
		for (unsigned int j = 1; j < N; j++) {
			double rho_j = rho_min + j * h;
			v_d(j - 1) += o * o * rho_j * rho_j + 1./rho_j;
		}

		Jacobi J(v_a, v_d);
		J.run();

		// Henta ut fyrste eigenvektor og spar honom
		eigmat.col(i) = J.vectors(true).col(0);
	}

	// Skriv tala til fil
	ostringstream os;
	os << "gstate_" << N << "_" << fixed << rho_max << ".csv";
	eigmat.save(os.str(), csv_ascii);

	return EXIT_SUCCESS;
}

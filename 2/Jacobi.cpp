#include <cmath>

#include "Jacobi.hpp"

using namespace std;
using namespace arma;

Jacobi::Jacobi(mat &A) : A(A) {
	if (!A.is_square()) {
		throw invalid_argument("matrisa lyt vera kvadratisk");
	}

	n = A.n_rows;
}

Jacobi::Jacobi(double a, double d, unsigned int n) {
	if (n < 2) {
		throw invalid_argument("n lyt vera minst 2");
	}

	this->n = n;
	A = mat(n, n, fill::zeros);
	A.diag() += d;
	A.diag(1) += a;
	A.diag(-1) += a;
}

Jacobi::Jacobi(vec a, vec d) {
	if (a.size() != d.size() - 1) {
		throw invalid_argument("d lyt ha eitt og berre eitt tal meir enn a");
	} else if (a.size() < 2) {
		throw invalid_argument("n lyt vera minst 2");
	}

	n = d.size();
	A = mat(n, n, fill::zeros);
	A.diag() = d;
	A.diag(1) = a;
	A.diag(-1) = a;
}

void Jacobi::initmax() {
	mi = Col<unsigned int>(n, fill::zeros);
	mx = vec(n, fill::zeros);

	for (unsigned int row = 0; row < n; row++) {
		double m = 0;
		for (unsigned int col = 0; col < n; col++) {
			if (abs(A(row, col)) > m && row != col) {
				m = abs(A(row, col));
				mi(row) = col;
				mx(row) = m;
			}
		}
	}

	// Hindra at diagonalelementet i fyrste rad er valt til max om rada elles er tom
	if (mi(0) == 0) {
		mi(0) = 1;
		mx(0) = abs(A(0, 1));
	}
}

// Finn det ikkje-diagonale elementet med høgst talverd
double Jacobi::max(unsigned int &k, unsigned int &l) {
	double m = 0;

	for (unsigned int row = 0; row < n; row++) {
		if (abs(A(row, mi(row))) > m) {
			k = row;
			l = mi(k);
			m = mx(k);
		}
	}

	return m;
}

void Jacobi::newmax(unsigned int k, unsigned int l) {
	double m;

	for (unsigned int row = 0; row < n; row++) {
		if (row != k && row != l) {
			// Ikkje rad k eller l, kan henda slepp me full gjennomgang
			if (mi(row) != k && mi(row) != l) {
				// Korkje kolonne k eller l hadde det største talet sist, sjå om dei har fått det no
				m = mx(row);
				if (abs(A(row, k)) > m) {
					mi(row) = k;
					mx(row) = abs(A(row, k));
				}
				if (abs(A(row, l)) > mx(row)) {
					mi(row) = l;
					mx(row) = abs(A(row, l));
				}

				// Veit kvar det største talet er, haldt fram med neste rad
				continue;
			}

			// Kolonne k eller l hadde det største talet sist, haldt fram om det enno er slik
			if (abs(A(row, mi(row))) == mx(row)) {
				continue;
			}
		}

		// Sjå gjennom heile rada å nyo (gjeld m.a. rad k og l)
		m = 0;
		for (unsigned int col = 0; col < n; col++) {
			if (abs(A(row, col)) > m && row != col) {
				m = abs(A(row, col));
				mi(row) = col;
				mx(row) = m;
			}
		}

		// Vel fyrste kolonne om heile rada er null (fyrste rad: andre kolonne i staden)
		if (m == 0) {
			if (row) {
				 mi(row) = mx(row) = 0;
			} else {
				mi(row) = 1;
				mx(row) = abs(A(row, 1));
			}
		}
	}
}

bool Jacobi::run() {
	// Utrekninga kan berre køyrast ein gong for kvar instans
	if (iter) return false;

	unsigned int k, l;
	R = mat(n, n, fill::eye);

	// Set opp max-vektoren
	initmax();

	// Køyr gjennom til alle ikkje-diagonale element er ovsmåe (< 10⁻⁸) eller bryt av etter n³ gonger
	for (iter = 0; iter < n*n*n && max(k, l) > eps; iter++) {
		rotate(k, l);
		newmax(k, l);
	}

	return true;
}


void Jacobi::rotate(unsigned int k, unsigned int l) {
	double s, c, t, tau;
	double A_kk, A_ll, A_ik, A_il, R_ik, R_il;

	// Finn cos(θ) og sin(θ)
	if (A(k, l) != 0) {
		tau = ( A(l, l) - A(k, k) ) / (2. * A(k, l));

		// Vel den minste rota
		if (tau > 0) {
			t = 1. / (tau + sqrt(1 + tau*tau)); // == -τ + sqrt(1 + τ²)
		} else {
			t = -1. / (-tau + sqrt(1 + tau*tau)); // == -τ - sqrt(1 + τ²)
		}

		c = 1. / sqrt(1 + t*t);
		s = c * t;
	} else {
		// Nemnaren i τ vert null, bruk kjend grenseverd i staden
		c = 1;
		s = 0;
	}

	/* Rekna ut sjølve rotasjonen */

	// Rekna ut diagonalelementa A(k, k) og A(l, l)
	A_kk = A(k, k);
	A_ll = A(l, l);
	A(k, k) = c*c * A_kk - 2 * c * s * A(k, l) + s*s * A_ll;
	A(l, l) = s*s * A_kk + 2 * c * s * A(k, l) + c*c * A_ll;

	// Set A(k, l) og A(l, k) til null for å sleppa avrundingsfråvik og reknekostnader
	A(k, l) = A(l, k) = 0;

	// Skjeringane mellom rad l og k og kolonne k og l er rekna ut, ta resten av rad/kolonne k og l
	for (unsigned int i = 0; i < n; i++) {
		if (i != k && i != l) {
			A_ik = A(i, k);
			A_il = A(i, l);
			A(i, k) = A(k, i) = c * A_ik - s * A_il;
			A(i, l) = A(l, i) = c * A_il + s * A_ik;
		}

		// Rekna ut nye eigenvektorar
		R_ik = R(i, k);
		R_il = R(i, l);
		R(i, k) = c * R_ik - s * R_il;
		R(i, l) = c * R_il + s * R_ik;
	}
}

mat Jacobi::input() {
	if (iter) {
		throw logic_error("Bad om matrisa A etter at run() er kalla");
	}

	return A;
}

vec Jacobi::values(bool sorted) {
	if (!iter) {
		throw logic_error("Bad om eigenverd før run() er kalla");
	}

	return sorted ? static_cast<vec>(sort(A.diag())) : A.diag();
}

mat Jacobi::vectors(bool sorted) {
	if (!iter) {
		throw logic_error("Bad om eigenvektorar før run() er kalla");
	}

	return sorted ? R.cols(sort_index(values())) : R;
}

unsigned int Jacobi::iterations() {
	if (!iter) {
		throw logic_error("Bad om talet på iterasjonar før run() er kalla");
	}

	return iter;
}

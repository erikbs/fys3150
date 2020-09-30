#include <cmath>
#include <iostream>

#include <armadillo>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#define private public
#include "Jacobi.hpp"

using namespace std;
using namespace arma;

TEST_CASE ("Finn .max() høgste ikkje-diagonale element?") {
	SECTION("Diagonalisert matrise") {
		vec d = {5, 5, 5, 5, 5};
		vec a = {3, 2, 1, 0};

		Jacobi J(a, d);
		J.initmax();
		unsigned int k, l;
		double m = J.max(k, l);

		REQUIRE(m == 3);
		REQUIRE(k == 0);
		REQUIRE(l == 1);
	}

	SECTION("Sjølvvald matrise") {
		mat A = mat(5, 5, fill::zeros);
		A(4, 4) = 10;
		A(1, 3) = -4;

		Jacobi J(A);
		J.initmax();
		unsigned int k, l;
		double m = J.max(k, l);

		REQUIRE(m == 4);
		REQUIRE(k == 1);
		REQUIRE(l == 3);
	}
}

TEST_CASE("Er eigenverda rette?") {
	unsigned int N = 7;
	int a = -1, d = 2;

	Jacobi J(a, d, N);
	J.run();
	vec eigval = J.values(true);

	vec exact = vec(N - 1);
	for (unsigned int j = 1; j < N; j++) {
		exact(j - 1) = d + 2 * a * cos(j * M_PI / N);
	}

	// Krev at differansen er innanfor ein toleranse på om lag 10⁻⁸ i kvart element
	REQUIRE(sum(abs(eigval - exact)) / (N - 1) < 10e-8);
}

TEST_CASE("Er ortogonaliteten halden ved lag?") {
	unsigned int N = 10;
	int a = -1, d = 2;

	Jacobi J(a, d, N);
	J.run();
	mat eigvec = J.vectors();

	// Sidan R er unitær, skal R*R' = I
	REQUIRE(accu(abs(eigvec.t() * eigvec - eye(N-1, N-1))) / ((N-1) * (N-1)) < 10e-8);
}

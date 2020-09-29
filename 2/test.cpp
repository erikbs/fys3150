#include <iostream>

#include <armadillo>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#define private public
#include "Jacobi.hpp"

using namespace std;
using namespace arma;

TEST_CASE ("Max", "[max]") {
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

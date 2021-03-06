#include "Tridiagonal.hpp"
#include <iostream>
#include <cmath>

using namespace std;
using namespace arma;

void Tridiagonal::init(unsigned int n, vec a, vec b, vec c) {
	#ifdef DEBUG
	cout << "DEBUG: init() kalla med argument:" << endl;
	cout << "n=" << n << endl;
	a.print("a=");
	b.print("b=");
	c.print("c=");
	#endif

	if (a.n_elem != n - 1 || b.n_elem != n || c.n_elem != n - 1) {
		throw invalid_argument("Ugild lengd på minst éin inngangsvektor");
	}

	this->n = n;
	this->a = a;
	this->b = b;
	this->c = c;
}

// Generell utgåve av løysingsalgoritmen
vec Tridiagonal::solve_gn(vec y) {
	#ifdef DEBUG
	cout << "DEBUG: solve_gn() kalla med argument:" << endl;
	y.print("y=");
	#endif

	if (y.n_elem != n) {
		throw invalid_argument("Ugild lengd på inngangsvektoren");
	}

	// Set opp vektorane b̃, ỹ og v
	vec bn = vec(n);
	vec yn = vec(n);
	vec v  = vec(n);

	// Set opp klokka
	clock_t tick, tock;
	tick = clock();

	// Framlenges-steget
	bn(0) = b(0);
	yn(0) = y(0);
	for (unsigned int i = 1; i < n; i++) {
		bn(i) = b(i) - (a(i - 1) * c(i - 1)) / bn(i - 1);
		yn(i) = y(i) - (a(i - 1) * yn(i - 1)) / bn(i - 1);
	}

	// Baklenges-steget
	v(n - 1) = yn(n - 1) / bn(n - 1);
	for (unsigned int i = n - 1; i; i--) {
		v(i - 1) = (yn(i - 1) - c(i - 1) * v(i)) / bn(i - 1);
	}

	// Meld om tida
	tock = clock();
	cout << "solve_gn() løyste likninga på " << (double)(tock - tick) / CLOCKS_PER_SEC << " sekund" << endl;

	#ifdef DEBUG
	cout << "DEBUG: solve_gn() gav løysinga:" << endl;
	v.print("v=");
	#endif

	return v;
}

// Spesialisert utgåve av løysingsalgoritmen
vec Tridiagonal::solve_sp(vec y) {
	#ifdef DEBUG
	cout << "DEBUG: solve_sp() kalla med argument:" << endl;
	y.print("y=");
	#endif

	if (y.n_elem != n) {
		throw invalid_argument("Ugild lengd på inngangsvektoren");
	}

	// Set opp vektorane b̃, ỹ og v
	vec bn = vec(n);
	vec yn = vec(n);
	vec v  = vec(n);

	// Set opp klokka
	clock_t tick, tock;
	tick = clock();

	// Framlenges-steget
	bn(0) = 2;
	yn(0) = y(0);
	for (unsigned int i = 1; i < n; i++) {
		bn(i) = (i + 2.f) / (i + 1);
		yn(i) = y(i) + yn(i - 1) / bn(i - 1);
	}

	// Baklenges-steget
	v(n - 1) = yn(n - 1) / bn(n - 1);
	for (unsigned int i = n - 1; i; i--) {
		v(i - 1) = (yn(i - 1) + v(i)) / bn(i - 1);
	}

	// Meld om tida
	tock = clock();
	cout << "solve_sp() løyste likninga på " << (double)(tock - tick) / CLOCKS_PER_SEC << " sekund" << endl;

	#ifdef DEBUG
	cout << "DEBUG: solve_sp() gav løysinga:" << endl;
	v.print("v=");
	#endif

	return v;
}

// Løys med LU-faktorisering
vec Tridiagonal::solve_lu(vec y) {
	#ifdef DEBUG
	cout << "DEBUG: solve_lu() kalla med argument:" << endl;
	y.print("y=");
	#endif

	if (y.n_elem != n) {
		throw invalid_argument("Ugild lengd på inngangsvektoren");
	}

	// Set opp A
	mat A = mat(n, n);
	A.zeros();
	A.diag()   = b;
	A.diag(-1) = a;
	A.diag(1)  = c;

	// Set opp klokka
	clock_t tick, tock;
	tick = clock();

	// Gjer LU-faktorisering (med P for å kunna nytta trimatu og trimatl, jf. Armadillo-dokumentasjonen)
	mat L, U, P;
	lu(L, U, P, A);
	vec v = solve(trimatu(U), solve(trimatl(L), P*y));

	// Meld om tida
	tock = clock();
	cout << "solve_lu() løyste likninga på " << (double)(tock - tick) / CLOCKS_PER_SEC << " sekund" << endl;

	#ifdef DEBUG
	cout << "DEBUG: solve_lu() gav løysinga:" << endl;
	v.print("v=");
	#endif

	return v;
}

void Tridiagonal::print() {
	#ifdef DEBUG
	cout << "DEBUG: print() kalla" << endl;
	#endif
}

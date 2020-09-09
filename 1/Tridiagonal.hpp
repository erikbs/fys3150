#ifndef TRIDIAGONAL_HPP
#define TRIDIAGONAL_HPP

#include <armadillo>

class Tridiagonal {
private:
	unsigned int n;
	arma::vec a, b, c;

public:
	void init(unsigned int n, arma::vec a, arma::vec b, arma::vec c);
	arma::vec solve_gen(arma::vec y);
	arma::vec solve_spec(arma::vec y);
	arma::vec solve_lu(arma::vec y);
	void print();
};

#endif

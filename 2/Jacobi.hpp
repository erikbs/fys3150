#ifndef JACOBI_HPP
#define JACOBI_HPP

#include <armadillo>

class Jacobi {
private:
	double eps = 1.0e-8;
	unsigned int n;
	arma::mat A;
	arma::mat R;
	arma::Col<unsigned int> mi;
	arma::vec mx;

	unsigned int iter = 0;

	void initmax();
	double max(unsigned int &k, unsigned int &l);
	void newmax(unsigned int k, unsigned int l);

	void rotate(unsigned int k, unsigned int l);
public:
	Jacobi(arma::mat &A);
	Jacobi(double a, double d, unsigned int n);
	Jacobi(arma::vec a, arma::vec d);

	/**
	 * Gjer utrekninga av eigenverd og -vektorar etter Jacobi-metoden
	 *
	 * @returns true om alt gjekk vel, false om utrekninga alt er køyrd
	 */
	bool run();

	/**
	 * Få ein kopi av inngangsmatrisa (berre før utrekninga er køyrd)
	 *
	 * @returns kopi av inngangsmatrisa
	 * @throws logic_error om run() er kalla og inngangsmatrisa med di ikkje finst meir
	 */
	arma::mat input();

	/**
	 * Henta ut dei utrekna eigenverda i utrekna rekkjefylgje (berre etter utrekninga)
	 *
	 * @param sorted	tala vert sorterte om denne er sann
	 * @returns vektor med eigenverda
	 * @throws logic_error om run() ikkje er kalla fyrst
	 */
	arma::vec values(bool sorted = false);

	/**
	 * Henta ut dei utrekna eigenvektorane i utrekna rekkjefylgje (berre etter utrekninga)
	 *
	 * @param sorted	kolonnene vert sorterte på eigenverd om denne er sann
	 * @returns matrise med eigenvektorane til kolonner
	 * @throws logic_error om run() ikkje er kalla fyrst
	 */
	arma::mat vectors(bool sorted = false);

	/**
	 * Få talet på iterasjonar i utrekninga (berre etter utrekninga)
	 *
	 * @returns talet på iterasjonar
	 * @throws logic_error om run() ikkje er kalla fyrst
	 */
	unsigned int iterations();
};

#endif

/*
 * Plasma2D.cpp
 *
 *  Created on: 23.11.2010
 *      Author: Martin Rueckl
 */

#include "Plasma2D.h"

#include <RandomNumberGenerator.h>
#include <Matrix2D.h>
#include <Matrix3D.h>
#include <MyException.h>
#include <vector>
//#include <boost/threadpool/pool.hpp>

#define EULER 0.57721566490153286060651209008240243104215933593992           	/*Eulers constant */

Plasma2D::Plasma2D(vector<Charge> _charges, double _base_length) :
	charges(_charges),
	charge_count(charges.size()),
	total_charge(0),
	base_width(_base_length),
	base_height(_base_length * sqrt(3.)),
	kSpace_base_width(2.* M_PI / base_width),
	kSpace_base_height(2. * M_PI / base_height),
	kSpace_cutoff(30), //es sollte gelten (k/k_min)^2<200 mit k_min=2pi/(sqrt(3)b)=kSpace_base_height damit folgt kSpace_cutoff=800pi^2/(3base_width^2)
	//für bisher benutzte parameter (base_width=10) ist die konvergenzbedinungn sogar weit besser erfüllt als nötig
	eta_squared(36. / (base_width * base_height)),
	realSpaceEnergyMatrix(charge_count, charge_count),
	kSpaceEnergyMatrix_real(kSpace_cutoff, kSpace_cutoff, charge_count),
	kSpaceEnergyMatrix_imag(kSpace_cutoff, kSpace_cutoff, charge_count),
	distance_matrix(charge_count, charge_count),
	rng()//seeds with time()
{
	for (unsigned u = 0; u < charge_count; u++) {
		total_charge += charges[u].charge;
		calculateParticleEnergy(u);
		updateDistanceMatrix(u);
	}
	energy = calculateEnergy();
}

Plasma2D::Plasma2D(unsigned _charge_count, double _base_length, double charge) :
	charge_count(_charge_count),
	total_charge(0),
	base_width(_base_length),
	base_height(_base_length * sqrt(3.)),
	kSpace_base_width(2. * M_PI / base_width),
	kSpace_base_height(2. * M_PI / base_height),
	kSpace_cutoff(30), //es sollte gelten (k/k_min)^2<200 mit k_min=2pi/(sqrt(3)b)=kSpace_base_height damit folgt kSpace_cutoff=800pi^2/(3base_width^2)
	//für bisher benutzte parameter (base_width=10) ist die konvergenzbedinungn sogar weit besser erfüllt als nötig
	eta_squared(36. / (base_width * base_height)),
	realSpaceEnergyMatrix(charge_count, charge_count),
	kSpaceEnergyMatrix_real(kSpace_cutoff, kSpace_cutoff, charge_count),
	kSpaceEnergyMatrix_imag(kSpace_cutoff, kSpace_cutoff, charge_count),
	distance_matrix(charge_count, charge_count),
	rng()//seeds with time()
{
	charges.reserve(charge_count);
	for (unsigned u = 0; u < charge_count; u++) {
		double x = (rng.getUniform() - 0.5) * base_width;
		double y = (rng.getUniform() - 0.5) * base_height;
		Charge c(x, y, charge);
		total_charge += c.charge;
		charges.push_back(c);
	}

	for (unsigned u = 0; u < charge_count; u++)
	{
		calculateParticleEnergy(u);
		updateDistanceMatrix(u);
	}

	energy = calculateEnergy();
}

void Plasma2D::modify()
{
	modifiedParticle = (unsigned)(rng.getUniform()*charge_count);
	oldLocation		 = charges[modifiedParticle].location;
	oldEnergy		 = energy;
	do{
		Point p = oldLocation + rng.getNormalDistributedPoint(base_width/10.,0) ;
//		p.x += (2 * rng.getUniform() - 1)*(base_width/2.);
//		p.y += (2 * rng.getUniform() - 1)*(base_heigth/2.);

		//keep all the particles in the box
		while (p.x > base_width / 2)	p.x = p.x - base_width;
		while (p.x < -base_width / 2)	p.x = p.x + base_width;
		while (p.y > base_height / 2)	p.y = p.y - base_height;
		while (p.y < -base_height / 2)	p.y = p.y + base_height;

		charges[modifiedParticle].location = p;
		updateDistanceMatrix(modifiedParticle);
	}while(collision(modifiedParticle,0.47));

	calculateParticleEnergy(modifiedParticle);
	energy = calculateEnergy();
}

void Plasma2D::undo()
{
	charges[modifiedParticle].location = oldLocation;
	calculateParticleEnergy(modifiedParticle);
	updateDistanceMatrix(modifiedParticle);
	energy = oldEnergy;//we dont need to sum over matrices because we stored old energy...
}

double Plasma2D::getEnergy() const
{
	return energy;
}

double Plasma2D::calculateEnergy()
{
	double realSpace_Energy = sumRealSpaceMatrix();
	double kSpace_Energy = sumKSpaceMatrices();
	/*
	 * This part is not needed for montecarlo simulations because it only adds a constant and therefore has no influence on energy differences between states.
	 *
	 * double background_Energy = total_charge * total_charge / 4 * (M_PI / (eta_squared * basis_length * basis_length) + (EULER + log(
	 * eta_quared / (M_PI * total_charge / (basis_length * basis_length)))) / charge_count);
	 */
	return 0.25 * realSpace_Energy + M_PI / (base_width * base_height) * kSpace_Energy;

}

double Plasma2D::sumRealSpaceMatrix() const
{
	double t = 0;
	for (unsigned i = 0; i < charge_count; i++)
		for (unsigned j = i; j < charge_count; j++)
			t += realSpaceEnergyMatrix(i, j);
	return t;
}

double Plasma2D::sumKSpaceMatrices() const
{
	double t = 0;
	for (int kx = 0; kx < kSpace_cutoff; kx++) {
		for (int ky = 0; ky < kSpace_cutoff; ky++) {
			double real = 0;
			double imag = 0;
			for (unsigned i = 0; i < charge_count; i++) {
				real += kSpaceEnergyMatrix_real(kx, ky, i);
				imag += kSpaceEnergyMatrix_imag(kx, ky, i);
			}
			t += real * real + imag * imag;
		}
	}
	return t;
}

void Plasma2D::calculateParticleEnergy(unsigned i)
{
	/* RealSpace :
	 * calculate only changed summands of double sum (i,j) and store them in energy matrix */

	for (unsigned j = 0; j <  charge_count ; j++)
	{
		if (i == j) continue;
		double qi = charges[i].charge;
		double qj = charges[j].charge;


		//since all the particles are forced to stay in the box by the cyclic boundary conditions from modify()
		//checking for eventually smaller distances to mirrored charges in the surounding boxes does not reproduce the formula provided
		//in the karch paper since they skipp the real space sum after n=(0,0)
		double d = (charges[i].location - charges[j].location).absPow2();
		double dl = (charges[i].location + Point(-base_width, 0) - charges[j].location).absPow2();
		double dr = (charges[i].location + Point(base_width, 0) - charges[j].location).absPow2();
		double dt = (charges[i].location + Point(0, base_height) - charges[j].location).absPow2();
		double db = (charges[i].location + Point(0, -base_height) - charges[j].location).absPow2();
		double dbr = (charges[i].location + Point(base_width, -base_height) - charges[j].location).absPow2();
		double dbl = (charges[i].location + Point(-base_width, -base_height) - charges[j].location).absPow2();
		double dtr = (charges[i].location + Point(base_width, base_height) - charges[j].location).absPow2();
		double dtl = (charges[i].location + Point(-base_width, base_height) - charges[j].location).absPow2();

		double distance_squared = min(d, min(dl, min(dr, min(dt, min(db, min(dbr, min(dbl, min(dtr, dtl))))))));
		double energy_ij = qi * qj * expint1_NRIC(eta_squared * distance_squared);
		realSpaceEnergyMatrix(i, j) = energy_ij;
		realSpaceEnergyMatrix(j, i) = energy_ij;
	}
	//cout << realSpaceEnergyMatrix << endl;

	/* kSpace :
	 * calculate only changed summands of triple sum (kx,ky,i) and store them in energy matrix */
	double cutoff_half = kSpace_cutoff / 2;
	for (double kx = -kSpace_cutoff / 2; kx < kSpace_cutoff / 2; kx++) {
		for (double ky = -kSpace_cutoff / 2; ky < kSpace_cutoff / 2; ky++) {
			bool cutoff = kx * kx + ky * ky > cutoff_half * cutoff_half;
			bool center = (kx == 0) && (ky == 0); //sollte hier nicht (kx == 0) && (ky == 0) stehen??
			if (cutoff | center) continue;
			Point k = Point(kx * kSpace_base_width, ky * kSpace_base_height);
			double preFaktor = exp(-k.absPow2() / (8 * eta_squared)) / k.abs();
			Point p = charges[i].location;
			kSpaceEnergyMatrix_real(kx + cutoff_half, ky + cutoff_half, i) = preFaktor * charges[i].charge * cos(k * p);
			kSpaceEnergyMatrix_imag(kx + cutoff_half, ky + cutoff_half, i) = preFaktor * charges[i].charge * sin(k * p);
		}
	}
}

bool Plasma2D::collision(unsigned i, double min_distance)
{
	for(unsigned j=0;j<charge_count;j++)
		if(distance_matrix(i,j)<=min_distance && i!=j){
			//cout<<"collision: i1="<<i<<" i2="<<j<<" distance="<< distance_matrix(i,j)<<endl;
			return true;
		}
	return false;
}

void Plasma2D::updateDistanceMatrix(unsigned j)
{
	for(unsigned i=0;i<charge_count;i++){
		double d = (charges[i].location - charges[j].location).abs();
		double dl = (charges[i].location + Point(-base_width, 0) - charges[j].location).abs();
		double dr = (charges[i].location + Point(base_width, 0) - charges[j].location).abs();
		double dt = (charges[i].location + Point(0, base_height) - charges[j].location).abs();
		double db = (charges[i].location + Point(0, -base_height) - charges[j].location).abs();
		distance_matrix(i, j) = min(d, min(dl, min(dr, min(dt, db))));
		distance_matrix(j, i) = min(d, min(dl, min(dr, min(dt, db))));
	}
}

const Matrix2D& Plasma2D::getDistanceMatrix()const{
	return distance_matrix;
}

Matrix2D Plasma2D::getPairCorrelation(double stepsize, double min, double max) const
{
	Matrix2D m((max - min) / stepsize + 1, 3);
	for (double rstep = 0; rstep * stepsize + min < max; rstep++) {
		double inner = rstep * stepsize + min;
		double outer = inner + stepsize;
		m(rstep, 0) = inner;
		double x = 0;
		double x_squared = 0;
		for (unsigned i = 0; i < charge_count; i++) {
			double t = 0;
			for (unsigned j = 0; j < charge_count; j++) {
				double d = (charges[i].location - charges[j].location).abs();
				double dl = (charges[i].location + Point(-base_width, 0) - charges[j].location).abs();
				double dr = (charges[i].location + Point(base_width, 0) - charges[j].location).abs();
				double dt = (charges[i].location + Point(0, base_height) - charges[j].location).abs();
				double db = (charges[i].location + Point(0, -base_height) - charges[j].location).abs();
				double dbr = (charges[i].location + Point(base_width, -base_height) - charges[j].location).abs();
				double dbl = (charges[i].location + Point(-base_width, -base_height) - charges[j].location).abs();
				double dtr = (charges[i].location + Point(base_width, base_height) - charges[j].location).abs();
				double dtl = (charges[i].location + Point(-base_width, base_height) - charges[j].location).abs();

				if (inner < d && d <= outer) t++;
				if (inner < dl && dl <= outer) t++;
				if (inner < dr && dr <= outer) t++;
				if (inner < dt && dt <= outer) t++;
				if (inner < db && db <= outer) t++;
				if (inner < dtr && dtr <= outer) t++;
				if (inner < dtl && dtl <= outer) t++;
				if (inner < dbr && dbr <= outer) t++;
				if (inner < dbl && dbl <= outer) t++;
			}
			x += t;
			x_squared += t * t;
		}
		m(rstep, 1) = x;
		m(rstep, 1) = m(rstep, 1) / (charge_count * 2);
		m(rstep, 1) = m(rstep, 1) / (inner * M_PI * stepsize);
		m(rstep, 1) = m(rstep, 1) / (charge_count / (base_width * base_height));
		m(rstep, 2) = sqrt(1. / (charge_count - 1) * (x_squared - x * x / charge_count));//standard deviation
	}
	return m;
}

ostream& operator <<(ostream &os, const Plasma2D &obj)
{
	Matrix2D m(obj.charge_count,3);
	for (unsigned c = 0; c < obj.charge_count; c++){
		m(c,0)=obj.charges[c].location.x;
		m(c,1)=obj.charges[c].location.y;
		m(c,2)=obj.charges[c].charge;
	}
	os << m;
	return os;
}

/*
 * reimplementation of the Algorithm to calculate the
 * Exponential integral E_1(z) from the numerical recipes
 * This is no complete reimplementation since only the integrl
 * E1 is needed. It is reimplemented as double valued function
 * by T.Kampf
 */

#define FPMIN 1.0e-100                	/*Close to smallest representable floating-point number.*/
#define EPS 1.0e-7                   	/*Desired relative error, not smaller than the machine precision */
double Plasma2D::expint1_NRIC(double x) const
{
	if (x <= 0) return 0;
	static int max_it = 1000; /*Maximum allowed number of iterations.*/
	if (x > 1.0) { /* using the Lentz's algorithm cf. Num Recipes in C */
		double b = x + 1.0;
		double c = 1.0 / FPMIN;
		double d = 1.0 / b;
		double h = d;
		for (int i = 1; i <= max_it; i++) {
			double a = -i * (i);
			b += 2.0;
			d = 1.0 / (a * d + b); /* Denominators cannot be zero.*/
			c = b + a / c;
			double del = c * d;
			h *= del;
			if (fabs(del - 1.0) < EPS) return h * exp(-x);
		}
		throw MyException("continued fraction failed in expint1_NRIC! x = %lf\n");

	} else { /* Evaluate series. */
		double ans = (-log(x) - EULER); /* Set fïrst term.*/
		double fact = 1.0;

		for (int i = 1; i <= max_it; i++) {
			fact *= -x / i;
			double del = -fact / (i);
			ans += del;
			if (fabs(del) < fabs(ans) * EPS) return ans;
		}
		throw MyException("series failed in expint1_NRIC! x = %lf\n");
	}
}

Plasma2D::~Plasma2D()
{
	// TODO Auto-generated destructor stub
}

/*
 * Exec.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: gstu0908
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <sstream>

#include <Matrix2D.h>

#include <util/boundaries/OuterBoundary.h>
#include <util/boundaries/CircleBoundary.h>
#include <util/boundaries/RectangularBoundary.h>
#include <util/boundaries/MixedBoundary.h>
#include <Simulator.h>
#include <SimulationEvaluator.h>
#include <SpinEchoSequence.h>
//#include <PeriodicMSEData.h>

#include <util/MPICommunicator.h>
#include <util/Capillary2D.h>
#include <util/CapillaryConfiguration.h>
#include <util/Field2D.h>
#include <util/PrecalculatedField2D.h>
#include <util/LiveField2D.h>
#include <util.h>
#include <util/IOManager.h>
#include <RandomNumberGenerator.h>
#include <PerfClock.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/filesystem.hpp>
#include <time.h>

using namespace std;

void rndPerformance();
void christian_ortsaufgeloest();
void test_matrix_io();
void diffusionsstatistik();

int main(int argc, char* argv[])
{
	ifstream file("test.txt");


	cout<<c<<(c=='\n'?"newline":"not newline...")<<endl;
	file.seekg(0,ios::end);
	file.get(c);
	cout<<c<<(c=='\n'?"newline":"not newline...")<<endl;
	file.seekg(1,ios::end);
	file.get(c);
	cout<<c<<(c=='\n'?"newline":"not newline...")<<endl;
	return 1;
	cout<<atof("asd")<<endl;
	cout<<atof("0.asd")<<endl;
	cout<<atof("asd.123")<<endl;
	cout<<atof("asd123")<<endl;
	cout<<atof("123asd")<<endl;
	cout<<atof("  123  asd")<<endl;
	//rndPerformance();
	//christian_ortsaufgeloest();
	//diffusionsstatistik();
	//test_matrix_io();
}

void diffusionsstatistik()
{
//	double Ra=8.3E-6;
//	double D=1E-9;
//	unsigned samples = 750000;
//	vector<PrecalculatedField2D> statistics(100);
//	ProgressMonitor pmon;
//
//	for(double Ri=1E-6;Ri<Ra;Ri+=Ra/10)
//	{
//		Capillary2D cap(Point(0,0),Ri);
//		CapillaryConfiguration capillaries(cap);
//		LiveField2D field(capillaries);
//		CircleBoundary boundary(Ra);
//		Simulator<CircleBoundary,LiveField2D> sim(boundary,field,0.01,D,1000,1);
//		sim.setCapillaryConfiguration(capillaries);
//		unsigned percent = sim.getTimeSteps()/100;
//		stringstream path1;
//		path1<<"Ri="<<fixed<<setprecision(2)<<Ri*1E6;
//		if (!boost::filesystem::is_directory(path1.str()))
//			boost::filesystem::create_directory(path1.str());
//		for(double r=Ri;r<Ra;r+=Ra/20)
//		{
//			stringstream path2;
//			path2<<path1.str()<<"/r="<<fixed<<setprecision(3)<<r*1E6;
//			if (boost::filesystem::is_directory(path2.str()))
//				continue;
//			else
//				boost::filesystem::create_directory(path2.str());
//
//			ProgressBar *progress = new ProgressBar("Sampling: Ri="+util::num2str(Ri)+" r="+util::num2str(r),samples);
//			pmon.addProgressBar(progress);
//			for(unsigned u=0;u<100;u++)
//				statistics[u]=PrecalculatedField2D(Point(-Ra,-Ra),Point(Ra,Ra),0.1E-6);
//			#pragma omp parallel for shared(pmon, sim, progress, statistics)
//			for(unsigned u=0;u<samples;u++)
//			{
//				Point p(0,r);
//				for(unsigned step=0;step<sim.getTimeSteps();step++)
//				{
//					p=sim.getNextStepLocation(p);
//					if(step%percent==0) (statistics[step/percent])(p)++;
//				}
//				progress->progress();
//			}
//
//			for(unsigned u=0;u<100;u++)
//			{
//				stringstream filename;
//				filename<<path2.str()<<"/t="<<fixed<<setprecision(4)<<(u*sim.getTimeIntervall()/100.)<<".csv";
//				ofstream file(filename.str().c_str());
//				file << statistics[u];
//				file.close();
//			}
//			progress->finish();
//			//delete progress;
//		}
//	}

}

void test_matrix_io(){
//	Matrix2D m(20,20);
//	Matrix2D m2(20,20);
//	for(unsigned u=0;u<20;u++)
//		for(unsigned n=0;n<20;n++)
//			m(u,n)=sqrt(M_PI*u*n)*exp(20.);
//	fstream file("original.csv", ios::out);
//	file << m;
//	file.close();
	string filename="\\\\MAMBA\\q\\home\\mrueckl\\data_wuerzburg\\tests\\geometry_tests\\cap.csv";
	if (!boost::filesystem::exists(filename)) throw MyException("The File " + filename + " does not exist.");
	Matrix2D m(1,1);
	fstream file(filename.c_str(), ios::in);
	file >> m;
	file.close();
	cout << m << endl;
}


void christian_ortsaufgeloest(){
//	double Ra	=9.4884*1E-6;//\ follows eta=0.084
//	double Ri	=2.75*1E-6;	 ///
//	double D	=1*1E-9;
//	double T	=0.2;
//	double ts	=10000;
//	double Bz	=1.12888;//follows dw=151 rad/s
//	double gamma=2.67522E8;
//	double chi	=1E-6;
//	double alpha=M_PI / 2;
//
//	CapillaryConfiguration capillaries(Capillary2D(Point(0, 0), Ri));
//	Field2D *field = new LiveField2D(capillaries);
//	OuterBoundary *boundary = new MixedBoundary(Ra);
//	Simulator simulator(T, D, ts, Bz,alpha, chi, gamma);
//	simulator.setCapillaryConfiguration(capillaries);
//	simulator.setField(field);
//	simulator.setBoundary(boundary);
//
//	fstream sample_trajectory_file("sample_trajectories.csv", ios::out);
//	sample_trajectory_file << simulator.getSampleTrajectory(simulator.getTimeSteps(), 10);
//
//	SimulationEvaluator evaluator(simulator, 100);
//	evaluator.getMagnetizationDistribution(50000000, 50, 10, Ra);
}

void rndPerformance(){
	double sigma = 0.141421E-6;
	double  mean = 0.1;
	boost::normal_distribution<double> norm(mean, sigma);
	boost::mt19937 gen;
	boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > norm_distrib_rnd(gen, norm);

	boost::rand48 gen2;
	boost::variate_generator<boost::rand48&, boost::normal_distribution<double> > norm_distrib_rnd2(gen2, norm);
	double d;
	RandomNumberGenerator rnd3;
	double t=0;
	d = clock();
	for (double u = 0; u < 500000000; u++) {
		t += rnd3.getUniform();
	}
	d = clock() - d;
	cout << "Own RNG: " << 500000000./(d / CLOCKS_PER_SEC) << "randoms/sec\t" << t<<endl;

	d = clock();
	t=0;
	for (double u = 0; u < 500000000; u++) {
		t += gen2();
	}
	d = clock() - d;
	cout << "boost_rand48 RNG: " << 500000000./(d / CLOCKS_PER_SEC) << "randoms/sec\t" << t<<endl;

	d = clock();
	t=0;
	for (double u = 0; u < 500000000; u++) {
		t += gen();
	}
	d = clock() - d;
	cout << "boost_mt19937 RNG: " << 500000000./(d / CLOCKS_PER_SEC) << "randoms/sec\t" << t<<endl;

	d = clock();	t=0;
	for (double u = 0; u < 5000000; u++) {
		t += norm_distrib_rnd2();
	}
	d = clock() - d;
	cout << "boost_norm distribution with rand48: " << 5000000./(d / CLOCKS_PER_SEC)<< "randoms/sec\t" << t<<endl;

	d = clock();	t=0;
	for (double u = 0; u < 5000000; u++) {
		t += norm_distrib_rnd();
	}
	d = clock() - d;
	cout << "boost_norm distribution with mt19937: " << 5000000./(d / CLOCKS_PER_SEC) << "randoms/sec\t" << t<<endl;

	d = clock();	t=0;
	for (double u = 0; u < 2500000; u++) {
		double dr;
		double rnd1, rnd2;
		do {
			rnd1 = 2.0 * rnd3.getUniform() - 1;
			rnd2 = 2.0 * rnd3.getUniform() - 1;
			dr = rnd1 * rnd1 + rnd2 * rnd2;
		} while (dr >= 1.);
		t += sqrt(-2.0 * log(dr) / dr) * rnd1 * sigma + mean;
		t += sqrt(-2.0 * log(dr) / dr) * rnd2 * sigma + mean;
	}
	d = clock() - d;
	cout << "Own norm distribution: " << 2*2500000*(d / CLOCKS_PER_SEC)<< "randoms/sec\t" << t<<endl;
}


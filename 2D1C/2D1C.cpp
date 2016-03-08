//============================================================================
// Name        : 2D1C.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Plasma2D.h"
#include <util.h>
#include <fstream>
#include <RandomNumberGenerator.h>
#include <Progress.h>
#include <Matrix2D.h>
using namespace std;

//parameters changeable by command line
unsigned equilibrate = 50000;
double T = 2;
double T_reduce = 1;
double T_reduce_steps = 0;
unsigned intermediate_steps = 10000;
double width = 20;
string plasma_file = "";
string prefix = "plasma";

string flg_plasma_file = "-plasma";
string flg_width = "-width";
string flg_equilibrate = "-equilibrate";
string flg_temperature = "-T";
string flg_prefix = "-prefix";
string flg_temperature_reduce = "-T-reduce";
string flg_temperature_reduce_steps = "-reduce-steps";
string flg_temperature_intermediate = "-intermediate-steps";

const string flg_help1 = "--help", flg_help2 = "-h", flg_help3 = "--h", flg_help4 = "-help";

void printHelp();
void writePlasmaFile(string filename, Plasma2D &plasma);
void writeMatrixFile(string str, const Matrix2D &distanceMatrix);
bool collision(const Plasma2D &plasma,double min_distance);

int main(int argc, char* argv[]) {
	//iterate over params and search for field params
	for (int a = 0; a < argc; a++) {
		if ((flg_help1.compare(argv[a]) == 0) | (flg_help2.compare(argv[a]) == 0) | (flg_help3.compare(argv[a]) == 0) | (flg_help4.compare(
				argv[a]) == 0)) {
			printHelp();
			return 0;
		} else if (flg_plasma_file.compare(argv[a]) == 0 && a + 1 < argc) {
			plasma_file = argv[a + 1];
		} else if (flg_width.compare(argv[a]) == 0 && a + 1 < argc) {
			width = atof(argv[a + 1]);
		} else if (flg_equilibrate.compare(argv[a]) == 0 && a + 1 < argc) {
			equilibrate = atoi(argv[a + 1]);
		} else if (flg_temperature.compare(argv[a]) == 0 && a + 1 < argc) {
			T = atof(argv[a + 1]);
		} else if (flg_prefix.compare(argv[a]) == 0 && a + 1 < argc) {
			prefix = argv[a + 1];
		} else if (flg_temperature_reduce.compare(argv[a]) == 0 && a + 1 < argc) {
			T_reduce = atof(argv[a + 1]);
		} else if (flg_temperature_reduce_steps.compare(argv[a]) == 0 && a + 1 < argc) {
			T_reduce_steps = atoi(argv[a + 1]);
		} else if (flg_temperature_intermediate.compare(argv[a]) == 0 && a + 1 < argc) {
			intermediate_steps = atoi(argv[a + 1]);
		}
	}
	ifstream file(plasma_file.c_str());
	Matrix2D charge_matrix(1,1);
	file>>charge_matrix;
	file.close();

	vector<Charge> charges;
	for (unsigned u = 0; u < charge_matrix.lines; u++) {
		charges.push_back(Charge(charge_matrix(u,0),charge_matrix(u,1),charge_matrix(u,2)));
	}

	cout << "Prefix: " << prefix << endl;
	cout << "Plasma file: " << plasma_file << endl;
	cout << "Charges: " << charges.size() << endl;
	cout << "width = " << width << endl;
	cout << "equilibrate = " << equilibrate << endl;
	cout << "T = " << T << endl;
	cout << "T-reduce = " << T_reduce << endl;
	cout << "intermediate-steps = " << intermediate_steps << endl;
	cout << "reduce-steps = " << T_reduce_steps << endl;

	Plasma2D plasma(charges, width);
	double last_energy = plasma.getEnergy();
	double current_energy = 0;
	unsigned accepts = 0;
	writePlasmaFile(prefix + ".csv", plasma);

	ProgressMonitor pmon;
	ProgressBar *pbar=new ProgressBar("Metropolis", equilibrate + T_reduce_steps * intermediate_steps);
	pmon.addProgressBar(pbar);
	RandomNumberGenerator rng;

	Matrix2D energy_flow(equilibrate + T_reduce_steps * intermediate_steps,2);

	//equilibrate
	for (unsigned i = 0; i < equilibrate; i++) {
		plasma.modify();
		current_energy = plasma.getEnergy();
		//energy difference:
		//lower energy->negative
		//higher energy->positive
		double dE = current_energy - last_energy;
		//acceptance probability
		//dE > 0 -> exp(-dE/T) < 1
		//dE < 0 -> exp(-dE/T) > 1 -> set p to 1
		double p = min(1., exp(-dE / T));
		//get random value in [0,1]
		double rnd = rng.getUniform();

		if (rnd <= p) {	//accept
			last_energy = current_energy;
			accepts = 1;
		} else { 		//discard
			plasma.undo();
			current_energy = last_energy;
			accepts = 0;
		}
		pbar->progress();
		energy_flow(i,0)=current_energy;
		energy_flow(i,1)=accepts;
	}

	writeMatrixFile(prefix + "_distances.csv", plasma.getDistanceMatrix());
	writeMatrixFile(prefix + util::num2str(equilibrate) + "_pair_correlation.csv", plasma.getPairCorrelation(0.01, 0, 10));
	writePlasmaFile(prefix + util::num2str(equilibrate) + ".csv", plasma);

	//run a simulated annealing
	for (unsigned u = 0; u < T_reduce_steps; u++) {
		for (unsigned i = 0; i < intermediate_steps; i++) {
			plasma.modify();
			current_energy = plasma.getEnergy();
			double dE = current_energy - last_energy;
			double p = min(1., exp(-dE / T));
			double rnd = rng.getUniform();
			if (rnd <= p) {	//accept
				last_energy = current_energy;
				accepts = 1;
			} else { 		//discard
				plasma.undo();
				current_energy = last_energy;
				accepts = 0;
			}
			energy_flow(u*i+equilibrate,0)=current_energy;
			energy_flow(u*i+equilibrate,1)=accepts;

			pbar->progress();
		}

		T *= T_reduce;
//		writeMatrixFile(prefix + util::num2str((u + 1) * intermediate_steps + equilibrate) + "_distances.csv", plasma.getDistanceMatrix());
		writeMatrixFile(prefix + util::num2str((u + 1) * intermediate_steps + equilibrate) + "_pair_correlation.csv", plasma.getPairCorrelation(0.01, 0, 10));
		writePlasmaFile(prefix + util::num2str((u + 1) * intermediate_steps + equilibrate) + ".csv", plasma);

	}
	ofstream energy_file((prefix + "_energy.csv").c_str());
	energy_file << energy_flow << endl;
	energy_file.close();
	return 0;
}

bool collision(const Plasma2D &plasma, double min_distance)
{
	const Matrix2D &distances = plasma.getDistanceMatrix();
	for(unsigned i=0;i<plasma.getChargeCount();i++)
		for(unsigned j=i+1;j<plasma.getChargeCount();j++)
			if(distances(i,j)<=min_distance){
				cout<<"collision: i1="<<i<<" i2="<<j<<" distance="<< distances(i,j)<<endl;
				return true;
			}
	return false;
}

void writePlasmaFile(string filename, Plasma2D &plasma) {
	ofstream plasma_file(filename.c_str());
	plasma_file << plasma;
	plasma_file.close();
}

void writeMatrixFile(string str, const Matrix2D &distanceMatrix) {
	ofstream file(str.c_str());
	file << distanceMatrix;
	file.close();
}

void printHelp() {
	cout << "============================================================================" << endl;
	cout << "Tool for simulating two dimensional plasma with log-potential via metropolis" << endl;
	cout << "algorithm. Implemented according to Physica A 369 (2006) 599–611: \"A Gibbs " << endl;
	cout << "point field model for the spatial pattern of coronary capillaries.\"" << endl;
	cout << "by Martin Rueckl" << endl;
	cout << "========================================================================" << endl;
	cout << "  Usage:" << endl;
	cout << "\t" << flg_plasma_file << " <file>\t\t: Input Plasma file in csv-format: x,y,q" << endl;
	cout << "\t" << flg_width << " <a>\t\t: Border length of unit cell." << endl;
	cout << "\t" << flg_temperature << " <T>\t\t\t: Initial temperature." << endl;
	cout << "\t" << flg_equilibrate << " <e>\t: Steps to perform before temperature decrease starts." << endl;
	cout << "\t" << flg_temperature_intermediate << " <i>\t: Steps between reducing temperature." << endl;
	cout << "\t" << flg_temperature_reduce << " <r>\t\t: Factor by which temperature gets reduced." << endl;
	cout << "\t" << flg_temperature_reduce_steps << " <s>\t: Times temperature gets reduced." << endl;
	cout << "\t" << flg_prefix << "\t\t\t: Prefix for output files. (\"" << prefix << "\")" << endl;
}

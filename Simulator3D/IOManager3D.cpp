/*
 * IOManager.cpp
 *
 *  Created on: 18.10.2010
 *      Author: Martin Rueckl
 */

#ifndef IOMANAGER3D_CPP_
#define IOMANAGER3D_CPP_

#include "IOManager3D.h"
#include <algorithm>
#include <util.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <Matrix2D.h>
#include <Matrix3D.h>


template<typename B,typename F>
IOManager3D<B,F>::IOManager3D(const SimulationEvaluator3D<B,F> &_evaluator, string _path):
	evaluator(_evaluator),
	simulator(_evaluator.getSimulator())
{
	unsigned u = 0;
	if (boost::filesystem::is_directory(_path)) {
		while (boost::filesystem::is_directory(_path + util::num2str(u)))
			u++;
		cerr << "Directory " << _path << " already exists, creating " << _path + util::num2str(u) << endl;
		_path = _path + util::num2str(u);
	}
	boost::filesystem::create_directory(_path);
	this->path = _path + "/";
	printInfo(cout);
	writeConfigurationFile		(path + "config.cfg");
	writeCapillaryConfiguration	(path + "capillary_configuration.csv");
	writeSequenceConfiguration	(path + "sequence_configuration.csv");
	writeField					(path + "field_configuration.csv");
	writeSampleTrajectory		(path + "sample_trajectories.csv");
	writeCollisionTree			(path + "collision_tree.csv");
}

template<typename B,typename F>
IOManager<B,F>::~IOManager()
{
	// TODO Auto-generated destructor stub
}

template<typename B,typename F>
void IOManager<B,F>::checkFile(string filename)
{
	if (!boost::filesystem::exists(filename)) throw MyException("The File " + filename + " does not exist.");
}

template<typename B,typename F>
void IOManager<B,F>::writeSequenceData(unsigned sequence, string filename, unsigned fract)
{
	filename = path + filename;
	cout << "writing sequence data to " << filename << "..." << endl;
	fstream real_signal_file((filename+"_real").c_str(), ios::out);
	fstream imag_signal_file((filename+"_imag").c_str(), ios::out);
	fstream error_file((filename + "_error").c_str(), ios::out);

	const Matrix2D &real_signal = evaluator.getRealSignal();
	const Matrix2D &imag_signal = evaluator.getImagSignal();
	const Matrix2D &error = evaluator.getError();

	real_signal_file << real_signal;
	imag_signal_file << imag_signal;
	error_file << error;

	real_signal_file.close();
	imag_signal_file.close();
	error_file.close();
}

template<typename B,typename F>
void IOManager<B,F>::writeConfigurationFile(string filename)
{
	cout << "Writing Configuration file " << filename << endl;
	fstream file(filename.c_str(), ios::out);
	printInfo(file);
}

template<typename B,typename F>
void IOManager<B,F>::printInfo(ostream &stream)
{
	stream << "Evaluator-Configuration:" << endl;
	stream << "\tSamples         =" << setw(12) << setfill(' ') << evaluator.getSamples() << endl;
	stream << "\tSequenceCount   =" << setw(12) << setfill(' ') << evaluator.getSequenceCount() << endl;
	stream << endl;

	stream << "Simulator-Configuration:" << endl;
	stream << "\tRandomSeed      =" << setw(12) << setfill(' ') << simulator.rng.getInitialSeed() << endl;
	stream << "\tStepSize        =" << setw(12) << setfill(' ') << simulator.getTimeStepSize() * 1E3 << " ms" << endl;
	stream << "\tSteps           =" << setw(12) << setfill(' ') << simulator.getTimeSteps() << endl;
	stream << "\tIntervall       =" << setw(12) << setfill(' ') << simulator.getTimeIntervall() * 1E3 << " ms" << endl;
	stream << "\tD               =" << setw(12) << setfill(' ') << simulator.getDiffusionConsant() * 1E9 << " um^2/ms" << endl;
	stream << "\tomega0          =" << setw(12) << setfill(' ') << simulator.getOmega0() << " rad/s" << endl;
	stream << "\tsigma           =" << setw(12) << setfill(' ') << simulator.getMeanStepSize() * 1E6 << " um" << "(sqrt(2*D*StepSize), MeanStepSize)" << endl;
	stream << "Boundary			 =" << setw(12) << setfill(' ') << simulator.getBoundaryCondition().getDescription() << endl;
	stream << endl;

	CapillaryConfiguration const &c = simulator.getCapillaryConfiguration();
	stream << "Capillary-Configuration:" << endl;
	stream << "\tmin(radius)     =" << setw(12) << setfill(' ') << c.getMinimalRadius() * 1E6 << " um" << endl;
	stream << "\tmax(radius)     =" << setw(12) << setfill(' ') << c.getMaximalRadius() * 1E6 << " um" << endl;
	stream << "\tnumber          =" << setw(12) << setfill(' ') << c.getCapillaryCount() << endl;
	stream << "\tA               =" << setw(12) << setfill(' ') << c.getCapillaryArea() * 1E12 << " um^2" << endl;
	stream << endl;

	stream << "Field-Configuration:" << endl;
	stream << "\tDescription:     " << simulator.getField().getDescription() << endl << endl;
}

template<typename B,typename F>
void IOManager<B,F>::writeCapillaryConfiguration(string filename)
{
	cout << "Writing Capillaries to " << filename << endl;
	fstream file(filename.c_str(), ios::out);
	file << simulator.getCapillaryConfiguration();
	file.close();
}

template<typename B,typename F>
void IOManager<B,F>::writeField(string filename, unsigned res)
{
	cout << "Writing Field2D-File " << filename << " ...";
	fstream file(filename.c_str(), ios::out);
	file << simulator.getField();
	file.close();
	cout << "Done" << endl;
}

template<typename B,typename F>
vector<vector<double> > IOManager<B,F>::readSequences(string sequenceFile)
{
	cout << "Reading Sequence-File " << sequenceFile << " ...";
	checkFile(sequenceFile);
	ifstream istream(sequenceFile.c_str());
	vector<vector<double> > ret;
	while (!istream.eof()) {
		vector<double> seq;
		while (istream.peek() != '\n' && !istream.eof()) {
			double d;
			istream >> d;
			istream.ignore();
			seq.push_back(d);
		}
		istream.ignore();
		ret.push_back(seq);
	}
	istream.close();
	cout << "Done" << endl;
	return ret;
}

/*
 * schreibt eine beispieltrajektorie in ein file mit gegebenem dateinamen
 */
template<typename B,typename F>
void IOManager<B,F>::writeSampleTrajectory(string filename)
{
	cout << "Writing Sample Trajectory " << filename << " ...";
	fstream sample_trajectory_file(filename.c_str(), ios::out);
	sample_trajectory_file << simulator.runSpacialResolvedTrajectory();
	sample_trajectory_file.close();
	cout << "Done" << endl;
}

template<typename B,typename F>
void IOManager<B,F>::writeCollisionTree(string filename)
{
	cout << "Writing Collision Tree " << filename << " ...";
	fstream tree_file(filename.c_str(), ios::out);
	tree_file << simulator.getCollisionTree();
	tree_file.close();
	cout << "Done" << endl;
}

template<typename B,typename F>
void IOManager<B,F>::writeSequenceConfiguration(string filename)
{
	cout << "Writing Sequences to " << filename << endl;
	fstream file(filename.c_str(), ios::out);
	file.setf(ios::showpos);
	const Matrix2D &sequences = evaluator.getSequences();
	file << sequences;
	file.close();
}

#endif /* IOMANAGER_CPP_ */

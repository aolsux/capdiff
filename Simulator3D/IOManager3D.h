/*
 * IOManager.h
 *
 *  Created on: 18.10.2010
 *      Author: Martin Rueckl
 */

#ifndef IOMANAGER_H_
#define IOMANAGER_H_

#include <string>
#include <fstream>
#include <iostream>
#include <time.h>
#include <boost/filesystem.hpp>

#include "Simulator3D.h"
#include "SimulationEvaluator3D.h"

using namespace std;

template<typename Simulator_Type>
class IOManager3D
{

public:
	IOManager3D(const SimulationEvaluator3D<Simulator_Type> &evaluator, string path);
	virtual ~IOManager3D();

	void printInfo(ostream &stream);

	void writeConfigurationFile(string filename);
	void writeCapillaryConfiguration(string filename);
	void writeSequenceConfiguration(string filename);
	void writeField(string filename, unsigned res = 1);
	void writeSampleTrajectory(string filename);
	void writeCollisionTree(string filename);

	void writeSequenceData(unsigned sequence, string filename, unsigned detailReduceLevel = 1);
	void writeSpatialData(string filename);

	static vector<vector<double> > readSequences(string sequenceFile);

private:

	static void checkFile(string filename);
	const SimulationEvaluator3D<Simulator_Type>   &evaluator;
	const Simulator_Type 						  &simulator;
	string 										  path;
};

#include "IOManager3D.cpp"

#endif /* IOMANAGER_H_ */

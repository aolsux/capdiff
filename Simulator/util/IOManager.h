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

#include "../Simulator.h"
#include "../SimulationEvaluator.h"
#include "PrecalculatedField2D.h"

using namespace std;

template<typename Boundary_Type,typename Field_Type>
class IOManager
{

public:
	IOManager(const SimulationEvaluator<Boundary_Type,Field_Type> &evaluator, string path);
	virtual ~IOManager();

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
	const SimulationEvaluator<Boundary_Type,Field_Type> &evaluator;
	const Simulator<Boundary_Type,Field_Type> 			&simulator;
	string 												path;
};

#include "IOManager.cpp"

#endif /* IOMANAGER_H_ */

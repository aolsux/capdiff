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
#include <MyException.h>
#include <boost/filesystem.hpp>

#include <util/boundaries/OuterBoundary.h>
#include <util/boundaries/CircleBoundary.h>
#include <util/boundaries/RectangularBoundary.h>
#include <util/boundaries/MixedBoundary.h>
#include <util/boundaries/ReflectiveRectangularBoundary.h>
#include <util/boundaries/CyclicRectangularBoundary.h>
#include <Simulator.h>
#include <SimulationEvaluator.h>
#include <SpinEchoSequence.h>

#include <util/MPICommunicator.h>
#include <util/Capillary2D.h>
#include <util/CapillaryConfiguration.h>
#include <util/Field2D.h>
#include <util/PrecalculatedField2D.h>
#include <util/LiveField2D.h>
#include <util.h>
#include <util/IOManager.h>

bool parseCommandLine(int argc, char* argv[]);
void printHelp();
void checkSettings();
void prepare();

template<typename Boundary_Type, typename Field_Type>
void work(Boundary_Type boundary, Field_Type field, CapillaryConfiguration capillaries);

const double GyromagneticRatio_Proton = 2.67522E8;/* 1 / (s T) */

const string flg_help1 = "--help", flg_help2 = "-h", flg_help3 = "--h", flg_help4 = "-help";

//parameters changeable by commandline
const string flg_Diffusion = "-D";//in um^2/ms
double diffusion = 1;// 1um^2/ms bis 5um^2/ms
const string flg_Samples = "-samples";
unsigned trajectorySamples = 25000;
const string flg_Steps = "-steps";
unsigned timeSteps = 0;

const string flg_Time = "-time";
double timeLength = 0;//20-40ms

//Off resonance settings
bool omega0set=false;
bool offresonanceset=false;
const string flg_Omega0 = "-omega";
double omega0 = 1000;
const string flg_Bz = "-B";
double Bz = 1;//7T
const string flg_DeltaChi = "-deltaChi";
double deltaChi = 1E-6;//8ppm
const string flg_GyromagRatio = "-gyromagRatio";
double gyromagneticRatio = GyromagneticRatio_Proton;//wasser

//const string flg_TiltingAngle = "-tiltingangle";
//double tiltingAngle = M_PI / 2;//90°

const string flg_MakeErrorMap = "-errormap";
bool errorMap = false;
double subsamples = 0.005*1E-6;
const string flg_Spatial = "-spatial";
bool spatial=false;
double spatial_res=0;
unsigned time_res=0;

//Sequence settings
const string flg_SequenceFile = "-sequencefile";
bool customSequences = false;
string sequenceFile;

//boundary settings
const string flg_CircleBoundary = "-radius";
const string flg_RectangularWidth = "-width";
const string flg_RectangularHeight = "-height";
const string flg_ReflectiveBoundary = "-reflective";
bool circleBoundary = false;
bool reflectiveBoundary = false;
double radius = 0;
double width = 0;//nächster nachbar 2*17um
double height = 0;//übernächster nachbar 2*17um*sqrt(3)

//capillaries
const string flg_CapillaryFile = "-capillaries";
bool customCapillaries = false;
string capillaryFile;

//field settings
const string flg_FieldFile = "-field";
bool customField = false;
string fieldFile;

//output settings
const string flg_OutputDirectory = "-out";
string outputDirectory = "DefaultOutput";

bool parseCommandLine(int argc, char* argv[])
{
	for (int a = 0; a < argc; a++) {
		if ((flg_help1.compare(argv[a]) == 0) | (flg_help2.compare(argv[a]) == 0) | (flg_help3.compare(argv[a]) == 0) | (flg_help4.compare(argv[a]) == 0)) {
			printHelp();
			return false;
		} else if (flg_Diffusion.compare(argv[a]) == 0 && a + 1 < argc) {
			diffusion = atof(argv[a + 1]) * 1E-9;
		} else if (flg_Samples.compare(argv[a]) == 0 && a + 1 < argc) {
			trajectorySamples = atoi(argv[a + 1]);
		} else if (flg_Steps.compare(argv[a]) == 0 && a + 1 < argc) {
			timeSteps = atoi(argv[a + 1]);
//		} else if (flg_TiltingAngle.compare(argv[a]) == 0 && a + 1 < argc) {
//			tiltingAngle = M_PI * atof(argv[a + 1]) / 180.;
		} else if (flg_Time.compare(argv[a]) == 0 && a + 1 < argc) {
			timeLength = atof(argv[a + 1]) * 1E-3;
		} else if (flg_Bz.compare(argv[a]) == 0 && a + 1 < argc) {
			Bz = atof(argv[a + 1]);
			offresonanceset=true;
		} else if (flg_DeltaChi.compare(argv[a]) == 0 && a + 1 < argc) {
			deltaChi = atof(argv[a + 1]) * 1E-6;
			offresonanceset=true;
		} else if (flg_GyromagRatio.compare(argv[a]) == 0 && a + 1 < argc) {
			gyromagneticRatio = atof(argv[a + 1]);
			offresonanceset=true;
		} else if (flg_Omega0.compare(argv[a]) == 0 && a + 1 < argc) {
			omega0 = atof(argv[a + 1]);
			omega0set=true;
			offresonanceset=true;
		}

		/*GEOMETRY*/
		else if (flg_CircleBoundary.compare(argv[a]) == 0 && a + 1 < argc) {
			radius = atof(argv[a + 1]) * 1E-6;
			circleBoundary = true;
			reflectiveBoundary = true;
		} else if (flg_RectangularHeight.compare(argv[a]) == 0 && a + 1 < argc) {
			height = atof(argv[a + 1]) * 1E-6;
		} else if (flg_RectangularWidth.compare(argv[a]) == 0 && a + 1 < argc) {
			width = atof(argv[a + 1]) * 1E-6;
		} else if (flg_ReflectiveBoundary.compare(argv[a]) == 0) {
			reflectiveBoundary = true;
		}

		/*SEQUENCE*/
		else if (flg_SequenceFile.compare(argv[a]) == 0 && a + 1 < argc) {
			sequenceFile = argv[a + 1];
			customSequences = true;
		}
		/*FIELD*/
		else if (flg_FieldFile.compare(argv[a]) == 0 && a + 1 < argc) {
			fieldFile = argv[a + 1];
			customField = true;
		}/*CAPILLARY*/
		else if (flg_CapillaryFile.compare(argv[a]) == 0 && a + 1 < argc) {
			capillaryFile = argv[a + 1];
			customCapillaries = true;
		}
		/*ERRORMAP*/
		else if (flg_MakeErrorMap.compare(argv[a]) == 0 && a + 1 < argc) {
			errorMap = true;
			subsamples = atof(argv[a + 1])*1E-6;
		}
		/*SPATIAL*/
		else if (flg_Spatial.compare(argv[a]) == 0 && a + 2 < argc) {
			spatial = true;
			spatial_res= atof(argv[a + 1])*1E-6;//abstand zwischen gitterpunkten in um
			time_res = atoi(argv[a + 2]);//samplerate der zeitschritte
		}
		/*OUTPUT*/
		else if (flg_OutputDirectory.compare(argv[a]) == 0 && a + 1 < argc) {
			outputDirectory = string("./") + argv[a + 1];
		}
	}
	return true;
}

void checkSettings()
{
	if (timeSteps == 0) 								cout << "WARNING: No value was set for the time step  size, it will be automatically adjusted to match an average step length of r_min/4." << endl;
	if (diffusion == 1) 								cout << "WARNING: No value was set for the diffusion constant. Using D="<<(diffusion=1E-9)<<"." << endl;
	if (trajectorySamples == 25000) 					cout << "WARNING: No value was set for the amount of trajectory samples. Using N=25000." << endl;
	if (outputDirectory.size() <= 2)					cout << "WARNING: Invalid output Directory, reset to " << (outputDirectory = "./DefaultOutput") << "!" << endl;
	//if (tiltingAngle == M_PI / 2)						cout << "WARNING: No value was set for the tilting angle. Using PHI=90°." << endl;
	if (Bz == 1 && !omega0set)							cout << "WARNING: No value was set for Bz. Using Bz=1T." << endl;
	if (deltaChi == 1 && !omega0set)					cout << "WARNING: No value was set for deltaChi. Using deltaChi=1." << endl;
	//if (gyromagneticRatio == GyromagneticRatio_Proton)cout << "WARNING: No value was set for the gyromagneticRatio. Using Protons gyromagnetic ratio." << endl;
	if (!offresonanceset)								cout << "WARNING: No value were set for off resonance, using omega0="+util::num2str(omega0)+"."<<endl;
	if (!omega0set)					omega0=0.5*GyromagneticRatio_Proton*deltaChi*Bz;

	if (timeLength == 0) 								throw MyException("No value was set for simulation time. Aborting.");
	if (circleBoundary && !reflectiveBoundary)			throw MyException("Invalid boundary condition: CyclicCircle. Aborting ");
	if (circleBoundary && radius <= 0) 					throw MyException("Invalid radius: "  + util::num2str(radius)+ ". Aborting.");
	if (!circleBoundary && width <= 0) 					throw MyException("Invalid width: " + util::num2str(width) + ". Aborting." );
	if (!circleBoundary && height <= 0) 				throw MyException("Invalid height: " + util::num2str(height) + ". Aborting.");
	if (customCapillaries && !boost::filesystem::exists(capillaryFile))	throw MyException("Invalid capillary file: " + capillaryFile + ". Aborting.");
	if (customField && !boost::filesystem::exists(fieldFile)) 			throw MyException("Invalid field file: " + fieldFile + ". Aborting.");
	if (customSequences && !boost::filesystem::exists(sequenceFile)) 	throw MyException("Invalid sequence file: " + sequenceFile + ". Aborting." );
	if (spatial && (time_res==0 || spatial_res==0)) 	throw MyException("Invalid spatial settings. Aborting." );
}

int main(int argc, char* argv[])
{
	MPICommunicator bla(argc, argv);
	if (!parseCommandLine(argc, argv)) return 0;
	checkSettings();
	prepare();
	return 0;
}

void prepare()
{
	CapillaryConfiguration capillaries;
	if(customCapillaries)
		capillaries = CapillaryConfiguration::readFile(capillaryFile);
	else
		capillaries.addCapillary(Capillary2D(Point(0, 0), 1E-6));

	if (customField) {
		PrecalculatedField2D field=PrecalculatedField2D::readFile(fieldFile, Point(0,0),Point(width, height));
		cout << "Field scaling errors:" << endl;
		cout << "\t" << field.getWidth() - ((double) field.columns - 1.) * field.getXResolution() << endl;
		cout << "\t" << field.getHeight() - ((double) field.lines - 1.) * field.getYResolution() << endl;
		if (circleBoundary){
			CircleBoundary boundary(radius);
			work(boundary, field, capillaries);
		}else{
			if(reflectiveBoundary){
				ReflectiveRectangularBoundary boundary(Point(0,0), Point(width, height));
				work(boundary, field, capillaries);
			}else{
				CyclicRectangularBoundary boundary(Point(0,0), Point(width, height));
				work(boundary, field, capillaries);
			}
		}
	}else{
		LiveField2D field(capillaries);
		if (circleBoundary){
			CircleBoundary boundary(radius);
			work(boundary, field, capillaries);
		}else{
			if(reflectiveBoundary){
				ReflectiveRectangularBoundary boundary(Point(0,0), Point(width, height));
				work(boundary, field, capillaries);
			}else{
				CyclicRectangularBoundary boundary(Point(0,0), Point(width, height));
				work(boundary, field, capillaries);
			}
		}
	}

}

template<typename Boundary_Type, typename Field_Type>
void work(Boundary_Type boundary, Field_Type field, CapillaryConfiguration capillaries)
{
	Simulator<Boundary_Type, Field_Type> simulator(boundary, field, timeLength, diffusion, timeSteps, omega0);
	simulator.setCapillaryConfiguration(capillaries);
	SimulationEvaluator<Boundary_Type,Field_Type> evaluator(simulator, trajectorySamples);
	//if (customSequences) evaluator.addSpinEchoSequences(IOManager<Boundary_Type,Field_Type>::readSequences(sequenceFile));

	IOManager<Boundary_Type,Field_Type> iomanager(evaluator, outputDirectory);

	if (errorMap)
	{
		cout << "Generating error map!" << endl;
		Matrix2D m_errorMap = simulator.getErrorMap(trajectorySamples, subsamples);
		fstream file("errorMap.csv", ios::out);
		file << m_errorMap;
		file.close();
	}
	else
	{
		if(!spatial)
		{
			evaluator.evaluateTrajectories();
			for (unsigned s = 0; s < evaluator.getSequenceCount(); s++)
				iomanager.writeSequenceData(s, "sequence" + util::num2str(s) + ".csv");
		}
		else
		{
			evaluator.getMagnetizationDistribution(spatial_res,time_res);
			iomanager.writeSpatialData("spatial.csv");
		}
	}
}

void printHelp()
{

	cout << "========================================================================" << endl;
	cout << "Tool for simulating spin dephasing." << endl;
	cout << "========================================================================" << endl;
	cout << "  General Parameters:" << endl;

	cout << "\t" << flg_Diffusion <<  "\t\t: Diffusion constant (" << diffusion << ")[um^2/ms]" << endl;
	cout << "\t" << flg_Samples <<      "\t: Amount of sampled trajectories (" << trajectorySamples << ")" << endl;
	cout << "\t" << flg_Steps <<      "\t\t: Amount of time steps in one trajectory (" << timeSteps << ", this means the value gets optimized internally)" << endl;
	cout << "\t" << flg_Time <<       "\t\t: Simulated length of time.[ms]" << endl;
	cout << "\t" << flg_Omega0 <<     "\t\t: Set off resonance frequency directly.[radiant/s]" << endl;
	cout << "\t" << flg_Bz << 		  "\t\t: Magnetic Field Strength in z-Direction . (" << Bz << ")[T]" << endl;
	cout << "\t" << flg_DeltaChi <<     "\t: Susceptibility difference from capillary to surrounding tissue. (" << deltaChi << ")[ppm]" << endl;
	cout << "\t" << flg_GyromagRatio << "\t: Gyromagnetic ratio for nuclii. (Hydrogen=" << GyromagneticRatio_Proton << ")[1/(sT)]" << endl;
//	cout << "\t" << flg_TiltingAngle << "\t: Tilting angle of capillary towards Bz. (90)[°]" << endl;
	cout << "\t" << flg_MakeErrorMap << "\t: Only Generate a error map for the phase errors per time step for the given settings. (" << subsamples << ")" << endl;
	cout << "\t" << flg_Spatial		 << " s t\t: Generate spatial and time resolved magnetization. spatial resolution given by s[um] with one time snapshot in t time steps." << endl;

	cout << "  Sequence Settings:" << endl;
	cout << "\t" << flg_SequenceFile << "\t: File for echo sequences to be simulated. (FID-signal)" << endl;

	cout << "  Boundary Settings:" << endl;
	cout << "\t" << flg_CircleBoundary << "\t\t: Radius for a circle shaped environment.[um]" << endl;
	cout << "\t" << flg_RectangularWidth << "\t\t: Width for a rectangular environment.[um]" << endl;
	cout << "\t" << flg_RectangularHeight << "\t\t: Height for a rectangular environment.[um]" << endl;
	cout << "\t" << flg_ReflectiveBoundary << "\t\t: Set Boundary Conditions to reflective for rectangular environments. (Default is cyclic)" << endl;

	cout << "  Capillary Settings:" << endl;
	cout << "\t" << flg_CapillaryFile << "\t: Capillary File. (Single Capillary with center (0/0) and radius 1um)" << endl;

	cout << "  Field Settings:" << endl;
	cout << "\t" << flg_FieldFile << "\t\t: Use a previously calculated Field from the file instead of a in time calculated one." << endl;

	cout << "  Output Settings:" << endl;
	cout << "\t" << flg_OutputDirectory << "\t\t: Redirect all output to this Folder. (\"DefaultOutput\")" << endl;

}


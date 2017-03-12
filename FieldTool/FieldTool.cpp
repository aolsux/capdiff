#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <limits>

#include <util.h>
#include <util/MPICommunicator.h>
#include <util/Capillary2D.h>
#include <util/CapillaryConfiguration.h>
#include <util/PrecalculatedField2D.h>
#include <util/Compartement.h>
#include <util/IOManager.h>

#include <Progress.h>
#include <PerfClock.h>

using namespace std;


PrecalculatedField2D calculateField(double cutoff, unsigned fieldRes, double width, double height);
void				 maskField(double cutOff);
PrecalculatedField2D shrinkField(unsigned subsampleresolution);
Matrix2D 			 makeHistogram(double omega_max, unsigned bins);
CapillaryConfiguration mirrorCapillaries(CapillaryConfiguration caps, Rectangle rect, double cut_off);


void printHelp();
void printFieldStats(const PrecalculatedField2D &field);


//omp
unsigned nthreads = 1;

//command line
string flg_cutoff = "-cutoff", flg_res = "-resolution", flg_width = "-width", flg_height = "-height";
string flg_inputCapillaryFile = "-load-capillaries", flg_inputFieldFile = "-load-field";

string flg_calculate = "-calculate";
string flg_scale = "-scale";
string flg_mask = "-mask";
string flg_shrink = "-shrink";
string flg_histo = "-histo";
string flg_centered = "-centered";
string flg_forceWrite = "-write";

string flg_maskThreshold = "-mask-cutoff";


const string flg_help1 = "--help", flg_help2 = "-h", flg_help3 = "--h", flg_help4 = "-help";

ProgressMonitor pmon;

//Global Variables
PrecalculatedField2D field, masked_field;
CapillaryConfiguration capillaries;
int zeros=-1;

//parameters changeable by comandline
double		width				= 20;
double		height				= 20;
double		cutoff				= 0;
double		field_res			= 100E6;
double		mask_threshold		= 0;
double		histo_omega_cutoff	= numeric_limits<double>::max();
unsigned	histo_bins			= 10000;
unsigned	shrink_factor		= 1;
double		scaling				= 1;
bool        force_write         = false;
string capillaryFileName = "";
string fieldFileName = "";

int main(int argc, char* argv[])
{
	MPICommunicator mpicomm(argc, argv);
	bool calculate			= false;
	bool shrink				= false;
	bool mask				= false;
	bool scale				= false;
	bool histo				= false;
	bool centered			= false;
	bool loadField			= false;
	bool loadCapillaries	= false;
	bool mask_threshold_set = false;
	bool width_set			= false;
	bool height_set			= false;


	//iterate over params and search for field params
	for (int a = 0; a < argc; a++) {
		if ((flg_help1.compare(argv[a]) == 0) | (flg_help2.compare(argv[a]) == 0) | (flg_help3.compare(argv[a]) == 0) | (flg_help4.compare(argv[a]) == 0)) {
			printHelp();
			return 0;
		} else if (flg_calculate.compare(argv[a]) == 0 ) {
			calculate = true;
		} else if (flg_cutoff.compare(argv[a]) == 0 && a + 1 < argc) {
			cutoff = atof(argv[a + 1]) * 1E-6;
		} else if (flg_res.compare(argv[a]) == 0 && a + 1 < argc) {
			field_res = atoi(argv[a + 1]) * 1E6;
		} else if (flg_width.compare(argv[a]) == 0 && a + 1 < argc) {
			width = atof(argv[a + 1]) * 1E-6;
			width_set=true;
		} else if (flg_height.compare(argv[a]) == 0 && a + 1 < argc) {
			height = atof(argv[a + 1]) * 1E-6;
			height_set=true;
		} else if (flg_inputCapillaryFile.compare(argv[a]) == 0 && a + 1 < argc) {
			capillaryFileName = argv[a + 1];
			loadCapillaries = true;
		} else if (flg_inputFieldFile.compare(argv[a]) == 0 && a + 1 < argc) {
			fieldFileName = argv[a + 1];
			loadField = true;
		} else if (flg_shrink.compare(argv[a]) == 0 && a + 1 < argc) {
			shrink_factor = atoi(argv[a + 1]);
			shrink = true;
		} else if (flg_mask.compare(argv[a]) == 0) {
			mask = true;
		} else if (flg_maskThreshold.compare(argv[a]) == 0 && a + 1 < argc) {
			mask_threshold = atof(argv[a + 1]);
			mask_threshold_set = true;
		} else if (flg_scale.compare(argv[a]) == 0 && a + 1 < argc) {
			scaling = atof(argv[a + 1]);
			scale = true;
		} else if (flg_centered.compare(argv[a]) == 0) {
			centered=true;
		}else if (flg_forceWrite.compare(argv[a]) == 0) {
			force_write=true;
		} else if (flg_histo.compare(argv[a]) == 0 ) {
			histo = true;
			if(a + 1 < argc)
			{
				histo_bins = (atoi(argv[a+1])==0 ? 10000 : atoi(argv[a+1]));
				if((histo_bins!=0) && (a + 2 < argc))
					histo_omega_cutoff = (atof(argv[a+2])==0 ? numeric_limits<double>::max() : atof(argv[a+2]));
			}
		}

	}
#pragma omp parallel
	{
		if (omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
	}

	if (!calculate && !loadField)					throw MyException("Cannot do anything without field.");
	if (!loadCapillaries && !loadField) 			throw MyException("No field or capillary file given.");
	if (!mask && mask_threshold_set) 				throw MyException("Cannot use -mask-cutoff without masking command.");
	if (calculate && !loadCapillaries)				throw MyException("Cannot calculate without capillary file.");
	if (calculate && (!width_set || !height_set))	throw MyException("Cannot calculate without width/height.");
	if (mask && (!width_set || !height_set))		throw MyException("Cannot mask without width/height.");
	if (mask && (!loadField || !loadCapillaries)) 	throw MyException("Cannot mask field without both, capillary file and field file.");
	if (histo && (!width_set || !height_set))						throw MyException("Cannot make a histogram without width/height.");
	if (histo && (!(loadField||calculate) || !loadCapillaries)) 	throw MyException("Cannot make a histogram without both, capillary file and field file.");


	if(loadCapillaries)
		capillaries = CapillaryConfiguration::readFile(capillaryFileName);

	//either generate field from capillary configuration or load it from file
	if(calculate){
		cout << "Calculating Field for "<< capillaryFileName <<
				" (width=" << width <<
				", height=" << height <<
				", resolution=" << field_res <<
				", cutoff=" << cutoff << ")"<<endl;
		field = calculateField(cutoff, field_res, width, height);
		string fieldFileName = capillaryFileName.substr(0,capillaryFileName.find('.')) +
				"_width=" + std::to_string(width)+
				"_height=" +  std::to_string(height)+
				"_res=" +  std::to_string(field_res)+
				"_cor=" +  std::to_string(cutoff)+".csv";
		cout << "Writing Field-File " << fieldFileName << " ...";
		fstream out_file(fieldFileName.c_str(), ios::out);
		out_file << field;
		out_file.close();
		cout << "Done" << endl;
	} else {
		if(centered)
			field = PrecalculatedField2D::readFile(fieldFileName, Point(-width/2,-height/2), Point(width/2, height/2));
		else
			field = PrecalculatedField2D::readFile(fieldFileName, Point(0,0), Point(width, height));
	}

	//make lower resolution field
	if (shrink){
		cout << "Shrinking " << fieldFileName << " by " << shrink_factor << endl;
		PrecalculatedField2D lowres_field = shrinkField(shrink_factor);
		cout << "Writing lower resolution Field-File lowres_" << fieldFileName << " ...";
		fstream lowres_file(("lowres_"+fieldFileName).c_str(), ios::out);
		lowres_file << lowres_field;
		lowres_file.close();
		cout << "Done" << endl;
	}

	//scale field
	if(scale){
		cout<< "Scaling Field by " << scaling << endl;
		field*=scaling;
		if(!histo | force_write)
		{
			cout << "Writing scaled Field-File scaled_field.csv ...";
			fstream scaled_file("scaled_field.csv", ios::out);
			scaled_file << field;
			scaled_file.close();
		}
		cout << "Done" << endl;
	}

	//mask field
	if (mask){
		cout << "Masking " << fieldFileName << " with " << capillaryFileName << endl;
		maskField(mask_threshold);
		cout << "Writing masked Field-File masked_field.csv ...";
		fstream masked_file("masked_field.csv", ios::out);
		masked_file << masked_field;
		masked_file.close();
		cout << "Done" << endl;
	}

	//make histogram
	if(histo){
		cout<< "Making Histogram for " << fieldFileName << endl;
		Matrix2D histo = makeHistogram(histo_omega_cutoff, histo_bins);
		cout << "Writing Histogram  histo.csv...";
		fstream histo_file("histo.csv", ios::out);
		histo_file << histo;
		histo_file.close();
		cout << "Done" << endl;
	}

	return 0;
}

Matrix2D makeHistogram(double omega_cutoff, unsigned histo_bins)
{
	/*Mirror capillaries with small cutoff radius in order to mask only "half circles" on the edge of the simulation environment.*/
	CapillaryConfiguration c_conf = mirrorCapillaries(capillaries, field, (field.getUpperRightCorner()-field.getLowerLeftCorner()).abs()*1.1);
	fstream capillary_file("mirrored_caps.csv", ios::out);
	capillary_file << c_conf;
	capillary_file.close();

	vector<const CollisionObject *> collisionObjects;
	for (unsigned int i = 0; i < c_conf.capillaries.size(); i++)
		collisionObjects.push_back(&(c_conf.capillaries[i]));
	Compartement *compartement = new Compartement(c_conf.getSystemSize(), collisionObjects);

	fstream collision_file("collision_tree.csv", ios::out);
	collision_file << *compartement;
	collision_file.close();

	ProgressBar *pbar=new ProgressBar("Making histogram", field.columns * field.lines);
	pmon.addProgressBar(pbar);
	omp_lock_t my_lock;
	omp_init_lock(&my_lock);
	double temp = 0;
#pragma omp parallel for shared(temp,compartement,field,my_lock)
	for (unsigned x = 0; x < field.lines; x++) {
		for (unsigned y = 0; y < field.columns; y++) {
			Point t = field.transform(x, y);
			if (compartement->getCollision(t) == NULL)
			{
				omp_set_lock(&my_lock);
				temp = max(temp,abs(field(x,y)));
				omp_unset_lock(&my_lock);
			}
		}
	}
	double		omega_max	= min(temp, omega_cutoff); //get max amplitude or cutoff
	double		domega		= 2.*omega_max/((double)histo_bins-1.);
	double collisions		= 0;
	Matrix2D histo(histo_bins,2);

	for (unsigned b = 0; b < histo_bins; b++) histo(b,0)=(double)b*domega-omega_max;

	Point shift(field.getXResolution()/2, field.getYResolution()/2);
	#pragma omp parallel shared(pbar,compartement,field,histo,my_lock)
	{
		Matrix2D thread_private_histo(histo_bins,2);
		double thread_private_collisions=0;
		#pragma omp for
		for (unsigned x = 0; x < field.lines-1; x++) //iterate over area elements, not supporting points
		{
			for (unsigned y = 0; y < field.columns-1; y++) //iterate over area elements, not supporting points
			{
				Point t = field.transform(x, y) + shift;
				if(compartement->getCollision(t) != NULL)
					thread_private_collisions++;
				else
				{
					if ((abs(field.getValue(t)) < omega_max))
					{
						unsigned bin = (field.getValue(t) + omega_max)/domega;
						thread_private_histo(bin,1)++;
					}
					else if(field.getValue(t) <= -omega_max)
						thread_private_histo(0,1)++;
					else if(field.getValue(t) >=  omega_max)
						thread_private_histo(histo_bins-1,1)++;
				}

				pbar->progress();
			}
		}
		/*Accumulate the private thread histograms...*/
		omp_set_lock(&my_lock);
		histo		+= thread_private_histo;
		collisions	+= thread_private_collisions;
		omp_unset_lock(&my_lock);
	}

	cout << "Excluded "<<collisions<<"("<<(double)collisions/((double)((field.lines-1)*(field.columns-1)))<<"%)"<<" points from histogram generation"<<endl;
	double		area_elem	= field.getXResolution()*field.getYResolution();
	double		area		= field.getArea()-collisions*area_elem;
	for (unsigned b = 0; b < histo_bins; b++) histo(b,1)=histo(b,1)*area_elem/(area*domega);
	pbar->finish();
	return histo;
}

void maskField(double cutOff)
{
	if(zeros!=-1)return;
	/*Mirror capillaries with small cutoff radius in order to mask only "half circles" on the edge of the simulation environment.*/
	CapillaryConfiguration c_conf = mirrorCapillaries(capillaries, field, (field.getUpperRightCorner()-field.getLowerLeftCorner()).abs()*1.1);
	zeros		 = 0;
	masked_field = field;
	ProgressBar *pbar=new ProgressBar("Masking Field", masked_field.columns * masked_field.lines);
	pmon.addProgressBar(pbar);

	vector<const CollisionObject *> collisionObjects;
	for (unsigned int i = 0; i < c_conf.capillaries.size(); i++)
		collisionObjects.push_back(&(c_conf.capillaries[i]));

	Compartement *compartement = new Compartement(c_conf.getSystemSize(), collisionObjects);

#pragma omp parallel for shared(pbar,masked_field) reduction(+:zeros)
	for (unsigned x = 0; x < masked_field.lines; x++) {
		for (unsigned y = 0; y < masked_field.columns; y++) {
			Point t = masked_field.transform(x, y);
			if ((compartement->getCollision(t) != NULL) | (abs(masked_field(x, y)) > cutOff && cutOff != 0))
			{
				masked_field(x, y) = 0;
				zeros++;
			}
			pbar->progress();
		}
	}
	pbar->finish();
	cout << "set: " << zeros << " point to zero!" << endl;
}

PrecalculatedField2D shrinkField(unsigned subsample)
{
	unsigned lines=field.getLines()/subsample;
	unsigned cols=field.getColumns()/subsample;
	PrecalculatedField2D small(Matrix2D(lines,cols), field.getLowerLeftCorner(), field.getUpperRightCorner());
	for (unsigned x = 0; x < small.lines; x++) {
		for (unsigned y = 0; y < small.columns; y++) {
			Point t = small.transform(x, y);
			small(x,y)=field(x*subsample,y*subsample);
		}
	}
	return small;
}

void printHelp()
{
	cout << "========================================================================" << endl;
	cout << "Tool for calculating and refactoring capillary fields. (c) Martin Rueckl" << endl;
	cout << "========================================================================" << endl;
	cout << "  Options:" << endl;
	cout << "\t" << flg_res <<  				" n\t\t: Target-Resolution of calculated field. [points per um]" << endl;
	cout << "\t" << flg_width << 				" w\t\t: Width (x-dimensional size) of calculated field. [um]" << endl;
	cout << "\t" << flg_height <<   			" h\t\t: Height (y-dimensional size) of calculated field. [um]" << endl;
	cout << "\t" << flg_cutoff << 				" c\t\t: Cutoff-Radius c for re-sampling of the given capillaries. [um]" << endl;
	cout << "\t" << flg_inputCapillaryFile <<  " path\t: Filename/Path of capillary-file." << endl;
	cout << "\t" << flg_inputFieldFile << 	 " path\t\t: Filename for field input." << endl;
	cout << "\t" << flg_centered <<				"  \t\t: Assume the given field and capillaries are centered around the origin."<<endl;
	cout << "\t" << flg_forceWrite <<			"    \t: Force write of temporarily generated fields and other data."<<endl;

	cout << "  Commands:" << endl;

	cout << "\t" << flg_calculate << " s\t\t: Calculate the field for a given capillary configuration." << endl;
	cout << "\t\t\t\t  (Needs: " << flg_inputCapillaryFile << ", "
								 << flg_width << ", "
								 << flg_height << ", "
								 << flg_res<< ", "
								 << flg_cutoff << ")" << endl;

	cout << "\t" << flg_shrink << " s\t\t: Create a file with smaller resolution. (Only take one point in s)" << endl;

	cout << "\t" << flg_mask << "\t\t\t: Mask the Field with a given capillary configuration." << endl;
	cout << "\t\t\t\t  (Needs: " << flg_inputCapillaryFile << ", "
								 << flg_inputFieldFile << ", "
								 << flg_width << ", "
								 << flg_height << ")" << endl;

	cout << "\t" << flg_histo << " b c\t\t\t: Make a histogram of the fields amplitudes, containing b bins and max frequency c. If no c and b are given, c=10000, no freuquency cap. The capillary area is excluded." << endl;
	cout << "\t\t\t\t  (Needs: " << flg_inputCapillaryFile << ", "
								 << flg_inputFieldFile << ", "
								 << flg_width << ", "
								 << flg_height << ")" << endl;

	cout << "\t\t\t\t  (Optional: " << flg_maskThreshold << " for a cutoff of field amplitudes)" << endl;

	cout << "\t" << flg_scale << " s\t\t: Scale fields amplitude by s." << endl;
}



PrecalculatedField2D calculateField(double cutOff, unsigned fieldRes, double width, double height)
{
	CapillaryConfiguration conf = mirrorCapillaries(capillaries,Rectangle(Point(0,0),Point(width, height)), cutOff);
	if (force_write)
	{
		std::ofstream mconf("mirrorer_capillaries.csv");
		mconf << conf;
	}
	PrecalculatedField2D field(Point(0,0),Point(width, height), 1. / fieldRes);
	ProgressBar *pbar=new ProgressBar("Calculating field (" + util::num2str(nthreads) + " CPUs)", field.columns * field.lines);
	pmon.addProgressBar(pbar);
	#pragma omp parallel for shared(pbar,field)
	for (unsigned ydot = 0; ydot < field.columns; ydot++) {
		for (unsigned xdot = 0; xdot < field.lines; xdot++) {
			Point t = field.transform(xdot, ydot);
			double val = 0;
			for (unsigned c = 0; c < conf.capillaries.size(); c++) {
				double r_c = conf.capillaries[c].getRadius();
				Point loc = conf.capillaries[c].getLocation();
				Point dif = loc - t;
				Point angle = conf.capillaries[c]._getAngle();
				//((x - y) (x + y) Cos[2 a] - 2 x y Sin[2 a])/(x^2 + y^2)^2
				val += r_c * r_c * ((dif.x - dif.y) * (dif.x + dif.y) * (angle.x * angle.x - angle.y * angle.y) - 2 * dif.x * dif.y * 2 * angle.x * angle.y)
						/ (dif.absPow2() * dif.absPow2());

			}
			field(xdot, ydot) = (std::isnan(val) | std::isinf(val))? 0.0 : val;
			pbar->progress();

		}
	}
	pbar->finish();
	return field;
}

void printFieldStats(const PrecalculatedField2D &field)
{
//	cout << "width=" << field.width << endl;
//	cout << "height=" << field.height << endl;
//	cout << "Xdots=" << field.lines << endl;
//	cout << "Ydots=" << field.columns << endl;
//	cout << "Xres=" << field.resX << endl;
//	cout << "Yres=" << field.resY << endl;
}


//
//void calculateField2()
//{
//	CapillaryConfiguration work = IOManager::readCapillaryConfiguration(capillaryFileName);
//	field = new PrecalculatedField2D(width, height, 1. / fieldRes);
//	c_conf.addCapillary(Capillary2D(Point(0, 0), 1E-6));//create capillary in center with 1um radius
//	c_conf.mirrorCapillaries(width, height, cutOff);
//
//	PerfClock::tick();
//	//calculate single capillary field
//	ProgressBar pbar("Calculating 1Capillary periodic field (" + util::num2str(work.capillaries.size()) + " Capillaries)", field->columns * field->lines
//			* work.capillaries.size());
//	//calculateFieldSlice(pbar, 0, field->columns);
//
//	//add field for each capillary
//	PrecalculatedField2D *single = field;
//	field = new PrecalculatedField2D(width, height, 1. / fieldRes);
//	pbar.set(0);
//#pragma omp parallel for shared(pbar, field, single, work)
//	for (unsigned u = 0; u < work.capillaries.size(); u++) {
//		const Point &loc = work.capillaries[u].getLocation();
//		const double r_squared = work.capillaries[u].getRadius() * work.capillaries[u].getRadius() / (1E-6 * 1E-6);
//		for (unsigned x = 0; x < field->lines; x++)
//			for (unsigned y = 0; y < field->columns; y++) {
//				Point p(x * field->resX - field->width / 2, y * field->resY - field->height / 2);
//				/* Vector from capillary to actual field point (x,y).
//				 This vector might point to a location not included in field.
//				 Hence periodic continuation is used.*/
//				Point dif = p - loc;
//				if (dif.x >= field->width / 2) dif.x -= field->width;
//				else if (dif.x <= -field->width / 2) dif.x += field->width;
//				if (dif.y >= field->height / 2) dif.y -= field->height;
//				else if (dif.y <= -field->height / 2) dif.y += field->height;
//				(*field)(x, y) += r_squared * single->getValue(dif);
//				pbar.progress();
//			}
//	}
//	cout << "duration:" << PerfClock::tack() << endl;
//	std::stringstream ss;
//	ss << capillaryFileName << "_width=" << width << "_height=" << height << "_res=" << fieldRes << "_cor=" << cutOff << "_superposition.csv";
//	IOManager::writeField(field, ss.str(), 1);
//	delete field;
//}



CapillaryConfiguration mirrorCapillaries(CapillaryConfiguration caps, Rectangle rect, double cut_off)
{
	if (cut_off <= 0) 					throw MyException("Invalid cutoff radius.");
	CapillaryConfiguration mirrored;
//	int x_low	= -1 * ceil(cut_off / rect.getWidth())/2 - 3;
//	int x_high	= 	   ceil(cut_off / rect.getWidth())/2 + 3;
//	int y_low	= -1 * ceil(cut_off / rect.getHeight())/2 - 3;
//	int y_high	= 	   ceil(cut_off / rect.getHeight())/2 + 3;

	double low  =  -1 * ceil(cut_off / min(rect.getWidth(),rect.getHeight())) - 1;
	double high =       ceil(cut_off / min(rect.getWidth(),rect.getHeight())) + 2;
	for (double xs = low ; xs <= high; xs++) {
		for (double ys = low; ys <= high; ys++) {
			for (unsigned i = 0; i < caps.getCapillaryCount(); i++)
			{
				Capillary2D cap = caps.capillaries[i];
				Point shiftetLocation = Point(xs * rect.getWidth(), ys * rect.getHeight()) + cap.getLocation();
				if ((shiftetLocation-rect.getCenter()).abs() <= cut_off)
					mirrored.addCapillary(Capillary2D(shiftetLocation, cap.getRadius(), cap.getAngle()));
			}
		}
	}
//	fstream cap_out_file(("mirrored_" + capillaryFileName.substr(0,capillaryFileName.find('.')) + "_cor=" + util::num2str(cut_off)+".csv").c_str(), ios::out);
//	cap_out_file << mirrored;
//	cap_out_file.close();
	return mirrored;
}


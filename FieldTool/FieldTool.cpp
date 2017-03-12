#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <omp.h>
#include <limits>

#include <util.h>
#include <util/MPICommunicator.h>
#include <util/Capillary2D.h>
#include <util/CapillaryConfiguration.h>
#include <util/PrecalculatedField2D.h>
#include <util/Compartement.h>
#include <util/IOManager.h>

#include <boost/format.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <Progress.h>
#include <PerfClock.h>

using namespace std;
namespace po = boost::program_options;


PrecalculatedField2D calculateField(const CapillaryConfiguration& caps, double cutoff, unsigned fieldRes, double width, double height, bool debug);
void				 maskField(double cutOff);
PrecalculatedField2D shrinkField(unsigned subsampleresolution);
Matrix2D 			 makeHistogram(double omega_max, unsigned bins);
CapillaryConfiguration mirrorCapillaries(CapillaryConfiguration caps, Rectangle rect, double cut_off);
void printFieldStats(const PrecalculatedField2D &field);

int main(int argc, char* argv[])
{
	MPICommunicator mpicomm(argc, argv);	
	
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("width", po::value<double>(), "Width (x-dimension size) of calculated field. [µm]")
		("height",po::value<double>(), "Height (y-dimensional size) of calculated field. [µm]")
		("cutoff",po::value<double>(), "Cutoff-Radius for re-sampling of capillaries around the field. [µm]")
		("resolution",po::value<double>(), "Target-Resolution of calculated field. [points per µm]")
		("debug",po::bool_switch()->default_value(false), "generate various additional outputs (also some files)")
		("threads",po::value<int>()->default_value(std::thread::hardware_concurrency()), "Number of threads to use for calculation\n"
																						 "(defaults to # of cpus of machine)")
		("calculate", po::value<std::string>(),"Generate a field from the given capillary configuration.\n"
											   "Filename to the capillary file to calculate the field for.\n"
											   "  NOTE: the cappillary file needs to follow this format:\n"
											   "     x,y,radius,angle\n"
											   "  where x,y, and radius are given in METERS\n")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    
	
	if(vm.count("width"))
		std::cout << "given width: " << vm["width"].as<double>() << std::endl;

	if (vm.count("help")) {
		cout << "========================================================================" << endl;
		cout << "Tool for calculating and refactoring capillary fields. (c) Martin Rueckl" << endl;
		cout << "========================================================================" << endl;
		cout << desc << "\n";
		return 1;
	}

	omp_set_dynamic(0);     // Explicitly disable dynamic teams
	omp_set_num_threads(vm["threads"].as<int>()); // Use n threads for all consecutive parallel regions
	
	// generate field from capillary configuration
	if(vm.count("calculate"))
	{
		if (!vm.count("width") || !vm.count("height") || !vm.count("resolution") || !vm.count("cutoff"))
		{
			std::cerr << "Error: the --calculate option requires --width, --height, --cutoff, and --resolution" << std::endl;
			std::cerr << "Usage:" << std::endl << desc << std::endl;			
			return 1;
		}
		
		double width = vm["width"].as<double>();
		double height = vm["height"].as<double>();
		double resolution = vm["resolution"].as<double>();
		double cutoff = vm["cutoff"].as<double>();
		std::string capfile = vm["calculate"].as<std::string>();
		
		CapillaryConfiguration capillaries = CapillaryConfiguration::readFile(capfile);
		
		std::cout << boost::format("Calculating Field for %s (width=%f, height=%f, resolution=%f, cutoff=%f)") 
						% capfile % width % height % resolution % cutoff << std::endl;
		auto field = calculateField(capillaries, cutoff, resolution, width, height, vm["debug"].as<bool>());
		std::string fieldFileName = capfile.substr(0,capfile.find('.')) +
				"_width=" + std::to_string(width)+
				"_height=" + std::to_string(height)+
				"_res=" +  std::to_string(resolution)+
				"_cor=" +  std::to_string(cutoff)+".csv";
		std::cout << "Writing Field-File " << fieldFileName << " ...";
		fstream out_file(fieldFileName.c_str(), ios::out);
		out_file << field;
		out_file.close();
		cout << "Done" << endl;
		return 0;
	} 
	else
	{
		std::cerr << "Dont know what to do... give either --calculate, or others (todo implement)" << std::endl;
		std::cerr << "Usage:\n" << desc << std::endl;
		return 1;
	}
	
	/*
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
	*/
	return 0;
}

/*
Matrix2D makeHistogram(double omega_cutoff, unsigned histo_bins)
{
	//Mirror capillaries with small cutoff radius in order to mask only "half circles" on the edge of the simulation environment.
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
		//Accumulate the private thread histograms...
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
*/

/*
void maskField(double cutOff)
{
	if(zeros!=-1)return;
	//Mirror capillaries with small cutoff radius in order to mask only "half circles" on the edge of the simulation environment.
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
}*/

/*
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
*/

PrecalculatedField2D calculateField(const CapillaryConfiguration& capillaries, double cutOff, unsigned fieldRes, double width, double height, bool force_write)
{
	CapillaryConfiguration conf = mirrorCapillaries(capillaries,Rectangle(Point(0,0),Point(width, height)), cutOff);
	if (force_write)
	{
		std::ofstream mconf("mirrorer_capillaries.csv");
		mconf << conf;
	}
	PrecalculatedField2D field(Point(0,0),Point(width, height), 1. / fieldRes);
	ProgressMonitor pmon;
	ProgressBar *pbar=new ProgressBar("Calculating field (" + std::to_string(std::thread::hardware_concurrency()) + " CPUs)", field.columns * field.lines);
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


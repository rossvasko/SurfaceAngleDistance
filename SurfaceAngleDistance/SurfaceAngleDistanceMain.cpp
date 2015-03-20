// SurfaceAngleDistance.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include "ijk.txx"
#include "ijkIO.txx"
#include "SurfaceAngleDistance.h"
#include "ANGLE_DIST_PARAM.h"

/////////////Linux includes
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath> 
#include <string.h>
#include <stdio.h>
#include <stdio.h>
#include <float.h>
#include <stdlib.h>   
////////////////////////////

using namespace std;

string file1 = "";
string file2 = "";

string prefix = "";

string extra_data_testing = "";
string run_compare = "";

bool epsilonSet = false;
bool deltaLengthSet = false;


bool extNeighborhood = false;
bool numSimplicesHistogram = false;
bool areaSimplicesHistogram = false;
double areaBelowCutoff = -1;
double out_tri_angle = -1;
double binary_coloring_angle = -1;

double edge_width = -1;
double epsilon = 0;
double deltaLength = 0;

int sub_resolution = 1;

bool ignore_surface_orientation = true;
bool out_tri = false;
bool testing_script_mode = false;
bool testing_script_mode_nums = false;
bool back_coloring = false;
bool adjust_coloring_range = false;
bool percent_done = false;
bool binary_coloring = false;
bool out_param = false;
bool show_time = false;
bool info = false;



void usage_error(){
	cout << "Usage: SurfaceAngleDistance [OPTIONS] <surface1.off> <surface2.off>" << endl;
	cout << "Required Args: " << endl;
	cout << "\t-e <x> Sets epsilon neighborhod to x" << endl;
	//cout << "\t\t-n <1/0> sets extended neighborhood to t/f (1/0)" << endl;
	cout << "\t-d <x> sets the delta edge distance to x" << endl;

	cout << "Optional Args: " << endl;
	cout << "\t-nh Generates histogram of number of simplices" << endl;
	cout << "\t-ah Generates histogram of normalized area" << endl;
	cout << "\t-ca <x> Prints proportion of area under the angle difference x" << endl;
	cout << "\t-w <x> Sets the box width of the box width intersection to x" << endl;
	cout << "\t-sr <x> Sets the dimensions of the sub-regions of the grid to x" << endl;
	cout << "\t-prefix <x> Sets prefix of file output to x" << endl;
	cout << "\t-out_tri <x> Outputs triangles above angle cutoff x" << endl;
	cout << "\t-orient Calculates angle differences and accounts for orientation" << endl;
	cout << "\t-back_color Outputs back color to off file" << endl;
	cout << "\t-adjust_coloring_range Outputs color in a 0-1 range" << endl;
	cout << "\t-percent_done Outputs percent done" << endl;
	cout << "\t-binary_coloring <x> changes coloring to be two colors with cutoff of x" << endl;
	cout << "\t-out_param outputs the set parameters" << endl;
	cout << "\t-time outputs timing of the program" << endl;
	cout << "\t-info outputs additional info" << endl;

	exit(0);
}

void compute_color(double * colors, int angle, double custom_color){

	double max_angle = 90;
	if (!ignore_surface_orientation){
		max_angle = 180;
	}

	if (angle < max_angle / 4){
		colors[0] = 0;
		colors[1] = custom_color * (angle - 0 * max_angle / 4) / (max_angle / 4);
		colors[2] = custom_color;
	}
	else if (angle < 2 * max_angle / 4){
		colors[0] = 0;
		colors[1] = custom_color;
		colors[2] = custom_color - custom_color * (angle - 1 * max_angle / 4) / (max_angle / 4);
	}
	else if (angle < 3 * max_angle / 4){
		colors[0] = custom_color * (angle - 2 * max_angle / 4) / (max_angle / 4);;
		colors[1] = custom_color;
		colors[2] = 0;
	}
	else{
		colors[0] = custom_color;
		colors[1] = custom_color - custom_color * (angle - 3 * max_angle / 4) / (max_angle / 4);
		colors[2] = 0;
	}
}

void parse_command_line(int argc, char * argv[]){

	if (argc >= 7){
		//Reads in two surface files
		file1 = argv[argc - 2];
		file2 = argv[argc - 1];

		for (int x = 1; x < argc - 2; x++){
			//Setting the epsilon neighborhood
			if (!strcmp(argv[x], "-e")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> epsilon;
				x++;
				epsilonSet = true;
			}
			//Setting the extended neighborhood boolean
			else if (!strcmp(argv[x], "-n")){
				extNeighborhood = strcmp(argv[x + 1], "0");
				x++;
			}
			//Setting the delta edge length
			else if (!strcmp(argv[x], "-d")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> deltaLength;
				x++;
				deltaLengthSet = true;
			}
			//Setting number of simplices histogram to true
			else if (!strcmp(argv[x], "-nh")){
				numSimplicesHistogram = true;
			}
			//Setting area of simplices histogram to true
			else if (!strcmp(argv[x], "-ah")){
				areaSimplicesHistogram = true;
			}
			//Setting variable to print the area of simplices whose max min angle is under
			else if (!strcmp(argv[x], "-ca")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> areaBelowCutoff;
				x++;
			}
			//Setting variable to of box width
			else if (!strcmp(argv[x], "-w")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> edge_width;
				x++;
			}
			//Setting variable for file prefix
			else if (!strcmp(argv[x], "-prefix")){
				prefix = argv[x + 1];
				x++;
			}
			//Setting variable to consider orientation
			else if (!strcmp(argv[x], "-orient")){
				ignore_surface_orientation = false;
			}
			//Setting variable for testing script
			else if (!strcmp(argv[x], "-tsm")){
				extra_data_testing = argv[x + 1];
				x++;
				run_compare = argv[x + 1];
				x++;
				testing_script_mode = true;
			}
			//Setting variable for outputting tri above angle cutoff
			else if (!strcmp(argv[x], "-out_tri")){
				out_tri = true;
				istringstream s;
				s.str(argv[x + 1]);
				s >> out_tri_angle;
				x++;
			}
			//Setting variable for outputting back color
			else if (!strcmp(argv[x], "-back_color")){
				back_coloring = true;
			}
			//Setting variable for outputting 0-1 color range
			else if (!strcmp(argv[x], "-adjust_coloring_range")){
				adjust_coloring_range = true;
			}
			//Setting variable for outputting percent done
			else if (!strcmp(argv[x], "-percent_done")){
				percent_done = true;
			}
			//Setting variable for binary coloring
			else if (!strcmp(argv[x], "-binary_coloring")){
				binary_coloring = true;
				istringstream s;
				s.str(argv[x + 1]);
				s >> binary_coloring_angle;
				x++;
			}
			//Setting variable for param output
			else if (!strcmp(argv[x], "-out_param")){
				out_param = true;
			}
			//Setting variable for param output
			else if (!strcmp(argv[x], "-time")){
				show_time = true;
			}
			//Setting variable for param output
			else if (!strcmp(argv[x], "-info")){
				info = true;
			}
			//Setting variable for testing script
			else if (!strcmp(argv[x], "-tsm_nums")){
				testing_script_mode_nums = true;
			}
			//Setting variable for testing script
			else if (!strcmp(argv[x], "-sr")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> sub_resolution;
				x++;
			}
			//An unknown argument was sent
			else{
				cout << "'" << argv[x] << "' is not a recognized argument." << endl;
				usage_error();
			}
		}
	}
	else{
		//There are not the required amount of arguments
		cout << "Not enough arguments were passed" << endl;
		usage_error();
	}

	if (!(deltaLengthSet && epsilonSet)){
		//The required parameters were not set
		cout << "Delta and epsilon were not set" << endl;
		usage_error();
	}

	if (!testing_script_mode && out_param){
		cout << "Arguments: " << endl;
		cout << "File name 1: " << file1 << endl;
		cout << "File name 2: " << file2 << endl;
		//cout << "Extended neighborhood: " << extNeighborhood << endl;
		cout << "Epsilon neighborhood: " << epsilon << endl;
		cout << "Delta edge length: " << deltaLength << endl;
		cout << "Sub-resolution length: " << sub_resolution << endl;
		cout << "Creating histogram of number of simplices: " << numSimplicesHistogram << endl;
		cout << "Creating histogram of area of simplices: " << areaSimplicesHistogram << endl;
		cout << "Printing area of simplices under cutoff: " << areaBelowCutoff << endl;
		cout << "File name prefix: " << prefix << endl;
		cout << "Coloring range: " << adjust_coloring_range << endl;
		cout << "Back color: " << back_coloring << endl;
		cout << "Percent done: " << percent_done << endl;
	}

}

void memory_exhaustion()
{
	cerr << "Error: Out of memory.  Terminating program." << endl;
	exit(10);
}

int main(int argc, char** argv)
{
	std::set_new_handler(memory_exhaustion);

	//Parsing command line
	parse_command_line(argc, argv);

	//Variables to read in surfaces
	vector<double> s1_points;
	vector<int> s1_simplices;

	vector<double> s2_points;
	vector<int> s2_simplices;

	int dim = 3;
	int dim_surface = 2;

	ifstream surface_file1, surface_file2;
	surface_file1.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	surface_file2.exceptions(std::ifstream::failbit | std::ifstream::badbit);

	try
	{
		surface_file1.open(file1.c_str());
	}
	catch (std::ifstream::failure e)
	{
		cerr << "Exiting" << endl;
		cerr << "Error opening file: " << file1 << endl;
		exit(10);
	}
	try
	{
		surface_file2.open(file2.c_str());
	}
	catch (std::ifstream::failure e)
	{
		cerr << "Error opening file: " << file2 << endl;
		cerr << "Exiting" << endl;
		exit(20);
	}

	//Reading surfaces
	try{
		IJK::ijkinOFF(surface_file1, dim, dim_surface, s1_points, s1_simplices);
	}
	catch (IJK::PROCEDURE_ERROR error){
		cerr << "Error while reading surface of " << file1 << endl;
		if (error.NumMessages() == 0) {
			cerr << "Unknown error." << endl;
		}
		else {
			error.Print(cerr);
		}
		cerr << "Exiting." << endl;
		exit(30);
	}
	catch (...){
		cerr << "Unknown error while reading surface of " << file1 << endl;
		exit(50);
	}

	if (!testing_script_mode && !testing_script_mode_nums){
		cout << "Size of simplices 1: " << s1_simplices.size() / 3 << endl;
	}
	try{
		IJK::ijkinOFF(surface_file2, dim, dim_surface, s2_points, s2_simplices);
	}
	catch (IJK::PROCEDURE_ERROR error){
		cerr << "Error while reading surface of " << file2 << endl;
		if (error.NumMessages() == 0) {
			cerr << "Unknown error." << endl;
		}
		else {
			error.Print(cerr);
		}
		cerr << "Exiting." << endl;
		exit(30);
	}
	catch (...){
		cerr << "Unknown error while reading surface of " << file2 << endl;
		exit(50);
	}

	if (!testing_script_mode && !testing_script_mode_nums){
		cout << "Size of simplices 2: " << s2_simplices.size() / 3 << endl;
	}

	//s1_simplices.erase(s1_simplices.begin() + 300, s1_simplices.end());
	//s2_simplices.erase(s2_simplices.begin(), s2_simplices.end() - 300);

	vector<double> s1_angle_distances(s1_simplices.size() / 3);
	vector<double> s2_angle_distances(s2_simplices.size() / 3);


	ANGLE_DIST_PARAM params;
	params.epsilon = epsilon;
	params.delta = deltaLength;
	params.width = edge_width;
	params.s1_angle_distances = &s1_angle_distances;
	params.s2_angle_distances = &s2_angle_distances;
	params.num_simplices_histogram = numSimplicesHistogram;
	params.area_simplices_histogram = areaSimplicesHistogram;
	params.print_area_angle_cutoff = areaBelowCutoff;
	params.prefix = prefix;
	params.orient = ignore_surface_orientation;
	params.testing_script_mode = testing_script_mode;
	params.testing_script_mode_nums = testing_script_mode_nums;
	params.percent_done = percent_done;
	params.time = show_time;
	params.info = info;
	params.sub_resolution = sub_resolution;

	double max_of_min_angles = surface_angle_distance::surface_angle_distance_extended(s1_points, s1_simplices, s2_points, s2_simplices, params);
	if (!testing_script_mode){
		//cout << "Maximum of minimum angles: " << max_of_min_angles << endl;
	}
	else{
		cout << extra_data_testing << "," << run_compare << endl;
	}


	if (out_tri){
		cout << "Simplices in surface 1 above angle cutoff: " << endl;
		for (int n = 0; n < s1_angle_distances.size(); n++){
			if (s1_angle_distances[n] > out_tri_angle){
				cout << "\tSimplex: " << n << " (";
				for (int m = 0; m < 3; m++){
					cout << s1_simplices[3 * n + m];
					if (m < 2){
						cout << ", ";
					}
					else{
						cout << ") ";
					}
				}
				cout << "Angle: " << s1_angle_distances[n];
				cout << " Vertex coords: ";
				for (int m = 0; m < 3; m++){
					cout << "(";
					for (int p = 0; p < 3; p++){
						cout << s1_points[3 * s1_simplices[3 * n + m] + p];
						if (p < 2){
							cout << " ";
						}
					}
					cout << ") ";
				}
				cout << endl;
			}
		}

		cout << "Simplices in surface 2 above angle cutoff: " << endl;
		for (int n = 0; n < s2_angle_distances.size(); n++){
			if (s2_angle_distances[n] > out_tri_angle){
				cout << "\tSimplex: " << n << " (";
				for (int m = 0; m < 3; m++){
					cout << s2_simplices[3 * n + m];
					if (m < 2){
						cout << ", ";
					}
					else{
						cout << ") ";
					}
				}
				cout << "Angle: " << s2_angle_distances[n];
				cout << " Vertex coords: ";
				for (int m = 0; m < 3; m++){
					cout << "(";
					for (int p = 0; p < 3; p++){
						cout << s2_points[3 * s2_simplices[3 * n + m] + p];
						if (p < 2){
							cout << " ";
						}
					}
					cout << ") ";
				}
				cout << endl;
			}
		}
	}

	string off_file = prefix + "surface_coloring_1.off";
	if (testing_script_mode){
		off_file = prefix + "perfect_coloring.off";
	}
	ofstream color_output(off_file.c_str());

	int numv_per_simplex = 3;

	double * coord = &s1_points[0];
	int numv = s1_points.size() / 3;

	int * simplex_vert = &s1_simplices[0];
	int nums = s1_simplices.size() / 3;

	std::vector<double> front_color;
	std::vector<double> back_color;


	double max_angle_calcs = 90;

	double custom_color = 255;
	if (adjust_coloring_range){
		custom_color = 1;
	}

	if (max_of_min_angles != 0){
		if (!binary_coloring){
			for (int x = 0; x < s1_angle_distances.size(); x++){

				/*
				int r = (255 - 255 * (double) (max_of_min_angles - s1_angle_distances[x]) / max_of_min_angles);
				int g = 0;
				int b = 255 * (max_of_min_angles - s1_angle_distances[x]) / max_of_min_angles;
				int alpha = 1;
				*/
				double r = (custom_color - custom_color * (double)(max_angle_calcs - s1_angle_distances[x]) / max_angle_calcs);
				double g = 0;
				double b = custom_color * (max_angle_calcs - s1_angle_distances[x]) / max_angle_calcs;
				double alpha = 1;

				double * colors = new double[3];
				compute_color(colors, s1_angle_distances[x], custom_color);
				r = colors[0];
				g = colors[1];
				b = colors[2];

				front_color.push_back(r);
				front_color.push_back(g);
				front_color.push_back(b);
				front_color.push_back(1);

				back_color.push_back(r);
				back_color.push_back(g);
				back_color.push_back(b);
				back_color.push_back(1);
			}
		}
		else{
			for (int x = 0; x < s1_angle_distances.size(); x++){

				double r = 0;
				double g = 0;
				double b = custom_color;
				double alpha = 1;

				if (s1_angle_distances[x] >= binary_coloring_angle){
					r = custom_color;
					b = 0;
				}

				front_color.push_back(r);
				front_color.push_back(g);
				front_color.push_back(b);
				front_color.push_back(1);

				back_color.push_back(r);
				back_color.push_back(g);
				back_color.push_back(b);
				back_color.push_back(1);
			}
		}
	}
	else{
		for (int x = 0; x < s1_angle_distances.size(); x++){
			double r = 0;
			double g = 0;
			double b = custom_color;
			double alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			front_color.push_back(alpha);

			back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			back_color.push_back(alpha);
		}
	}

	try{
		if (adjust_coloring_range){
			if (back_coloring){
				IJK::ijkoutColorFacesOFF(color_output, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, back_color);
			}
			else{
				IJK::ijkoutColorFacesOFF(color_output, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, std::vector<double>());
			}
		}
		else{
			std::vector<int> front_color_int;
			std::vector<int> back_color_int;
			for (int n = 0; n < front_color.size(); n++){
				front_color_int.push_back((int)front_color[n]);
				back_color_int.push_back((int)back_color[n]);
			}
			if (back_coloring){
				IJK::ijkoutColorFacesOFF(color_output, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color_int, back_color_int);
			}
			else{
				IJK::ijkoutColorFacesOFF(color_output, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color_int, std::vector<int>());
			}
		}
	}
	catch (IJK::PROCEDURE_ERROR error){
		cerr << "Error while writing surface of " << off_file << endl;
		if (error.NumMessages() == 0) {
			cerr << "Unknown error." << endl;
		}
		else {
			error.Print(cerr);
		}
		cerr << "Exiting." << endl;
		exit(30);
	}
	catch (...){
		cerr << "Unknown Error while writing surface of " << off_file << endl;
		exit(50);
	}

	off_file = prefix + "surface_coloring_2.off";
	if (testing_script_mode){
		off_file = prefix + "religrad_coloring.off";
	}
	ofstream color_output2(off_file.c_str());

	numv_per_simplex = 3;

	coord = &s2_points[0];
	numv = s2_points.size() / 3;

	simplex_vert = &s2_simplices[0];
	nums = s2_simplices.size() / 3;

	front_color.clear();
	back_color.clear();



	if (max_of_min_angles != 0){
		if (!binary_coloring){
			for (int x = 0; x < s2_angle_distances.size(); x++){
				/*
				int r = (int)(255 - 255 * (double)(max_of_min_angles - s2_angle_distances[x]) / max_of_min_angles);
				int g = 0;
				int b = 255 * (max_of_min_angles - s2_angle_distances[x]) / max_of_min_angles;
				int alpha = 1;
				*/

				double r = (int)(custom_color - custom_color * (double)(max_angle_calcs - s2_angle_distances[x]) / max_angle_calcs);
				double g = 0;
				double b = 255 * (max_angle_calcs - s2_angle_distances[x]) / max_angle_calcs;
				double alpha = 1;

				double * colors = new double[3];
				compute_color(colors, s2_angle_distances[x], custom_color);
				r = colors[0];
				g = colors[1];
				b = colors[2];

				front_color.push_back(r);
				front_color.push_back(g);
				front_color.push_back(b);
				front_color.push_back(1);

				back_color.push_back(r);
				back_color.push_back(g);
				back_color.push_back(b);
				back_color.push_back(1);
			}
		}
		else{
			for (int x = 0; x < s2_angle_distances.size(); x++){

				double r = 0;
				double g = 0;
				double b = custom_color;
				double alpha = 1;

				if (s2_angle_distances[x] >= binary_coloring_angle){
					r = custom_color;
					b = 0;
				}

				front_color.push_back(r);
				front_color.push_back(g);
				front_color.push_back(b);
				front_color.push_back(1);

				back_color.push_back(r);
				back_color.push_back(g);
				back_color.push_back(b);
				back_color.push_back(1);
			}
		}
	}
	else{
		for (int x = 0; x < s2_angle_distances.size(); x++){
			double r = 0;
			double g = 0;
			double b = custom_color;
			double alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			front_color.push_back(alpha);

			back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			back_color.push_back(alpha);
		}
	}

	try{
		if (adjust_coloring_range){
			if (back_coloring){
				IJK::ijkoutColorFacesOFF(color_output2, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, back_color);
			}
			else{
				IJK::ijkoutColorFacesOFF(color_output2, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, std::vector<double>());
			}
		}
		else{
			std::vector<int> front_color_int;
			std::vector<int> back_color_int;
			for (int n = 0; n < front_color.size(); n++){
				front_color_int.push_back((int)front_color[n]);
				back_color_int.push_back((int)back_color[n]);
			}
			if (back_coloring){
				IJK::ijkoutColorFacesOFF(color_output2, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color_int, back_color_int);
			}
			else{
				IJK::ijkoutColorFacesOFF(color_output2, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color_int, std::vector<int>());
			}
		}
	}
	catch (IJK::PROCEDURE_ERROR error){
		cerr << "Error while writing surface of " << off_file << endl;
		if (error.NumMessages() == 0) {
			cerr << "Unknown error." << endl;
		}
		else {
			error.Print(cerr);
		}
		cerr << "Exiting." << endl;
		exit(30);
	}
	catch (...){
		cerr << "Unknown error while writing surface of " << off_file << endl;
		exit(50);
	}

	return 0;
}




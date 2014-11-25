// SurfaceAngleDistance.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "ijk.txx"
#include "ijkIO.txx"
#include "SurfaceAngleDistance.h"

using namespace std;

string file1 = "";
string file2 = "";

string prefix = "";

bool epsilonSet = false;
bool deltaLengthSet = false;


bool extNeighborhood = false;
bool numSimplicesHistogram = false;
bool areaSimplicesHistogram = false;
double areaBelowCutoff = -1;

double edge_width = -1;
double epsilon = 0;
double deltaLength = 0;

void usage_error(){
	cout << "Usage: SurfaceAngleDistance <surface1.off> <surface2.off>" << endl;
	cout << "Required Args: " << "-e <x> Sets epsilon neighborhod to x" << endl;
	//cout << "\t\t-n <1/0> sets extended neighborhood to t/f (1/0)" << endl;
	cout << "\t\t-d <x> sets the delta edge distance to x" << endl;

	exit(0);
}

void parse_command_line(int argc, char * argv[]){
	
	if (argc >= 7){
		//Reads in two surface files
		file1 = argv[1];
		file2 = argv[2];

		for (int x = 3; x < argc; x++){
			//Setting the epsilon neighborhood
			if (!strcmp(argv[x], "-e")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> epsilon;

				epsilonSet = true;
			}
			//Setting the extended neighborhood boolean
			else if (!strcmp(argv[x], "-n")){
				extNeighborhood = strcmp(argv[x + 1], "0");
			}
			//Setting the delta edge length
			else if (!strcmp(argv[x], "-d")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> deltaLength;

				deltaLengthSet = true;
			}
			//Setting number of simplices histogram to true
			else if(!strcmp(argv[x], "-nh")){
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
			}
			//Setting variable to print the area of simplices whose max min angle is under
			else if (!strcmp(argv[x], "-w")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> edge_width;
			}
			//Setting variable for file prefix
			else if (!strcmp(argv[x], "-prefix")){
				prefix = argv[x + 1];
			}
		}
	}
	else{
		//There are not the required amount of arguments
		usage_error();
	}

	if (!(deltaLengthSet && epsilonSet)){
		//The required parameters were not set
		usage_error();
	}

	cout << "Arguments: " << endl;
	cout << "File name 1: " << file1 << endl;
	cout << "File name 2: " << file2 << endl;
	//cout << "Extended neighborhood: " << extNeighborhood << endl;
	cout << "Epsilon neighborhood: " << epsilon << endl;
	cout << "Delta edge length: " << deltaLength << endl;
	cout << "Creating histogram of number of simplices: " << numSimplicesHistogram << endl;
	cout << "Creating histogram of area of simplices: " << areaSimplicesHistogram << endl;
	cout << "Printing area of simplices under cutoff: " << areaBelowCutoff << endl;
	cout << "File name prefix: " << prefix << endl;

}

int main(int argc, char** argv)
{

	//Parsing command line
	parse_command_line(argc, argv);

	//Variables to read in surfaces
	vector<double> s1_points;
	vector<int> s1_simplices;

	vector<double> s2_points;
	vector<int> s2_simplices;

	int dim = 3;
	int dim_surface = 2;

	ifstream surface_file1(file1);
	ifstream surface_file2(file2);

	//Reading surfaces
	IJK::ijkinOFF(surface_file1, dim, dim_surface, s1_points, s1_simplices);
	cout << "Size of simplices 1: " << s1_simplices.size() / 3 << endl;

	IJK::ijkinOFF(surface_file2, dim, dim_surface, s2_points, s2_simplices);	
	cout << "Size of simplices 2: " << s2_simplices.size() / 3 << endl;

	//s1_simplices.erase(s1_simplices.begin() + 300, s1_simplices.end());
	//s2_simplices.erase(s2_simplices.begin(), s2_simplices.end() - 300);

	vector<double> s1_angle_distances;
	vector<double> s2_angle_distances;

	double max_of_min_angles = surface_angle_distance::surface_angle_distance_extended(s1_points, s1_simplices, s2_points, s2_simplices, epsilon, deltaLength, edge_width, s1_angle_distances, s2_angle_distances, numSimplicesHistogram, areaSimplicesHistogram, areaBelowCutoff, prefix);
	cout << "Maximum of minimum angles: " << max_of_min_angles << endl;


	string off_file = prefix + "surface_coloring_1.off";
	ofstream color_output(off_file);

	int numv_per_simplex = 3;

	double * coord = &s1_points[0];
	int numv = s1_points.size() / 3;

	int * simplex_vert = &s1_simplices[0];
	int nums = s1_simplices.size() / 3;

	std::vector<double> front_color;
	std::vector<double> back_color;

	if (max_of_min_angles != 0){
		for (int x = 0; x < s1_angle_distances.size(); x++){
			double r = (1 - 1 * (double) (max_of_min_angles - s1_angle_distances[x]) / max_of_min_angles);
			double g = 0;
			double b = (double) (max_of_min_angles - s1_angle_distances[x]) / max_of_min_angles;
			double alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			//front_color.push_back(1);

			back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			//back_color.push_back(1);
		}
	}
	else{
		for (int x = 0; x < s1_angle_distances.size(); x++){
			int r = 1;
			int g = 0;
			int b = 0;
			int alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			//front_color.push_back(alpha);

			/*back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			back_color.push_back(alpha);*/
		}
	}



	IJK::ijkoutColorFacesOFF(color_output, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, back_color);

	off_file = prefix + "surface_coloring_2.off";
	ofstream color_output2(off_file);

	numv_per_simplex = 3;

	coord = &s2_points[0];
	numv = s2_points.size() / 3;

	simplex_vert = &s2_simplices[0];
	nums = s2_simplices.size() / 3;

	front_color.clear();
	back_color.clear();

	if (max_of_min_angles != 0){
		for (int x = 0; x < s2_angle_distances.size(); x++){
			double r = (1 - 1 * (double)(max_of_min_angles - s2_angle_distances[x]) / max_of_min_angles);
			double g = 0;
			double b = (double)(max_of_min_angles - s2_angle_distances[x]) / max_of_min_angles;
			double alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			//front_color.push_back(1);

			back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			//back_color.push_back(1);
		}
	}
	else{
		for (int x = 0; x < s2_angle_distances.size(); x++){
			int r = 1;
			int g = 0;
			int b = 0;
			int alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			//front_color.push_back(alpha);

			/*back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			back_color.push_back(alpha);*/
		}
	}

	IJK::ijkoutColorFacesOFF(color_output2, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, back_color);

	return 0;
}




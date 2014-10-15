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

double epsilon = 0;
bool extNeighborhood = false;

double deltaLength = 0;

void usage_error(){
	cout << "Usage: SurfaceAngleDistance <surface1.off> <surface2.off>" << endl;
	cout << "Required Args: " << "-e <x> Sets epsilon neighborhod to x" << endl;
	cout << "\t\t-n <1/0> sets extended neighborhood to t/f (1/0)" << endl;
	cout << "\t\t-d <x> sets the delta edge distance to x" << endl;

	exit(0);
}

void parse_command_line(int argc, char * argv[]){
	
	if (argc == 9){
		//Reads in two surface files
		file1 = argv[1];
		file2 = argv[2];

		for (int x = 3; x < argc; x++){
			//Setting the epsilon neighborhood
			if (!strcmp(argv[x], "-e")){
				istringstream s;
				s.str(argv[x + 1]);
				s >> epsilon;
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
			}
		}
	}
	else{
		//There are not the required amount of arguments
		usage_error();
	}

	cout << "Arguments: " << endl;
	cout << "File name 1: " << file1 << endl;
	cout << "File name 2: " << file2 << endl;
	cout << "Extended neighborhood: " << extNeighborhood << endl;
	cout << "Epsilon neighborhood: " << epsilon << endl;
	cout << "Delta edge length: " << deltaLength << endl;
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
	IJK::ijkinOFF(surface_file2, dim, dim_surface, s2_points, s2_simplices);	

	cout << "Size of simplices 1: " << s1_simplices.size() / 3 << endl;
	cout << "Size of simplices 2: " << s2_simplices.size() / 3 << endl;

	vector<double> s1_angle_distances;
	vector<double> s2_angle_distances;

	double max_of_min_angles = surface_angle_distance::surface_angle_distance_extended(s1_points, s1_simplices, s2_points, s2_simplices, epsilon, deltaLength, s1_angle_distances, s2_angle_distances);
	cout << "Maximum of minimum angles: " << max_of_min_angles << endl;


	string off_file = "surface_coloring_1.off";
	ofstream color_output(off_file);

	int numv_per_simplex = 3;

	double * coord = &s1_points[0];
	int numv = s1_points.size() / 3;

	int * simplex_vert = &s1_simplices[0];
	int nums = s1_simplices.size() / 3;

	std::vector<double> front_color;
	std::vector<double> back_color;

	cout << nums << endl;

	if (max_of_min_angles != 0){
		for (int x = 0; x < s1_angle_distances.size(); x++){
			double r = (1 - 1 * (double) (max_of_min_angles - s1_angle_distances[x]) / max_of_min_angles);
			double g = 0;
			double b = (double) (max_of_min_angles - s1_angle_distances[x]) / max_of_min_angles;
			double alpha = 1;

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
			int r = 255;
			int g = 0;
			int b = 0;
			int alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			front_color.push_back(alpha);

			/*back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			back_color.push_back(alpha);*/
		}
	}



	IJK::ijkoutColorFacesOFF(color_output, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, back_color);

	off_file = "surface_coloring_2.off";
	ofstream color_output2(off_file);

	numv_per_simplex = 3;

	coord = &s2_points[0];
	numv = s2_points.size() / 3;

	simplex_vert = &s2_simplices[0];
	nums = s2_simplices.size() / 3;

	front_color.clear();
	back_color.clear();

	cout << nums << endl;

	if (max_of_min_angles != 0){
		for (int x = 0; x < s2_angle_distances.size(); x++){
			double r = (1 - 1 * (double)(max_of_min_angles - s2_angle_distances[x]) / max_of_min_angles);
			double g = 0;
			double b = (double)(max_of_min_angles - s2_angle_distances[x]) / max_of_min_angles;
			double alpha = 1;

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
			int r = 1;
			int g = 0;
			int b = 0;
			int alpha = 1;

			front_color.push_back(r);
			front_color.push_back(g);
			front_color.push_back(b);
			front_color.push_back(alpha);

			/*back_color.push_back(r);
			back_color.push_back(g);
			back_color.push_back(b);
			back_color.push_back(alpha);*/
		}
	}

	IJK::ijkoutColorFacesOFF(color_output2, dim, numv_per_simplex, coord, numv, simplex_vert, nums, front_color, back_color);

	return 0;
}




#ifndef ANGLE_DIST_PARAM_H
#define ANGLE_DIST_PARAM_H

#include <vector>
#include <stdlib.h>   

class ANGLE_DIST_PARAM{

public:
	double epsilon;
	double delta;
	double width;
	std::vector<double> * s1_angle_distances;
	std::vector<double> * s2_angle_distances;
	bool num_simplices_histogram;
	bool area_simplices_histogram;
	double print_area_angle_cutoff;
	std::string prefix;
	bool orient;
	bool testing_script_mode;
	bool testing_script_mode_nums;
	bool percent_done;
	bool time;
	bool info;

};

#endif
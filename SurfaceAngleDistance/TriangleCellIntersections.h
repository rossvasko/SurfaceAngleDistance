#ifndef TriangleCellIntersection_h
#define TriangleCellIntersection_h

#include <vector>

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

class SubGrid{

public:
	SubGrid(){

	}

	bool non_empty = false;
	void *sub_grid;

};

class TriangleCellIntersection{

public:
	TriangleCellIntersection(const std::vector<double> & points, const std::vector<int> & simplices, double width, int sub_resolution, double & average_boxes, double & intersection_per_box);
	double get_width();
	void find_simplices_in_shell_distance_from_cube(const double point[], int shell_level, std::vector<int> & simplices);
	void free_memory();

private:

	int lengths[3];
	int original_lengths[3];
	double mins[3];
	double maxes[3];
	double width;
	int sub_resolution;

	SubGrid *list_of_intersecting_simplices;

	void find_mins_and_maxes(const std::vector<double> & points);
	int find_box_index(int x, int y, int z);
	void find_boxes_triangle_intersects(const double triangle_points[], std::vector<int> & boxes_intersected);
	bool vector_normal_to_triangle(const double triangle[], double normal_vector[]);
	bool region_test(const double triangle[], int a, int b, int c);
	bool triangle_square_intersection(const double triangle2D[], int a, int b);
	bool line_square_intersection(const double line[], const double point[], int a, int b);
	int find_box_index_sub(int x, int y, int z);

};

#endif
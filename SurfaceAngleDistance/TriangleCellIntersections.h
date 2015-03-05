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

class TriangleCellIntersection{

public:
	TriangleCellIntersection(const std::vector<double> & points, const std::vector<int> & simplices, double width, double & average_boxes, double & intersection_per_box);
	double get_width();
	void find_simplices_in_shell_distance_from_cube(const std::vector<double> & point, int shell_level, std::vector<int> & simplices);

private:

	int lengths[3];
	double mins[3];
	double maxes[3];
	double width;

	std::vector<std::vector<int> > list_of_intersecting_simplices;
	
	void find_mins_and_maxes(const std::vector<double> & points);
	int find_box_index(int x, int y, int z);
	void find_boxes_triangle_intersects(const std::vector<double> & triangle_points, std::vector<int> & boxes_intersected);
	bool vector_normal_to_triangle(const std::vector<double> & triangle, std::vector<double> & normal_vector);
	bool region_test(const std::vector<double> triangle, int a, int b, int c);
	bool triangle_square_intersection(const std::vector<double> & triangle2D, int a, int b);
	bool line_square_intersection(const std::vector<double> & line, const std::vector<double> & point, int a, int b);

};

#endif
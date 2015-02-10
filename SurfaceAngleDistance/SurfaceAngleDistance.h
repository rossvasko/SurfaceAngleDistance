#include "TriangleCellIntersections.h"
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

class segmentedTriangle{

public:
	segmentedTriangle(){
	
	}
	segmentedTriangle::segmentedTriangle(int size);

	std::vector<double> points;
	std::vector<int> simplices;

	std::vector<double> boundary_points;

	std::vector<double> normal;
	std::vector<double> normal_ab;
	std::vector<double> normal_ac;
	std::vector<double> normal_bc;

	bool valid;

	std::vector<double> angle_diff_of_closest;
	std::vector<double> closest_distance;
	std::vector<double> smallest_distance;
	std::vector<int> closest;
	std::vector<double> closest_point;

	bool checked;

};

namespace surface_angle_distance
{
	double surface_angle_distance_extended(const std::vector<double> & s1_points, const std::vector<int> & s1_simplices, const std::vector<double> & s2_points, const std::vector<int> & s2_simplices, double epsilon, double delta, double width, std::vector<double> & s1_angle_distances, std::vector<double> & s2_angle_distances, bool num_simplices_histogram, bool area_simplices_histogram, double print_area_angle_cutoff, std::string prefix, bool orient, bool testing_script_mode, bool percent_done, bool time, bool info);
	void segment_triangle(std::vector<double> & triangle, double delta);
	double point_to_point_distance(const std::vector<double> & p1, const std::vector<double> & p2);
	double edge_length(const std::vector<double> & triangle);
	void create_segmented_triangles(const std::vector<double> & s_points, const std::vector<int> & s_simplices, std::vector<segmentedTriangle> & segmented_triangles, double delta);
	bool normal_unit_vector_to_triangle(const std::vector<double> & triangle, std::vector<double> & normal_vector);
	void compact_simplices(std::vector<int> & simplices, std::vector<double> & points);
	void surface_to_surface_extended(std::vector<segmentedTriangle> & triangles1, std::vector<segmentedTriangle> & triangles2, TriangleCellIntersection * grid, double epsilon, int surface_num, bool percent_done);
	void center_of_triangle(const std::vector<double> & points, std::vector<double> & point);
	void point_to_surface_nearest_point(const std::vector<double> & point, std::vector<segmentedTriangle> & surface, const std::vector<double> & normal_to_triangle, double & distance, double & angle, std::vector<double> & temp_point, TriangleCellIntersection * grid, const std::vector<double> & current_normal, int & minimum_angle_index, double epsilon, std::vector<double> & closest_point, double & distance_angle, double & smallest_distance);
	void point_to_triangle(const std::vector<double> & point, const segmentedTriangle & triangle, double & distance, std::vector<double> & closest_point);
	bool normal_cross_product(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & normal_vector);
	void vector_point_to_point(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & vector_ab);
	void create_vectors_from_triangle(const std::vector<double> & triangle, std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3);
	double dot_product(const std::vector<double> & a, const std::vector<double> & b);
	void point_to_line(const std::vector<double> & point, const std::vector<double> & line_a, const std::vector<double> & line_b, std::vector<double> & closest_point);
	void vector_projection(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & projection);
	double vector_magnitude(const std::vector<double> & a);
	double angle_of_vectors(const std::vector<double> & a, const std::vector<double> & b);
	void compute_angle_difference(const std::vector<double> & point, const std::vector<segmentedTriangle> & surface, const std::vector<double> & normal_to_triangle, double & distance, int & minimum_angle_index, double & minimum_angle_diff, std::vector<double> & minimum_angle_point, const std::vector<double> & point_to_compute_distance_from, double epsilon, TriangleCellIntersection * grid);
	void create_num_simplices_histogram(const std::vector<double> & angles1, const std::vector<double> & angles2, std::string prefix, bool testing_script_mode);
	void create_area_simplices_histogram(const std::vector<segmentedTriangle> & segmented_triangles_s1, const std::vector<segmentedTriangle> & segmented_triangles_s2, std::string prefix, bool testing_script_mode);
	double magnitude_cross_product(const std::vector<double> & a, const std::vector<double> & b);
	void print_info_under_cutoff(const std::vector<segmentedTriangle> & segmented_triangles_s1, const std::vector<segmentedTriangle> & segmented_triangles_s2, const std::vector<double> & angles1, const std::vector<double> & angles2, double cutoff, bool testing_script_mode);
	double smallest_edge_length(const std::vector<double> & triangle);
	int edge_length_vertex(const std::vector<double> & triangle);
}
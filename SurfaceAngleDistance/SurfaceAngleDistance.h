#include "TriangleCellIntersections.h"
#include "ANGLE_DIST_PARAM.h"
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

	void set_size(int size);

	int simplices_size;

	double *points;
	int *simplices;

	double boundary_points[9];

	double normal[3];
	double normal_ab[3];
	double normal_ac[3];
	double normal_bc[3];

	double a_to_b[3];
	double b_to_c[3];
	double a_to_c[3];
	double c_to_b[3];

	double ab_mag;
	double bc_mag;
	double ac_mag;

	double a_dot;
	double b_dot;
	double c_dot;

	double sphere_center[3];
	double sphere_radius;

	bool valid;

	double *angle_diff_of_closest;
	double *smallest_distance;
	int *closest;
	double *closest_point;

	bool checked;

	double area;
	double edge_length;

	int triangle_case;

};

namespace surface_angle_distance
{
	double surface_angle_distance_extended(const std::vector<double> & s1_points, const std::vector<int> & s1_simplices, const std::vector<double> & s2_points, const std::vector<int> & s2_simplices, ANGLE_DIST_PARAM & params);
	void segment_triangle(std::vector<double> & triangle, double delta);
	double point_to_point_distance(const double p1[], const double p2[]);
	double edge_length(const std::vector<double> & triangle);
	void create_segmented_triangles(const std::vector<double> & s_points, const std::vector<int> & s_simplices, std::vector<segmentedTriangle> & segmented_triangles, double delta);
	bool normal_unit_vector_to_triangle(const std::vector<double> & triangle, std::vector<double> & normal_vector);
	void compact_simplices(std::vector<int> & simplices, std::vector<double> & points);
	void surface_to_surface_extended(std::vector<segmentedTriangle> & triangles1, std::vector<segmentedTriangle> & triangles2, TriangleCellIntersection * grid, double epsilon, int surface_num, bool percent_done);
	void center_of_triangle(const double triangle[], double point[]);
	void point_to_surface_nearest_point(const double point[], std::vector<segmentedTriangle> & surface, const double normal_to_triangle[], double & angle, double closest_point[], TriangleCellIntersection * grid, const double current_normal[], int & minimum_angle_index, double epsilon, double & smallest_distance);
	void point_to_triangle(const double point[], const segmentedTriangle & triangle, double & distance, double closest_point[]);
	bool normal_cross_product(const double a[], const double b[], double normal_vector[]);
	void vector_point_to_point(const double a[], const double b[], double ab[]);
	void create_vectors_from_triangle(const double triangle[], double vec1[], double vec2[], double vec3[]);
	double dot_product(const double a[], const double b[]);
	void point_to_line(const double point[], const double line_a[], const double line_b[], double closest_point[], const double line_ab[], double magnitude);
	//void point_to_line(const double point[], const double line_a[], const double line_b[], double closest_point[]);
	void vector_projection(const double a[], const double b[], double projection[], double magnitude);
	double vector_magnitude(const double a[]);
	double angle_of_vectors(const double a[], const double b[]);
	void compute_angle_difference(const std::vector<double> & point, const std::vector<segmentedTriangle> & surface, const std::vector<double> & normal_to_triangle, double & distance, int & minimum_angle_index, double & minimum_angle_diff, std::vector<double> & minimum_angle_point, const std::vector<double> & point_to_compute_distance_from, double epsilon, TriangleCellIntersection * grid);
	void create_num_simplices_histogram(const std::vector<double> & angles1, const std::vector<double> & angles2, std::string prefix, bool testing_script_mode);
	void create_area_simplices_histogram(const std::vector<segmentedTriangle> & segmented_triangles_s1, const std::vector<segmentedTriangle> & segmented_triangles_s2, std::string prefix, bool testing_script_mode);
	double magnitude_cross_product(const double a[], const double b[]);
	void print_info_under_cutoff(const std::vector<segmentedTriangle> & segmented_triangles_s1, const std::vector<segmentedTriangle> & segmented_triangles_s2, const std::vector<double> & angles1, const std::vector<double> & angles2, double cutoff, bool testing_script_mode);
	double smallest_edge_length(const std::vector<double> triangle);
	int edge_length_vertex(const std::vector<double> triangle);
	double cross_product(const double a[], const double b[], double normal_vector[]);
	void vector_projection_normal(const double a[], const double b[], double projection[]);
	double point_to_point_distance_squared(const double p1[], const double p2[]);
	double vector_magnitude_squared(const double a[]);
	void free_memory(std::vector<segmentedTriangle> & triangles, TriangleCellIntersection & cells);
}
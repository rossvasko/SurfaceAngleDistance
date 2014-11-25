#include "TriangleCellIntersections.h"
#include <vector>

class segmentedTriangle{

public:
	segmentedTriangle(){
	
	}

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
	std::vector<int> closest;
	std::vector<double> closest_point;

};

namespace surface_angle_distance
{
	double surface_angle_distance_extended(std::vector<double> s1_points, std::vector<int> s1_simplices, std::vector<double> s2_points, std::vector<int> s2_simplices, double epsilon, double delta, double width, std::vector<double> & s1_angle_distances, std::vector<double> & s2_angle_distances, bool num_simplices_histogram, bool area_simplices_histogram, double print_area_angle_cutoff, std::string prefix);
	void segment_triangle(std::vector<double> & triangle, double delta);
	double point_to_point_distance(std::vector<double> p1, std::vector<double> p2);
	double edge_length(std::vector<double> triangle);
	void create_segmented_triangles(std::vector<double> s_points, std::vector<int> s_simplices, std::vector<segmentedTriangle> & segmented_triangles, double delta);
	bool normal_unit_vector_to_triangle(std::vector<double> triangle, std::vector<double> & normal_vector);
	void compact_simplices(std::vector<int> & simplices, std::vector<double> & points);
	void surface_to_surface_extended(std::vector<segmentedTriangle> & triangles1, std::vector<segmentedTriangle> triangles2, TriangleCellIntersection * grid, double epsilon, int surface_num);
	void center_of_triangle(std::vector<double> points, std::vector<double> & point);
	void point_to_surface_nearest_point(std::vector<double> point, std::vector<segmentedTriangle> surface, std::vector<double> normal_to_triangle, double & distance, std::vector<double> & closest_point, TriangleCellIntersection * grid);
	void point_to_triangle(std::vector<double> point, segmentedTriangle triangle, double & distance, std::vector<double> & closest_point);
	bool normal_cross_product(std::vector<double> a, std::vector<double> b, std::vector<double> & normal_vector);
	void vector_point_to_point(std::vector<double> a, std::vector<double> b, std::vector<double> & vector_ab);
	void create_vectors_from_triangle(std::vector<double> triangle, std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3);
	double dot_product(std::vector<double> a, std::vector<double> b);
	void point_to_line(std::vector<double> point, std::vector<double> line_a, std::vector<double> line_b, std::vector<double> & closest_point);
	void vector_projection(std::vector<double> a, std::vector<double> b, std::vector<double> & projection);
	double vector_magnitude(std::vector<double> a);
	double angle_of_vectors(std::vector<double> a, std::vector<double> b);
	void compute_angle_difference(std::vector<double> point, std::vector<segmentedTriangle> surface, std::vector<double> normal_to_triangle, double & distance, int & minimum_angle_index, double & minimum_angle_diff, std::vector<double> & minimum_angle_point, std::vector<double> point_to_compute_distance_from, double epsilon, TriangleCellIntersection * grid);
	void create_num_simplices_histogram(std::vector<double> angles1, std::vector<double> angles2, std::string prefix);
	void create_area_simplices_histogram(std::vector<segmentedTriangle> segmented_triangles_s1, std::vector<segmentedTriangle> segmented_triangles_s2, std::string prefix);
	double magnitude_cross_product(std::vector<double> a, std::vector<double> b);
	void print_area_simplices_under_cutoff(std::vector<segmentedTriangle> segmented_triangles_s1, std::vector<segmentedTriangle> segmented_triangles_s2, double cutoff);
}
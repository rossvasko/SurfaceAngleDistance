#ifndef TriangleCellIntersection_h
#define TriangleCellIntersection_h

#include <vector>

class TriangleCellIntersection{

public:
	TriangleCellIntersection(std::vector<double> points, std::vector<int> simplices, double width, double & average_boxes);
	double get_width();
	void find_simplices_in_shell_distance_from_cube(std::vector<double> point, int shell_level, std::vector<int> & simplices);

private:

	int lengths[3];
	double mins[3];
	double maxes[3];
	double width;

	std::vector<std::vector<int>> list_of_intersecting_simplices;
	
	void find_mins_and_maxes(std::vector<double> points);
	int find_box_index(int x, int y, int z);
	void find_boxes_triangle_intersects(std::vector<double> triangle_points, std::vector<int> & boxes_intersected);
	bool vector_normal_to_triangle(std::vector<double> triangle, std::vector<double> & normal_vector);

};

#endif
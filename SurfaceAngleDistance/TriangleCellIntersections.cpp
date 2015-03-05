#include "TriangleCellIntersections.h"
#include <math.h>
#include <iostream>

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

const int DIM3 = 3;

TriangleCellIntersection::TriangleCellIntersection(const std::vector<double> & points, const std::vector<int> & simplices, double box_width, double & average_boxes, double & intersection_per_box){

	average_boxes = 0;

	//Finding the mins and maxes
	find_mins_and_maxes(points);

	width = box_width;

	//Assigning box lengths
	for (int n = 0; n < 3; n++){
		lengths[n] = (maxes[n] - mins[n]) / (width) + 1;
	}

	//Preparing triangle intersection vector
	int num_boxes = lengths[0] * lengths[1] * lengths[2];
	

	try
	{
		list_of_intersecting_simplices.resize(num_boxes);

	//Find the boxes the simplices intersect
	for (int n = 0; n < simplices.size(); n += 3){
		//Storing the points of a triangle
		std::vector<double> triangle_points(DIM3 * 3);

		//Looping through the three simplices
		for (int m = 0; m < 3; m++){
			//Looping throught the three points of each simplice
			for (int p = 0; p < 3; p++){
				triangle_points[3 * m + p] = (points[3 * simplices[n + m] + p]);
			}
		}
		//Find the boxes the triangles intersects
		std::vector<int> boxes_intersected;
		find_boxes_triangle_intersects(triangle_points, boxes_intersected);

		average_boxes += boxes_intersected.size();

		//Storing the simplex id to the corresponding list
		for (int m = 0; m < boxes_intersected.size(); m++){
			list_of_intersecting_simplices[boxes_intersected[m]].push_back(n / 3);
		}
	}

	intersection_per_box = average_boxes;
	intersection_per_box /= num_boxes;
	average_boxes /= (simplices.size() / 3);

	}
	catch (...)
	{
		std::cerr << "Errors creating intersection grid." << endl;
		std::cerr << "Use a larger -w parameter" << endl;
		std::cerr << "Exiting" << endl;
		exit(10);
	}
}

void TriangleCellIntersection::find_boxes_triangle_intersects(const std::vector<double> & triangle_points, std::vector<int> & boxes_intersected){

	boxes_intersected.clear();

	double triangle_mins[3];
	double triangle_maxes[3];
	std::vector<double> normal_vector(DIM3);
	
	if (vector_normal_to_triangle(triangle_points, normal_vector)){

		//Finding the constant in the plane equation containing the triangle
		double constant = 0;
		for (int n = 0; n < 3; n++){
			constant -= triangle_points[n] * normal_vector[n];
		}

		//Initializing mins and maxes
		for (int n = 0; n < 3; n++){
			triangle_maxes[n] = triangle_points[n];
			triangle_mins[n] = triangle_points[n];
		}



		//Finding mins and maxes
		//Looping through rest of points
		for (int n = 3; n < triangle_points.size(); n += 3){

			//Looping through x, y, and z
			for (int m = 0; m < 3; m++){
				//If the current point is less than the min, update
				if (triangle_points[n + m] < triangle_mins[m]){
					triangle_mins[m] = triangle_points[m + n];
				}
				//If the current point is greater than the max, update
				if (triangle_points[n + m] > triangle_maxes[m]){
					triangle_maxes[m] = triangle_points[m + n];
				}
			}
		}

		//Subtracting mins
		for (int n = 0; n < 3; n++){
			triangle_mins[n] -= mins[n];
			triangle_maxes[n] -= mins[n];
		}

		//The lower and upper boxes for x, y, and z
		int lower[3];
		int upper[3];

		for (int n = 0; n < 3; n++){
			lower[n] = (int)(triangle_mins[n] / width);
			upper[n] = (int)(triangle_maxes[n] / width);
		}

		std::vector<double> adjusted_triangle_points(triangle_points.size());

		for (int n = 0; n < triangle_points.size(); n++){
			adjusted_triangle_points[n] = ((triangle_points[n] - mins[n % 3]) / width);
		}
		
		//Recording intersected boxes
		for (int a = lower[0]; a <= upper[0] && a < lengths[0]; a++){
			for (int b = lower[1]; b <= upper[1] && b < lengths[1]; b++){
				for (int c = lower[2]; c <= upper[2] && c < lengths[2]; c++){

					//Records the value of the dot product of the vector normal to the plane and the vector between cube vertices and a triangle point
					std::vector<double> plane_values_vertices(8);
					for (int x = 0; x <= 1; x++){
						for (int y = 0; y <= 1; y++){
							for (int z = 0; z <= 1; z++){
								double value = constant;
								value += (normal_vector[0] * ((a + x) * width + mins[0]));
								value += (normal_vector[1] * ((b + y) * width + mins[1]));
								value += (normal_vector[2] * ((c + z) * width + mins[2]));
								plane_values_vertices[4 * x + 2 * y + z] = (value);
							}
						}
					}

					//If diagonal vertices have a different value, the cube is in the plane
					bool add_box = false;
					for (int x = 0; x < 4 && !add_box; x++){
						add_box = (plane_values_vertices[x] * plane_values_vertices[7 - x] <= 0);
					}

					//Record intersected boxes
					if (add_box){
						if (region_test(adjusted_triangle_points, a, b, c)){
							boxes_intersected.push_back(find_box_index(a, b, c));
						}
					}
				}
			}
		}


	}
}

bool TriangleCellIntersection::region_test(const std::vector<double> triangle, int a, int b, int c){
	
	std::vector<double> xyTriangle(6);
	int index = 0;
	for (int n = 0; n < triangle.size() - 1; n++){
		if (n % 3 == 2){
			n++;
		}
		xyTriangle[index++] = (triangle[n]);
	}

	std::vector<double> xzTriangle(6);
	index = 0;
	for (int n = 0; n < triangle.size(); n++){
		if (n % 3 == 1){
			n++;
		}
		xzTriangle[index++] = (triangle[n]);
	}

	std::vector<double> yzTriangle(6);
	index = 0;
	for (int n = 0; n < triangle.size(); n++){
		if (n % 3 == 0){
			n++;
		}
		yzTriangle[index++] = (triangle[n]);
	}

	return (triangle_square_intersection(xyTriangle, a, b) && triangle_square_intersection(xzTriangle, a, c) && triangle_square_intersection(yzTriangle, b, c));
}

bool TriangleCellIntersection::triangle_square_intersection(const std::vector<double> & triangle2D, int a, int b){

	std::vector<double> line3(3);
	line3[0] = (triangle2D[3] - triangle2D[1]);
	line3[1] = (triangle2D[0] - triangle2D[2]);
	line3[2] = (triangle2D[1] * triangle2D[2] - triangle2D[0] * triangle2D[3]);
	std::vector<double> point3(2);
	point3[0] = (triangle2D[4]);
	point3[1] = (triangle2D[5]);

	std::vector<double> line2(3);
	line2[0] = (triangle2D[1] - triangle2D[5]);
	line2[1] = (triangle2D[4] - triangle2D[0]);
	line2[2] = (triangle2D[5] * triangle2D[0] - triangle2D[4] * triangle2D[1]);
	std::vector<double> point2(2);
	point2[0] = (triangle2D[2]);
	point2[1] = (triangle2D[3]);

	std::vector<double> line1(3);
	line1[0] = (triangle2D[5] - triangle2D[3]);
	line1[1] = (triangle2D[2] - triangle2D[4]);
	line1[2] = (triangle2D[3] * triangle2D[4] - triangle2D[2] * triangle2D[5]);
	std::vector<double> point1(2);
	point1[0] = (triangle2D[0]);
	point1[1] = (triangle2D[1]);


	return line_square_intersection(line1, point1, a, b) && line_square_intersection(line2, point2, a, b) && line_square_intersection(line3, point3, a, b);
}

bool TriangleCellIntersection::line_square_intersection(const std::vector<double> & line, const std::vector<double> & point, int a, int b){

	double epsilon_value = .001;

	double point_value = line[0] * point[0] + line[1] * point[1] + line[2];

	bool return_value = false;

	for (int x = 0; x <= 1 && !return_value; x++){
		for (int y = 0; y <= 1 && !return_value; y++){
			double current_value = line[0] * (a + x) + line[1] * (b + y) + line[2];
			if (point_value > 0){
				return_value = current_value >= -epsilon_value;
			}else{
				return_value = current_value <= epsilon_value;
			}
		}
	}

	return return_value;
}

void TriangleCellIntersection::find_simplices_in_shell_distance_from_cube(const std::vector<double> & point, int shell_level, std::vector<int> & simplices){

	simplices.clear();

	int box_coords[3];

	//Finding the coordinates of the box the point is in
	for (int n = 0; n < 3; n++){
		box_coords[n] = (point[n] - mins[n]) / width;
	}

	//Looping through shell of cube
	for (int a = fmax(0, box_coords[0] - shell_level); a < lengths[0] && a <= box_coords[0] + shell_level; a++){
		for (int b = fmax(0, box_coords[1] - shell_level); b < lengths[1] && b <= box_coords[1] + shell_level; b++){
			int c = fmax(0, box_coords[2] - shell_level);

			bool c_ran = false;
			while (c < lengths[2] && c <= box_coords[2] + shell_level){
				c_ran = true;
				//Finding box index of a, b, c
				int current_box_index = find_box_index(a, b, c);
				//Finding all simplices that intersect that box and adds them to the list of simplices

				for (int n = 0; n < list_of_intersecting_simplices[current_box_index].size(); n++){
					int current_simplex = list_of_intersecting_simplices[current_box_index][n];
					//if (std::find(simplices.begin(), simplices.end(), current_simplex) == simplices.end()){
						simplices.push_back(current_simplex);
					//}
				}
				//If at starting or ending face of cube, increment z normally
				if (a == box_coords[0] - shell_level || a == box_coords[0] + shell_level){
					c++;
				}
				//If at "bottom" or "top" face increment normally
				else if (b == box_coords[1] - shell_level || b == box_coords[1] + shell_level){
					c++;
				}
				//If in the middle of "side" faces, skip the interior
				else{
					if (!c_ran){
						c = fmin(lengths[2], box_coords[2] + shell_level);
					}else{
						c++;
					}
				}
			}
		}
	}	
}

double TriangleCellIntersection::get_width(){
	return width;
}

void TriangleCellIntersection::find_mins_and_maxes(const std::vector<double> & points){
	
	if (points.size() >= 3){
		//Initializing
		for (int n = 0; n < 3; n++){
			mins[n] = points[n];
			maxes[n] = points[n];
		}

		//Looping through rest of points
		for (int n = 3; n < points.size(); n+=3){
			
			//Looping through x, y, and z
			for (int m = 0; m < 3; m++){
				//If the current point is less than the min, update
				if (points[n + m] < mins[m]){
					mins[m] = points[m + n];
				}
				//If the current point is greater than the max, update
				if (points[n + m] > maxes[m]){
					maxes[m] = points[m + n];
				}
			}
		}
	}
}

/*
Find the normal cross product
*/
bool TriangleCellIntersection::vector_normal_to_triangle(const std::vector<double> & triangle, std::vector<double> & normal_vector){

	std::vector<double> a(DIM3);
	std::vector<double> b(DIM3);

	//Finding vectors of edges of triangles
	for (int n = 0; n < 3; n++){
		a[n] = (triangle[n + 3] - triangle[n + 0]);
	}
	for (int n = 0; n < 3; n++){
		b[n] = (triangle[n + 6] - triangle[n + 0]);
	}

	//Solving cross product
	double term1 = a[1] * b[2] - a[2] * b[1];
	double term2 = -(a[0] * b[2] - a[2] * b[0]);
	double term3 = a[0] * b[1] - a[1] * b[0];

	double magnitude = sqrt(term1*term1 + term2*term2 + term3*term3);

	if (magnitude == 0){
		//Returning false because the magnitude of the cross product is zero
		return false;
	}
	else{
		//Adding to vector
		normal_vector[0] = (term1);
		normal_vector[1] = (term2);
		normal_vector[2] = (term3);
	}

	//Returning success
	return true;
}

int TriangleCellIntersection::find_box_index(int x, int y, int z){

	int index = x * lengths[1] * lengths[2] + y * lengths[2] + z;
	return index;

}
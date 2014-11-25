#include "stdafx.h"
#include "TriangleCellIntersections.h"
#include <math.h>
#include <iostream>

using namespace std;

TriangleCellIntersection::TriangleCellIntersection(std::vector<double> points, std::vector<int> simplices, double box_width, double & average_boxes){

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
	for (int n = 0; n < num_boxes; n++){
		std::vector<int> temp_vector;
		list_of_intersecting_simplices.push_back(temp_vector);
	}

	//Find the boxes the simplices intersect
	for (int n = 0; n < simplices.size(); n += 3){
		//Storing the points of a triangle
		std::vector<double> triangle_points;

		//Looping through the three simplices
		for (int m = 0; m < 3; m++){
			//Looping throught the three points of each simplice
			for (int p = 0; p < 3; p++){
				triangle_points.push_back(points[3 * simplices[n + m] + p]);
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

	average_boxes /= (simplices.size() / 3);
}

void TriangleCellIntersection::find_boxes_triangle_intersects(std::vector<double> triangle_points, std::vector<int> & boxes_intersected){

	boxes_intersected.clear();

	double triangle_mins[3];
	double triangle_maxes[3];
	std::vector<double> normal_vector;
	
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


		
		//Recording intersected boxes
		for (int a = lower[0]; a <= upper[0] && a < lengths[0]; a++){
			for (int b = lower[1]; b <= upper[1] && b < lengths[1]; b++){
				for (int c = lower[2]; c <= upper[2] && c < lengths[2]; c++){

					//Records the value of the dot product of the vector normal to the plane and the vector between cube vertices and a triangle point
					std::vector<double> plane_values_vertices;
					for (int x = 0; x <= 1; x++){
						for (int y = 0; y <= 1; y++){
							for (int z = 0; z <= 1; z++){
								double value = constant;
								value += (normal_vector[0] * ((a + x) * width + mins[0]));
								value += (normal_vector[1] * ((b + y) * width + mins[1]));
								value += (normal_vector[2] * ((c + z) * width + mins[2]));
								plane_values_vertices.push_back(value);
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
						boxes_intersected.push_back(find_box_index(a, b, c));
					}
				}
			}
		}
	}
}

void TriangleCellIntersection::find_simplices_in_shell_distance_from_cube(std::vector<double> point, int shell_level, std::vector<int> & simplices){

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
					if (std::find(simplices.begin(), simplices.end(), current_simplex) == simplices.end()){
						simplices.push_back(current_simplex);
					}
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

void TriangleCellIntersection::find_mins_and_maxes(std::vector<double> points){
	
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
bool TriangleCellIntersection::vector_normal_to_triangle(std::vector<double> triangle, std::vector<double> & normal_vector){

	std::vector<double> a;
	std::vector<double> b;

	//Finding vectors of edges of triangles
	for (int n = 0; n < 3; n++){
		a.push_back(triangle[n + 3] - triangle[n + 0]);
	}
	for (int n = 0; n < 3; n++){
		b.push_back(triangle[n + 6] - triangle[n + 0]);
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
		normal_vector.clear();
		normal_vector.push_back(term1);
		normal_vector.push_back(term2);
		normal_vector.push_back(term3);
	}

	//Returning success
	return true;
}

int TriangleCellIntersection::find_box_index(int x, int y, int z){

	int index = x * lengths[1] * lengths[2] + y * lengths[2] + z;
	return index;

}
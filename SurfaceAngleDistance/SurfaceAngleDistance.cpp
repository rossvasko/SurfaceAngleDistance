#include "stdafx.h"
#include <vector>
#include <iostream>
#include "SurfaceAngleDistance.h"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>


using namespace std;


	/*
	Computes the minimum angle difference between two surfaces
	*/
double surface_angle_distance::surface_angle_distance_extended(std::vector<double> s1_points, std::vector<int> s1_simplices, std::vector<double> s2_points, std::vector<int> s2_simplices, double epsilon, double delta, double width, std::vector<double> & s1_angle_distances, std::vector<double> & s2_angle_distances, bool num_simplices_histogram, bool area_simplices_histogram, double print_area_angle_cutoff, std::string prefix){
		
		s1_angle_distances.clear();
		s2_angle_distances.clear();

		std::vector<segmentedTriangle> segmented_triangles_s1;
		std::vector<segmentedTriangle> segmented_triangles_s2;
		
		create_segmented_triangles(s1_points, s1_simplices, segmented_triangles_s1, delta);
		create_segmented_triangles(s2_points, s2_simplices, segmented_triangles_s2, delta);
		
		double total_area_s1 = 0;
		double total_edge_length_s1 = 0;
		int number_of_segmented_triangles_s1 = 0;

		for (int n = 0; n < segmented_triangles_s1.size(); n++){
			if (segmented_triangles_s1[n].valid){
				std::vector<double> vec1;
				std::vector<double> vec2;
				std::vector<double> vec3;

				create_vectors_from_triangle(segmented_triangles_s1[n].boundary_points, vec1, vec2, vec3);
				total_area_s1 += .5 * magnitude_cross_product(vec1, vec2);
				number_of_segmented_triangles_s1 += segmented_triangles_s1[n].simplices.size() / 3;

				for (int m = 0; m < segmented_triangles_s1[n].simplices.size(); m += 3){
					std::vector<double> temp_triangle;
					for (int p = 0; p < 3; p++){
						temp_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + p] + 0]);
						temp_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + p] + 1]);
						temp_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + p] + 2]);
					}
					total_edge_length_s1 += edge_length(temp_triangle);
				}
			}
		}
		double average_area_s1 = total_area_s1 / number_of_segmented_triangles_s1;
		double average_edge_length_s1 = total_edge_length_s1 / number_of_segmented_triangles_s1;

		double total_area_s2 = 0;
		double total_edge_length_s2 = 0;
		double number_of_segmented_triangles_s2 = 0;

		for (int n = 0; n < segmented_triangles_s2.size(); n++){
			if (segmented_triangles_s2[n].valid){
				std::vector<double> vec1;
				std::vector<double> vec2;
				std::vector<double> vec3;
				create_vectors_from_triangle(segmented_triangles_s2[n].boundary_points, vec1, vec2, vec3);
				total_area_s2 += .5 * magnitude_cross_product(vec1, vec2);
				number_of_segmented_triangles_s2 += segmented_triangles_s2[n].simplices.size() / 3;


				for (int m = 0; m < segmented_triangles_s2[n].simplices.size(); m += 3){
					std::vector<double> temp_triangle;
					for (int p = 0; p < 3; p++){
						temp_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + p] + 0]);
						temp_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + p] + 1]);
						temp_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + p] + 2]);
					}
					total_edge_length_s2 += edge_length(temp_triangle);
				}
			}
		}
		double average_area_s2 = total_area_s2 / number_of_segmented_triangles_s2;
		double average_edge_length_s2 = total_edge_length_s2 / number_of_segmented_triangles_s2;

		double width_s1 = 0;
		double width_s2 = 0;

		if (width > 0){
			width_s1 = width;
			width_s2 = width;
		}else{
			width_s1 = 3 * sqrt(average_area_s1 * 4 / sqrt(3));
			width_s2 = 3 * sqrt(average_area_s2 * 4 / sqrt(3));
		}

		std::cout << std::endl;
		std::cout << "Average area of the segmented triangles in surface 1 is: " << average_area_s1 << std::endl;
		std::cout << "Average edge length of the segmented triangles in surface 1 is: " << average_edge_length_s1 << std::endl;
		std::cout << "Number of segmented triangles in surface 1 is: " << number_of_segmented_triangles_s1 << std::endl;
		std::cout << "Edge length of surface 1 regular grid: " << width_s1 << std::endl;

		double average_number_of_boxes = 0;
		TriangleCellIntersection s1_grid(s1_points, s1_simplices, width_s1, average_number_of_boxes);
		std::cout << "Average number of boxes intersected in surface 1: " << average_number_of_boxes << std::endl << std::endl;

		std::cout << "Average area of the segmented triangles in surface 2 is: " << average_area_s2 << std::endl;
		std::cout << "Average edge length of the segmented triangles in surface 2 is: " << average_edge_length_s2 << std::endl;
		std::cout << "Number of segmented triangles in surface 2 is: " << number_of_segmented_triangles_s2 << std::endl;
		std::cout << "Edge length of surface 2 regular grid: " << width_s2 << std::endl;

		TriangleCellIntersection s2_grid(s2_points, s2_simplices, width_s2, average_number_of_boxes);
		std::cout << "Average number of boxes intersected in surface 2: " << average_number_of_boxes << std::endl << std::endl;

		//std::cout << "0%" << std::endl; 
		surface_to_surface_extended(segmented_triangles_s1, segmented_triangles_s2, &s2_grid, epsilon, 1);
		//std::cout << "50%" << std::endl;
		surface_to_surface_extended(segmented_triangles_s2, segmented_triangles_s1, &s1_grid, epsilon, 2);
		//std::cout << "100%" << std::endl;
		
		double max_angle = -1;

		//Finds the maximum angle
		for (int x = 0; x < segmented_triangles_s1.size(); x++){
			double max_of_simplex = -1;
			if (segmented_triangles_s1[x].valid){
				for (int y = 0; y < segmented_triangles_s1[x].angle_diff_of_closest.size(); y++){
					if (max_angle < segmented_triangles_s1[x].angle_diff_of_closest[y]){
						max_angle = segmented_triangles_s1[x].angle_diff_of_closest[y];
					}
					if (max_of_simplex < segmented_triangles_s1[x].angle_diff_of_closest[y]){
						max_of_simplex = segmented_triangles_s1[x].angle_diff_of_closest[y];
					}
				}
				s1_angle_distances.push_back(max_of_simplex);
			}
		}
		for (int x = 0; x < segmented_triangles_s2.size(); x++){
			double max_of_simplex = -1;
			if (segmented_triangles_s2[x].valid){
				for (int y = 0; y < segmented_triangles_s2[x].angle_diff_of_closest.size(); y++){
					if (max_angle < segmented_triangles_s2[x].angle_diff_of_closest[y]){
						max_angle = segmented_triangles_s2[x].angle_diff_of_closest[y];
					}
					if (max_of_simplex < segmented_triangles_s2[x].angle_diff_of_closest[y]){
						max_of_simplex = segmented_triangles_s2[x].angle_diff_of_closest[y];
					}
				}
				s2_angle_distances.push_back(max_of_simplex);
			}
		}

		if (num_simplices_histogram){
			create_num_simplices_histogram(s1_angle_distances, s2_angle_distances, prefix);
		}

		if (area_simplices_histogram){
			create_area_simplices_histogram(segmented_triangles_s1, segmented_triangles_s2, prefix);
		}

		if (print_area_angle_cutoff >= 0){
			print_area_simplices_under_cutoff(segmented_triangles_s1, segmented_triangles_s2, print_area_angle_cutoff);
		}
		
		//Returns the max angle
		return max_angle;
	}


	void surface_angle_distance::print_area_simplices_under_cutoff(std::vector<segmentedTriangle> segmented_triangles_s1, std::vector<segmentedTriangle> segmented_triangles_s2, double cutoff){

		double total_area = 0;
		double area_below_cutoff = 0;

		for (int n = 0; n < segmented_triangles_s1.size(); n++){

			if (segmented_triangles_s1[n].valid){
				std::vector<double> vec1;
				std::vector<double> vec2;
				std::vector<double> vec3;
				create_vectors_from_triangle(segmented_triangles_s1[n].boundary_points, vec1, vec2, vec3);
				total_area += .5 * magnitude_cross_product(vec1, vec2);

				for (int m = 0; m < segmented_triangles_s1[n].simplices.size(); m += 3){
					if (segmented_triangles_s1[n].angle_diff_of_closest[m/3] <= cutoff){
						std::vector<double> vec1_seg;
						std::vector<double> vec2_seg;
						std::vector<double> vec3_seg;

						std::vector<double> new_triangle;
						for (int x = 0; x < 3; x++){
							new_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 0]);
							new_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 1]);
							new_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 2]);
						}

						create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
						double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

						area_below_cutoff += increment_area;
					}
				}
			}
		}

		std::cout << area_below_cutoff / total_area << " of the area (" << area_below_cutoff << "/" << total_area << ") is below the cutoff of " << cutoff << ".\n";


		total_area = 0;
		area_below_cutoff = 0;

		for (int n = 0; n < segmented_triangles_s2.size(); n++){
			if (segmented_triangles_s2[n].valid){
				std::vector<double> vec1;
				std::vector<double> vec2;
				std::vector<double> vec3;
				create_vectors_from_triangle(segmented_triangles_s2[n].boundary_points, vec1, vec2, vec3);
				total_area += .5 * magnitude_cross_product(vec1, vec2);

				for (int m = 0; m < segmented_triangles_s2[n].simplices.size(); m+=3){
					if (segmented_triangles_s2[n].angle_diff_of_closest[m / 3] <= cutoff){

						std::vector<double> vec1_seg;
						std::vector<double> vec2_seg;
						std::vector<double> vec3_seg;

						std::vector<double> new_triangle;
						for (int x = 0; x < 3; x++){
							new_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 0]);
							new_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 1]);
							new_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 2]);
						}

						create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
						double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

						area_below_cutoff += increment_area;
					}
				}
			}
		}

		std::cout << area_below_cutoff / total_area << " of the area (" << area_below_cutoff << "/" << total_area << ") is below the cutoff of " << cutoff << ".\n";

	}

	/*
	Creates a histogram of the area of simplices in range of angles
	*/
	void surface_angle_distance::create_area_simplices_histogram(std::vector<segmentedTriangle> segmented_triangles_s1, std::vector<segmentedTriangle> segmented_triangles_s2, std::string prefix){
	
		std::string histo1_name = prefix + "area_simplices_histogram_1.dat";
		std::ofstream histo1(histo1_name);

		double bins[91];

		for (int n = 0; n < 91; n++){
			bins[n] = 0;
		}

		double total_area = 0;

		for (int n = 0; n < segmented_triangles_s1.size(); n++){
			if (segmented_triangles_s1[n].valid){
				std::vector<double> vec1;
				std::vector<double> vec2;
				std::vector<double> vec3;
				create_vectors_from_triangle(segmented_triangles_s1[n].boundary_points, vec1, vec2, vec3);
				total_area += .5 * magnitude_cross_product(vec1, vec2);
				
				for (int m = 0; m < segmented_triangles_s1[n].simplices.size(); m+=3){

					int index = (int)floor(segmented_triangles_s1[n].angle_diff_of_closest[m/3]);
					std::vector<double> vec1_seg;
					std::vector<double> vec2_seg;
					std::vector<double> vec3_seg;

					std::vector<double> new_triangle;
					for (int x = 0; x < 3; x++){
						new_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 0]);
						new_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 1]);
						new_triangle.push_back(segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 2]);
					}

					create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
					double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

					bins[index] += increment_area;
				}
			}
		}

		bins[89] += bins[90];

		histo1 << "#Angle\Normalized area of simplices\n";
		for (int n = 0; n < 90; n++){
			histo1 << (n + 1) << "\t" << (bins[n] / total_area) << "\n";
		}

		histo1.close();


		std::string histo2_name = prefix + "area_simplices_histogram_2.dat";
		std::ofstream histo2(histo2_name);

		bins[91];

		for (int n = 0; n < 91; n++){
			bins[n] = 0;
		}

		total_area = 0;

		for (int n = 0; n < segmented_triangles_s2.size(); n++){
			if (segmented_triangles_s2[n].valid){
				std::vector<double> vec1;
				std::vector<double> vec2;
				std::vector<double> vec3;
				create_vectors_from_triangle(segmented_triangles_s2[n].boundary_points, vec1, vec2, vec3);
				total_area += .5 * magnitude_cross_product(vec1, vec2);

				for (int m = 0; m < segmented_triangles_s2[n].simplices.size(); m+=3){

					int index = (int)floor(segmented_triangles_s2[n].angle_diff_of_closest[m/3]);
					std::vector<double> vec1_seg;
					std::vector<double> vec2_seg;
					std::vector<double> vec3_seg;

					std::vector<double> new_triangle;
					for (int x = 0; x < 3; x++){
						new_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 0]);
						new_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 1]);
						new_triangle.push_back(segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 2]);
					}

					create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
					double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

					bins[index] += increment_area;
				}
			}
		}

		bins[89] += bins[90];
		
		histo2 << "#Angle\Normalized area of simplices\n";
		for (int n = 0; n < 90; n++){
			histo2 << (n + 1) << "\t" << (bins[n] / total_area) << "\n";
		}

		histo2.close();
	}

	/*
	Creates a histogram of the number of simplices in ranges of angles
	*/
	void surface_angle_distance::create_num_simplices_histogram(std::vector<double> angles1, std::vector<double> angles2, std::string prefix){

		std::string histo1_name = prefix + "num_simplices_histogram_1.dat";
		std::ofstream histo1(histo1_name);

		int bins[91];
		for (int n = 0; n < 91; n++){
			bins[n] = 0;
		}
		for (int n = 0; n < angles1.size(); n++){
			int index = floor(angles1[n]);
			bins[index]++;
		}
		bins[89] += bins[90];

		histo1 << "#Angle\tNumber of simplices\n";
		for (int n = 0; n < 90; n++){
			histo1 << (n + 1) << "\t" << bins[n] << "\n";
		}

		histo1.close();

		std::string histo2_name = prefix + "num_simplices_histogram_2.dat";
		std::ofstream histo2(histo2_name);

		for (int n = 0; n < 91; n++){
			bins[n] = 0;
		}

		for (int n = 0; n < angles2.size(); n++){
			int index = floor(angles2[n]);
			bins[index]++;
		}
		bins[89] += bins[90];

		histo2 << "#Angle\tNumber of simplices\n";
		for (int n = 0; n < 90; n++){
			histo2 << (n + 1) << "\t" << bins[n] << "\n";
		}

		histo2.close();

	}

	/*
	Calculates the extended surface to surface information from triangles1 to triangles2
	*/
	void surface_angle_distance::surface_to_surface_extended(std::vector<segmentedTriangle> & triangles1, std::vector<segmentedTriangle> triangles2, TriangleCellIntersection * grid, double epsilon, int num_surface){
		
		//Loops through all triangles in surface 1
		for (int x = 0; x < triangles1.size(); x++){
			//std::cout << x << " " << triangles1.size() << "\n";

			//std::cout << (int)((100.0 * (x) / triangles1.size()) / 10) << " " << (int)((100.0 * (x - 1) / triangles1.size())) << std::endl;

			if (num_surface == 1){
				if ((int)((100.0 * (x) / triangles1.size()) / 10) >(int)((100.0 * (x - 1) / triangles1.size() / 10))){
					//std::cout << (int)(100.0 * (x) / triangles1.size() / 10) * 5 << "%" << std::endl;
				}
			}
			else{
				if ((int)((100.0 * (x) / triangles1.size()) / 10) >(int)((100.0 * (x - 1) / triangles1.size() / 10))){
					//std::cout << 50 + (int)(100.0 * (x) / triangles1.size() / 10) * 5 << "%" << std::endl;
				}
			}

			segmentedTriangle * current_triangle = &triangles1[x];

			//Compute info if the current triangle is valid
			if (current_triangle->valid){

				//Loops through all segmented trianlges in the current triangle
				for (int y = 0; y < current_triangle->simplices.size(); y += 3){
					double angle_diff = 360;
					double distance = 0;
					int closest = -1;
					std::vector<double> closest_point;
					std::vector<double> temporary_point;
					std::vector<double> center_point;
					std::vector<double> triangle;

					//Prepares a triangle
					for (int z = y; z < y + 3; z++){
						triangle.push_back(current_triangle->points[3 * current_triangle->simplices[z] + 0]);
						triangle.push_back(current_triangle->points[3 * current_triangle->simplices[z] + 1]);
						triangle.push_back(current_triangle->points[3 * current_triangle->simplices[z] + 2]);
					}

					//Computes center of a triangle
					center_of_triangle(triangle, center_point);

					double nearest_distance = DBL_MAX;

					//Finds the nearest point to the center of the triangle
					point_to_surface_nearest_point(center_point, triangles2, current_triangle->normal, nearest_distance, temporary_point, grid);

					//Computes the minimum angle difference
					compute_angle_difference(center_point, triangles2, current_triangle->normal, distance, closest, angle_diff, closest_point, temporary_point, epsilon, grid);

					//Stores information to the current triangle
					current_triangle->angle_diff_of_closest.push_back(angle_diff);
					current_triangle->closest_distance.push_back(distance);
					current_triangle->closest.push_back(closest);
					current_triangle->closest_point.push_back(closest_point[0]);
					current_triangle->closest_point.push_back(closest_point[1]);
					current_triangle->closest_point.push_back(closest_point[2]);
				}				
			}
		}
	}

	/*
	Computes the minimum angle difference between a point with a predefined normal to another surface within the distance epsilon from another point
	*/
	void surface_angle_distance::compute_angle_difference(std::vector<double> point, std::vector<segmentedTriangle> surface, std::vector<double> normal_to_triangle, double & distance, int & minimum_angle_index, double & minimum_angle_diff, std::vector<double> & minimum_angle_point, std::vector<double> point_to_compute_distance_from, double epsilon, TriangleCellIntersection * grid){

		//Simplices to check
		std::vector<int> simplices_to_check;
		int current_shell = 0;

		//Make a vector of all simplices possibly in the epsilon range
		while (epsilon >= grid->get_width()*(current_shell - 1)){
			std::vector<int> temp_vector;
			grid->find_simplices_in_shell_distance_from_cube(point_to_compute_distance_from, current_shell, temp_vector);
			
			for (int n = 0; n < temp_vector.size(); n++){
				if (std::find(simplices_to_check.begin(), simplices_to_check.end(), temp_vector[n]) == simplices_to_check.end()){
					simplices_to_check.push_back(temp_vector[n]);
				}
			}
			current_shell++;
		}

		//Loop through all simplices that could possibly meet the epsilon distance
		for (int n = 0; n < simplices_to_check.size(); n++){
		
			segmentedTriangle current_triangle = surface[simplices_to_check[n]];

			//If the surface is valid continue
			if (current_triangle.valid){
				double distance_from_point = 0;
				std::vector<double> temp_vector;

				//Computes the distance from a point to the triangle
				point_to_triangle(point_to_compute_distance_from, current_triangle, distance_from_point, temp_vector);

				//If the distance is less than the specified epsilon, continue
				if (distance_from_point < epsilon){
					//Calculates the angle between the normal and the normal of the current triangle
					double currentAngle = angle_of_vectors(normal_to_triangle, current_triangle.normal);
					//If the current angle is less than the minimum, update the new information
					if (currentAngle < minimum_angle_diff){
						minimum_angle_diff = currentAngle;
						minimum_angle_index = simplices_to_check[n];
						point_to_triangle(point, current_triangle, distance, minimum_angle_point);
					}
				}
			}
		}

		/* Old code


		//Loops through all of the triangles that make up a surface
		for (int x = 0; x < surface.size(); x++){
			
			//If the surface is valid continue
			if (surface[x].valid){
				double distance_from_point = 0;
				std::vector<double> temp_vector;

				//Computes the distance from a point to the triangle
				point_to_triangle(point_to_compute_distance_from, surface[x], distance_from_point, temp_vector);

				//If the distance is less than the specified epsilon, continue
				if (distance_from_point < epsilon){
					//Calculates the angle between the normal and the normal of the current triangle
					double currentAngle = angle_of_vectors(normal_to_triangle, surface[x].normal);
					//If the current angle is less than the minimum, update the new information
					if (currentAngle < minimum_angle_diff){
						minimum_angle_diff = currentAngle;
						minimum_angle_index = x;
						point_to_triangle(point, surface[x], distance, minimum_angle_point);
					}
				}
			}
		}
		*/
	}

	/*
	Computes information for point to surface
	*/
	void surface_angle_distance::point_to_surface_nearest_point(std::vector<double> point, std::vector<segmentedTriangle> surface, std::vector<double> normal_to_triangle, double & distance, std::vector<double> & closest_point, TriangleCellIntersection * grid){
	
		distance = DBL_MAX;
		int current_shell = 0;
		std::vector<int> checked;

		//Loop until it is not possible to find a smaller distance
		while (distance >= grid->get_width() * (current_shell - 1)){
			std::vector<int> simplices_in_shell;
			//Find simplices in current shell

			grid->find_simplices_in_shell_distance_from_cube(point, current_shell, simplices_in_shell);

			//Compute distance of simplices in shell
			for (int n = 0; n < simplices_in_shell.size(); n++){
				if (std::find(checked.begin(), checked.end(), simplices_in_shell[n]) == checked.end()){
					segmentedTriangle current_triangle = surface[simplices_in_shell[n]];
					if (current_triangle.valid){

						double temp_distance = 0;
						std::vector<double> temp_vector;

						//Calculate the distance from the point to the trianlge
						point_to_triangle(point, current_triangle, temp_distance, temp_vector);

						//If the current distance is less than the current minimum, update the information
						if (temp_distance < distance){
							distance = temp_distance;
							closest_point.clear();
							for (int x = 0; x < 3; x++){
								closest_point.push_back(temp_vector[x]);
							}
						}
					}
					checked.push_back(simplices_in_shell[n]);
				}
			}
			//std::cout << current_shell << " " << distance << std::endl;
			current_shell++;
		}

		/* Old code

		distance = DBL_MAX;

		//Loops through the triangles that make up the surface
		for (int x = 0; x < surface.size(); x++){
		
			//If the triangle is valid, continue
			if (surface[x].valid){

				double temp_distance = 0;
				std::vector<double> temp_vector;

				//Calculate the distance from the point to the trianlge
				point_to_triangle(point, surface[x], temp_distance, temp_vector);

				//If the current distance is less than the current minimum, update the information
				if (temp_distance < distance){
					distance = temp_distance;
					closest_point.clear();
					for (int x = 0; x < 3; x++){
						closest_point.push_back(temp_vector[x]);
					}
				}
			}
		}
		*/
	}

	void surface_angle_distance::point_to_triangle(std::vector<double> point, segmentedTriangle triangle, double & distance, std::vector<double> & closest_point){
		
		std::vector<double> vec_a;
		std::vector<double> vec_b;
		std::vector<double> vec_c;

		//Creates vectors of the triangle points

		for (int x = 0; x < 3; x++){
			vec_a.push_back(triangle.boundary_points[x]);
		}
		for (int x = 3; x < 6; x++){
			vec_b.push_back(triangle.boundary_points[x]);
		}
		for (int x = 6; x < 9; x++){
			vec_c.push_back(triangle.boundary_points[x]);
		}

		std::vector<double> projection_vector;

		//Find a vector from the point to a point on the plane of the triangle
		std::vector<double> point_to_plane;
		vector_point_to_point(point, vec_a, point_to_plane);

		//Projects the point on to the normal vector
		vector_projection(point_to_plane, triangle.normal, projection_vector);

		//Calculates the distance of the projection vector
		distance = vector_magnitude(projection_vector);

		//Finds a point projected onto the plane
		std::vector<double> projected_point;
		for (int x = 0; x < 3; x++){
			projected_point.push_back(point[x] + projection_vector[x]);
		}

		/*
		Follows the algorithm described in mesh paper
		*/
		if (dot_product(projected_point, triangle.normal_ab) > dot_product(vec_a, triangle.normal_ab)){
		//Region 1

			//Calculates the distance from the project point to line ab
			point_to_line(projected_point, vec_a, vec_b, closest_point);

			//Calculates the distance from the point on the line to the original point
			distance = point_to_point_distance(point, closest_point);
		}
		else if (dot_product(projected_point, triangle.normal_bc) > dot_product(vec_b, triangle.normal_bc)){
		//Region 2 or 3

			std::vector<double> vector_cb;
			std::vector<double> vector_cp;

			//Creates a vector from c to b
			vector_point_to_point(vec_c, vec_b, vector_cb);

			//Creates a vector from c to the original point
			vector_point_to_point(vec_c, projected_point, vector_cp);

			//Comparing the dot product of cb and cp to 0
			if (dot_product(vector_cb, vector_cp) >= 0){
				//If true, region 2

				//Calculates the distance from the project point to line bc
				point_to_line(projected_point, vec_b, vec_c, closest_point);
			}else{
				//If false, region 3

				//Calculates the distance from the project point to line ac
				point_to_line(projected_point, vec_a, vec_c, closest_point);
			}

			//Calculates the distance from the point on the line to the original point
			distance = point_to_point_distance(point, closest_point);
		}
		else if (dot_product(projected_point, triangle.normal_ac) > dot_product(vec_c, triangle.normal_ac)){
			//Region 4

			//Calculates the distance from the project point to line ac
			point_to_line(projected_point, vec_a, vec_c, closest_point);

			//Calculates the distance from the point on the line to the original point
			distance = point_to_point_distance(point, closest_point);
		}
		else{
			//The point is in the triangle and the projected point is the closest point
			closest_point.clear();
			
			for (int x = 0; x < 3; x++){
				closest_point.push_back(projected_point[x]);
			}
		}
	}

	void surface_angle_distance::point_to_line(std::vector<double> point, std::vector<double> line_a, std::vector<double> line_b, std::vector<double> & closest_point){
		

		closest_point.clear();
		//Creates vectors for point a to the point and from point a to point b
		std::vector<double> a_to_point;
		std::vector<double> line_ab;

		vector_point_to_point(line_a, point, a_to_point);
		vector_point_to_point(line_a, line_b, line_ab);

		if (dot_product(a_to_point, line_ab) < 0){

			//Point a is the closest point
			for (int x = 0; x < 3; x++){
				closest_point.push_back(line_a[x]);
			}
		}
		else if (pow(vector_magnitude(line_ab), 2) < dot_product(a_to_point, line_ab)){

			//Point b is the closest point
			for (int x = 0; x < 3; x++){
				closest_point.push_back(line_b[x]);
			}
		}
		else{

			//A point on the line is the closest point
			std::vector<double> projected_vector;

			//Projects the vector a_to_point onto line_ab
			vector_projection(a_to_point, line_ab, projected_vector);

			//Solves for the vector from the point on the line to the original point
			std::vector<double> line_to_point;
			for (int x = 0; x < 3; x++){
				line_to_point.push_back(a_to_point[x] - projected_vector[x]);
			}

			//Stores in the closest point
			for (int x = 0; x < 3; x++){
				closest_point.push_back(point[x] - line_to_point[x]);
			}
		}
	}

	void surface_angle_distance::vector_point_to_point(std::vector<double> a, std::vector<double> b, std::vector<double> & ab){
		
		ab.clear();

		//Calculates the differences of the terms of the vectors
		for (int x = 0; x < a.size(); x++){
			ab.push_back(b[x] - a[x]);
		}
	}

	/*
	Takes triangles and creates a vector of the triangles segmented
	*/
	void surface_angle_distance::create_segmented_triangles(std::vector<double> s_points, std::vector<int> s_simplices, std::vector<segmentedTriangle> & segmented_triangles, double delta){
		
		segmented_triangles.clear();

		//Loops through all simplices
		for (int n = 0; n < s_simplices.size(); n += 3){
			std::vector<double> triangle_points;

			//Creates a triangle of the points in the simplices
			for (int m = n; m < n + 3; m++){
				triangle_points.push_back(s_points[3 * s_simplices[m] + 0]);
				triangle_points.push_back(s_points[3 * s_simplices[m] + 1]);
				triangle_points.push_back(s_points[3 * s_simplices[m] + 2]);
			}

			//Vector normal to the triangle
			std::vector<double> normal_unit_vector;

			//Vector a to b
			std::vector<double> vec1;
			//Vector c to a
			std::vector<double> vec2;
			//Vector b to c
			std::vector<double> vec3;

			//The following vectors are also in the plane
			//Vector normal to ab
			std::vector<double> normal_vec1;
			//Vector normal to ca
			std::vector<double> normal_vec2;
			//Vector normal to bc
			std::vector<double> normal_vec3;

			//Creates the vectors
			create_vectors_from_triangle(triangle_points, vec1, vec2, vec3);
			
			//Variable to keep track of the triangle being added
			bool added = false;

			//Computes normal 
			if (normal_unit_vector_to_triangle(triangle_points, normal_unit_vector)){
				if (normal_cross_product(normal_unit_vector, vec1, normal_vec1)){
					if (normal_cross_product(normal_unit_vector, vec2, normal_vec2)){
						if (normal_cross_product(normal_unit_vector, vec3, normal_vec3)){

							std::vector<double> triangle_points_copy(triangle_points);

							//Segments triangle to the edge length delta
							segment_triangle(triangle_points_copy, delta);

							//Prepares the segmented simplices
							std::vector<int> segmented_simplices;

							for (int a = 0; a < triangle_points_copy.size() / 3; a++){
								segmented_simplices.push_back(a);
							}

							//Compacts the simplices so that there is no "redundant" data
							compact_simplices(segmented_simplices, triangle_points_copy);

							segmentedTriangle segmented_triangle;

							//Sets the triangle validity to true because
							segmented_triangle.valid = true;

							//Sets points to the newly computed segmentations
							segmented_triangle.points = triangle_points_copy;
							segmented_triangle.simplices = segmented_simplices;

							//Sets the normal vector to the normal unit vector computed
							segmented_triangle.normal = normal_unit_vector;

							//This information is needed for later algorithms
							//If any of these are true, the triangle is obtuse and the angle at point c is obtuse

							if (dot_product(vec1, vec2) > 0){

								/*
								New assignment of normals
									normal_ab = normal_bc'
									normal_ac = normal_ab'
									normal_bc = normal_ac'
								*/

								segmented_triangle.normal_ab = normal_vec3;
								segmented_triangle.normal_ac = normal_vec1;
								segmented_triangle.normal_bc = normal_vec2;

								/*
								New assignment of points
									a = b'
									b = c'
									c = a'
								*/

								for (int x = 3; x < 9; x++){
									segmented_triangle.boundary_points.push_back(triangle_points[x]);
								}
								for (int x = 0; x < 3; x++){
									segmented_triangle.boundary_points.push_back(triangle_points[x]);
								}


							}
							else if (dot_product(vec1, vec3) > 0){

								/*
								New assignment of normals
								normal_ab = normal_ac'
								normal_ac = normal_bc'
								normal_bc = normal_ab'
								*/

								segmented_triangle.normal_ab = normal_vec2;
								segmented_triangle.normal_ac = normal_vec3;
								segmented_triangle.normal_bc = normal_vec1;

								/*
								New assignment of points
								a = c'
								b = a'
								c = b'
								*/

								for (int x = 6; x < 9; x++){
									segmented_triangle.boundary_points.push_back(triangle_points[x]);
								}
								for (int x = 0; x < 3; x++){
									segmented_triangle.boundary_points.push_back(triangle_points[x]);
								}
								for (int x = 3; x < 6; x++){
									segmented_triangle.boundary_points.push_back(triangle_points[x]);
								}


							}
							else{
								/*
								New assignment of normals
								normal_ab = normal_ab'
								normal_ac = normal_ac'
								normal_bc = normal_bc'
								*/

								segmented_triangle.normal_ab = normal_vec1;
								segmented_triangle.normal_ac = normal_vec2;
								segmented_triangle.normal_bc = normal_vec3;

								/*
								New assignment of points
								a = a'
								b = b'
								c = c'
								*/

								for (int x = 0; x < 9; x++){
									segmented_triangle.boundary_points.push_back(triangle_points[x]);
								}


							}
							//Add to the list of segmented triangles
							segmented_triangles.push_back(segmented_triangle);
							added = true;
						}
					}
				}
			}

			if(!added){
				segmentedTriangle segmented_triangle;
				segmented_triangle.valid = false;
				segmented_triangles.push_back(segmented_triangle);
			}

		}
	}

	/*
	Compacts the triangle storage to be less repetitive and use less data
	*/
	void surface_angle_distance::compact_simplices(std::vector<int> & simplices, std::vector<double> & points){
		
		//Loops through all simplices
		for (int x = 0; x < simplices.size(); x++){

			//Stores the current x, y, and z
			double current_x = points[3 * simplices[x] + 0];
			double current_y = points[3 * simplices[x] + 1];
			double current_z = points[3 * simplices[x] + 2];
			
			//Looks for repeats in the rest of the vector
			for (int y = x + 1; y < simplices.size(); y++){
				//If the two simplices aren't pointing to the same point, continue
				if (simplices[x] != simplices[y]){

					//Finds a new x, y, and z
					double new_x = points[3 * simplices[y] + 0];
					double new_y = points[3 * simplices[y] + 1];
					double new_z = points[3 * simplices[y] + 2];

					if (new_x == current_x && new_y == current_y && new_z == current_z){
						//If all of the points equal each other, the information is held twice and can be deleted
						points.erase(points.begin() + 3 * simplices[y]);
						points.erase(points.begin() + 3 * simplices[y]);
						points.erase(points.begin() + 3 * simplices[y]);

						//Shifts remaining elements in the list
						for (int z = y + 1; z < simplices.size(); z++){
							if (simplices[z] > simplices[y]){
								simplices[z] = simplices[z] - 1;
							}
						}

						//Have the yth index equal the xth index
						simplices[y] = simplices[x];
					}
				}
			}
		}
	}

	/*
	Creates a normal unit vector to triangle
	*/
	bool surface_angle_distance::normal_unit_vector_to_triangle(std::vector<double> triangle, std::vector<double> & normal_vector){
		
		//Creating vectors that correspond to the sides of the triangle
		std::vector<double> vec1;
		std::vector<double> vec2;

		//Making two vectors of the sides of the triangles
		for (int n = 0; n < 3; n++){
			vec1.push_back(triangle[n + 6] - triangle[n + 0]);
			vec2.push_back(triangle[n + 3] - triangle[n + 0]);
		}
		
		//Solving cross product
		double term1 = vec1[1] * vec2[2] - vec1[2] * vec2[1];
		double term2 = -(vec1[0] * vec2[2] - vec1[2] * vec2[0]);
		double term3 = vec1[0] * vec2[1] - vec1[1] * vec2[0];

		double magnitude = sqrt(term1*term1 + term2*term2 + term3*term3);


		if (magnitude == 0){
			//If magnitude is zero return false
			return false;
		}else{
			//Normalizes vector
			term1 /= magnitude;
			term2 /= magnitude;
			term3 /= magnitude;

			//Adding to vector
			normal_vector.clear();
			normal_vector.push_back(term1);
			normal_vector.push_back(term2);
			normal_vector.push_back(term3);
		}

		return true;

	}

	/*
	Segments a trianlge into multiple trianlges with edge lengths less than delta
	*/
	void surface_angle_distance::segment_triangle(std::vector<double> & triangle, double delta){
	
		if (edge_length(triangle) > delta){
			//Creating three new triangle
			std::vector<double> t1;
			std::vector<double> t2;
			std::vector<double> t3;
			std::vector<double> t4;

			//Finding midpoints of current edges
			std::vector<double> midpoint01;
			std::vector<double> midpoint02;
			std::vector<double> midpoint12;

			for (int n = 0; n < 3; n++){
				midpoint01.push_back((triangle[n + 0] + triangle[n + 3]) / 2.0);
				midpoint02.push_back((triangle[n + 0] + triangle[n + 6]) / 2.0);
				midpoint12.push_back((triangle[n + 3] + triangle[n + 6]) / 2.0);
			}

			//Creating the 4 new triangles
			for (int n = 0; n < 3; n++){
				t1.push_back(triangle[n + 0]);
				t2.push_back(triangle[n + 3]);
				t3.push_back(triangle[n + 6]);
				t4.push_back(midpoint01[n]);
			}

			for (int n = 0; n < 3; n++){
				t1.push_back(midpoint01[n]);
				t2.push_back(midpoint12[n]);
				t3.push_back(midpoint02[n]);
				t4.push_back(midpoint12[n]);
			}

			for (int n = 0; n < 3; n++){
				t1.push_back(midpoint02[n]);
				t2.push_back(midpoint01[n]);
				t3.push_back(midpoint12[n]);
				t4.push_back(midpoint02[n]);
			}

			//Recursively segmenting the newly created triangles
			segment_triangle(t1, delta);
			segment_triangle(t2, delta);
			segment_triangle(t3, delta);
			segment_triangle(t4, delta);

			//Replacing triangle with the new segmented triangles
			triangle.clear();
			triangle.insert(triangle.end(), t1.begin(), t1.end());
			triangle.insert(triangle.end(), t2.begin(), t2.end());
			triangle.insert(triangle.end(), t3.begin(), t3.end());
			triangle.insert(triangle.end(), t4.begin(), t4.end());

		}

	}

	/*
	Finds distance between two points
	*/
	double surface_angle_distance::point_to_point_distance(std::vector<double> p1, std::vector<double> p2){
		
		double distance = 0;

		//Sums to squares of the differences of the points
		for (int n = 0; n < 3; n++){
			distance += (p1[n] - p2[n])*(p1[n] - p2[n]);
		}

		//Calculates the square root of the sum of squares
		distance = sqrt(distance);

		return distance;

	}

	/*
	Find the center of triangle
	*/
	void surface_angle_distance::center_of_triangle(std::vector<double> triangle, std::vector<double> & point){
		
		point.clear();
		
		for (int x = 0; x < 3; x++){
			//Finds the average of the points
			double average = 0;
			average += triangle[x + 0];
			average += triangle[x + 3];
			average += triangle[x + 6];
			average /= 3.0;

			point.push_back(average);
		}
	}

	/*
	Finds the edge length of a triangle
	*/
	double surface_angle_distance::edge_length(std::vector<double> triangle){
		
		double edge_length = 0;

		//Points that define the triangle
		std::vector<double> p1(triangle.begin() + 0, triangle.begin() + 3);
		std::vector<double> p2(triangle.begin() + 3, triangle.begin() + 6);
		std::vector<double> p3(triangle.begin() + 6, triangle.begin() + 9);

		double max_length = fmax(point_to_point_distance(p1, p2), point_to_point_distance(p1, p3));
		max_length = fmax(max_length, point_to_point_distance(p2, p3));

		//Sum of distance between points
		//edge_length += point_to_point_distance(p1, p2);
		//edge_length += point_to_point_distance(p1, p3);
		//edge_length += point_to_point_distance(p2, p3);

		return max_length;
	}

	void surface_angle_distance::create_vectors_from_triangle(std::vector<double> triangle, std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3){
	
		//Difference from a to b
		for (int x = 0; x < 3; x++){
			vec1.push_back(triangle[x + 3] - triangle[x + 0]);
		}

		//Difference from c to a
		for (int x = 0; x < 3; x++){
			vec2.push_back(triangle[x + 0] - triangle[x + 6]);
		}

		//Difference from b to c
		for (int x = 0; x < 3; x++){
			vec3.push_back(triangle[x + 6] - triangle[x + 3]);
		}

	}

	/*
	Computes the dot product of two vectors
	*/
	double surface_angle_distance::dot_product(std::vector<double> a, std::vector<double> b){
	
		double value = 0;

		//Calculating each term of the dot product
		for (int x = 0; x < a.size(); x++){
			value += a[x] * b[x];
		}

		return value;

	}

	/*
	Calculates the project of vector a on to vector b
	*/
	void surface_angle_distance::vector_projection(std::vector<double> a, std::vector<double> b, std::vector<double> & projection){

		projection.clear();

		//Copying vector b
		for (int x = 0; x < 3; x++){
			projection.push_back(b[x]);
		}
	
		//Dot product of a and b
		double dot = dot_product(a, b);

		//Magnitude of vector b
		double magnitude = vector_magnitude(b);

		//Storing each term of the projection vector proj_b(a) = (a * b)/(|b||b|) * b
		if (magnitude != 0){
			for (int x = 0; x < a.size(); x++){
				projection[x] = (dot * projection[x]) / (magnitude * magnitude);
			}
		}
		else{
			for (int x = 0; x < a.size(); x++){
				projection[x] = 0;
			}
		}
	}

	/*
	Calculates the magnitude of a vector
	*/
	double surface_angle_distance::vector_magnitude(std::vector<double> a){
		
		double magnitude = 0;

		//Summing the squares of each term
		for (int x = 0; x < a.size(); x++){
			magnitude += a[x] * a[x];
		}

		//Taking the square root of the sum
		magnitude = sqrt(magnitude);

		return magnitude;

	}

	/*
	Find the normal cross product
	*/
	bool surface_angle_distance::normal_cross_product(std::vector<double> a, std::vector<double> b, std::vector<double> & normal_vector){

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
			//Normalizing
			term1 /= magnitude;
			term2 /= magnitude;
			term3 /= magnitude;

			//Adding to vector
			normal_vector.clear();
			normal_vector.push_back(term1);
			normal_vector.push_back(term2);
			normal_vector.push_back(term3);
		}

		//Returning success
		return true;
	}

	/*
	Find the cross product
	*/
	double surface_angle_distance::magnitude_cross_product(std::vector<double> a, std::vector<double> b){

		//Solving cross product
		double term1 = a[1] * b[2] - a[2] * b[1];
		double term2 = -(a[0] * b[2] - a[2] * b[0]);
		double term3 = a[0] * b[1] - a[1] * b[0];

		double magnitude = sqrt(term1*term1 + term2*term2 + term3*term3);

		return magnitude;
	}

	/*
	Computes the angle between vectors
	*/
	double surface_angle_distance::angle_of_vectors(std::vector<double> a, std::vector<double> b){
	
		//Dot product of vectors a and b
		double dot = dot_product(a, b);

		//Magnitude of vector a times magnitude of vector b
		double magnitudes = vector_magnitude(a) * vector_magnitude(b);

		if (magnitudes == 0){
			//Returning 90 because a vector is the zero vector
			return 90;
		}
		else{
			//a * b = |a||b|cos(x), acos((a * b)/(|a||b|)) = x
			double x = dot / magnitudes;
			double angle = 0;
			angle = acos(x) * 180 / M_PI;

			if (angle != angle)
				angle = 0;

			//Adjusting the range to 0 - 90
			angle = 90 - abs(angle - 90);

			return angle;
		}
	}


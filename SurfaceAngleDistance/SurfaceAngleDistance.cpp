#include <vector>
#include <iostream>
#include "SurfaceAngleDistance.h"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

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

///////Timing//////
#ifdef _WIN32
#include <time.h>
#include <windows.h>
#endif

#ifdef __unix
#include <time.h>
#include <sys/time.h>

uint64_t getTimeMs(void)
{
	struct timeval tv;

	gettimeofday(&tv, 0);
	return uint64_t(tv.tv_sec) * 1000 + tv.tv_usec / 1000;
}
#endif
///////////////////

using namespace std;
bool ignore_orientation = true;
const int DIM3 = 3;

segmentedTriangle::segmentedTriangle(int size) : angle_diff_of_closest(size), closest_distance(size), closest(size), smallest_distance(size), closest_point(size * 3), boundary_points(DIM3 * 3)
{

}

	/*
	Computes the minimum angle difference between two surfaces
	*/
double surface_angle_distance::surface_angle_distance_extended(const std::vector<double> & s1_points, const std::vector<int> & s1_simplices, const std::vector<double> & s2_points, const std::vector<int> & s2_simplices, double epsilon, double delta, double width, std::vector<double> & s1_angle_distances, std::vector<double> & s2_angle_distances, bool num_simplices_histogram, bool area_simplices_histogram, double print_area_angle_cutoff, std::string prefix, bool orient, bool testing_script_mode, bool percent_done, bool time, bool info){
		

		#ifdef _WIN32
			long start = GetTickCount();
		#elif __unix
			uint64_t start = getTimeMs();
		#endif

		ignore_orientation = orient;

		std::vector<segmentedTriangle> segmented_triangles_s1(s1_simplices.size() / 3);
		std::vector<segmentedTriangle> segmented_triangles_s2(s2_simplices.size() / 3);
		
		#ifdef _WIN32
			long segment_start = GetTickCount();
		#elif __unix
			uint64_t segment_start = getTimeMs();
		#endif

		create_segmented_triangles(s1_points, s1_simplices, segmented_triangles_s1, delta);
		create_segmented_triangles(s2_points, s2_simplices, segmented_triangles_s2, delta);

		#ifdef _WIN32
			long segment_triangles_time = GetTickCount() - segment_start;
			segment_start = GetTickCount();
		#elif __unix
			uint64_t segment_triangles_time = getTimeMs() - segment_start;
			segment_start = getTimeMs();
		#endif

		double total_area_s1 = 0;
		double total_edge_length_s1 = 0;
		int number_of_segmented_triangles_s1 = 0;

		for (int n = 0; n < segmented_triangles_s1.size(); n++){
			if (segmented_triangles_s1[n].valid){
				std::vector<double> vec1(DIM3);
				std::vector<double> vec2(DIM3);
				std::vector<double> vec3(DIM3);

				create_vectors_from_triangle(segmented_triangles_s1[n].boundary_points, vec1, vec2, vec3);
				total_area_s1 += .5 * magnitude_cross_product(vec1, vec2);
				number_of_segmented_triangles_s1 += segmented_triangles_s1[n].simplices.size() / 3;
				total_edge_length_s1 += edge_length(segmented_triangles_s1[n].boundary_points);
			}
		}
		//double average_area_s1 = total_area_s1 / number_of_segmented_triangles_s1;
		double average_area_s1 = total_area_s1 / s1_simplices.size();
		double average_edge_length_s1 = total_edge_length_s1 / s1_simplices.size();

		double total_area_s2 = 0;
		double total_edge_length_s2 = 0;
		double number_of_segmented_triangles_s2 = 0;


		for (int n = 0; n < segmented_triangles_s2.size(); n++){
			if (segmented_triangles_s2[n].valid){
				std::vector<double> vec1(DIM3);
				std::vector<double> vec2(DIM3);
				std::vector<double> vec3(DIM3);
				create_vectors_from_triangle(segmented_triangles_s2[n].boundary_points, vec1, vec2, vec3);
				total_area_s2 += .5 * magnitude_cross_product(vec1, vec2);
				number_of_segmented_triangles_s2 += segmented_triangles_s2[n].simplices.size() / 3;
				total_edge_length_s2 += edge_length(segmented_triangles_s2[n].boundary_points);
			}
		}
		double average_area_s2 = total_area_s2 / s2_simplices.size();
		double average_edge_length_s2 = total_edge_length_s2 / s2_simplices.size();

		double width_s1 = 0;
		double width_s2 = 0;

		if (width > 0){
			width_s1 = width;
			width_s2 = width;
		}else{
			width_s1 = 3 * sqrt(average_area_s1 * 4 / sqrt(3));
			width_s2 = 3 * sqrt(average_area_s2 * 4 / sqrt(3));
		}

		if (!testing_script_mode){
			if (info){
				std::cout << std::endl;
				std::cout << "Average area of the triangles in surface 1 is: " << average_area_s1 << std::endl;
				std::cout << "Average edge length of the triangles in surface 1 is: " << average_edge_length_s1 << std::endl;
				std::cout << "Number of segmented triangles in surface 1 is: " << number_of_segmented_triangles_s1 << std::endl;
				std::cout << "Edge length of surface 1 regular grid: " << width_s1 << std::endl;
			}
			else{
				std::cout << "Number of segmented triangles in surface 1 is: " << number_of_segmented_triangles_s1 << std::endl;
			}
		}

		double average_number_of_boxes = 0;
		TriangleCellIntersection s1_grid(s1_points, s1_simplices, width_s1, average_number_of_boxes);
		if (!testing_script_mode){
			if (info){
				std::cout << "Average number of boxes intersected in surface 1: " << average_number_of_boxes << std::endl << std::endl;

				std::cout << "Average area of the triangles in surface 2 is: " << average_area_s2 << std::endl;
				std::cout << "Average edge length of the triangles in surface 2 is: " << average_edge_length_s2 << std::endl;
				std::cout << "Number of segmented triangles in surface 2 is: " << number_of_segmented_triangles_s2 << std::endl;
				std::cout << "Edge length of surface 2 regular grid: " << width_s2 << std::endl;
			}
			else{
				std::cout << "Number of segmented triangles in surface 2 is: " << number_of_segmented_triangles_s2 << std::endl;
			}
		}

		TriangleCellIntersection s2_grid(s2_points, s2_simplices, width_s2, average_number_of_boxes);
		if (!testing_script_mode){
			if (info){
				std::cout << "Average number of boxes intersected in surface 2: " << average_number_of_boxes << std::endl << std::endl;
			}
		}

		#ifdef _WIN32
				long cell_intersection_time = GetTickCount() - segment_start;
				segment_start = GetTickCount();
		#elif __unix
				uint64_t cell_intersection_time = getTimeMs() - segment_start
				segment_start = getTimeMs();
		#endif

		if (percent_done){
			std::cout << "0%" << std::endl;
		}
		surface_to_surface_extended(segmented_triangles_s1, segmented_triangles_s2, &s2_grid, epsilon, 1, percent_done);
		if (percent_done){
			std::cout << "50%" << std::endl;
		}
		surface_to_surface_extended(segmented_triangles_s2, segmented_triangles_s1, &s1_grid, epsilon, 2, percent_done);
		if (percent_done){
			std::cout << "100%" << std::endl;
		}
		
		double max_angle_1 = -1;
		double max_angle_2 = -1;

		double max_dist_1 = -1;
		double max_dist_2 = -1;

		//Finds the maximum angle and distance
		for (int x = 0; x < segmented_triangles_s1.size(); x++){
			double max_of_simplex = -1;
			if (segmented_triangles_s1[x].valid){
				for (int y = 0; y < segmented_triangles_s1[x].angle_diff_of_closest.size(); y++){
					if (max_angle_1 < segmented_triangles_s1[x].angle_diff_of_closest[y]){
						max_angle_1 = segmented_triangles_s1[x].angle_diff_of_closest[y];
					}
					if (max_of_simplex < segmented_triangles_s1[x].angle_diff_of_closest[y]){
						max_of_simplex = segmented_triangles_s1[x].angle_diff_of_closest[y];
					}
					if (max_dist_1 < segmented_triangles_s1[x].smallest_distance[y]){
						max_dist_1 = segmented_triangles_s1[x].smallest_distance[y];
					}
				}
				s1_angle_distances[x] = (max_of_simplex);
			}
			else{
				s1_angle_distances[x] = (-1);
			}
		}
		for (int x = 0; x < segmented_triangles_s2.size(); x++){
			double max_of_simplex = -1;
			if (segmented_triangles_s2[x].valid){
				for (int y = 0; y < segmented_triangles_s2[x].angle_diff_of_closest.size(); y++){
					if (max_angle_2 < segmented_triangles_s2[x].angle_diff_of_closest[y]){
						max_angle_2 = segmented_triangles_s2[x].angle_diff_of_closest[y];
					}
					if (max_of_simplex < segmented_triangles_s2[x].angle_diff_of_closest[y]){
						max_of_simplex = segmented_triangles_s2[x].angle_diff_of_closest[y];
					}
					if (max_dist_2 < segmented_triangles_s2[x].smallest_distance[y]){
						max_dist_2 = segmented_triangles_s2[x].smallest_distance[y];
					}
				}
				s2_angle_distances[x] = (max_of_simplex);
			}
			else{
				s2_angle_distances[x] = (-1);
			}
		}

		if (testing_script_mode){
			std::cout << max_angle_1 << "," << max_angle_2;
		}

		if (num_simplices_histogram){
			create_num_simplices_histogram(s1_angle_distances, s2_angle_distances, prefix, testing_script_mode);
		}

		if (area_simplices_histogram){
			create_area_simplices_histogram(segmented_triangles_s1, segmented_triangles_s2, prefix, testing_script_mode);
		}

		if (print_area_angle_cutoff >= 0){
			print_info_under_cutoff(segmented_triangles_s1, segmented_triangles_s2, s1_angle_distances, s2_angle_distances, print_area_angle_cutoff, testing_script_mode);
		}
		
		if (time){
#ifdef _WIN32
			cout << "Triangle segment time: " << segment_triangles_time << endl;
			cout << "Cell intersection time: " << cell_intersection_time << endl;
			cout << "Angle comp time: " << (GetTickCount() - segment_start) << endl;
			cout << "Total time: " << (GetTickCount() - start) << endl;
#elif __unix
			cout << "Triangle segment time: " << segment_triangles_time << endl;
			cout << "Cell intersection time: " << cell_intersection_time << endl;
			cout << "Angle comp time: " << (getTimeMs() - segment_start) << endl;
			cout << "Total time: " << (getTimeMs() - start) << endl;
#endif

		}

		cout << "Max of min angle from surface 1 to surface 2: " << max_angle_1 << endl;
		cout << "Max of min angle from surface 2 to surface 1: " << max_angle_2 << endl;
		cout << "Max of min distance from surface 1 to surface 2: " << max_dist_1 << endl;
		cout << "Max of min distance from surface 2 to surface 1: " << max_dist_2 << endl;

		//Returns the max angle
		return fmax(max_angle_1, max_angle_2);
	}


	void surface_angle_distance::print_info_under_cutoff(const std::vector<segmentedTriangle> & segmented_triangles_s1, const std::vector<segmentedTriangle> & segmented_triangles_s2, const std::vector<double> & angles1, const std::vector<double> & angles2, double cutoff, bool testing_script_mode){

		double total_area = 0;
		double area_below_cutoff = 0;

		int total_simplices = 0;
		int simplices_below = 0;

		int num_invalid = 0;

		for (int n = 0; n < angles1.size(); n++){
			if (angles1[n] != -1){
				total_simplices++;
				if (angles1[n] <= cutoff){
					simplices_below++;
				}
			}else{
				num_invalid++;
			}
		}

		for (int n = 0; n < segmented_triangles_s1.size(); n++){

			if (segmented_triangles_s1[n].valid){

				std::vector<double> vec1(DIM3);
				std::vector<double> vec2(DIM3);
				std::vector<double> vec3(DIM3);
				create_vectors_from_triangle(segmented_triangles_s1[n].boundary_points, vec1, vec2, vec3);

				for (int m = 0; m < segmented_triangles_s1[n].simplices.size(); m += 3){

					std::vector<double> vec1_seg(DIM3);
					std::vector<double> vec2_seg(DIM3);
					std::vector<double> vec3_seg(DIM3);

					std::vector<double> new_triangle(DIM3 * 3);
					for (int x = 0; x < 3; x++){
						new_triangle[3 * x + 0] = (segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 0]);
						new_triangle[3 * x + 1] = (segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 1]);
						new_triangle[3 * x + 2] = (segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 2]);
					}

					create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
					double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

					total_area += increment_area;

					if (segmented_triangles_s1[n].angle_diff_of_closest[m / 3] <= cutoff){
						area_below_cutoff += increment_area;
					}
				}
			}
		}

		if (testing_script_mode){
			std::cout << "," << cutoff << "," << ((total_area - area_below_cutoff) / total_area) << "," << (total_area - area_below_cutoff) << "," << total_area;
			std::cout << "," << ((double)(total_simplices - simplices_below) / (double)total_simplices) << "," << (total_simplices - simplices_below) << "," << total_simplices;
			std::cout << "," << num_invalid;
		}else {
			std::cout << area_below_cutoff / total_area << " of the area (" << area_below_cutoff << "/" << total_area << ") is below the cutoff of " << cutoff << "." << endl;
		}

		total_area = 0;
		area_below_cutoff = 0;

		total_simplices = 0;
		simplices_below = 0;
		
		num_invalid = 0;

		for (int n = 0; n < angles2.size(); n++){
			if (angles2[n] != -1){
				total_simplices++;
				if (angles2[n] <= cutoff){
					simplices_below++;
				}
			}
			else{
				num_invalid++;
			}
		}

		for (int n = 0; n < segmented_triangles_s2.size(); n++){
			if (segmented_triangles_s2[n].valid){
				std::vector<double> vec1(DIM3);
				std::vector<double> vec2(DIM3);
				std::vector<double> vec3(DIM3);
				create_vectors_from_triangle(segmented_triangles_s2[n].boundary_points, vec1, vec2, vec3);

				for (int m = 0; m < segmented_triangles_s2[n].simplices.size(); m+=3){
					

					std::vector<double> vec1_seg(DIM3);
					std::vector<double> vec2_seg(DIM3);
					std::vector<double> vec3_seg(DIM3);

					std::vector<double> new_triangle(3 * DIM3);
					for (int x = 0; x < 3; x++){
						new_triangle[3 * x + 0] = (segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 0]);
						new_triangle[3 * x + 1] = (segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 1]);
						new_triangle[3 * x + 2] = (segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 2]);
					}

					create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
					double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

					total_area += increment_area;

					if (segmented_triangles_s2[n].angle_diff_of_closest[m / 3] <= cutoff){
						area_below_cutoff += increment_area;
					}
				}
			}
		}

		if (testing_script_mode){
			std::cout << "," << ((total_area - area_below_cutoff) / total_area) << "," << (total_area - area_below_cutoff) << "," << total_area;
			std::cout << "," << ((double)(total_simplices - simplices_below) / (double)total_simplices) << "," << (total_simplices - simplices_below) << "," << total_simplices;
			std::cout << "," << num_invalid << ",";
		}
		else {
			std::cout << area_below_cutoff / total_area << " of the area (" << area_below_cutoff << "/" << total_area << ") is below the cutoff of " << cutoff << "." << endl;
		}

	}

	/*
	Creates a histogram of the area of simplices in range of angles
	*/
	void surface_angle_distance::create_area_simplices_histogram(const std::vector<segmentedTriangle> & segmented_triangles_s1, const std::vector<segmentedTriangle> & segmented_triangles_s2, std::string prefix, bool testing_script_mode){
	
		std::string histo1_name = prefix + "area_simplices_histogram_1.dat";
		if (testing_script_mode){
			histo1_name = prefix + "area_simplices_histogram_perfect.dat";
		}
		std::ofstream histo1(histo1_name.c_str());

		int size = 181;

		if (ignore_orientation){
			size = 91;
		}

		double * bins = new double[size];

		for (int n = 0; n < size; n++){
			bins[n] = 0;
		}

		double total_area = 0;

		for (int n = 0; n < segmented_triangles_s1.size(); n++){
			if (segmented_triangles_s1[n].valid){
				std::vector<double> vec1(DIM3);
				std::vector<double> vec2(DIM3);
				std::vector<double> vec3(DIM3);
				create_vectors_from_triangle(segmented_triangles_s1[n].boundary_points, vec1, vec2, vec3);
				total_area += .5 * magnitude_cross_product(vec1, vec2);
				
				for (int m = 0; m < segmented_triangles_s1[n].simplices.size(); m+=3){

					int index = (int)floor(segmented_triangles_s1[n].angle_diff_of_closest[m/3]);
					std::vector<double> vec1_seg(DIM3);
					std::vector<double> vec2_seg(DIM3);
					std::vector<double> vec3_seg(DIM3);

					std::vector<double> new_triangle(DIM3 * 3);
					for (int x = 0; x < 3; x++){
						new_triangle[3 * x + 0] = (segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 0]);
						new_triangle[3 * x + 1] = (segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 1]);
						new_triangle[3 * x + 2] = (segmented_triangles_s1[n].points[3 * segmented_triangles_s1[n].simplices[m + x] + 2]);
					}

					create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
					double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

					bins[index] += increment_area;
				}
			}
		}

		bins[size - 2] += bins[size - 1];

		histo1 << "#Angle\tNormalized area of simplices" << endl;
		for (int n = 0; n < size - 1; n++){
			histo1 << (n + 1) << "\t" << (bins[n] / total_area) << endl;
		}

		histo1.close();


		std::string histo2_name = prefix + "area_simplices_histogram_2.dat";
		if (testing_script_mode){
			histo2_name = prefix + "area_simplices_histogram_religrad.dat";
		}
		std::ofstream histo2(histo2_name.c_str());

		//bins[91];

		for (int n = 0; n < size; n++){
			bins[n] = 0;
		}

		total_area = 0;

		for (int n = 0; n < segmented_triangles_s2.size(); n++){
			if (segmented_triangles_s2[n].valid){
				std::vector<double> vec1(DIM3);
				std::vector<double> vec2(DIM3);
				std::vector<double> vec3(DIM3);
				create_vectors_from_triangle(segmented_triangles_s2[n].boundary_points, vec1, vec2, vec3);
				total_area += .5 * magnitude_cross_product(vec1, vec2);

				for (int m = 0; m < segmented_triangles_s2[n].simplices.size(); m+=3){

					int index = (int)floor(segmented_triangles_s2[n].angle_diff_of_closest[m/3]);
					std::vector<double> vec1_seg(DIM3);
					std::vector<double> vec2_seg(DIM3);
					std::vector<double> vec3_seg(DIM3);

					std::vector<double> new_triangle(3 * DIM3);
					for (int x = 0; x < 3; x++){
						new_triangle[3 * x + 0] = (segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 0]);
						new_triangle[3 * x + 1] = (segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 1]);
						new_triangle[3 * x + 2] = (segmented_triangles_s2[n].points[3 * segmented_triangles_s2[n].simplices[m + x] + 2]);
					}

					create_vectors_from_triangle(new_triangle, vec1_seg, vec2_seg, vec3_seg);
					double increment_area = .5 * magnitude_cross_product(vec1_seg, vec2_seg);

					bins[index] += increment_area;
				}
			}
		}

		bins[size - 2] += bins[size - 1];
		
		histo2 << "#Angle\tNormalized area of simplices" << endl;
		for (int n = 0; n < size - 1; n++){
			histo2 << (n + 1) << "\t" << (bins[n] / total_area) << endl;
		}

		histo2.close();
	}

	/*
	Creates a histogram of the number of simplices in ranges of angles
	*/
	void surface_angle_distance::create_num_simplices_histogram(const std::vector<double> & angles1, const std::vector<double> & angles2, std::string prefix, bool testing_script_mode){

		std::string histo1_name = prefix + "num_simplices_histogram_1.dat";
		std::ofstream histo1(histo1_name.c_str());


		int size = 181;

		if (ignore_orientation){
			size = 91;
		}

		int  * bins = new int[size];
		for (int n = 0; n < size; n++){
			bins[n] = 0;
		}
		for (int n = 0; n < angles1.size(); n++){
			int index = floor(angles1[n]);
			bins[index]++;
		}
		bins[size - 2] += bins[size - 1];

		histo1 << "#Angle\tNumber of simplices" << endl;
		for (int n = 0; n < size - 1; n++){
			histo1 << (n + 1) << "\t" << bins[n] << endl;
		}

		histo1.close();

		std::string histo2_name = prefix + "num_simplices_histogram_2.dat";
		std::ofstream histo2(histo2_name.c_str());

		for (int n = 0; n < size; n++){
			bins[n] = 0;
		}

		for (int n = 0; n < angles2.size(); n++){
			int index = floor(angles2[n]);
			bins[index]++;
		}
		bins[size - 2] += bins[size - 1];

		histo2 << "#Angle\tNumber of simplices" << endl;
		for (int n = 0; n < size - 1; n++){
			histo2 << (n + 1) << "\t" << bins[n] << endl;
		}

		histo2.close();

	}

	/*
	Calculates the extended surface to surface information from triangles1 to triangles2
	*/
	void surface_angle_distance::surface_to_surface_extended(std::vector<segmentedTriangle> & triangles1, std::vector<segmentedTriangle> & triangles2, TriangleCellIntersection * grid, double epsilon, int num_surface, bool percent_done){
		
		//Loops through all triangles in surface 1
		for (int x = 0; x < triangles1.size(); x++){
			//std::cout << x << " " << triangles1.size() << "\n";

			//std::cout << (int)((100.0 * (x) / triangles1.size()) / 10) << " " << (int)((100.0 * (x - 1) / triangles1.size())) << std::endl;

			if (num_surface == 1){
				if ((int)((100.0 * (x) / triangles1.size()) / 10) >(int)((100.0 * (x - 1) / triangles1.size() / 10))){
					if (percent_done){
						std::cout << (int)(100.0 * (x) / triangles1.size() / 10) * 5 << "%" << std::endl;
					}
				}
			}
			else{
				if ((int)((100.0 * (x) / triangles1.size()) / 10) >(int)((100.0 * (x - 1) / triangles1.size() / 10))){
					if (percent_done){
						std::cout << 50 + (int)(100.0 * (x) / triangles1.size() / 10) * 5 << "%" << std::endl;
					}
				}
			}

			segmentedTriangle * current_triangle = &triangles1[x];

			//Compute info if the current triangle is valid
			if (current_triangle->valid){

				/*current_triangle->angle_diff_of_closest.reserve(current_triangle->simplices.size() / 3);
				current_triangle->closest_distance.reserve(current_triangle->simplices.size() / 3);
				current_triangle->closest.reserve(current_triangle->simplices.size() / 3);
				current_triangle->smallest_distance.reserve(current_triangle->simplices.size() / 3);
				current_triangle->closest_point.reserve(current_triangle->simplices.size());*/

				//Loops through all segmented trianlges in the current triangle
				for (int y = 0; y < current_triangle->simplices.size(); y += 3){
					double angle_diff = 360;
					double distance = 0;
					int closest = -1;
					std::vector<double> closest_point(DIM3);
					std::vector<double> temporary_point(DIM3);
					std::vector<double> center_point(DIM3);
					std::vector<double> triangle(DIM3 * 3);

					//Prepares a triangle
					for (int z = y; z < y + 3; z++){
						triangle[3 * (z - y) + 0] = (current_triangle->points[3 * current_triangle->simplices[z] + 0]);
						triangle[3 * (z - y) + 1] = (current_triangle->points[3 * current_triangle->simplices[z] + 1]);
						triangle[3 * (z - y) + 2] = (current_triangle->points[3 * current_triangle->simplices[z] + 2]);
					}

					//Computes center of a triangle
					center_of_triangle(triangle, center_point);

					double nearest_distance = DBL_MAX;
					double distance_angle = DBL_MAX;
					double smallest_distance = DBL_MAX;

					//Finds the nearest point to the center of the triangle
					point_to_surface_nearest_point(center_point, triangles2, current_triangle->normal, nearest_distance, angle_diff, temporary_point, grid, current_triangle->normal, closest, epsilon, closest_point, distance_angle, smallest_distance);

					//Computes the minimum angle difference
					//compute_angle_difference(center_point, triangles2, current_triangle->normal, distance, closest, angle_diff, closest_point, temporary_point, epsilon, grid);
					

					//Stores information to the current triangle
					current_triangle->angle_diff_of_closest[y / 3] = (angle_diff);
					current_triangle->closest_distance[y / 3] = (distance_angle);
					current_triangle->closest[y / 3] = (closest);
					current_triangle->smallest_distance[y / 3] = (smallest_distance);
					current_triangle->closest_point[y / 3 + 0] = (closest_point[0]);
					current_triangle->closest_point[y / 3 + 1] = (closest_point[1]);
					current_triangle->closest_point[y / 3 + 2] = (closest_point[2]);

				}				
			}
		}
	}

	/*
	Computes the minimum angle difference between a point with a predefined normal to another surface within the distance epsilon from another point
	*/
	void surface_angle_distance::compute_angle_difference(const std::vector<double> & point, const std::vector<segmentedTriangle> & surface, const std::vector<double> & normal_to_triangle, double & distance, int & minimum_angle_index, double & minimum_angle_diff, std::vector<double> & minimum_angle_point, const std::vector<double> & point_to_compute_distance_from, double epsilon, TriangleCellIntersection * grid){

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
				std::vector<double> temp_vector(DIM3);

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
	void surface_angle_distance::point_to_surface_nearest_point(const std::vector<double> & point, std::vector<segmentedTriangle> & surface, const std::vector<double> & normal_to_triangle, double & distance, double & angle, std::vector<double> & closest_point, TriangleCellIntersection * grid, const std::vector<double> & current_normal, int & minimum_angle_index, double epsilon, std::vector<double> & minimum_angle_point, double & distance_angle, double & smallest_distance){

		distance = DBL_MAX;
		distance_angle = DBL_MAX;
		angle = DBL_MAX;
		int current_shell = 0;
		std::vector<int> checked;

		std::vector<int> simplex_index;
		std::vector<double> distances;
		std::vector<double> angle_differences;
		std::vector<double> closest_points;

		//Loop until it is not possible to find a smaller distance
		while (distance >= grid->get_width() * (current_shell - 1)){
			std::vector<int> simplices_in_shell;
			//Find simplices in current shell

			grid->find_simplices_in_shell_distance_from_cube(point, current_shell, simplices_in_shell);

			//Compute distance of simplices in shell
			for (int n = 0; n < simplices_in_shell.size(); n++){
				//if (std::find(checked.begin(), checked.end(), simplices_in_shell[n]) == checked.end()){
				if (!surface[simplices_in_shell[n]].checked){
					segmentedTriangle current_triangle = surface[simplices_in_shell[n]];
					if (current_triangle.valid){

						double temp_distance = 0;
						std::vector<double> temp_vector(DIM3);

						//Calculate the distance from the point to the trianlge
						point_to_triangle(point, current_triangle, temp_distance, temp_vector);

						//If the current distance is less than the current minimum, update the information
						if (temp_distance < distance){
							distance = temp_distance;
							for (int x = 0; x < 3; x++){
								closest_point[x] = (temp_vector[x]);
							}
						}

						//Calculate the angles between normals
						double current_angle = angle_of_vectors(current_normal, current_triangle.normal);

						distances.push_back(temp_distance);
						angle_differences.push_back(current_angle);
						simplex_index.push_back(simplices_in_shell[n]);
						for (int x = 0; x < 3; x++){
							closest_points.push_back(temp_vector[x]);
						}

					}

					surface[simplices_in_shell[n]].checked = true;
					checked.push_back(simplices_in_shell[n]);
				}
			}
			//std::cout << current_shell << " " << distance << std::endl;
			current_shell++;
		}

		smallest_distance = distance;

		

		for (int n = 0; n < distances.size(); n++){
			if (distances[n] <= (distance + epsilon)){
				if (angle > angle_differences[n]){
					angle = angle_differences[n];
					minimum_angle_index = simplex_index[n];
					distance_angle = distances[n];

					minimum_angle_point[0] = (closest_points[3 * n + 0]);
					minimum_angle_point[1] = (closest_points[3 * n + 1]);
					minimum_angle_point[2] = (closest_points[3 * n + 2]);
				}
			}
		}


		/////////Added//////////////////
		//Make a vector of all simplices possibly in the epsilon range

		//Simplices to check
		std::vector<int> simplices_to_check;

		while (epsilon + distance >= grid->get_width()*(current_shell - 1)){
			std::vector<int> temp_vector;
			grid->find_simplices_in_shell_distance_from_cube(point, current_shell, temp_vector);

			for (int n = 0; n < temp_vector.size(); n++){
				//if (std::find(simplices_to_check.begin(), simplices_to_check.end(), temp_vector[n]) == simplices_to_check.end() && std::find(checked.begin(), checked.end(), temp_vector[n]) == checked.end()){
				if (!surface[temp_vector[n]].checked){
					simplices_to_check.push_back(temp_vector[n]);
					surface[temp_vector[n]].checked = true;
					checked.push_back(temp_vector[n]);
				}
			}
			current_shell++;
		}

		//Loop through all simplices that could possibly meet the epsilon distance
		for (int n = 0; n < simplices_to_check.size(); n++){

			segmentedTriangle current_triangle = surface[simplices_to_check[n]];

			//If the surface is valid continue
			if (current_triangle.valid){

				//Calculates the angle between the normal and the normal of the current triangle
				double currentAngle = angle_of_vectors(normal_to_triangle, current_triangle.normal);

				//If the distance is less than the specified epsilon, continue
				if (currentAngle < angle){

					double distance_from_point = 0;
					std::vector<double> temp_vector(DIM3);

					//Computes the distance from a point to the triangle
					point_to_triangle(point, current_triangle, distance_from_point, temp_vector);

					//If the current angle is less than the minimum, update the new information
					if (distance_from_point < epsilon + distance){
						angle = currentAngle;
						minimum_angle_index = simplices_to_check[n];

						distance_angle = distance_from_point;

						minimum_angle_point[0] = (temp_vector[0]);
						minimum_angle_point[1] = (temp_vector[1]);
						minimum_angle_point[2] = (temp_vector[2]);
					}
				}
			}
		}

		for (int n = 0; n < checked.size(); n++){
			surface[checked[n]].checked = false;
		}
	}

	void surface_angle_distance::point_to_triangle(const std::vector<double> & point, const segmentedTriangle & triangle, double & distance, std::vector<double> & closest_point){
		
		std::vector<double> vec_a(DIM3);
		std::vector<double> vec_b(DIM3);
		std::vector<double> vec_c(DIM3);

		//Creates vectors of the triangle points

		for (int x = DIM3 * 0; x < DIM3 * 1; x++){
			vec_a[x - DIM3 * 0] = triangle.boundary_points[x];
		}
		for (int x = DIM3 * 1; x < DIM3 * 2; x++){
			vec_b[x - DIM3 * 1] = triangle.boundary_points[x];
		}
		for (int x = DIM3 * 2; x < DIM3 * 3; x++){
			vec_c[x - DIM3 * 2] = triangle.boundary_points[x];
		}

		std::vector<double> projection_vector(DIM3);

		//Find a vector from the point to a point on the plane of the triangle
		std::vector<double> point_to_plane(DIM3);
		vector_point_to_point(point, vec_a, point_to_plane);

		//Projects the point on to the normal vector
		vector_projection(point_to_plane, triangle.normal, projection_vector);

		//Calculates the distance of the projection vector
		distance = vector_magnitude(projection_vector);

		//Finds a point projected onto the plane
		std::vector<double> projected_point(DIM3);
		for (int x = 0; x < 3; x++){
			projected_point[x] = (point[x] + projection_vector[x]);
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

			std::vector<double> vector_cb(DIM3);
			std::vector<double> vector_cp(DIM3);

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
			
			for (int x = 0; x < 3; x++){
				closest_point[x] = (projected_point[x]);
			}
		}
	}

	void surface_angle_distance::point_to_line(const std::vector<double> & point, const std::vector<double> & line_a, const std::vector<double> & line_b, std::vector<double> & closest_point){
		
		//Creates vectors for point a to the point and from point a to point b
		std::vector<double> a_to_point(DIM3);
		std::vector<double> line_ab(DIM3);

		vector_point_to_point(line_a, point, a_to_point);
		vector_point_to_point(line_a, line_b, line_ab);

		if (dot_product(a_to_point, line_ab) < 0){

			//Point a is the closest point
			for (int x = 0; x < DIM3; x++){
				closest_point[x] = (line_a[x]);
			}
		}
		else if (pow(vector_magnitude(line_ab), 2) < dot_product(a_to_point, line_ab)){

			//Point b is the closest point
			for (int x = 0; x < DIM3; x++){
				closest_point[x] = (line_b[x]);
			}
		}
		else{

			//A point on the line is the closest point
			std::vector<double> projected_vector(DIM3);

			//Projects the vector a_to_point onto line_ab
			vector_projection(a_to_point, line_ab, projected_vector);

			//Solves for the vector from the point on the line to the original point
			std::vector<double> line_to_point(DIM3);
			for (int x = 0; x < DIM3; x++){
				line_to_point[x] = (a_to_point[x] - projected_vector[x]);
			}

			//Stores in the closest point
			for (int x = 0; x < DIM3; x++){
				closest_point[x] = (point[x] - line_to_point[x]);
			}
		}
	}

	void surface_angle_distance::vector_point_to_point(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & ab){
		
		//Calculates the differences of the terms of the vectors
		for (int x = 0; x < DIM3; x++){
			ab[x] = (b[x] - a[x]);
		}
	}

	/*
	Takes triangles and creates a vector of the triangles segmented
	*/
	void surface_angle_distance::create_segmented_triangles(const std::vector<double> & s_points, const std::vector<int> & s_simplices, std::vector<segmentedTriangle> & segmented_triangles, double delta){

		//Loops through all simplices
		for (int n = 0; n < s_simplices.size(); n += 3){
			std::vector<double> triangle_points(DIM3 * 3);

			//Creates a triangle of the points in the simplices
			for (int m = n; m < n + 3; m++){
				triangle_points[DIM3 * (m - n) + 0] = (s_points[3 * s_simplices[m] + 0]);
				triangle_points[DIM3 * (m - n) + 1] = (s_points[3 * s_simplices[m] + 1]);
				triangle_points[DIM3 * (m - n) + 2] = (s_points[3 * s_simplices[m] + 2]);
			}

			//Vector normal to the triangle
			std::vector<double> normal_unit_vector(DIM3);

			//Vector a to b
			std::vector<double> vec1(DIM3);
			//Vector c to a
			std::vector<double> vec2(DIM3);
			//Vector b to c
			std::vector<double> vec3(DIM3);

			//The following vectors are also in the plane
			//Vector normal to ab
			std::vector<double> normal_vec1(DIM3);
			//Vector normal to ca
			std::vector<double> normal_vec2(DIM3);
			//Vector normal to bc
			std::vector<double> normal_vec3(DIM3);

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
							std::vector<int> segmented_simplices(triangle_points_copy.size() / 3);

							for (int a = 0; a < triangle_points_copy.size() / 3; a++){
								segmented_simplices[a] = a;
							}

							//Compacts the simplices so that there is no "redundant" data
							compact_simplices(segmented_simplices, triangle_points_copy);

							segmentedTriangle segmented_triangle(segmented_simplices.size() / 3);

							//Sets the triangle validity to true because
							segmented_triangle.valid = true;

							//Sets checked to false
							segmented_triangle.checked = false;

							//Sets points to the newly computed segmentations
							segmented_triangle.points = triangle_points_copy;
							segmented_triangle.simplices = segmented_simplices;
							/*
							ofstream off_file;
							off_file.open("triangle_test.off");

							off_file << "DOFF" << endl;
							off_file << triangle_points_copy.size() / 3 << " " << segmented_simplices.size() / 3 << " 0" << endl;
							for (int q = 0; q < triangle_points_copy.size(); q++){
								off_file << triangle_points_copy[q];
								if (q % 3 == 2){
									off_file << endl;
								}
								else {
									off_file << " ";
								}
							}
							off_file << endl;
							for (int q = 0; q < segmented_simplices.size(); q++){
								if (q % 3 == 0){
									off_file << "3 ";
								}
								off_file << segmented_simplices[q];
								if (q % 3 == 2){
									off_file << endl;
								}
								else {
									off_file << " ";
								}
							}


							off_file.close();*/

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
									segmented_triangle.boundary_points[x - 3] = (triangle_points[x]);
								}
								for (int x = 0; x < 3; x++){
									segmented_triangle.boundary_points[x + 6] = (triangle_points[x]);
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
									segmented_triangle.boundary_points[x - 6] = (triangle_points[x]);
								}
								for (int x = 0; x < 3; x++){
									segmented_triangle.boundary_points[x + 3] = (triangle_points[x]);
								}
								for (int x = 3; x < 6; x++){
									segmented_triangle.boundary_points[x + 3] = (triangle_points[x]);
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
									segmented_triangle.boundary_points[x] = (triangle_points[x]);
								}


							}
							//Add to the list of segmented triangles
							segmented_triangles[n / 3] = (segmented_triangle);
							added = true;
						}
					}
				}
			}

			if(!added){
				segmentedTriangle segmented_triangle;
				segmented_triangle.valid = false;
				segmented_triangles[n / 3] = (segmented_triangle);
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
	bool surface_angle_distance::normal_unit_vector_to_triangle(const std::vector<double> & triangle, std::vector<double> & normal_vector){
		
		//Creating vectors that correspond to the sides of the triangle
		std::vector<double> vec1(DIM3);
		std::vector<double> vec2(DIM3);

		//Making two vectors of the sides of the triangles
		for (int n = 0; n < 3; n++){
			vec2[n] = (triangle[n + 3] - triangle[n + 0]);
			vec1[n] = (triangle[n + 6] - triangle[n + 0]);
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
			normal_vector[0] = (term1);
			normal_vector[1] = (term2);
			normal_vector[2] = (term3);
		}

		return true;

	}

	/*
	Segments a triangle into multiple triangles with edge lengths less than delta
	*/
	void surface_angle_distance::segment_triangle(std::vector<double> & triangle, double delta){
	
		if (edge_length(triangle) > delta){

			if (edge_length(triangle) > 1.25 * smallest_edge_length(triangle)){
			
				int vertex_index = edge_length_vertex(triangle);

				//Midpoint
				std::vector<double> midpoint(DIM3);

				//Shifted points
				std::vector<double> triangle_shifted(DIM3 * 3);

				//Creating new triangles
				std::vector<double> t1(DIM3 * 3);
				std::vector<double> t2(DIM3 * 3);

				if (vertex_index == 1){
					for (int n = 0; n < 3; n++){
						midpoint[n] = ((triangle[n + 3] + triangle[n + 6]) / 2.0);
					}
					for (int n = 0; n < 9; n++){
						triangle_shifted[n] = (triangle[n]);
					}
				}
				else if (vertex_index == 2){
					for (int n = 0; n < 3; n++){
						midpoint[n] = ((triangle[n + 0] + triangle[n + 6]) / 2.0);
					}
					for (int n = 3; n < 9; n++){
						triangle_shifted[n - 3] = (triangle[n]);
					}
					for (int n = 0; n < 3; n++){
						triangle_shifted[n + 6] = (triangle[n]);
					}
				}
				else{
					for (int n = 0; n < 3; n++){
						midpoint[n] = ((triangle[n + 0] + triangle[n + 3]) / 2.0);
					}
					for (int n = 6; n < 9; n++){
						triangle_shifted[n - 6] = (triangle[n]);
					}
					for (int n = 0; n < 6; n++){
						triangle_shifted[n + 3] = (triangle[n]);
					}
				}

				//Creating the 2 new triangles
				for (int n = 0; n < 3; n++){
					t1[0 + n] = (triangle_shifted[n + 0]);
					t2[0 + n] = (triangle_shifted[n + 0]);
				}
				for (int n = 0; n < 3; n++){
					t1[3 + n] = (triangle_shifted[n + 3]);
					t2[3 + n] = (midpoint[n]);
				}
				for (int n = 0; n < 3; n++){
					t1[6 + n] = (midpoint[n]);
					t2[6 + n] = (triangle_shifted[n + 6]);
				}

				

				//Recursively segmenting the newly created triangles
				segment_triangle(t1, delta);
				segment_triangle(t2, delta);

				//Replacing triangle with the new segmented triangles
				triangle.clear();
				triangle.insert(triangle.end(), t1.begin(), t1.end());
				triangle.insert(triangle.end(), t2.begin(), t2.end());

			}
			else{

				//Creating new triangles
				std::vector<double> t1(DIM3 * 3);
				std::vector<double> t2(DIM3 * 3);
				std::vector<double> t3(DIM3 * 3);
				std::vector<double> t4(DIM3 * 3);

				//Finding midpoints of current edges
				std::vector<double> midpoint01(DIM3);
				std::vector<double> midpoint02(DIM3);
				std::vector<double> midpoint12(DIM3);


				for (int n = 0; n < 3; n++){
					midpoint01[n] = ((triangle[n + 0] + triangle[n + 3]) / 2.0);
					midpoint02[n] = ((triangle[n + 0] + triangle[n + 6]) / 2.0);
					midpoint12[n] = ((triangle[n + 3] + triangle[n + 6]) / 2.0);
				}

				//Creating the 4 new triangles
				for (int n = 0; n < 3; n++){
					t1[n] = (triangle[n + 0]);
					t2[n] = (triangle[n + 3]);
					t3[n] = (triangle[n + 6]);
					t4[n] = (midpoint01[n]);
				}

				for (int n = 0; n < 3; n++){
					t1[n + 3] = (midpoint01[n]);
					t2[n + 3] = (midpoint12[n]);
					t3[n + 3] = (midpoint02[n]);
					t4[n + 3] = (midpoint12[n]);
				}

				for (int n = 0; n < 3; n++){
					t1[n + 6] = (midpoint02[n]);
					t2[n + 6] = (midpoint01[n]);
					t3[n + 6] = (midpoint12[n]);
					t4[n + 6] = (midpoint02[n]);
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

	}

	/*
	Finds distance between two points
	*/
	double surface_angle_distance::point_to_point_distance(const std::vector<double> & p1, const std::vector<double> & p2){
		
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
	void surface_angle_distance::center_of_triangle(const std::vector<double> & triangle, std::vector<double> & point){
		
		for (int x = 0; x < 3; x++){
			//Finds the average of the points
			double average = 0;
			average += triangle[x + 0];
			average += triangle[x + 3];
			average += triangle[x + 6];
			average /= 3.0;

			point[x] = (average);
		}
	}

	/*
	Finds the edge length of a triangle
	*/
	double surface_angle_distance::edge_length(const std::vector<double> & triangle){
		
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

	/*
	Finds the vertex opposite of the longest edge
	*/
	int surface_angle_distance::edge_length_vertex(const std::vector<double> & triangle){

		int vertex_index = 0;

		//Points that define the triangle
		std::vector<double> p1(triangle.begin() + 0, triangle.begin() + 3);
		std::vector<double> p2(triangle.begin() + 3, triangle.begin() + 6);
		std::vector<double> p3(triangle.begin() + 6, triangle.begin() + 9);

		if (point_to_point_distance(p1, p2) > point_to_point_distance(p1, p3)){

			if (point_to_point_distance(p1, p2) > point_to_point_distance(p2, p3)){
				vertex_index = 3;
			}
			else{
				vertex_index = 1;
			}
		}else{

			if (point_to_point_distance(p1, p3) > point_to_point_distance(p2, p3)){
				vertex_index = 2;
			}
			else{
				vertex_index = 1;
			}
		}

		return vertex_index;
	}

	/*
	Finds the smallest edge length of a triangle
	*/
	double surface_angle_distance::smallest_edge_length(const std::vector<double> & triangle){

		double edge_length = 0;

		//Points that define the triangle
		std::vector<double> p1(triangle.begin() + 0, triangle.begin() + 3);
		std::vector<double> p2(triangle.begin() + 3, triangle.begin() + 6);
		std::vector<double> p3(triangle.begin() + 6, triangle.begin() + 9);

		double max_length = fmin(point_to_point_distance(p1, p2), point_to_point_distance(p1, p3));
		max_length = fmin(max_length, point_to_point_distance(p2, p3));

		//Sum of distance between points
		//edge_length += point_to_point_distance(p1, p2);
		//edge_length += point_to_point_distance(p1, p3);
		//edge_length += point_to_point_distance(p2, p3);

		return max_length;
	}

	void surface_angle_distance::create_vectors_from_triangle(const std::vector<double> & triangle, std::vector<double> & vec1, std::vector<double> & vec2, std::vector<double> & vec3){
		

		//Difference from a to b
		for (int x = 0; x < DIM3; x++){
			vec1[x] = (triangle[x + DIM3 * 1] - triangle[x + DIM3 * 0]);
		}

		//Difference from c to a
		for (int x = 0; x < DIM3; x++){
			vec2[x] = (triangle[x + DIM3 * 0] - triangle[x + DIM3 * 2]);
		}

		//Difference from b to c
		for (int x = 0; x < DIM3; x++){
			vec3[x] = (triangle[x + DIM3 * 2] - triangle[x + DIM3 * 1]);
		}

	}

	/*
	Computes the dot product of two vectors
	*/
	double surface_angle_distance::dot_product(const std::vector<double> & a, const std::vector<double> & b){
	
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
	void surface_angle_distance::vector_projection(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & projection){

		//Copying vector b
		for (int x = 0; x < 3; x++){
			projection[x] = (b[x]);
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
	double surface_angle_distance::vector_magnitude(const std::vector<double> & a){
		
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
	bool surface_angle_distance::normal_cross_product(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & normal_vector){

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
			normal_vector[0] = term1;
			normal_vector[1] = term2;
			normal_vector[2] = term3;
		}

		//Returning success
		return true;
	}

	/*
	Find the cross product
	*/
	double surface_angle_distance::magnitude_cross_product(const std::vector<double> & a, const std::vector<double> & b){

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
	double surface_angle_distance::angle_of_vectors(const std::vector<double> & a, const std::vector<double> & b){
	
		//Dot product of vectors a and b
		double dot = dot_product(a, b);

		//Magnitude of vector a times magnitude of vector b
		double magnitudes = vector_magnitude(a) * vector_magnitude(b);

		if (magnitudes == 0){
			//Returning 180 or 90 because a vector is the zero vector
			if (ignore_orientation)
				return 90;
			else
				return 180;
		}
		else{
			//a * b = |a||b|cos(x), acos((a * b)/(|a||b|)) = x
			double x = dot / magnitudes;
			double angle = 0;
			angle = acos(x) * 180 / M_PI;

			if (angle != angle)
				angle = 0;

			if (ignore_orientation){
				//Adjusting the range to 0 - 90
				angle = 90 - abs(angle - 90);
			}

			return angle;
		}
	}


/*
 * sampling.h
 *
 *  Created on: 02.03.2019
 *      Author: willy
 */

#include <vector>
#include <assert.h>

#include <random>
#include <algorithm>
#include "mathFunctions.h"



using namespace std;

double MAXVALUE = 10000000000;

unsigned int minDistance(const std::vector< std::vector<double> > &points,
							const std::vector<double> d_i){
	/*
	 * From the points returns the maximal minimal distance of all other points to the point in ind.
	 */


}

void drawOnePoint(	const std::vector< std::vector<double> > &points,
					std::vector< std::vector<double> > &selected_points,
					vector<unsigned int> &selected_points_indices,
					unsigned int ind,
					std::vector< std::vector<double> > &d){
	/*
	 * adds the point at index ind from points to selected_points
	 * deletes the point from points
	 * sets the distances in d to MAXVALUE for the row ind and column ind
	 */

	selected_points[0].push_back(points[0][ind]);
	selected_points[1].push_back(points[1][ind]);
	selected_points[2].push_back(points[2][ind]);

	selected_points_indices.push_back(ind);

	for(unsigned int i = 0; i < d.size(); i++){
		d[i][ind] = MAXVALUE;
//		d[ind][i] = MAXVALUE;
	}
}

vector<unsigned int> euclidean_farthest_point_sampling( const std::vector< std::vector<double> > &points,
										std::vector< std::vector<double> > &selected_points,
										const unsigned int K){
	/*
	 * as described in memoli:
	 * K ... number of points to sample.
	 * randomly choose a starting point.
	 * Then take the point with the maximal distance.
	 * Subsequent points are selected maximizing the minimal distance of the
	 * already selected points.
	 *
	 * In memoli 4000 points were selected in a first step
	 */

	assert(points[0].size() >= K);

	// calculate the distance matrix
	std::vector< std::vector<double> > d(points[0].size(), vector<double>(points[0].size(),0));
	dist_matrix(d,points);

	// we subsequently sample points from points_rest
//	std::vector< std::vector<double> > points_rest = points;

	srand(time(NULL));

	// get a random start point
	unsigned int random_start_ind = rand() % points[0].size();

	selected_points[0].resize(0);
	selected_points[1].resize(0);
	selected_points[2].resize(0);



	vector<unsigned int> selected_points_indices;

	// contains the index to the point with the minimal distance in the rest of the points
	// meaning in the points that have not been selected yet
	// given an i \in N
	// min_distance_indices[i] ... 	is the index in points that has the min distance to
	//								selected_points[i]
	// min_distances[i] ...			is the the min distance to
	//								selected_points[i]
	vector<unsigned int> min_distance_indices;
	vector<double> min_distances;

//	printMatrix(d);


	drawOnePoint(points, selected_points,selected_points_indices,random_start_ind,d);

//	printMatrix(d);
	auto mm = std::min_element(d[random_start_ind].begin(),d[random_start_ind].end());
	unsigned int max_min_ind = std::distance(d[random_start_ind].begin(), mm);

//	cout << "mm: " << *mm << endl;
//	cout << "random_start_ind: " << random_start_ind << endl;
//	cout << "max_min_ind: " << max_min_ind << endl;

	min_distances.push_back(d[random_start_ind][max_min_ind]);
	min_distance_indices.push_back(max_min_ind);


	// sn ... number of currently sampled points
	for(unsigned int sn = 1; sn < K; sn++){
//		printPoints(selected_points);

		// get the maximal distance of the minimal distances
		auto mm = std::max_element(min_distances.begin(),min_distances.end());
		unsigned int min_distances_ind = std::distance(min_distances.begin(), mm);
		unsigned int selected_point_ind = min_distance_indices[min_distances_ind];


//		cout << "sn = " << sn << endl;
//		cout << "selected_points_indices: \n" << selected_points_indices << endl;
//		cout << "minDistances: \n" << min_distances << endl;
//		cout << "min_distance_indices: \n" << min_distance_indices << endl;
//		cout << "mm: " << *mm << endl;
//		cout << "selected_point_ind: " << selected_point_ind << endl;
//
//		cout << "----------------------------" << endl;


		drawOnePoint(points, selected_points, selected_points_indices, selected_point_ind, d);

		// update the minimal distances
		for(unsigned int i = 0; i < sn; i++){
			// if we selected the point that was previously the point that had the minimal distance
			// to the given point we have to recalculate the distance
			if(min_distance_indices[i] == selected_point_ind){

//				double min_value = *std::min_element(d[selected_points_indices[i]].begin(),d[selected_points_indices[i]].end());

				auto mm2 = std::min_element(d[selected_points_indices[i]].begin(),d[selected_points_indices[i]].end());
				unsigned int min_ind = std::distance(d[selected_points_indices[i]].begin(), mm2);

				min_distances[i] = *mm2;
				min_distance_indices[i] = min_ind;
			}
		}


		auto mm2 = std::min_element(d[max_min_ind].begin(),d[max_min_ind].end());
		unsigned int min_ind = std::distance(d[max_min_ind].begin(), mm2);

		min_distances.push_back(*mm2);
		min_distance_indices.push_back(min_ind);



	}

	return(selected_points_indices);
}





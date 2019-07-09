/*
 * sampling.h
 *
 *  Created on: 02.03.2019
 *      Author: willy
 */

#include <vector>


#ifndef FIRSTLOWERBOUNDRELATIONOFPOSANDNEG_SRC_SAMPLING_H_
#define FIRSTLOWERBOUNDRELATIONOFPOSANDNEG_SRC_SAMPLING_H_


unsigned int maxMinDistance(const std::vector< std::vector<double> > &points,
							const unsigned int ind);

void drawOnePoint(	const std::vector< std::vector<double> > &points,
					std::vector< std::vector<double> > &selected_points,
					vector<unsigned int> &selected_points_indices,
					unsigned int ind,
					std::vector< std::vector<double> > &d);

vector<unsigned int> euclidean_farthest_point_sampling( const std::vector< std::vector<double> > &points,
										std::vector< std::vector<double> > &selected_points,
										const unsigned int K);




#endif /* FIRSTLOWERBOUNDRELATIONOFPOSANDNEG_SRC_SAMPLING_H_ */

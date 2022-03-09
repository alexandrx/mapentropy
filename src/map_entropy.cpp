/*
 * Copyright 2015-2019 Autoware Foundation. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 PointCloud map entropy calculation and filtering

 Alexander Carballo, 2019/08/20
 */

/*
 * This code is originally written by David Droeschel & Jan Razlaw, and is based on their paper:
 * Razlaw, J., Droeschel, D., Holz, D., & Behnke, S. "Evaluation of registration methods for sparse 3D laser scans." 
 * In 2015 IEEE European Conference on Mobile Robots (ECMR), September 2015.
 * http://www.ais.uni-bonn.de/papers/ECMR_2015_Razlaw.pdf
 */

#include <omp.h>

#include <iostream>
#include <unistd.h> //for sync
#include <string>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/geometry.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>
#include <pcl/search/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/filters/passthrough.h>

#include <map_entropy/entropy_point_types.h>
#include <map_entropy/map_entropy.h>

namespace map_entropy {

MapEntropy::MapEntropy():
	stepSize_(1),
	radius_(0.3),
	minNeighbors_(15),
	minLimit_(-std::numeric_limits<double>::infinity()),
	maxLimit_(std::numeric_limits<double>::infinity()),
	filterField_(std::string("")),
	punishSolitaryPoints_(true),
	withPlaneVariance_(true),
	saveASCII_(false),
	passFilter_(false),
	hasXField_(false),
	hasYField_(false),
	hasZField_(false),
	hasIntensityField_(false),
	hasEntropyField_(false),
	hasIntensityEntropyField_(false),
	hasPlaneVarianceField_(false),
	entropySum_(0.f),
	planeVarianceSum_(0.f),
	lonelyPoints_(0),
	meanMapEntropy_(0.f),
	meanPlaneVariance_(0.f),
	cloud_xyz_(new pcl::PointCloud<pcl::PointXYZ>),
	cloud_entropy_(new pcl::PointCloud< PointXYZWithEntropy >),
	inputCloud_(new pcl::PCLPointCloud2),
	outputCloud_(new pcl::PCLPointCloud2),
	grid_divide_(false),
	divide_only_(false),
	grid_size_(0.f),
    output_dir_("./"),
    name_prefix_("grid_")
{
	fileIndices_.clear();
	fileNames_.clear();
	ptcldGrid_.clear();
}

double MapEntropy::computeEntropy( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud )
{
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covarianceMatrixNormalized = Eigen::Matrix3f::Identity();;

	// estimate the XYZ centroid and the normalized covariance matrix
	pcl::compute3DCentroid (*cloud, centroid);
	pcl::computeCovarianceMatrixNormalized (*cloud, centroid, covarianceMatrixNormalized);

	// compute the determinant and return the entropy
	double determinant = static_cast<double>((( 2 * M_PI * M_E) * covarianceMatrixNormalized).determinant());

	return 0.5f*log(determinant);
}

double MapEntropy::computeEntropy( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud, const std::vector<int>& indices )
{
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covarianceMatrixNormalized = Eigen::Matrix3f::Identity();;

	// estimate the XYZ centroid and the normalized covariance matrix
	pcl::compute3DCentroid (*cloud, centroid);
	pcl::computeCovarianceMatrixNormalized (*cloud, indices, centroid, covarianceMatrixNormalized);

	// compute the determinant and return the entropy
	double determinant = static_cast<double>((( 2 * M_PI * M_E) * covarianceMatrixNormalized).determinant());

	return 0.5f*log(determinant);
}

// double computeEntropy( pcl::PointCloud< pcl::PointXYZI >::Ptr cloud, double &intensity_entropy)
//{
// 	Eigen::Vector4f centroid;
// 	Eigen::Matrix3f covarianceMatrixNormalized = Eigen::Matrix3f::Identity();
// 	if (cloud.size() < 3) {
// 		return std::numeric_limits<double>::infinity();
// 	}
// 	pcl::PointCloud<pcl::PointXYZ> cloud2;
// 	pcl::copyPointCloud(cloud, cloud2);
// 	// estimate the XYZ centroid and the normalized covariance matrix
// 	pcl::compute3DCentroid (*cloud2, centroid);
// 	pcl::computeCovarianceMatrixNormalized (*cloud2, centroid, covarianceMatrixNormalized);
// 	//Have to create a histogram to compute the probabilities and then compute entropy
// 	//// estimate mean intensity
// 	// double mean_intensity = 0.0;
// 	// for (auto p : cloud.points) {
// 	// 	mean_intensity += p.intensity;
// 	// }
// 	// mean_intensity /= cloud.size();
// 	// double cov_intensity = 0.0;
// 	// for (auto p : cloud.points) {
// 	// 	cov_intensity += (p.intensity - mean_intensity)*(p.intensity - mean_intensity);
// 	// }
// 	// cov_intensity /= cloud.size();
// 	// compute the determinant and return the entropy
// 	double determinant = static_cast<double>((( 2 * M_PI * M_E) * covarianceMatrixNormalized).determinant());
// 	//intensity_entropy = 0.5f*log(( 2 * M_PI * M_E) * cov_intensity);
// 	return 0.5f*log(determinant);
// }

double MapEntropy::computePlaneVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud )
{
	double meanDistTopQuarter = 0;

	if (cloud->points.size() < 3) {
		return std::numeric_limits<double>::infinity();
	}

	std::vector<double> sortedDistances;

	// fit plane using RANSAC
	pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
	pcl::PointIndices::Ptr inliers (new pcl::PointIndices);

	pcl::SACSegmentation< pcl::PointXYZ > seg;
	seg.setOptimizeCoefficients (true);
	seg.setModelType (pcl::SACMODEL_PLANE);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setDistanceThreshold (0.01); //(0.005);
	seg.setInputCloud (cloud);
	seg.segment (*inliers, *coefficients);

	if( inliers->indices.size() < 3 ){
		//std::cout << "ERROR: Could not estimate a planar model for the given subset of points." << std::endl;
		meanDistTopQuarter = std::numeric_limits<double>::infinity();
	}else{
		// compute the distances of the points to the plane
		for( size_t i = 0; i < cloud->points.size(); ++i ){
			double distancePointToPlane = (cloud->points[i].x * coefficients->values[0]) + (cloud->points[i].y * coefficients->values[1]) +( cloud->points[i].z * coefficients->values[2] ) + coefficients->values[3];
			sortedDistances.push_back(fabs(distancePointToPlane));
		}
		// sort distances
		std::sort(sortedDistances.begin(), sortedDistances.end());

		// compute mean of quartile that contains the largest distances
		int quarterOfArray = sortedDistances.size() / 4;
		for( size_t i = quarterOfArray * 3; i < sortedDistances.size(); i++ ){
			meanDistTopQuarter += sortedDistances[i];
		}
		meanDistTopQuarter /= static_cast<double> (quarterOfArray);
	}

	return meanDistTopQuarter;
}

double MapEntropy::computePlaneVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud, const std::vector<int>& indices )
{
	double meanDistTopQuarter = 0;

	std::vector<double> sortedDistances;

	// fit plane using RANSAC
	pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
	pcl::PointIndices::Ptr inliers (new pcl::PointIndices);

	pcl::SACSegmentation< pcl::PointXYZ > seg;
	seg.setOptimizeCoefficients (true);
	seg.setModelType (pcl::SACMODEL_PLANE);
	seg.setMethodType (pcl::SAC_RANSAC);
	seg.setDistanceThreshold (0.01); //(0.005);
	seg.setInputCloud (cloud);
	pcl::PointIndices::Ptr pIndices (new pcl::PointIndices);
	pIndices->indices = indices;
	seg.setIndices (pIndices);
	seg.segment (*inliers, *coefficients);

	if( inliers->indices.size() < 3 ){
		//std::cout << "Could not estimate a planar model for the given subset of points." << std::endl;
		meanDistTopQuarter = std::numeric_limits<double>::infinity();
	}else{
		// compute the distances of the points to the plane
		for( size_t i = 0; i < cloud->points.size(); ++i ){
			double distancePointToPlane = (cloud->points[i].x * coefficients->values[0]) + (cloud->points[i].y * coefficients->values[1]) +( cloud->points[i].z * coefficients->values[2] ) + coefficients->values[3];
			sortedDistances.push_back(fabs(distancePointToPlane));
		}
		// sort distances
		std::sort(sortedDistances.begin(), sortedDistances.end());

		// compute mean of quartile that contains the largest distances
		int quarterOfArray = sortedDistances.size() / 4;
		for( size_t i = quarterOfArray * 3; i < sortedDistances.size(); i++ ){
			meanDistTopQuarter += sortedDistances[i];
		}
		meanDistTopQuarter /= static_cast<double> (quarterOfArray);
	}

	return meanDistTopQuarter;
}

void MapEntropy::computeEntroyAndVariance(const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud )
{
	pcl::KdTreeFLANN< pcl::PointXYZ > kdtree;
	kdtree.setInputCloud (cloud);
	int me_lonelyPointsSum = 0;
	double me_entropySum = 0.0;
	double me_planeVarianceSum = 0.0;
	cloud_entropy_.reset(new pcl::PointCloud< PointXYZWithEntropy >);

	#pragma omp parallel reduction (+:me_entropySum, me_planeVarianceSum, me_lonelyPointsSum)
	{
		#pragma omp for schedule(dynamic)
		for (size_t i = 0; i < cloud->points.size(); i += stepSize_ ) {

			// print status
			if( i % (cloud->points.size()/20) == 0 ){
				int percent = i * 100 / cloud->points.size();
				std::cout << percent << " %" << std::endl;
			}

			// search for neighbors in radius
			std::vector<int> pointIdxRadiusSearch;
			std::vector<float> pointRadiusSquaredDistance;
			int numberOfNeighbors = kdtree.radiusSearch (cloud->points[i], radius_, pointIdxRadiusSearch, pointRadiusSquaredDistance);

			// compute values if enough neighbors found
			double localEntropy = 0;
			double localPlaneVariance = 0;
			if( numberOfNeighbors > minNeighbors_ || !punishSolitaryPoints_ ){

				// save neighbors in localCloud
				pcl::PointCloud< pcl::PointXYZ >::Ptr localCloud (new pcl::PointCloud< pcl::PointXYZ >);

				for( size_t iz = 0; iz < pointIdxRadiusSearch.size(); ++iz ){
					localCloud->points.push_back(cloud->points[ pointIdxRadiusSearch[iz] ] );
				}

				// compute entropy and plane variance
				localEntropy = computeEntropy(localCloud);
				localPlaneVariance = computePlaneVariance(localCloud);
			}else{
				localEntropy = std::numeric_limits<double>::infinity();
				localPlaneVariance = std::numeric_limits<double>::infinity();
				me_lonelyPointsSum++;
			}

			// save values in new point
			PointXYZWithEntropy p;
			p.x = cloud->points[i].x;
			p.y = cloud->points[i].y;
			p.z = cloud->points[i].z;

			if (std::isfinite(localPlaneVariance)){
				me_planeVarianceSum += localPlaneVariance;
				p.planeVariance = static_cast<float>(localPlaneVariance);
			}else{
				// handle cases where no value could be computed
				if( !punishSolitaryPoints_ ){
					p.planeVariance = 0;
				}else{
					me_planeVarianceSum += radius_;
					p.planeVariance = static_cast<float>(radius_);
				}
			}
			if (std::isfinite(localEntropy)){
				me_entropySum += localEntropy;
				p.entropy = static_cast<float>(localEntropy);
			}else{
				// handle cases where no value could be computed
				p.entropy = 0;
			}
			// p.planeEntropy = p.entropy * (1.0 - p.planeVariance); 

			// add new point to output cloud
			#pragma omp critical
			{
				cloud_entropy_->push_back( p );
			}
		}  // parallel for
	}  // pragma parallel

	// compute mean
	meanMapEntropy_ = me_entropySum / (static_cast<double>(cloud->points.size() / stepSize_));
	meanPlaneVariance_ = me_planeVarianceSum / (static_cast<double>(cloud->points.size() / stepSize_));

	std::cout << "--- " << std::endl;
	std::cout << "Mean Map Entropy is " << meanMapEntropy_ << std::endl;
	if (withPlaneVariance_)
		std::cout << "Mean Plane Variance is " << meanPlaneVariance_ << std::endl;

	//concatenate the pointcloud fields
	outputCloud_.reset (new pcl::PCLPointCloud2);
	pcl::toPCLPointCloud2 (*cloud_entropy_, *outputCloud_);

	pcl::PCLPointCloud2::Ptr outputCloudConcat (new pcl::PCLPointCloud2);
	pcl::concatenateFields(*inputCloud_, *outputCloud_, *outputCloudConcat);
	outputCloud_.reset (new pcl::PCLPointCloud2);
	pcl::copyPointCloud(*outputCloudConcat, *outputCloud_);

	int pointsActuallyUsed = (cloud->points.size() / stepSize_) - me_lonelyPointsSum;
	
	if( punishSolitaryPoints_ && (pointsActuallyUsed < me_lonelyPointsSum) ){
		std::cout << "Used more solitary than not-solitary points to compute the values. You should consider changing the parameters." << std::endl;
	}
}

/** \note: This code is under development, not tested. 
 *  The idea is to speed up the entropy calculation using a grid instead of the costly radius search.
 */
void MapEntropy::computeEntroyAndVarianceFromGrid()
{
	int me_lonelyPointsSum = 0;
	double me_entropySum = 0.0;
	double me_planeVarianceSum = 0.0;
	size_t points_count = 0;

	const auto maxCols = ptcldGrid_.divX();
	const auto maxRows = ptcldGrid_.divY();

	//#pragma omp parallel reduction (+:me_entropySum, me_planeVarianceSum, me_lonelyPointsSum)
	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic)
		for (int idx = 0; idx < ptcldGrid_.size(); idx++) {
			std::pair<size_t, size_t> rc = ptcldGrid_.getRowCol(idx);
			size_t row = rc.first;
			size_t col = rc.second;
			pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);

			// print status
			if( idx % (ptcldGrid_.size()/20) == 0 ){
				int percent = idx * 100 / ptcldGrid_.size();
				std::cout << percent << " %" << std::endl;
			}

			//if no points in this grid, skip it
			if (!ptcldGrid_.at(row,col).localCloud->points.size()) {
				//std::cout << "row: " << row << ", col: " << col << ", id: " << idx << " is empty!!! skipping this cell" << std::endl;
				continue;
			}
			std::cout << "row: " << row << ", col: " << col << ", id: " << idx << ", finding neighboring cells.. maxRows: " << maxRows << ", maxCols: " << maxCols << std::endl;

			//find neighboring grid cells
			if (row == 0 && col == 0) { // bottom left corner
				// self + 3 neighbors 
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col+1).localCloud;
			} else if (row == 0 && col == maxCols-1) { // bottom right corner
				// self + 3 neighbors 
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col).localCloud;
			} else if (row == maxRows-1 && col == 0) { // top left corner
				// self + 3 neighbors
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row-1,col).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col+1).localCloud;
			} else if (row == maxRows-1 && col == maxCols-1) { // top right corner
				// self + 3 neighbors
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row-1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
			} else if (row == 0 && col < maxCols-1) { // bottom border
				// self + 5 neighbors
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col+1).localCloud;
			} else if (row < maxRows-1 && col == 0) { //left most border
				// self + 5 neighbors
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row-1,col).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col+1).localCloud;
			} else if (row == maxRows-1 && col < maxCols-1) { // top border
				// self + 5 neighbors
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row-1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col+1).localCloud;
			} else if (row < maxRows-1 && col == maxCols-1) { //right most border
				// self + 5 neighbors
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row-1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col).localCloud;
			} else {
				// the general case
				//concatenate the neighbors
				*cloud += *ptcldGrid_.at(row-1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col).localCloud;
				*cloud += *ptcldGrid_.at(row-1,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row,col).localCloud;
				*cloud += *ptcldGrid_.at(row,col+1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col-1).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col).localCloud;
				*cloud += *ptcldGrid_.at(row+1,col+1).localCloud;
			}
			std::cout << "row: " << row << ", col: " << col << ", id: " << idx << ", localCloud size: " << ptcldGrid_.at(row,col).localCloud->points.size() << ", merged cloud size: " << cloud->points.size() << std::endl;

			//create the KdTree for the local neighborhood
			pcl::KdTreeFLANN< pcl::PointXYZ > kdtree;
			kdtree.setInputCloud (cloud);
			std::cout << "kdtree with " << cloud->points.size() << " points" << std::endl;

			//process the concatenated cloud
			for (size_t i = 0; i < ptcldGrid_.at(row,col).localCloud->points.size(); i += stepSize_ ) {
				// search for neighbors in radius
				std::vector<int> pointIdxRadiusSearch;
				std::vector<float> pointRadiusSquaredDistance;
				int numberOfNeighbors = kdtree.radiusSearch (ptcldGrid_[idx].localCloud->points[i], radius_, pointIdxRadiusSearch, pointRadiusSquaredDistance);

				// compute values if enough neighbors found
				double localEntropy = 0;
				double localPlaneVariance = 0;
				if( numberOfNeighbors > minNeighbors_ || !punishSolitaryPoints_ ){
					// save neighbors in localCloud
					pcl::PointCloud< pcl::PointXYZ >::Ptr localCloud (new pcl::PointCloud< pcl::PointXYZ >);

					for( size_t iz = 0; iz < pointIdxRadiusSearch.size(); ++iz ){
						localCloud->points.push_back(ptcldGrid_[idx].localCloud->points[ pointIdxRadiusSearch[iz] ] );
					}

					// compute entropy and plane variance
					localEntropy = computeEntropy(localCloud);
					localPlaneVariance = computePlaneVariance(localCloud);
				}else{
					localEntropy = std::numeric_limits<double>::infinity();
					localPlaneVariance = std::numeric_limits<double>::infinity();
					ptcldGrid_[idx].lonelyPointsSum++;
				}

				// save values in new point
				PointXYZWithEntropy p;
				p.x = ptcldGrid_[idx].localCloud->points[i].x;
				p.y = ptcldGrid_[idx].localCloud->points[i].y;
				p.z = ptcldGrid_[idx].localCloud->points[i].z;

				if (std::isfinite(localPlaneVariance)){
					ptcldGrid_[idx].planeVarianceSum += localPlaneVariance;
					p.planeVariance = static_cast<float>(localPlaneVariance);
				}else{
					// handle cases where no value could be computed
					if( !punishSolitaryPoints_ ){
						p.planeVariance = 0;
					}else{
						ptcldGrid_[idx].planeVarianceSum += radius_;
						p.planeVariance = static_cast<float>(radius_);
					}
				}
				if (std::isfinite(localEntropy)){
					ptcldGrid_[idx].entropySum += localEntropy;
					p.entropy = static_cast<float>(localEntropy);
				}else{
					// handle cases where no value could be computed
					p.entropy = 0;
				}

				// add new point to output cloud
				#pragma omp critical
				{
					ptcldGrid_[idx].cloud_entropy->push_back( p );
				}
			} // for local cloud points
			std::cout << "..... done!" << std::endl;
		}  // parallel for
	}  // pragma parallel

	for (int idx = 0; idx < ptcldGrid_.size(); idx++) {
		me_entropySum += ptcldGrid_[idx].entropySum;
		ptcldGrid_[idx].meanEntropy = ptcldGrid_[idx].entropySum / (static_cast<double>(ptcldGrid_[idx].cloud_entropy->points.size() / stepSize_)); 
		ptcldGrid_[idx].meanPlaneVariance = ptcldGrid_[idx].planeVarianceSum / (static_cast<double>(ptcldGrid_[idx].cloud_entropy->points.size() / stepSize_)); 
		me_planeVarianceSum += ptcldGrid_[idx].planeVarianceSum;
		me_lonelyPointsSum += ptcldGrid_[idx].lonelyPointsSum;
		points_count += ptcldGrid_[idx].localCloud->points.size();
	}

	// compute mean
	meanMapEntropy_ = me_entropySum / (static_cast<double>(points_count / stepSize_));
	meanPlaneVariance_ = me_planeVarianceSum / (static_cast<double>(points_count / stepSize_));

	std::cout << "--- " << std::endl;
	std::cout << "Mean Map Entropy is " << meanMapEntropy_ << std::endl;
	if (withPlaneVariance_)
		std::cout << "Mean Plane Variance is " << meanPlaneVariance_ << std::endl;

	for (int idx = 0; idx < ptcldGrid_.size(); idx++) {
		//concatenate the pointcloud fields
		pcl::toPCLPointCloud2 (*ptcldGrid_[idx].cloud_entropy, *ptcldGrid_[idx].outputCloud);
		pcl::PCLPointCloud2::Ptr outputCloudConcat (new pcl::PCLPointCloud2);
		pcl::concatenateFields(*ptcldGrid_[idx].inputCloud, *ptcldGrid_[idx].outputCloud, *outputCloudConcat);
		ptcldGrid_[idx].outputCloud.reset (new pcl::PCLPointCloud2);
		pcl::copyPointCloud(*outputCloudConcat, *ptcldGrid_[idx].outputCloud);

		int pointsActuallyUsed = (ptcldGrid_[idx].localCloud->points.size() / stepSize_) - me_lonelyPointsSum;
		if( punishSolitaryPoints_ && (pointsActuallyUsed < me_lonelyPointsSum) ){
			std::cout << "Used more solitary than not-solitary points to compute the values. You should consider changing the parameters." << std::endl;
		}
	}
}


int MapEntropy::readPointCloud()
{
	pcl::PCDReader pcdR;
	const size_t LARGE_POINTCLOUD_ = 500000; // experimental value

	if (!fileNames_.size() || !fileIndices_.size()) 
	{
		std::cout << "Couldn't read any file. Please provide the name of at least one PCD file as input." << std::endl;
		return 0;
	}
	
	inputCloud_.reset (new pcl::PCLPointCloud2);
	if (pcdR.read(fileNames_.at(0), *inputCloud_, origin_, orientation_, version_) < 0)
	{
		std::cout << "ERROR: Couldn't read file \"" << fileNames_.at(0) << "\"" << std::endl;
		return 0;
	}
	
	//Check input pointcloud fields
	hasXField_ = ((pcl::getFieldIndex(*inputCloud_, std::string("x"))) > -1);
	hasYField_ = ((pcl::getFieldIndex(*inputCloud_, std::string("y"))) > -1);
	hasZField_ = ((pcl::getFieldIndex(*inputCloud_, std::string("z"))) > -1);
	hasIntensityField_ = ((pcl::getFieldIndex(*inputCloud_, std::string("intensity"))) > -1);
	hasEntropyField_ = ((pcl::getFieldIndex(*inputCloud_, std::string("entropy"))) > -1) || ((pcl::getFieldIndex(*inputCloud_, std::string("range_entropy"))) > -1);
	hasIntensityEntropyField_ = ((pcl::getFieldIndex(*inputCloud_, std::string("intensity_entropy"))) > -1);
	hasPlaneVarianceField_ = ((pcl::getFieldIndex(*inputCloud_, std::string("planeVariance"))) > -1);

	std::cout << "Input Cloud fields: [";
	for (auto& f : inputCloud_->fields)
	{
		std::cout << f.name << ", ";
	}
	std::cout << "]  Input cloud width: " << inputCloud_->width << ", height: " << inputCloud_->height << std::endl;
	
	if (!(hasXField_ && hasYField_ && hasZField_))
	{
		std::cout << "ERROR: Incorrect format in file, does not have XYZ fields. Aborting!" << std::endl;
		return 0;
	}
   
    //for passFilter the pcd file needs to have entropy field
    if (passFilter_ && ((pcl::getFieldIndex(*inputCloud_, filterField_)) <= -1)) {
    	std::cout << "ERROR: \"-passFilter\" option requires an input pointcloud with \"" << filterField_ << "\" field, the provided PCD file doesn't have such field. Aborting!" << std::endl;
    	return 0;
    }

    //check the input pointcloud size, if too large advise to use -griddive param.
    if (inputCloud_->width * inputCloud_->height > LARGE_POINTCLOUD_ && !grid_divide_) {
    	std::cout << "Your pointcloud is larger than " << LARGE_POINTCLOUD_ << ", using the \"-griddivide\" (with corresponding \"-grid_size\") parameter is advised!" << std::endl;
    }

    if (grid_divide_ || divide_only_) {
    	pointCloudToGrid();
    }

    return 1;
}

int MapEntropy::writePointCloud()
{
	pcl::PCDWriter pcdW;

	if ((grid_divide_ || divide_only_) && ptcldGrid_.size()) {
	    for (int i = 0; i < ptcldGrid_.size(); i++) {
	        if (ptcldGrid_[i].inputCloud && ptcldGrid_[i].inputCloud->width * ptcldGrid_[i].inputCloud->height > 0) {
	        	//save each grid cells' pointcloud2 
	        	if (pcdW.write(ptcldGrid_[i].filename, *ptcldGrid_[i].inputCloud, origin_, orientation_, (!saveASCII_)) >= 0) {
	            	std::cout << "Wrote " << (ptcldGrid_[i].inputCloud->width * ptcldGrid_[i].inputCloud->height) << " points to file \"" << ptcldGrid_[i].filename << "\"" << std::endl;
	        	}
	        }
	        if (ptcldGrid_[i].outputCloud && ptcldGrid_[i].outputCloud->width * ptcldGrid_[i].outputCloud->height > 0) {
	            //and save the cell's entropy cloud
	            std::string saveDestination = ptcldGrid_[i].filename;
	            std::string postname = "_entropy__r_" + std::to_string(radius_) + "_s_" + std::to_string(stepSize_) + "_sol_" + std::to_string(punishSolitaryPoints_) + "_minN_" + std::to_string(minNeighbors_);
				postname = postname + "__meanEnt_" + std::to_string(ptcldGrid_[i].meanEntropy) + "__meanPlanVar_" + std::to_string(ptcldGrid_[i].meanPlaneVariance) + "__lonely_" + std::to_string(ptcldGrid_[i].lonelyPointsSum) + ".";
				saveDestination.replace(saveDestination.find_last_of("."),1,postname);
				if (pcdW.write(saveDestination, *ptcldGrid_[i].outputCloud, origin_, orientation_, (!saveASCII_)) >= 0) {
					//std::cout << "Output saved as \"" << saveDestination << "\"" << std::endl;
					std::cout << "Wrote grid cell " << i << " to file \"" << saveDestination << "\"" << std::endl;
				}
	        }
	    }

		sync();  //forces fs sync
    } else if (!grid_divide_ && !divide_only_) {
		// save output cloud in the directory of the input cloud
		std::string saveDestination = fileNames_[0];
		if (!passFilter_) {
			std::string postname = "_entropy__r_" + std::to_string(radius_) + "_s_" + std::to_string(stepSize_) + "_sol_" + std::to_string(punishSolitaryPoints_) + "_minN_" + std::to_string(minNeighbors_);
			postname = postname + "__meanEnt_" + std::to_string(meanMapEntropy_) + "__meanPlanVar_" + std::to_string(meanPlaneVariance_) + "__lonely_" + std::to_string(lonelyPoints_) + ".";
			saveDestination.replace(saveDestination.find_last_of("."),1,postname);
		} else {
			std::string postname = "_filtered__f_" + filterField_ + "_min_" + std::to_string(minLimit_) + "_max_" + std::to_string(maxLimit_) + ".";
			saveDestination.replace(saveDestination.find_last_of("."),1,postname);
		}
		if ( outputCloud_->data.size() > 0 )
		{
			std::cout << "Output Cloud fields: [";
			for (auto& f : outputCloud_->fields)
			{
				std::cout << f.name << ", ";
			}
			std::cout << "]" << std::endl;

			if (pcdW.write(saveDestination, *outputCloud_, origin_, orientation_, (!saveASCII_)) >= 0) {
				std::cout << "Output saved as \"" << saveDestination << "\"" << std::endl;
			}

			sync();  //forces fs sync
		} else {
			std::cout << "ERROR: Empty cloud. Saving error." << std::endl;
			return 0;
		}
	}
	return 1; 
}

/**
  * This code was inspired and partially based on the pcd_grid_divider implementation in Autoware
  * Original author: Yuki Kitsukawa (yuki.kitsukawa@tier4.jp)
  * This version can divide a pointcloud2 into grids, no need to know beforehand the pointcloud type
  */
void MapEntropy::pointCloudToGrid()
{
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();
    int fieldIdxX = -1;
    int fieldIdxY = -1;
    int fieldIdxZ = -1;

    fieldIdxX = pcl::getFieldIndex(*inputCloud_, "x");
    fieldIdxY = pcl::getFieldIndex(*inputCloud_, "y");
    fieldIdxZ = pcl::getFieldIndex(*inputCloud_, "z");
    const auto offX = inputCloud_->fields[fieldIdxX].offset;
    const auto offY = inputCloud_->fields[fieldIdxY].offset;
    const auto offZ = inputCloud_->fields[fieldIdxZ].offset;


    // Search minimum and maximum points along x and y axis.
    //for (sensor_msgs::PointCloud2ConstIterator<float> it(*inputCloud_, "x"); it != it.end(); ++it) { //<-- cannot use sensor_msgs
    // it[0] is x and it[1] is y
    const std::uint8_t* inptr = &inputCloud_->data[0];
    const auto incr = inputCloud_->point_step;
    for (size_t point = 0; point < inputCloud_->height*inputCloud_->width; point++, inptr += incr) {
    	const float* x = reinterpret_cast<const float*>(inptr + offX);
    	const float* y = reinterpret_cast<const float*>(inptr + offY);
    	//float* z = reinterpret_cast<float*>(inptr + offZ);
	    if (*x < min_x) {
	        min_x = *x;
	    }
	    if (*x > max_x) {
	        max_x = *x;
	    }
	    if (*y < min_y) {
	        min_y = *y;
	    }
	    if (*y > max_y) {
	        max_y = *y;
	    }
    } // for 

    // Find minimum and maximum boundary
    int min_x_b = grid_size_ * static_cast<int>(floor(min_x / grid_size_));
    int max_x_b = grid_size_ * static_cast<int>(floor(max_x / grid_size_) + 1);
    int min_y_b = grid_size_ * static_cast<int>(floor(min_y / grid_size_));
    int max_y_b = grid_size_ * static_cast<int>(floor(max_y / grid_size_) + 1);
    ptcldGrid_.cellSize(grid_size_);
    ptcldGrid_.setBounds(min_x, min_y, max_x, max_y);

    if (ptcldGrid_.cellCount() == 1) {
    	std::cout << "Grid division not possible: cell size of " << grid_size_ << " x " << grid_size_ << " is larger than the pointcloud bountaries [" << (max_x - min_x) << " x " << (max_y - min_y) << "]" << std::endl;
    	grid_divide_ = false;
    	return;
    }

    std::cout << "Dividing the input pointcloud into " << grid_size_ << " x " << grid_size_ << " cells, number of cells in the grid: " << ptcldGrid_.cellCount() << std::endl;

    // Define filename, lower/upper bound of every grid
    for (int row = 0; row < ptcldGrid_.divY(); row++) {
      for (int col = 0; col < ptcldGrid_.divX(); col++) {
        int id = ptcldGrid_.divX() * row + col;
        ptcldGrid_.at(row, col).grid_id = id;
        ptcldGrid_.at(row, col).grid_id_x = col;
        ptcldGrid_.at(row, col).grid_id_y = row;
        ptcldGrid_.at(row, col).lower_bound_x = ptcldGrid_.minX() + grid_size_ * col;
        ptcldGrid_.at(row, col).lower_bound_y = ptcldGrid_.minY() + grid_size_ * row;
        ptcldGrid_.at(row, col).upper_bound_x = ptcldGrid_.minX() + grid_size_ * (col + 1);
        ptcldGrid_.at(row, col).upper_bound_y = ptcldGrid_.minY() + grid_size_ * (row + 1);
        ptcldGrid_.at(row, col).filename = output_dir_ + name_prefix_ + "_" + std::to_string(grid_size_) + "_" +
        	std::to_string(id) + "_" + std::to_string(ptcldGrid_.at(row, col).lower_bound_x) + "_" + std::to_string(ptcldGrid_.at(row, col).lower_bound_y) + ".pcd";

 		//create the local XYZ pointcloud for computing entropy
        ptcldGrid_.at(row, col).localCloud.reset(new pcl::PointCloud<pcl::PointXYZ>);
  		ptcldGrid_.at(row, col).localCloud->header   = inputCloud_->header;
  		ptcldGrid_.at(row, col).localCloud->width    = 0;
  		ptcldGrid_.at(row, col).localCloud->height   = inputCloud_->height;
  		ptcldGrid_.at(row, col).localCloud->is_dense = inputCloud_->is_dense == 1;

  		//create the input pointcloud2 used to write the submaps later
        ptcldGrid_.at(row, col).inputCloud.reset (new pcl::PCLPointCloud2);
        ptcldGrid_.at(row, col).inputCloud->header = inputCloud_->header;
        ptcldGrid_.at(row, col).inputCloud->height = inputCloud_->height;
        ptcldGrid_.at(row, col).inputCloud->width = 0;  // value is updated later
        ptcldGrid_.at(row, col).inputCloud->is_bigendian = inputCloud_->is_bigendian;
        ptcldGrid_.at(row, col).inputCloud->point_step = inputCloud_->point_step;
        ptcldGrid_.at(row, col).inputCloud->row_step = 0;  // value is updated later
        ptcldGrid_.at(row, col).inputCloud->is_dense = inputCloud_->is_dense;
        //copy all the fields as they are
        ptcldGrid_.at(row, col).inputCloud->fields = inputCloud_->fields;

        //create the cloud for entropy information
        ptcldGrid_.at(row, col).cloud_entropy.reset(new pcl::PointCloud< PointXYZWithEntropy >);
      }  // for col
    }  // for row


    // Assign all points to appropriate grid according to their x/y value
    inptr = &inputCloud_->data[0];
    //for (sensor_msgs::PointCloud2ConstIterator<float> it(*inputCloud_, "x"); it != it.end(); ++it) {
    // it[0] is x and it[1] is y and it[2] is z
    for (size_t point = 0; point < inputCloud_->height*inputCloud_->width; point++) {
    	const float* x = reinterpret_cast<const float*>(inptr + offX);
    	const float* y = reinterpret_cast<const float*>(inptr + offY);
    	const float* z = reinterpret_cast<const float*>(inptr + offZ);
    	// int idx = static_cast<int>(floor((*x - static_cast<float>(min_x_b)) / grid_size_));
     //  	int idy = static_cast<int>(floor((*y - static_cast<float>(min_y_b)) / grid_size_));
     //  	int id = idy * div_x + idx;

      	// add one point to the XYZ cloud
      	pcl::PointXYZ tmp = pcl::PointXYZ(*x, *y, *z);
      	ptcldGrid_.atCoords(*x, *y).localCloud->points.push_back(tmp);
      	// and add one point to the PointCloud2 cloud
    	for (size_t i = 0; i < incr; i++) {
    		ptcldGrid_.atCoords(*x, *y).inputCloud->data.push_back(*inptr++);
    	}
    	ptcldGrid_.atCoords(*x, *y).inputCloud->width++;
    	ptcldGrid_.atCoords(*x, *y).inputCloud->row_step = ptcldGrid_.atCoords(*x, *y).inputCloud->width * ptcldGrid_.atCoords(*x, *y).inputCloud->point_step;
    }
}

void MapEntropy::filterPointCloud()
{
	std::cout << "Filtering by \"" << filterField_ << "\"  field."  << std::endl;
	pcl::PassThrough<pcl::PCLPointCloud2> entropyFilter;
	entropyFilter.setInputCloud(inputCloud_);
	entropyFilter.setKeepOrganized(false);
	entropyFilter.setFilterFieldName(filterField_);
	entropyFilter.setFilterLimits(minLimit_, maxLimit_);
	// outputCloud_.reset (new pcl::PCLPointCloud2);
	entropyFilter.filter(*outputCloud_);
	
}

void MapEntropy::computeEntropyOrFilter(const pcl::PCLPointCloud2::Ptr& inputCloud)
{
	inputCloud_.reset (new pcl::PCLPointCloud2);
	pcl::copyPointCloud(*inputCloud, *inputCloud_);
	pcl::fromPCLPointCloud2 (*inputCloud_, *cloud_xyz_);
	computeEntroyAndVariance(cloud_xyz_);
	if (passFilter_) {
		filterPointCloud();
	}
}

void MapEntropy::computeEntropyOrFilter()
{
	if (!passFilter_) {
		if (grid_divide_) {
			computeEntroyAndVarianceFromGrid();
		} else {
			pcl::fromPCLPointCloud2 (*inputCloud_, *cloud_xyz_);
			computeEntroyAndVariance(cloud_xyz_);
		}
	} else {
		filterPointCloud();
	}
}

}  // namespace

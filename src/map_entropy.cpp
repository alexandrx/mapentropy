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

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/common/geometry.h>
#include <pcl/common/centroid.h>
#include <pcl/common/transforms.h>
#include <pcl/common/time.h>
#include <pcl/console/parse.h>
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

double MapEntropy::computeEntropy( const pcl::PointCloud< pcl::PointXYZ >::Ptr cloud )
{
	Eigen::Vector4f centroid;
	Eigen::Matrix3f covarianceMatrixNormalized = Eigen::Matrix3f::Identity();

	if (cloud->size() < 3) {
		return std::numeric_limits<double>::infinity();
	}

	// estimate the XYZ centroid and the normalized covariance matrix
	pcl::compute3DCentroid (*cloud, centroid);
	pcl::computeCovarianceMatrixNormalized (*cloud, centroid, covarianceMatrixNormalized);

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

double MapEntropy::computePlaneVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr cloud )
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
	seg.segment (*inliers, *coefficients);

	if( inliers->indices.size() < 3 ){
		//PCL_ERROR ("Could not estimate a planar model for the given subset of points.");
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

// void MapEntropy::computeEntroyAndVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr cloud )
// {
// 	pcl::KdTreeFLANN< pcl::PointXYZ > kdtree;
// 	kdtree.setInputCloud (cloud);

// 	pcl::copyPointCloud(*cloud, *cloudi_entropy_); //only initializes XYZI fields, we initialize entropy and planevariance here

// 	#pragma omp parallel reduction (+:entropySum_, planeVarianceSum_, lonelyPoints_)
// 	{
// 		#pragma omp for schedule(dynamic)
// 		for (size_t i = 0; i < cloud->points.size(); i += stepSize_ ) {

// 			// print status
// 			if( i % (cloud->points.size()/20) == 0 ){
// 				int percent = i * 100 / cloud->points.size();
// 				std::cout << percent << " %" << std::endl;
// 			}

// 			// search for neighbors in radius
// 			std::vector<int> pointIdxRadiusSearch;
// 			std::vector<float> pointRadiusSquaredDistance;
// 			int numberOfNeighbors = kdtree.radiusSearch (cloud->points[i], radius_, pointIdxRadiusSearch, pointRadiusSquaredDistance);

// 			// compute values if enough neighbors found
// 			double localEntropy = 0;
// 			double localPlaneVariance = 0;
// 			if( numberOfNeighbors > minNeighbors_ || !punishSolitaryPoints_ ){
// 				// save neighbors in localCloud
// 				pcl::PointCloud< pcl::PointXYZ >::Ptr localCloud (new pcl::PointCloud< pcl::PointXYZ >);

// 				for( size_t iz = 0; iz < pointIdxRadiusSearch.size(); ++iz ){
// 					localCloud->points.push_back(cloud->points[ pointIdxRadiusSearch[iz] ] );
// 				}

// 				// compute entropy and plane variance
// 				localEntropy = computeEntropy(localCloud);
// 				if (withPlaneVariance_) {
// 					localPlaneVariance = computePlaneVariance(localCloud);
// 				}
// 			}else{
// 				localEntropy = std::numeric_limits<double>::infinity();
// 				if (withPlaneVariance_)
// 					localPlaneVariance = std::numeric_limits<double>::infinity();
// 				lonelyPoints_++;
// 			}

// 			float out_entropy = std::numeric_limits<float>::infinity();
// 			float out_planevariance = std::numeric_limits<float>::infinity();
// 			if (withPlaneVariance_) {
// 				if (std::isfinite(localPlaneVariance)){
// 					planeVarianceSum_ += localPlaneVariance;
// 					out_planevariance = static_cast<float>(localPlaneVariance);
// 				}else{
// 					// handle cases where no value could be computed
// 					if( !punishSolitaryPoints_ ){
// 						out_planevariance = 0;
// 					}else{
// 						planeVarianceSum_ += radius_;
// 						out_planevariance = static_cast<float>(radius_);
// 					}
// 				}
// 			}
// 			if (std::isfinite(localEntropy)){
// 				entropySum_ += localEntropy;
// 				out_entropy = static_cast<float>(localEntropy);
// 			} else if( !punishSolitaryPoints_ ){
// 				out_entropy = 0;
// 			} 

// 			// add new point to output cloud
// 			// #pragma omp critical
// 			// {
// 			// 	cloudi_entropy_->push_back( p );
// 			// }
// 			cloudi_entropy_->points[i].entropy = out_entropy;
// 			cloudi_entropy_->points[i].planeVariance = out_planevariance;
// 		}
// 	}

// 	// compute mean
// 	double meanMapEntropy = entropySum_ / (static_cast<double>(cloud->points.size() / stepSize_));
// 	double meanPlaneVariance = 0.0;
// 	if (withPlaneVariance_)
// 		meanPlaneVariance = planeVarianceSum_ / (static_cast<double>(cloud->points.size() / stepSize_));

// 	std::cout << "--- " << std::endl;
// 	std::cout << "Mean Map Entropy is " << meanMapEntropy << std::endl;
// 	if (withPlaneVariance_)
// 		std::cout << "Mean Plane Variance is " << meanPlaneVariance << std::endl;
		
// 	//concatenate the pointcloud fields
// 	outputCloud_.reset (new pcl::PCLPointCloud2);
// 	if (hasIntensityField_) {
// 		pcl::toPCLPointCloud2 (*cloudi_entropy_, *outputCloud_);
// 	} else {
// 		pcl::copyPointCloud(*cloudi_entropy_, *cloud_entropy_);
// 		pcl::toPCLPointCloud2 (*cloud_entropy_, *outputCloud_);
// 	}

// 	// for (auto& f : outputCloud_->fields)
// 	// {
// 	// 	PCL_INFO("field: %s\n", f.name.c_str());
// 	// }
// 	pcl::PCLPointCloud2::Ptr outputCloudConcat (new pcl::PCLPointCloud2);
// 	pcl::concatenateFields(*outputCloud_, *inputCloud_, *outputCloudConcat);
// 	outputCloud_.reset (new pcl::PCLPointCloud2);
// 	pcl::copyPointCloud(*outputCloudConcat, *outputCloud_);
	
// 	//std::cout << "Used " << entropyTimer.getTime() << " milliseconds to compute values for " << cloud_xyz->points.size() << " points." << std::endl;

// 	int pointsActuallyUsed = (cloud->points.size() / stepSize_) - lonelyPoints_;

// 	if( punishSolitaryPoints_ && (pointsActuallyUsed < lonelyPoints_) ){
// 		std::cout << "Used more solitary than not-solitary points to compute the values. You should consider changing the parameters." << std::endl;
// 	}
// }

void MapEntropy::computeEntroyAndVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr cloud )
{
	pcl::KdTreeFLANN< pcl::PointXYZ > kdtree;
	kdtree.setInputCloud (cloud);

	pcl::copyPointCloud(*cloud, *cloudi_entropy_); //only initializes XYZI fields, we initialize entropy and planevariance here

	int me_lonelyPointsSum = 0;
	float me_entropySum = 0.0;
	float me_planeVarianceSum = 0.0;

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
				if (withPlaneVariance_) {
					localPlaneVariance = computePlaneVariance(localCloud);
				}
			}else{
				localEntropy = std::numeric_limits<double>::infinity();
				if (withPlaneVariance_)
					localPlaneVariance = std::numeric_limits<double>::infinity();
				me_lonelyPointsSum++;
			}

			float out_entropy = std::numeric_limits<float>::infinity();
			float out_planevariance = std::numeric_limits<float>::infinity();
			if (withPlaneVariance_) {
				if (std::isfinite(localPlaneVariance)){
					me_planeVarianceSum += localPlaneVariance;
					out_planevariance = static_cast<float>(localPlaneVariance);
				}else{
					// handle cases where no value could be computed
					if( !punishSolitaryPoints_ ){
						out_planevariance = 0;
					}else{
						me_planeVarianceSum += radius_;
						out_planevariance = static_cast<float>(radius_);
					}
				}
			}
			if (std::isfinite(localEntropy)){
				me_entropySum += localEntropy;
				out_entropy = static_cast<float>(localEntropy);
			} else if( !punishSolitaryPoints_ ){
				out_entropy = 0;
			}
			cloudi_entropy_->points[i].entropy = out_entropy;
			cloudi_entropy_->points[i].planeVariance = out_planevariance;
		}
	}
	lonelyPoints_ = me_lonelyPointsSum;
        entropySum_ = me_entropySum;
        planeVarianceSum_ = me_planeVarianceSum;

	// compute mean
	double meanMapEntropy = entropySum_ / (static_cast<double>(cloud->points.size() / stepSize_));
	double meanPlaneVariance = 0.0;
	if (withPlaneVariance_)
		meanPlaneVariance = planeVarianceSum_ / (static_cast<double>(cloud->points.size() / stepSize_));

	std::cout << "--- " << std::endl;
	std::cout << "Mean Map Entropy is " << meanMapEntropy << std::endl;
	if (withPlaneVariance_)
		std::cout << "Mean Plane Variance is " << meanPlaneVariance << std::endl;
		
	//concatenate the pointcloud fields
	outputCloud_.reset (new pcl::PCLPointCloud2);
	if (hasIntensityField_) {
		pcl::toPCLPointCloud2 (*cloudi_entropy_, *outputCloud_);
	} else {
		pcl::copyPointCloud(*cloudi_entropy_, *cloud_entropy_);
		pcl::toPCLPointCloud2 (*cloud_entropy_, *outputCloud_);
	}

	// for (auto& f : outputCloud_->fields)
	// {
	// 	PCL_INFO("field: %s\n", f.name.c_str());
	// }
	pcl::PCLPointCloud2::Ptr outputCloudConcat (new pcl::PCLPointCloud2);
	pcl::concatenateFields(*outputCloud_, *inputCloud_, *outputCloudConcat);
	outputCloud_.reset (new pcl::PCLPointCloud2);
	pcl::copyPointCloud(*outputCloudConcat, *outputCloud_);
	
	//std::cout << "Used " << entropyTimer.getTime() << " milliseconds to compute values for " << cloud_xyz->points.size() << " points." << std::endl;

	int pointsActuallyUsed = (cloud->points.size() / stepSize_) - lonelyPoints_;

	if( punishSolitaryPoints_ && (pointsActuallyUsed < lonelyPoints_) ){
		std::cout << "Used more solitary than not-solitary points to compute the values. You should consider changing the parameters." << std::endl;
	}
}

MapEntropy::MapEntropy():
	stepSize_(1),
	radius_(0.3),
	minNeighbors_(15),
	minLimit_(std::numeric_limits<double>::lowest()),
	maxLimit_(std::numeric_limits<double>::infinity()),
	filterField_(std::string("entropy")),
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
	cloud_xyz_(new pcl::PointCloud<pcl::PointXYZ>),
	cloud_xyzi_(new pcl::PointCloud<pcl::PointXYZI>),
	cloud_entropy_(new pcl::PointCloud< PointXYZWithEntropy >),
	cloudi_entropy_(new pcl::PointCloud< PointXYZIWithEntropy >)
{

	fileIndices_.clear();
	fileNames_.clear();
}

MapEntropy::MapEntropy( int argc, char** argv ) :
	stepSize_(1),
	radius_(0.3),
	minNeighbors_(15),
	minLimit_(std::numeric_limits<double>::lowest()),
	maxLimit_(std::numeric_limits<double>::infinity()),
	filterField_(std::string("entropy")),
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
	cloud_xyz_(new pcl::PointCloud<pcl::PointXYZ>),
	cloud_xyzi_(new pcl::PointCloud<pcl::PointXYZI>),
	cloud_entropy_(new pcl::PointCloud< PointXYZWithEntropy >),
	cloudi_entropy_(new pcl::PointCloud< PointXYZIWithEntropy >)
{
	pcl::console::parse_argument (argc, argv, "-stepsize", stepSize_);
    pcl::console::parse_argument (argc, argv, "-radius", radius_);
	pcl::console::parse_argument (argc, argv, "-punishSolitaryPoints", punishSolitaryPoints_);
    pcl::console::parse_argument (argc, argv, "-minNeighbors", minNeighbors_);
    pcl::console::parse_argument (argc, argv, "-planevariance", withPlaneVariance_);
    pcl::console::parse_argument (argc, argv, "-ascii", saveASCII_);
    pcl::console::parse_argument (argc, argv, "-passfilter", passFilter_);
	pcl::console::parse_argument (argc, argv, "-filterField", filterField_);
    pcl::console::parse_argument (argc, argv, "-minLimit", minLimit_);
    pcl::console::parse_argument (argc, argv, "-maxLimit", maxLimit_);

	fileIndices_.clear();
	fileNames_.clear();
}

int MapEntropy::readPointCloud( int argc, char** argv )
{
	pcl::PCDReader pcdR;

	// get pointcloud
	fileIndices_ = pcl::console::parse_file_extension_argument (argc, argv, ".pcd");
	for (auto i : fileIndices_) {
		fileNames_.push_back( std::string( argv[i] ) );
	}

	if (!fileNames_.size()) 
	{
		PCL_ERROR ("Couldn't read any file. Please provide the name of at least one PCD file as input.\n");
		return 0;
	}
	
	inputCloud_.reset (new pcl::PCLPointCloud2);
	if (pcdR.read(argv[fileIndices_.at (0)], *inputCloud_, origin_, orientation_, version_) < 0)
	{
		PCL_ERROR ("Couldn't read file.\n");
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

	PCL_INFO("Input Cloud fields: [");
	for (auto& f : inputCloud_->fields)
	{
		PCL_INFO("%s, ", f.name.c_str());
	}
	PCL_INFO("]\n");
	//std::cout << "X field: " << hasXField_ << ", Y field: " << hasYField_ << ", Z field: " << hasZField_ << ", intensity field: " << hasIntensityField_ << ", entropy field: " << hasEntropyField_ << ", plane variance field: " << hasPlaneVarianceField_ << std::endl;
	
	if (!(hasXField_ && hasYField_ && hasZField_))
	{
		PCL_ERROR ("Incorrect format in file, does not have XYZ fields. Aborting!\n");
		return 0;
	}
   
    //for passFilter the pcd file needs to have entropy field
    if (passFilter_ && ((pcl::getFieldIndex(*inputCloud_, filterField_)) <= -1)) {
    	PCL_ERROR("\"-passFilter\" option requires an input pointcloud with \"%s\" field, the provided PCD file doesn't have such field. Aborting!\n", filterField_.c_str());
    	return 0;
    }

    return 1;
}

int MapEntropy::writePointCloud()
{
	pcl::PCDWriter pcdW;

	// save output cloud in the directory of the input cloud
	std::string saveDestination = fileNames_[0];
	if (!passFilter_) {
		std::string postname = "_entropy__r_" + std::to_string(radius_) + "_s_" + std::to_string(stepSize_) + "_sol_" + std::to_string(punishSolitaryPoints_) + "_minN_" + std::to_string(minNeighbors_)  + ".";
		saveDestination.replace(saveDestination.find_last_of("."),1,postname);
	} else {
		std::string postname = "_filtered__f_" + filterField_ + "_min_" + std::to_string(minLimit_) + "_max_" + std::to_string(maxLimit_) + ".";
		saveDestination.replace(saveDestination.find_last_of("."),1,postname);
	}
	if ( outputCloud_->data.size() > 0 )
	{
		PCL_INFO("Output Cloud fields: [");
		for (auto& f : outputCloud_->fields)
		{
			PCL_INFO("%s, ", f.name.c_str());
		}
		PCL_INFO("]\n");

		pcdW.write(saveDestination, *outputCloud_, origin_, orientation_, (!saveASCII_));
		std::cout << "Output saved as \"" << saveDestination << "\"" << std::endl;
	} else {
		PCL_ERROR ("Empty cloud. Saving error.\n");
		return 0;
	}
	return 1; 
}

void MapEntropy::filterPointCloud()
{
	PCL_WARN("Filtering by \"%s\" field.\n", filterField_.c_str());
	pcl::PassThrough<pcl::PCLPointCloud2> entropyFilter;
	entropyFilter.setInputCloud(inputCloud_);
	entropyFilter.setKeepOrganized(false);
	entropyFilter.setFilterFieldName(filterField_);
	entropyFilter.setFilterLimits(minLimit_, maxLimit_);
	outputCloud_.reset (new pcl::PCLPointCloud2);
	entropyFilter.filter(*outputCloud_);
	
}

void MapEntropy::computeEntropyOrFilter()
{
	pcl::StopWatch entropyTimer;
	entropyTimer.reset();

	if (!passFilter_) {
		// if (hasIntensityField_) {
		// 	pcl::fromPCLPointCloud2 (*inputCloud_, *cloud_xyzi_);
		// 	computeEntroyAndVariance(cloud_xyzi_);
		// } else {
			pcl::fromPCLPointCloud2 (*inputCloud_, *cloud_xyz_);
			//pcl::copyPointCloud(*cloud_xyz_, *cloud_xyzi_);
			// computeEntroyAndVariance(cloud_xyzi_);
			computeEntroyAndVariance(cloud_xyz_);
		// }
	} else {
		filterPointCloud();
	}

	std::cout << "Used " << entropyTimer.getTime() << " milliseconds to filter values. Input size " << inputCloud_->width * inputCloud_->height << " points, output size " << outputCloud_->width * outputCloud_->height << " points." << std::endl;
}

/*
int main( int argc, char** argv ) 
{
	pcl::PCDReader pcdR;
	pcl::PCDWriter pcdW;
	pcl::PCLPointCloud2::Ptr inputCloud;
	pcl::PCLPointCloud2::Ptr outputCloud;
	Eigen::Vector4f origin;
    Eigen::Quaternionf orientation;
    int version;

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_xyz (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_xyzi (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud< PointXYZWithEntropy >::Ptr cloud_entropy (new pcl::PointCloud< PointXYZWithEntropy >);
	pcl::PointCloud< PointXYZWithEntropy >::Ptr cloudi_entropy (new pcl::PointCloud< PointXYZIWithEntropy >);
	bool hasXField = false;
	bool hasYField = false;
	bool hasZField = false;
	bool hasIntensityField = false;
	bool hasEntropyField = false;
	bool hasIntensityEntropyField = false;
	bool hasPlaneVarianceField = false;

	double entropySum = 0.f;
	double planeVarianceSum = 0.f;
	int lonelyPoints = 0;
	
	// get parameters if given
	int stepSize = 1;
	double radius = 0.3;
	int minNeighbors = 15;
	double minLimit = std::numeric_limits<double>::lowest();
	double maxLimit = std::numeric_limits<double>::infinity();
	std::string filterField("entropy");

	pcl::console::parse_argument (argc, argv, "-stepsize", stepSize);
    pcl::console::parse_argument (argc, argv, "-radius", radius);
	bool punishSolitaryPoints = pcl::console::find_switch (argc, argv, "-punishSolitaryPoints");
    pcl::console::parse_argument (argc, argv, "-minNeighbors", minNeighbors);
    bool withPlaneVariance = pcl::console::find_switch (argc, argv, "-planevariance");
    bool saveASCII = pcl::console::find_switch (argc, argv, "-ascii");
    bool passFilter = pcl::console::find_switch (argc, argv, "-passfilter");
	pcl::console::parse_argument (argc, argv, "-filterField", filterField);
    pcl::console::parse_argument (argc, argv, "-minLimit", minLimit);
    pcl::console::parse_argument (argc, argv, "-maxLimit", maxLimit);

	// get pointcloud
	std::vector<int> fileIndices = pcl::console::parse_file_extension_argument (argc, argv, ".pcd");
	
	inputCloud.reset (new pcl::PCLPointCloud2);
	if (pcdR.read(argv[fileIndices.at (0)], *inputCloud, origin, orientation, version) < 0)
	{
		PCL_ERROR ("Couldn't read file.\n");
		return (-1);
	}
	
	//Check input pointcloud fields
	hasXField = ((pcl::getFieldIndex(*inputCloud, std::string("x"))) > -1);
	hasYField = ((pcl::getFieldIndex(*inputCloud, std::string("y"))) > -1);
	hasZField = ((pcl::getFieldIndex(*inputCloud, std::string("z"))) > -1);
	hasIntensityField = ((pcl::getFieldIndex(*inputCloud, std::string("intensity"))) > -1);
	hasEntropyField = ((pcl::getFieldIndex(*inputCloud, std::string("entropy"))) > -1) || ((pcl::getFieldIndex(*inputCloud, std::string("range_entropy"))) > -1);
	hasIntensityEntropyField = ((pcl::getFieldIndex(*inputCloud, std::string("intensity_entropy"))) > -1);
	hasPlaneVarianceField = ((pcl::getFieldIndex(*inputCloud, std::string("planeVariance"))) > -1);
	
	if (!(hasXField && hasYField && hasZField))
	{
		PCL_ERROR ("Incorrect format in file, does not have XYZ fields. Aborting!\n");
		return (-1);
	}
	
      
	// if (pcl::io::loadPCDFile< pcl::PointXYZ> (argv[fileIndices[0]], *inputCloud) == -1)
	// {
	// 	PCL_ERROR ("Couldn't read file.\n");
	// 	return (-1);
	// }

   
    //for passFilter the pcd file needs to have entropy field
    if (passFilter && ((pcl::getFieldIndex(*inputCloud, filterField)) <= -1)) {
    	PCL_ERROR("\"-passFilter\" option requires an input pointcloud with \"%s\" field, the provided PCD file doesn't have such field. Aborting!\n", filterField.c_str());
    	return (-1);
    }


	pcl::StopWatch entropyTimer;
	entropyTimer.reset();

	if (!passFilter)
	{
		if (hasEntropyField || hasPlaneVarianceField)
		{
			PCL_WARN("Input point cloud already has \"entropy\" and \"planeVariance\" fields and they will be overwritten.\nAbort now if this is not desired.\n");
		}
		
		std::cout << "Stepsize = " << stepSize << std::endl;
		std::cout << "Radius for neighborhood search = " << radius << std::endl;
		if( !punishSolitaryPoints ){
			std::cout << "Paper version" << std::endl;
		}else{
			std::cout << "Punishing solitary points. \nMinimal number of neighbors that have to be found = " << minNeighbors << std::endl;
		}
		std::cout << "--- " << std::endl;
		
	    //Copy input point cloud
	    if (hasIntensityField) {
	    	pcl::fromPCLPointCloud2 (*inputCloud, *cloudi_xyz);
	    } else {
    		pcl::fromPCLPointCloud2 (*inputCloud, *cloud_xyz);
    	}
	
		pcl::KdTreeFLANN< pcl::PointXYZ> kdtree;
		kdtree.setInputCloud (cloud_xyz);

		#pragma omp parallel reduction (+:entropySum, planeVarianceSum, lonelyPoints)
		{
			#pragma omp for schedule(dynamic)
			for (size_t i = 0; i < cloud_xyz->points.size(); i += stepSize ) {

				// print status
				if( i % (cloud_xyz->points.size()/20) == 0 ){
					int percent = i * 100 / cloud_xyz->points.size();
					std::cout << percent << " %" << std::endl;
				}

				// search for neighbors in radius
				std::vector<int> pointIdxRadiusSearch;
				std::vector<float> pointRadiusSquaredDistance;
				int numberOfNeighbors = kdtree.radiusSearch (cloud_xyz->points[i], radius, pointIdxRadiusSearch, pointRadiusSquaredDistance);

				// compute values if enough neighbors found
				double localEntropy = 0;
				double localPlaneVariance = 0;
				if( numberOfNeighbors > minNeighbors || !punishSolitaryPoints ){

					// save neighbors in localCloud
					pcl::PointCloud< pcl::PointXYZ>::Ptr localCloud (new pcl::PointCloud< pcl::PointXYZ>);

					for( size_t iz = 0; iz < pointIdxRadiusSearch.size(); ++iz ){
						localCloud->points.push_back(cloud_xyz->points[ pointIdxRadiusSearch[iz] ] );
					}

					// compute entropy and plane variance
					localEntropy = computeEntropy< pcl::PointXYZ >(localCloud);
					if (withPlaneVariance)
						localPlaneVariance = computePlaneVariance< pcl::PointXYZ >(localCloud);
				}else{
					localEntropy = std::numeric_limits<double>::infinity();
					if (withPlaneVariance)
						localPlaneVariance = std::numeric_limits<double>::infinity();
					lonelyPoints++;
				}

				// save values in new point
				PointXYZWithEntropy p;
				p.x = cloud_xyz->points[i].x;
				p.y = cloud_xyz->points[i].y;
				p.z = cloud_xyz->points[i].z;
				p.planeVariance = std::numeric_limits<float>::infinity();
				p.entropy = std::numeric_limits<float>::infinity();
			
				if (withPlaneVariance) {
					if (pcl_isfinite(localPlaneVariance)){
						planeVarianceSum += localPlaneVariance;
						p.planeVariance = static_cast<float>(localPlaneVariance);
					}else{
						// handle cases where no value could be computed
						if( !punishSolitaryPoints ){
							p.planeVariance = 0;
						}else{
							planeVarianceSum += radius;
							p.planeVariance = static_cast<float>(radius);
						}
					}
				}
				if (pcl_isfinite(localEntropy)){
					entropySum += localEntropy;
					p.entropy = static_cast<float>(localEntropy);
				} else if( !punishSolitaryPoints ){
					p.entropy = 0;
				}/*else{
					// handle cases where no value could be computed
					p.entropy = 0;
				}

				// add new point to output cloud
				#pragma omp critical
				{
					cloud_entropy->push_back( p );
				}
			}
		}

		// compute mean
		double meanMapEntropy = entropySum / (static_cast<double>(cloud_xyz->points.size() / stepSize));
		double meanPlaneVariance = 0.0;
		if (withPlaneVariance)
			meanPlaneVariance = planeVarianceSum / (static_cast<double>(cloud_xyz->points.size() / stepSize));

		std::cout << "--- " << std::endl;
		std::cout << "Mean Map Entropy is " << meanMapEntropy << std::endl;
		if (withPlaneVariance)
			std::cout << "Mean Plane Variance is " << meanPlaneVariance << std::endl;
			
		//concatenate the pointcloud fields
		PCL_INFO("checkpoint 1\n");
		outputCloud.reset (new pcl::PCLPointCloud2);
		pcl::toPCLPointCloud2 (*cloud_entropy, *outputCloud);
		for (auto& f : outputCloud->fields)
		{
			PCL_INFO("field: %s\n", f.name.c_str());
		}
		PCL_INFO("checkpoint 2a\n");
		pcl::PCLPointCloud2::Ptr outputCloudConcat (new pcl::PCLPointCloud2);
		pcl::concatenateFields(*outputCloud, *inputCloud, *outputCloudConcat);
		for (auto& f : outputCloudConcat->fields)
		{
			PCL_INFO("field: %s\n", f.name.c_str());
		}
		PCL_INFO("checkpoint 2b\n");
		
		std::cout << "Used " << entropyTimer.getTime() << " milliseconds to compute values for " << cloud_xyz->points.size() << " points." << std::endl;

		int pointsActuallyUsed = (cloud_xyz->points.size() / stepSize) - lonelyPoints;
	
		if( punishSolitaryPoints && (pointsActuallyUsed < lonelyPoints) ){
			std::cout << "Used more solitary than not-solitary points to compute the values. You should consider changing the parameters." << std::endl;
		}
	
		// save output cloud in the directory of the input cloud
		std::string saveDestination = argv[fileIndices[0]];
		saveDestination.replace(saveDestination.find_last_of("."),1,"_entropy.");
		if ( outputCloudConcat->data.size() > 0 )
		{
			PCL_INFO("checkpoint 3\n");
			pcdW.write(saveDestination, *outputCloudConcat, origin, orientation, (!saveASCII));
			PCL_INFO("checkpoint 4\n");
		}
		else
			PCL_ERROR ("Empty cloud. Saving error.\n");
	}
	else 
	{
		PCL_WARN("Filtering by \"%s\" field.\n", filterField.c_str());
    	pcl::PassThrough<pcl::PCLPointCloud2> entropyFilter;
    	entropyFilter.setInputCloud(inputCloud);
    	entropyFilter.setKeepOrganized(false);
    	entropyFilter.setFilterFieldName(filterField);
    	entropyFilter.setFilterLimits(minLimit, maxLimit);
    	outputCloud.reset (new pcl::PCLPointCloud2);
    	entropyFilter.filter(*outputCloud);
    	
		std::string saveDestination = argv[fileIndices[0]];
		saveDestination.replace(saveDestination.find_last_of("."),1,"_filtered.");
		if ( outputCloud->data.size() > 0 )
		{
			pcdW.write(saveDestination, *outputCloud, origin, orientation, (!saveASCII));
		}
		else
			PCL_ERROR ("Empty cloud. Saving error.\n");
			
		std::cout << "Used " << entropyTimer.getTime() << " milliseconds to filter values. Input size " << inputCloud->width * inputCloud->height << " points, output size " << outputCloud->width * outputCloud->height << " points. Output saved as \"" << saveDestination << "\"" << std::endl;
	}


	return 0;
}
*/

}  // namespace

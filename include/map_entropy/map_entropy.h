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
 Compute mean map entropy and add to pointcloud files
 @author Alexander Carballo (Nagoya University)
 */

#ifndef MAP_ENTROPY_H_
#define MAP_ENTROPY_H_

#include <string>
#include <vector>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/PCLPointCloud2.h>

#include <map_entropy/entropy_point_types.h>

namespace map_entropy {

class MapEntropy 
{
public:
    MapEntropy();
    MapEntropy( int argc, char** argv );

    double computeEntropy( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud );
    double computeEntropy( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud, const std::vector<int>& indices );

    double computePlaneVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud );
    double computePlaneVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud, const std::vector<int>& indices );

    void computeEntroyAndVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud );

    void filterPointCloud();
    void computeEntropyOrFilter();

    int readPointCloud( int argc, char** argv );
    int writePointCloud();

private:
    int stepSize_;
    double radius_;
    int minNeighbors_;
    double minLimit_;
    double maxLimit_;
    std::string filterField_;
    bool punishSolitaryPoints_;
    bool withPlaneVariance_;
    bool saveASCII_;
    bool passFilter_;
    bool hasXField_;
    bool hasYField_;
    bool hasZField_;
    bool hasIntensityField_;
    bool hasEntropyField_;
    bool hasIntensityEntropyField_;
    bool hasPlaneVarianceField_;
    double entropySum_;
    double planeVarianceSum_;
    int lonelyPoints_;
    double meanMapEntropy_;
    double meanPlaneVariance_;
    pcl::PCLPointCloud2::Ptr inputCloud_;
    pcl::PCLPointCloud2::Ptr outputCloud_;
    Eigen::Vector4f origin_;
    Eigen::Quaternionf orientation_;
    int version_;
    pcl::PointCloud< pcl::PointXYZ >::Ptr cloud_xyz_;
    pcl::PointCloud< pcl::PointXYZI >::Ptr cloud_xyzi_;
    pcl::PointCloud< PointXYZWithEntropy >::Ptr cloud_entropy_;
    pcl::PointCloud< PointXYZIWithEntropy >::Ptr cloudi_entropy_;
    std::vector<int> fileIndices_;
    std::vector< std::string > fileNames_;
};

}  // namespace

#endif // MAP_ENTROPY_H_
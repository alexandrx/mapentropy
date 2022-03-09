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
 PointCloud with entropy and plane variance fields
 @author Alexander Carballo (Nagoya University)
 */
#ifndef ENTROPY_POINT_TYPES_H_
#define ENTROPY_POINT_TYPES_H_

#include <pcl/point_types.h>
#include <pcl/common/common.h>

namespace map_entropy {

struct PointXYZWithEntropy
{
    PCL_ADD_POINT4D;                    // preferred way of adding a XYZ + padding
    float entropy;
    float planeVariance;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW     // make sure our new allocators are aligned
} EIGEN_ALIGN16;                            // enforce SSE padding for correct memory alignment

struct PointXYZIWithEntropy
{
    PCL_ADD_POINT4D;                    // preferred way of adding a XYZ + padding
    float intensity;                    // laser intensity re ading
    float entropy;                     // entropy of the point
    float planeVariance;             // plane variance of the point
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // make sure our new allocators are aligned
} EIGEN_ALIGN16;                            // enforce SSE padding for correct memory alignment

struct PointTypeWithIntensityEntropy
{
    PCL_ADD_POINT4D;                    // preferred way of adding a XYZ + padding
    float range_entropy;
    float intensity_entropy;
    float planeVariance;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW     // make sure our new allocators are aligned
} EIGEN_ALIGN16;                            // enforce SSE padding for correct memory alignment

}  // namespace

POINT_CLOUD_REGISTER_POINT_STRUCT (map_entropy::PointXYZWithEntropy,
        (float, x, x)
        (float, y, y)
        (float, z, z)
        (float, entropy, entropy)
        (float, planeVariance, planeVariance)
)

POINT_CLOUD_REGISTER_POINT_STRUCT (map_entropy::PointXYZIWithEntropy,
        (float, x, x)
        (float, y, y)
        (float, z, z)
        (float, intensity, intensity)
        (float, entropy, entropy)
        (float, planeVariance, planeVariance)
)

POINT_CLOUD_REGISTER_POINT_STRUCT (map_entropy::PointTypeWithIntensityEntropy,
        (float, x, x)
        (float, y, y)
        (float, z, z)
        (float, range_entropy, range_entropy)
        (float, intensity_entropy, intensity_entropy)
        (float, planeVariance, planeVariance)
)

#endif // ENTROPY_POINT_TYPES_H_
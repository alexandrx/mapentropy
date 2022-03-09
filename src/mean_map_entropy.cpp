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

 Alexander Carballo (Nagoya University)
 */

/*
 * This code is originally written by David Droeschel & Jan Razlaw, and is based on their paper:
 * Razlaw, J., Droeschel, D., Holz, D., & Behnke, S. "Evaluation of registration methods for sparse 3D laser scans." 
 * In 2015 IEEE European Conference on Mobile Robots (ECMR), September 2015.
 * http://www.ais.uni-bonn.de/papers/ECMR_2015_Razlaw.pdf
 */
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/common/time.h>
#include <pcl/console/parse.h>

#include <map_entropy/map_entropy.h>

class MeanMapEntropyPCL {
public:
    MeanMapEntropyPCL() {};
    MeanMapEntropyPCL( int argc, char** argv ); 
    void printUsage();
    void computeEntropy();

private:
    int _stepSize;
    double _radius;
    int _minNeighbors;
    double _minLimit;
    double _maxLimit;
    std::string _filterField;
    bool _punishSolitaryPoints;
    bool _withPlaneVariance;
    bool _saveASCII;
    bool _passFilter;
    bool _grid_divide;
    bool _divide_only;
    float _grid_size;
    std::string _output_dir;
    std::string _name_prefix;
    std::vector< std::string > _fileNames;
    std::vector<int> _fileIndices;
    map_entropy::MapEntropy _mmapEntropy;
};

MeanMapEntropyPCL::MeanMapEntropyPCL( int argc, char** argv ) :
    _stepSize(1),
    _radius(1.0),
    _minNeighbors(15),
    _minLimit(-std::numeric_limits<double>::infinity()),
    _maxLimit(std::numeric_limits<double>::infinity()),
    _filterField(std::string()),
    _punishSolitaryPoints(true),
    _withPlaneVariance(true),
    _saveASCII(false),
    _passFilter(false),
    _grid_divide(false),
    _divide_only(false),
    _output_dir(std::string()),
    _name_prefix(std::string()),
    _fileNames(std::vector< std::string >()),
    _fileIndices(std::vector<int>()),
    _mmapEntropy()
{
    //get the program arguments
    pcl::console::parse_argument (argc, argv, "-stepsize", _stepSize);
    pcl::console::parse_argument (argc, argv, "-radius", _radius);
    _punishSolitaryPoints = pcl::console::find_switch (argc, argv, "-punishSolitaryPoints");
    pcl::console::parse_argument (argc, argv, "-minNeighbors", _minNeighbors);
    _withPlaneVariance = pcl::console::find_switch (argc, argv, "-planevariance");
    _saveASCII = pcl::console::find_switch (argc, argv, "-ascii");
    _passFilter = pcl::console::find_switch (argc, argv, "-passfilter");
    pcl::console::parse_argument (argc, argv, "-filterField", _filterField);
    pcl::console::parse_argument (argc, argv, "-minLimit", _minLimit);
    pcl::console::parse_argument (argc, argv, "-maxLimit", _maxLimit);
    _grid_divide = pcl::console::find_switch (argc, argv, "-grid_divide");
    _divide_only = pcl::console::find_switch (argc, argv, "-divide_only");
    pcl::console::parse_argument (argc, argv, "-grid_size", _grid_size);
    pcl::console::parse_argument (argc, argv, "-grid_size", _grid_size);
    pcl::console::parse_argument (argc, argv, "-output_dir", _output_dir);
    pcl::console::parse_argument (argc, argv, "-name_prefix", _name_prefix);

    // get pointcloud file name(s)
    _fileIndices = pcl::console::parse_file_extension_argument (argc, argv, ".pcd");
    for (auto i : _fileIndices) {
        _fileNames.push_back( std::string( argv[i] ) );
    }

    if (!_fileNames.size()) 
    {
        PCL_ERROR ("Couldn't read any file. Please provide the name of at least one PCD file as input.\n");
        return;
    }

    _mmapEntropy.setStepSize(_stepSize);
    _mmapEntropy.setRadius(_radius);
    _mmapEntropy.seMinNeighbors(_minNeighbors);
    _mmapEntropy.setMinLimit(_minLimit);
    _mmapEntropy.setMaxLimit(_maxLimit);
    _mmapEntropy.setFilterField(_filterField);
    _mmapEntropy.setPunishSolitary(_punishSolitaryPoints);
    _mmapEntropy.setWithPlaneVariance(_withPlaneVariance);
    _mmapEntropy.setSaveAscii(_saveASCII);
    _mmapEntropy.setPassFilter(_passFilter);
    _mmapEntropy.setFileNames(_fileNames);
    _mmapEntropy.setFileIndices(_fileIndices);
    _mmapEntropy.setGridDivide(_grid_divide);
    _mmapEntropy.setDivideOnly(_divide_only);
    _mmapEntropy.setGridSize(_grid_size);
    _mmapEntropy.setOutputDir(_output_dir);
    _mmapEntropy.setNamePrefix(_name_prefix);

    PCL_INFO("This program is running with the following parameter values:\n");
    PCL_INFO("-stepsize= %s\n", std::to_string(_stepSize).c_str());
    PCL_INFO("-radius= %s\n", std::to_string(_radius).c_str());
    PCL_INFO("-minNeighbors= %s\n", std::to_string(_minNeighbors).c_str());
    PCL_INFO("-punishSolitaryPoints= %s\n", std::to_string(_punishSolitaryPoints).c_str());
    PCL_INFO("-planevariance= %s\n", std::to_string(_withPlaneVariance).c_str());
    PCL_INFO("-ascii= %s\n", std::to_string(_saveASCII).c_str());
    PCL_INFO("-passfilter= %s\n", std::to_string(_passFilter).c_str());
    PCL_INFO("-filterField= \"%s\"\n", _filterField.c_str());
    PCL_INFO("-minLimit= %s\n", std::to_string(_minLimit).c_str());
    PCL_INFO("-maxLimit= %s\n", std::to_string(_maxLimit).c_str());
    PCL_INFO("-grid_divide= %s\n", std::to_string(_grid_divide).c_str());
    PCL_INFO("-grid_size= %s\n", std::to_string(_grid_size).c_str());
    PCL_INFO("-output_dir= \"%s\"\n", _output_dir.c_str());
    PCL_INFO("-name_prefix= \"%s\"\n", _name_prefix.c_str());
}

void MeanMapEntropyPCL::printUsage()
{
    PCL_INFO("The arguments of this program are:\n");
    PCL_INFO("\"-stepsize\": step (increase) size to query the points, default %s\n", std::to_string(_stepSize).c_str());
    PCL_INFO("\"-radius\": size of the search radius (kdtree based), default value %s\n", std::to_string(_radius).c_str());
    PCL_INFO("\"-minNeighbors\": minimum number of neighbors in search radius, default value %s\n", std::to_string(_minNeighbors).c_str());
    PCL_INFO("\"-punishSolitaryPoints\": penalize points with very few neighbors, default value %s\n", std::to_string(_punishSolitaryPoints).c_str());
    PCL_INFO("\"-planevariance\": compute also the plane variance, default value %s\n", std::to_string(_withPlaneVariance).c_str());
    PCL_INFO("\"-ascii\": save output as ASCII pcd, default value %s\n", std::to_string(_saveASCII).c_str());
    PCL_INFO("\"-passfilter\": Do not compute entropy, instead bandpass filter points in the pointcloud, default value %s\n", std::to_string(_passFilter).c_str());
    PCL_INFO("\"-filterField\": Field name (ex., \"x\",\"y\",\"z\",\"intensity\", etc) to filter data, default value %s\n", _filterField.c_str());
    PCL_INFO("\"-minLimit\": minimum passing threshold, default value %s\n", std::to_string(_minLimit).c_str());
    PCL_INFO("\"-maxLimit\": maximum passing threshold, default value %s\n", std::to_string(_maxLimit).c_str());
    PCL_INFO("\"-grid_divide\": (specially for large pointclouds), divide the pointcloud into grids and process each grid and its 8 neighbors, default value %s\n", std::to_string(_grid_divide).c_str());
    PCL_INFO("\"-grid_size\": size of the grid cell, must be larger than or equal to the search radius, default value %s\n", std::to_string(_grid_size).c_str());
    PCL_INFO("\"-output_dir\": output directory name for the grid files, default value %s\n", _output_dir.c_str());
    PCL_INFO("\"-name_prefix\": prefix name for the grid files, default value %s\n", _name_prefix.c_str());
    PCL_INFO("The input PCD file(s)");
}

void MeanMapEntropyPCL::computeEntropy()
{
    pcl::StopWatch entropyTimer;
    entropyTimer.reset();
    if (_mmapEntropy.readPointCloud()) {
	if (!_divide_only) {
            _mmapEntropy.computeEntropyOrFilter();
	}
        _mmapEntropy.writePointCloud();
        PCL_INFO("Used %.6f seconds to filter values. Input size is %d points, output size is %d points.\n", entropyTimer.getTimeSeconds(), (_mmapEntropy.getInputPtCld()->width * _mmapEntropy.getInputPtCld()->height), (_mmapEntropy.getOutputPtCld()->width * _mmapEntropy.getOutputPtCld()->height));
    } else {
        printUsage();
    }
}

int main( int argc, char** argv ) 
{
    MeanMapEntropyPCL entropy(argc, argv);

    entropy.computeEntropy();

	return 0;
}

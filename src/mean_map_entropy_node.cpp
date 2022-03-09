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
#include<vector>
#include<string>

#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/PointField.h>
#include <sensor_msgs/point_field_conversion.h>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>

#include <map_entropy/map_entropy.h>

class MeanMapEntropyNode {
public:
    MeanMapEntropyNode() {};
    MeanMapEntropyNode( int argc, char** argv ); 
    void inputCallback(const sensor_msgs::PointCloud2::ConstPtr& inputCloud);

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
    std::string _input_topic;
    std::string _output_topic;
    std::string _output_frame;
    ros::Publisher _output_topic_pub;
    map_entropy::MapEntropy _mmapEntropy;
};

MeanMapEntropyNode::MeanMapEntropyNode( int argc, char** argv ) :
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
    _input_topic(std::string()),
    _output_topic(std::string()),
    _output_frame(std::string()),
    _mmapEntropy()
{
    ros::init(argc, argv, "mean_map_entropy");
    ros::NodeHandle nh;
    ros::NodeHandle private_nh("~");

    private_nh.getParam("stepsize", _stepSize);
    private_nh.getParam("radius", _radius);
    private_nh.getParam("punishSolitaryPoints", _punishSolitaryPoints);
    private_nh.getParam("minNeighbors", _minNeighbors);
    private_nh.getParam("planevariance", _withPlaneVariance);
    private_nh.getParam("passfilter", _passFilter);
    private_nh.getParam("filterField", _filterField);
    private_nh.getParam("minLimit", _minLimit);
    private_nh.getParam("maxLimit", _maxLimit);

    private_nh.getParam("input_topic", _input_topic);
    private_nh.getParam("output_topic", _output_topic);
    private_nh.getParam("output_frame", _output_frame);

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

    if (!_input_topic.size() || !_output_topic.size()) {
        ROS_ERROR("Must provide input and output topic names.");
        return;
    }

    //Subscribers and Publishers
    _output_topic_pub = nh.advertise<sensor_msgs::PointCloud2>(_output_topic, 10);
    ros::Subscriber input_topic_sub = nh.subscribe(_input_topic, 10, &MeanMapEntropyNode::inputCallback, this);

    ros::spin();
}

void MeanMapEntropyNode::inputCallback(const sensor_msgs::PointCloud2::ConstPtr& inputCloud) {
    pcl::PCLPointCloud2::Ptr in_pcl_pc2(new pcl::PCLPointCloud2);
    pcl_conversions::toPCL(*inputCloud, *in_pcl_pc2);
    
    _mmapEntropy.computeEntropyOrFilter(in_pcl_pc2);

    pcl::PCLPointCloud2::Ptr out_pcl_pc2(new pcl::PCLPointCloud2);
    sensor_msgs::PointCloud2::Ptr outputCloud_msg_ptr(new sensor_msgs::PointCloud2);
    out_pcl_pc2 = _mmapEntropy.getOutputPtCld();
    pcl_conversions::fromPCL(*out_pcl_pc2, *outputCloud_msg_ptr);

    //output the changed cloud
    if (_output_frame.size()) {
        outputCloud_msg_ptr->header.frame_id = _output_frame;
    }
    _output_topic_pub.publish(*outputCloud_msg_ptr);
}

int main( int argc, char** argv ) 
{
    MeanMapEntropyNode node(argc, argv);

	return 0;
}

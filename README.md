# Pointcloud Evaluation Tool
-----

Computes the *Mean Map Entropy* and the *Mean Plane Variance* of a point cloud.

Furthermore a pointcloud is generated that encodes the *Map Entropy* and *Plane Variance* in each point.

Filtering by field type was also included. ROS-mode and ROS-independent mode support was also added.

## Credits 
This code is based on the work by David Droeschel & Jan Razlaw:
```
Razlaw, J., Droeschel, D., Holz, D., & Behnke, S. "Evaluation of registration methods for sparse 3D laser scans." 
In 2015 IEEE European Conference on Mobile Robots (ECMR), September 2015.
```
Here are the links to their [Github](https://github.com/AIS-Bonn/pointcloud_evaluation_tool.git) and their [Paper](http://www.ais.uni-bonn.de/papers/ECMR_2015_Razlaw.pdf)

## Install
-------
### ROS mode
The following will also produce the ROS-independent build so you can run the program directly.
```
mkdir <ROS-WORKSPACE>/mean_map_entropy/src -p
cd <ROS-WORKSPACE>/mean_map_entropy/src
git clone <THIS-CODE-REPOSITORY>
cd ../
catkin build  # <- you can also use: colcon build --cmake-args -DCMAKE_BUILD_TYPE=Release
```

### ROS independent mode
```
mkdir <SOME-FOLDER>/mean_map_entropy -p
cd <SOME-FOLDER>/mean_map_entropy
git clone <THIS-CODE-REPOSITORY>
mkdir build
cd build
cmake ..
make
```

## Usage
-----
### ROS mode
```
cd <ROS-WORKSPACE>/mean_map_entropy
source install/setup.bash
roslaunch mean_map_entropy mean_map_entropy.launch input_topic:=<INPUT_CLOUD_TOPIC> output_topic:=<OUTPUT_CLOUD_TOPIC> output_frame:=<FRAME_NAME> 
```
You can use other arguments as defined in the launch file. And you can use RVIZ to visualize the output pointcloud topic.
 
### ROS independent mode
```
./mean_map_entropy_pcl path/to/pointcloud.pcd [-stepsize int] [-radius double] [-punishSolitaryPoints] [-planevariance] [-ascii] [-minNeighbors int] [-passfilter] [-filterField string] [-minLimit double] [-maxLimit double]
```
Available arguments are as follows:

* **stepsize:** stepsize used iterating through the pointcloud (default: 1)

* **radius:** radius used for search for neighboring points (default: 0.3)

* **punishSolitaryPoints:** (flag) punishes points with bad values where the number of neighbors is under minNeighbors 

* **planevariance:** (flag) whether to compute the plane Variance or not 

* **ascii:** (flag) whether to save the output point cloud file as ASCII PCD type

* **minNeighbors:** threshold (default: 15)

* **passfilter:** perform band pass filter on the values of the specified field

* **filterField:** the field name for filtering

* **minLimit:** minimum limit for filtering

* **minLimit:** maximum limit for filtering

Use the point cloud viewer of your choice to visualize the output (e.g. pcl_viewer).

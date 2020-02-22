# Pointcloud Evaluation Tool

This code based on the work by David Droeschel & Jan Razlaw:
https://github.com/AIS-Bonn/pointcloud_evaluation_tool.git

Computes the *Mean Map Entropy* and the *Mean Plane Variance* of a point cloud. 
Furthermore a pointcloud is generated that encodes the *Map Entropy* and *Plane Variance* in each point. 
Filtering by field type was also included. 

Use the point cloud viewer of your choice to visualize these measures (e.g. pcl_viewer). 

Install
-------
```
mkdir build 
cd build
cmake .. 
make
```
Usage
-----
```
./mean_map_entropy path/to/pointcloud.pcd [-stepsize int] [-radius double] [-punishSolitaryPoints] [-minNeighbors int] [-passfilter] [-filterField string] [-minLimit double] [-maxLimit double]
```
* **stepsize:** stepsize used iterating through the pointcloud (default: 1)
* **radius:** radius used for search for neighboring points (default: 0.3)
* **punishSolitaryPoints:** punishes points with bad values where the number of neighbors is under minNeighbors (default: false)
* **minNeighbors:** threshold (default: 15)
* **passfilter:** perform band pass filter on the values of the specified field
* **filterField:** the field name for filtering
* **minLimit:** minimum limit for filtering
* **minLimit:** maximum limit for filtering

---



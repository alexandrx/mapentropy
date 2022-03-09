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

#include <pcl/point_types.h>
#include <pcl/PCLPointCloud2.h>

#include <map_entropy/entropy_point_types.h>

namespace map_entropy {

class MapEntropy 
{
public:
    typedef struct ptcld_grid {
        size_t grid_id;
        size_t grid_id_x;
        size_t grid_id_y;
        int lower_bound_x;
        int lower_bound_y;
        int upper_bound_x;
        int upper_bound_y;
        std::string filename;
        pcl::PCLPointCloud2::Ptr inputCloud;
    } ptcld_grid;

    typedef struct ptcld_entropy_grid {
        size_t grid_id;
        size_t grid_id_x;
        size_t grid_id_y;
        int lower_bound_x;
        int lower_bound_y;
        int upper_bound_x;
        int upper_bound_y;
        int lonelyPointsSum;
        double entropySum ;
        double planeVarianceSum;
        double meanEntropy ;
        double meanPlaneVariance;
        std::string filename;
        pcl::PCLPointCloud2::Ptr inputCloud;
        pcl::PCLPointCloud2::Ptr outputCloud;
        pcl::PointCloud< pcl::PointXYZ >::Ptr localCloud;
        pcl::PointCloud< PointXYZWithEntropy >::Ptr cloud_entropy;
    } ptcld_entropy_grid;


    template <typename T>
    struct pointcloud_grid {
        pointcloud_grid() {
            clear();
        }   
        T& at(size_t index) { 
            return grid[index];
        }
        T& at(size_t row, size_t col) {
            size_t id = divisions_x * row + col;
            return grid[id];
        }
        T& atCoords(float x, float y) {
            int idx = static_cast<int>(floor((x - static_cast<float>(min_x_b)) / cell_size));
            int idy = static_cast<int>(floor((y - static_cast<float>(min_y_b)) / cell_size));
            int id = idy * divisions_x + idx;
            return at(id);
        }
        T& operator[] (size_t index) { 
            return at(index); 
        }
        std::pair<size_t, size_t> getRowCol(size_t index) {
            return std::pair<size_t, size_t>(grid[index].grid_id_y, grid[index].grid_id_x);
        }
        size_t getIndex(size_t row, size_t col) {
            return at(row, col).grid_id;
        }
        void setBounds(float minX, float minY, float maxX, float maxY) {
            if (cell_size > 0.f) {
                min_x_b = cell_size * static_cast<int>(floor(minX / cell_size));
                max_x_b = cell_size * static_cast<int>(floor(maxX / cell_size) + 1);
                min_y_b = cell_size * static_cast<int>(floor(minY / cell_size));
                max_y_b = cell_size * static_cast<int>(floor(maxY / cell_size) + 1);
                divisions_x = (max_x_b - min_x_b) / cell_size;
                divisions_y = (max_y_b - min_y_b) / cell_size;
                grid_num = divisions_x * divisions_y;
                grid.resize(grid_num);
            } else {
                min_x_b = 0;
                max_x_b = 0;
                min_y_b = 0;
                max_y_b = 0;
            }
        }
        void clear() {
            grid_num = 0;
            divisions_x = 0.f;
            divisions_y = 0.f;
            min_x_b = 0;
            max_x_b = 0;
            min_y_b = 0;
            max_y_b = 0;
            cell_size = 0;
            grid.clear();
        }
        void resize(size_t size) {
            grid.resize(size);
        }
        size_t size() {
            return grid.size();
        }
        float cellSize() {
            return cell_size;
        }
        void cellSize(float cellsize) {
            cell_size = cellsize;
        }
        size_t divX() {
            return divisions_x;
        }
        void divX(size_t divX) {
            divisions_x = divX;
        }
        size_t divY() {
            return divisions_y;
        }
        void divY(size_t divY) {
            divisions_y = divY;
        }
        float minX() {
            return min_x_b;
        }
        float minY() {
            return min_y_b;
        }
        float maxX() {
            return max_x_b;
        }
        float maxY() {
            return max_y_b;
        }
        size_t cellCount() {
            return grid_num;
        }
        //virtual void pointcloud2grid(const pcl::PCLPointCloud2::Ptr& cloud) = 0;

        size_t grid_num;
        size_t divisions_x;
        size_t divisions_y;
        int min_x_b;
        int max_x_b;
        int min_y_b;
        int max_y_b;
        float cell_size;
        std::vector< T > grid;
    };

    typedef struct pointcloud_grid<ptcld_entropy_grid> pointcloud_entropy_grid;

public:
    MapEntropy();

    double computeEntropy( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud );
    double computeEntropy( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud, const std::vector<int>& indices );

    double computePlaneVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud );
    double computePlaneVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud, const std::vector<int>& indices );

    void computeEntroyAndVariance( const pcl::PointCloud< pcl::PointXYZ >::Ptr& cloud );
    void computeEntroyAndVarianceFromGrid();

    void pointCloudToGrid();

    void filterPointCloud();
    void computeEntropyOrFilter();
    void computeEntropyOrFilter(const pcl::PCLPointCloud2::Ptr& inputCloud);

    int readPointCloud();
    int writePointCloud();

    //setters
    void setStepSize(int stepSize) { stepSize_ = stepSize; }
    void setRadius(double radius) { radius_ = radius; }
    void seMinNeighbors(int minNeighbors) { minNeighbors_ = minNeighbors; }
    void setMinLimit(double minLimit) { minLimit_ = minLimit; }
    void setMaxLimit(double maxLimit) { maxLimit_ = maxLimit; }
    void setFilterField(std::string filterField) { filterField_ = filterField; }
    void setPunishSolitary(bool punishSolitaryPoints) { punishSolitaryPoints_ = punishSolitaryPoints; }
    void setWithPlaneVariance(bool withPlaneVariance) { withPlaneVariance_ = withPlaneVariance; }
    void setSaveAscii(bool saveASCII) { saveASCII_ = saveASCII; }
    void setPassFilter(bool passFilter) { passFilter_ = passFilter; }
    void setFileNames(std::vector< std::string > fileNames) { fileNames_ = fileNames; }
    void setFileIndices(std::vector<int> fileIndices) { fileIndices_ = fileIndices; }
    void setGridDivide(bool grid_divide) { grid_divide_ = grid_divide; }
    void setDivideOnly(bool divide_only) { divide_only_ = divide_only; }
    void setGridSize(float grid_size) { grid_size_ = grid_size; }
    void setOutputDir(std::string output_dir) { output_dir_ = output_dir; }
    void setNamePrefix(std::string name_prefix) { name_prefix_ = name_prefix; }

    pcl::PCLPointCloud2::Ptr getOutputPtCld() { return outputCloud_; }
    pcl::PCLPointCloud2::Ptr getInputPtCld() { return inputCloud_; }

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
    pcl::PointCloud< PointXYZWithEntropy >::Ptr cloud_entropy_;
    std::vector<int> fileIndices_;
    std::vector< std::string > fileNames_;
    bool grid_divide_;
    bool divide_only_;
    float grid_size_;
    std::string output_dir_;
    std::string name_prefix_;
    pointcloud_entropy_grid ptcldGrid_;
};

}  // namespace

#endif // MAP_ENTROPY_H_

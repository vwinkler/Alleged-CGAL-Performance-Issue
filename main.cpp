#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <tbb/task_scheduler_init.h>
#define CGAL_LINKED_WITH_TBB
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<uint64_t, K> Vb;
typedef CGAL::Triangulation_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Parallel_tag> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> CT;
typedef CT::Lock_data_structure LockDS;

auto triangulate(const std::vector<std::pair<CT::Point, uint64_t>>& points,
                 const CGAL::Bbox_3& bounds) {
    const uint gridOccupancy = 100;
    
    LockDS lockingDS(bounds,
                     std::min(std::max(1u, gridOccupancy),
                              (uint) std::floor(std::cbrt(std::numeric_limits<int>::max()))));
    CT t(&lockingDS);
    t.insert(points.begin(), points.end());
    return t;
}

struct Input
{
    CGAL::Bbox_3 bounds;
    std::vector<std::pair<CT::Point, uint64_t>> points;
};

Input loadFromStream(std::istream& stream){
    std::vector<std::pair<CT::Point, uint64_t>> points;
    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double zmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymax = std::numeric_limits<double>::min();
    double zmax = std::numeric_limits<double>::min();

    std::string line;
    
    uint64_t count = 0;
    double x, y, z;
        
    while(stream >> x >> y >> z)
    {
        xmin = std::min(xmin, x);
        ymin = std::min(ymin, y);
        zmin = std::min(zmin, z);
        xmax = std::max(xmax, x);
        ymax = std::max(ymax, y);
        zmax = std::max(zmax, z);

        points.emplace_back(std::make_pair(CT::Point(x, y, z), count++));
    }

    return { {xmin, ymin, zmin, xmax, ymax, zmax}, std::move(points) };
}


int main(int /*argc*/, char** /*args*/) {
    size_t threads = 16;

    auto input = loadFromStream(std::cin);
    
    tbb::task_scheduler_init init(threads);
    
    auto start = std::chrono::steady_clock::now();
    auto dt = triangulate(input.points, input.bounds);
    auto end = std::chrono::steady_clock::now();
    

    auto numSimplices = dt.number_of_cells();
    
    printf("%lu cells : %lu ms\n",
           numSimplices,
           std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
}

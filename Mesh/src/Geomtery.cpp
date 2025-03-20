#include <gmsh.h>

#include "Geometry.h"
#include "Constants.h"

namespace Geometry
{
    int createDrop(double x, double y,
                   double z, double r,
                   int tag)
    {
        int drop_id = gmsh::model::occ::addSphere(x, y, z, r, -1, 0, Constants::PI_2);
        return drop_id;
    }
    int createElectrode(double x, double y, double z,
                        double length, double r, int tag, double angle)
    {
        int cylinder_id = gmsh::model::occ::addCylinder(x, y, z, 0, -length, 0, r, -1, Constants::PI);
        gmsh::model::occ::rotate({{3, cylinder_id}}, x, y, z, 0, -1, 0, Constants::PI_2);
        int sphere_id = gmsh::model::occ::addSphere(x, y, z, r, -1, 0, Constants::PI / 2, Constants::PI);
        gmsh::vectorpair electreodeDimTag;
        std::vector<gmsh::vectorpair> outDimTagsMap;
        gmsh::model::occ::fuse({{3, cylinder_id}}, {{3, sphere_id}}, electreodeDimTag, outDimTagsMap);
        int electrode_id = electreodeDimTag[0].second;
        gmsh::model::occ::rotate(electreodeDimTag, x, y, z, 0, 0, -1, angle);
        return electrode_id;
    }
    int createInsulator(double x, double y, double z, double h, double r1, double r2)
    {
        int bottom_ellipse_id = gmsh::model::occ::addEllipse(x, y, z, r1, r2);
        int bottom_id = gmsh::model::occ::addCurveLoop({bottom_ellipse_id});
        int top_ellipse_id = gmsh::model::occ::addEllipse(x, y, z + h, r1, r2);
        int top_id = gmsh::model::occ::addCurveLoop({top_ellipse_id});
        std::vector<std::pair<int, int>> insulator;
        gmsh::model::occ::addThruSections({bottom_id, top_id}, insulator);
        return insulator[0].second;
    }
    int createSurroundingSpace(double x, double y,
                               double z, double r,
                               double r1, int tag)
    {
        int sphere_id = gmsh::model::occ::addSphere(x, y, z, r, -1, 0, Constants::PI / 2);
        gmsh::model::occ::dilate({{3, sphere_id}}, x, y, z, 1, r1 / r, 1);
        return sphere_id;
    }
}
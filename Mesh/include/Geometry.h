#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace Geometry
{
    int createDrop(double x, double y,
                   double z, double r,
                   int tag = -1);
    int createElectrode(double x, double y, double z,
                        double length, double r, int tag = -1, double angle = 0);
    int createInsulator(double x, double y, double z, double h, double r1, double r2);
    int createSurroundingSpace(double x, double y,
                               double z, double r,
                               double r1, int tag = -1);
}

#endif
#include <fstream>
#include <iostream>
#include <set>
#include <gmsh.h>

#include "include/Constants.h"
#include "include/Geometry.h"
#include "include/Mesh.h"

int main(int argc, char *argv[])
{
    Mesh::init("Drop");
    int drop_id = Geometry::createDrop(Constants::X_DROP, Constants::Y_DROP, Constants::Z_DROP, Constants::R_DROP);

    int source_id = Geometry::createElectrode(Constants::X_DROP, Constants::Y_DROP - 3 * Constants::R_DROP, Constants::Z_DROP,
                                              50, Constants::R_DROP, -1);
    int drain_id = Geometry::createElectrode(Constants::X_DROP, Constants::Y_DROP + 3 * Constants::R_DROP, Constants::Z_DROP,
                                             50, Constants::R_DROP, -1, Constants::PI);
    int gate_id = Geometry::createElectrode(Constants::X_DROP + 3 * Constants::R_DROP, Constants::Y_DROP, Constants::Z_DROP,
                                            50, Constants::R_DROP, -1, -Constants::PI_2);
    int insulator_id = Geometry::createInsulator(0, 0, 0, Constants::Z_DROP, Constants::R_AIR, Constants::R_AIR_1);
    int surrounding_space_id = Geometry::createSurroundingSpace(0, 0, Constants::Z_DROP, Constants::R_AIR, Constants::R_AIR_1);

    Mesh::createMesh("Drop.msh");

    Mesh::addPhysicalGroup(3, {insulator_id}, -1, "insulator");
    Mesh::addPhysicalGroup(3, {surrounding_space_id}, -1, "Surrounding space");

    Mesh::addPhysicalGroup3DSurface(drop_id, "Drop");
    Mesh::addPhysicalGroup3DSurface(source_id, "Source");
    Mesh::addPhysicalGroup3DSurface(drain_id, "Drain");
    Mesh::addPhysicalGroup3DSurface(gate_id, "Gate");

    Mesh::visual(argc, argv);
    Mesh::finalize();
    return 0;
}
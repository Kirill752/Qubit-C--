#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>

namespace Mesh {
    void init(const std::string &name);
    void createMesh(const std::string &filename);
    void visual(int argc, char *argv[]);
    void finalize();
    int addPhysicalGroup(int dim, const std::vector<int> &tags, int tag = -1, const std::string &name = "");
    int addPhysicalGroup3DSurface(int tag, const std::string &name);
}

#endif
#include <string>
#include <set>
#include <gmsh.h>

namespace Mesh
{
    void init(const std::string &name)
    {
        gmsh::initialize();
        gmsh::model::add("Drop");
    }
    void createMesh(const std::string &filename)
    {
        gmsh::model::occ::removeAllDuplicates();
        gmsh::model::occ::synchronize();
        gmsh::model::mesh::generate(3);
        gmsh::write("Drop.msh");
    }
    void visual(int argc, char *argv[])
    {
        std::set<std::string> args(argv, argv + argc);
        if (!args.count("-nopopup"))
            gmsh::fltk::run();
    }
    void finalize()
    {
        gmsh::finalize();
    }
    int addPhysicalGroup(int dim, const std::vector<int> &tags, int tag, const std::string &name)
    {
        return gmsh::model::addPhysicalGroup(dim, tags, tag, name);
    }
    int addPhysicalGroup3DSurface(int tag, const std::string &name) {
        gmsh::vectorpair out;
        gmsh::model::getBoundary({{3, tag}}, out);
        std::vector<int> surface_ids;
        for (int i = 0; i < out.size(); i++)
        {
            surface_ids.push_back(out[i].second);
        }
        return gmsh::model::addPhysicalGroup(2, surface_ids, -1, name);
    }
}

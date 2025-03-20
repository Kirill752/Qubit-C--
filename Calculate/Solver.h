#include <string>
#include <mfem.hpp>

namespace Solve
{
    void PoisonEquation(std::string mesh_file,
                        int argc, char *argv[],
                        int order = 1,
                        bool static_cond = false,
                        bool pa = false,
                        bool fa = false,
                        const char *device_config = "cpu",
                        bool visualization = true,
                        bool algebraic_ceed = false);
    double Capacity(int num_attr, int order, int dim, mfem::Mesh &mesh, mfem::GridFunction &ugrad);
}

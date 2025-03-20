#include "Solver.h"
#include <eigen3/Eigen/Dense>

using namespace mfem;
using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace Solve
{
    void PoisonEquation(std::string mesh_file,
                        int argc, char *argv[],
                        int order,
                        bool static_cond,
                        bool pa,
                        bool fa,
                        const char *device_config,
                        bool visualization,
                        bool algebraic_ceed)
    {
        bool delete_fec;
        ConstantCoefficient zero(0.0);
        ConstantCoefficient one(1.0);
        ConstantCoefficient antione(-1.0);
        OptionsParser args(argc, argv);
        args.AddOption(&mesh_file, "-m", "--mesh",
                       "Mesh file to use.");
        args.AddOption(&order, "-o", "--order",
                       "Finite element order (polynomial degree) or -1 for"
                       " isoparametric space.");
        args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                       "--no-static-condensation", "Enable static condensation.");
        args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                       "--no-partial-assembly", "Enable Partial Assembly.");
        args.AddOption(&fa, "-fa", "--full-assembly", "-no-fa",
                       "--no-full-assembly", "Enable Full Assembly.");
        args.AddOption(&device_config, "-d", "--device",
                       "Device configuration string, see Device::Configure().");
        args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                       "--no-visualization",
                       "Enable or disable GLVis visualization.");
        args.Parse();
        if (!args.Good())
        {
            args.PrintUsage(std::cout);
            return;
        }
        args.PrintOptions(std::cout);

        Device device(device_config);
        device.Print();

        Mesh mesh(mesh_file, 1, 1);
        int dim = mesh.Dimension();

        {
            int ref_levels =
                (int)floor(log(5000. / mesh.GetNE()) / log(2.) / dim);
            for (int l = 0; l < ref_levels; l++)
            {
                mesh.UniformRefinement();
            }
        }

        {
            int ref_levels =
                (int)floor(log(50000. / mesh.GetNE()) / log(2.) / dim);
            for (int l = 0; l < ref_levels; l++)
            {
                mesh.UniformRefinement();
            }
        }

        // 1 - первая капля
        // 2 - вторая капля
        // 3, 4 - электроды
        // 5, 6 - затворы
        // cap[i][j]: i - строка; j - столбец
        double cap[6][6] = {};
        FiniteElementCollection *fec;
        if (order > 0)
        {
            fec = new H1_FECollection(order, dim);
            delete_fec = true;
        }
        else if (mesh.GetNodes())
        {
            fec = mesh.GetNodes()->OwnFEC();
            delete_fec = false;
        }
        else
        {
            fec = new H1_FECollection(order = 1, dim);
            delete_fec = true;
        }
        FiniteElementCollection *fecgrad = new ND_FECollection(order, dim);
        FiniteElementSpace fespace(&mesh, fec);
        FiniteElementSpace fespacegrad(&mesh, fecgrad);
        Array<int> ess_tdof_list;
        if (mesh.bdr_attributes.Size())
        {
            Array<int> ess_bdr(mesh.bdr_attributes.Max()); // определяем массив, длина которого равна количеству физических поверхностей
            ess_bdr = 1;                                   // определяем все физические поверхности как "существенные"
            fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
        }

        LinearForm b(&fespace);
        b.AddDomainIntegrator(new DomainLFIntegrator(zero));
        b.Assemble();
        GridFunction x(&fespace);
        BilinearForm a(&fespace);
        // Задаём диэлектрическую проницаемость среды.
        // В воздузе eps = 1
        // В изоляторе eps = 9.4
        Vector eps(mesh.attributes.Max());
        eps = 1.0;
        eps(9) = eps(10) * 9.4;
        PWConstCoefficient eps_func(eps);
        a.AddDomainIntegrator(new DiffusionIntegrator(eps_func));
        a.Assemble();
        a.Finalize();
        DiscreteLinearOperator grad(&fespace, &fespacegrad);
        grad.AddDomainInterpolator(new GradientInterpolator);
        grad.Assemble();
        grad.Finalize();
        GridFunction ugrad(&fespacegrad);
        Array<int> ess_bdr(mesh.bdr_attributes.Max());
        Coefficient *coeff[1];

        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                if (i != j)
                {
                    ess_bdr = 0; 
                    ess_bdr[j] = 1;
                    coeff[0] = &antione;
                }
                if (i == j)
                {
                    ess_bdr = 0; 
                    ess_bdr[0] = 1;
                    ess_bdr[1] = 1;
                    ess_bdr[2] = 1;
                    ess_bdr[3] = 1;
                    ess_bdr[4] = 1;
                    ess_bdr[5] = 1;
                    coeff[0] = &one;
                }
                x = 0; 
                x.ProjectBdrCoefficient(coeff, ess_bdr);
                OperatorPtr A;
                Vector B, X;
                a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

                if (!pa)
                {
#ifndef MFEM_USE_SUITESPARSE
                    GSSmoother M((SparseMatrix &)(*A));
                    PCG(*A, M, B, X, 1, 200, 1e-12, 0.0);
#else
                    UMFPackSolver umf_solver;
                    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
                    umf_solver.SetOperator(*A);
                    umf_solver.Mult(B, X);
#endif
                }
                else
                {
                    if (UsesTensorBasis(fespace))
                    {
                        if (algebraic_ceed)
                        {
                            ceed::AlgebraicSolver M(a, ess_tdof_list);
                            PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
                        }
                        else
                        {
                            OperatorJacobiSmoother M(a, ess_tdof_list);
                            PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
                        }
                    }
                    else
                    {
                        CG(*A, B, X, 1, 400, 1e-12, 0.0);
                    }
                }
                a.RecoverFEMSolution(X, b, x);
                grad.Mult(x, ugrad);
                ugrad *= -1.0;
                cap[i][j] = Capacity(i, order, dim, mesh, ugrad);
            }
        }

        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                cap[i][j] = (cap[i][j] + cap[j][i] / 2);
                cap[j][i] = cap[i][j];
            }
        }

        double beta[6][6];
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                if (i == j)
                {
                    beta[i][j] = cap[0][j] + cap[1][j] + cap[2][j] + cap[3][j] + cap[4][j] + cap[5][j];
                }
                else if (i != j)
                {
                    beta[i][j] = -cap[i][j];
                }
            }
        }

        MatrixXd cap_eigen(6, 6);
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                cap_eigen(i, j) = cap[i][j];
            }
        }

        std::cout << "Матрица ёмкостных коэффициентов: " << std::endl;
        std::cout << cap_eigen << std::endl;

        MatrixXd beta_eigen(6, 6);
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                beta_eigen(i, j) = beta[i][j];
            }
        }

        std::cout << "Матрица электростатической индукции: " << std::endl;
        std::cout << beta_eigen << std::endl;
        MatrixXd beta_dd(2, 2);
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                beta_dd(i, j) = beta[i][j];
            }
        }

        std::cout << "Матрица электростатической индукции для двух капель: " << std::endl;
        std::cout << beta_dd << std::endl;

        VectorXd phi_eigen(2);
        phi_eigen << 1, 1;

        // Вектор зарядов на каплях
        VectorXd q_dd(2);
        // Создаём цикл для вычисления зарядовой энегрии
        // n1 - количество электронов на первой капле
        // n2 - количество электронов на второй капле
        std::ofstream fout("charge_energy.txt");
        for (int n1 = -10; n1 < 11; ++n1)
        {
            for (int n2 = -10; n2 < 11; ++n2)
            {
                q_dd << n1, n2;
                fout << n1 << "   " << n2 << "   " << 0.5 * q_dd.transpose() * beta_dd.inverse() * q_dd << std::endl;
            }
        }
        fout.close();
        {
            std::ostringstream mesh_name, sol_name;
            mesh_name << "mesh." << std::setfill('0');
            sol_name << "sol." << std::setfill('0');

            std::ofstream mesh_ofs(mesh_name.str().c_str());
            mesh_ofs.precision(8);
            mesh.Print(mesh_ofs);

            std::ofstream sol_ofs(sol_name.str().c_str());
            sol_ofs.precision(8);
            x.Save(sol_ofs);
        }
        {
            std::ostringstream gsol_name;
            gsol_name << "gsol." << std::setfill('0');

            std::ofstream gsol_ofs(gsol_name.str().c_str());
            gsol_ofs.precision(8);
            ugrad.Save(gsol_ofs);
        }

        // // 16. Отпрвляем решение на GLVIS сервер.
        // if (visualization)
        // {
        //     char vishost[] = "localhost";
        //     int visport = 19916;
        //     socketstream sol_sock(vishost, visport);
        //     sol_sock.precision(8);
        //     sol_sock << "solution\n"
        //              << mesh << x << flush;
        // }
        // if (visualization)
        // {
        //     char vishost[] = "localhost";
        //     int visport = 19916;
        //     socketstream sol_sock(vishost, visport);
        //     sol_sock.precision(8);
        //     sol_sock << "grad\n"
        //              << mesh << ugrad << flush;
        // }

        delete fec;
        delete fecgrad;
    }
    double Capacity(int num_attr, int order, int dim, Mesh &mesh, GridFunction &ugrad)
    {
        // Рассчет емкости.
        // Создаём отделное пространство конечных элементов для расчета ёмкости.
        FiniteElementCollection *fec_capacity = new H1_FECollection(order, dim);
        FiniteElementSpace fespace_capacity(&mesh, fec_capacity);
        // Задаем линейную форму. Тут мы рассчитываем заряд по теореме Гаусса. Но так как потенциал на капле равен еденице,
        // то сразу получаем ёмкость.
        LinearForm capacity(&fespace_capacity);
        // Передаём значения вектора напряженности поля из сеточной функции в векторный кожффициент,
        // который будет играть роль вектора напряженности в интеграле.
        VectorGridFunctionCoefficient E_gridfunc(&ugrad);
        // Задаём существенные поверхности. Это по верхности, по которым будет вестись интегрирование.
        //  1 -- существенная поверхность
        //  0 -- несущественная поверхность
        Array<int> bdr_attr(mesh.bdr_attributes.Max());
        bdr_attr = 0;
        bdr_attr[num_attr] = 1;
        // Добавляем интегратор по поверхности вида ((Е, n), v), где v - базисная функция.
        //  bdr_atr - указатель на объект, по поверхности которого будет идти интеграрование.
        capacity.AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(E_gridfunc), bdr_attr);
        capacity.Assemble();
        // Задаем множитель перед интегралом
        GridFunction ones(&fespace_capacity);
        ones = -1. / (4 * 3.1415);
        // Возвращаем значение ёмкости.
        delete fec_capacity;
        return capacity(ones);
    }
}
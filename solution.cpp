#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

void solver(int i, int j, mfem::ParGridFunction& x, mfem::ParMesh& pmesh,  mfem::ParLinearForm& b, mfem::ParBilinearForm& a, mfem::Array<int> ess_tdof_list){

   // 1. Задаём начальные условия.
    {
      ConstantCoefficient zero(0.0);
      ConstantCoefficient one(1.0);
      ConstantCoefficient antione(-1.0);
      Array<int> ess_bdr(pmesh.bdr_attributes.Max());
      Coefficient* coeff[1]; // задаем массив значений, которые хотим присвоить как граничное условие Дирихле
      ess_bdr = 0; // делаем "несущественными" все физические поверхности
      //поиск собственных ёмкостей
      if (i == j){
      ess_bdr[0] = 1; 
      ess_bdr[1] = 1; 
      ess_bdr[2] = 1;
      ess_bdr[3] = 1;
      ess_bdr[4] = 1;
      ess_bdr[5] = 1;
      coeff[0]=&one;
      }
      else {
         ess_bdr[j] = 1;
         coeff[0] = &antione;
      }
      x.ProjectBdrCoefficient(coeff, ess_bdr); 
   } 

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   // 2.Решаем линейную систему A X = B.
   Solver *prec = NULL;
   prec = new HypreBoomerAMG;

   CGSolver cg(MPI_COMM_WORLD);
   cg.SetRelTol(1e-12);
   cg.SetMaxIter(2000);
   cg.SetPrintLevel(1);
   if (prec) { cg.SetPreconditioner(*prec); }
   cg.SetOperator(*A);
   cg.Mult(B, X);
   delete prec;

   // 3. Восстанавливем решение.
   a.RecoverFEMSolution(X, b, x);
};


double Capacity(int num_attr, int order, int dim, ParMesh& pmesh, ParGridFunction& ugrad){
  //Рассчет емкости.
  //Создаём отделное пространство конечных элементов для расчета ёмкости.
  FiniteElementCollection *fec_capacity = new H1_FECollection(order, dim); 
  FiniteElementSpace fespace_capacity(&pmesh, fec_capacity);
  //Задаем линейную форму. Тут мы рассчитываем заряд по теореме Гаусса. Но так как потенциал на капле равен еденице, 
  //то сразу получаем ёмкость.
  LinearForm capacity(&fespace_capacity);
  // Передаём значения вектора напряженности поля из сеточной функции в векторный кожффициент,
  //который будет играть роль вектора напряженности в интеграле.
  VectorGridFunctionCoefficient E_gridfunc(&ugrad);
  //Задаём существенные поверхности. Это по верхности, по которым будет вестись интегрирование.
  // 1 -- существенная поверхность
  // 0 -- несущественная поверхность
  Array<int> bdr_attr(pmesh.bdr_attributes.Max());
  bdr_attr = 0;
  bdr_attr[num_attr] = 1;
  //Добавляем интегратор по поверхности вида ((Е, n), v), где v - базисная функция.
  // bdr_atr - указатель на объект, по поверхности которого будет идти интеграрование.
  capacity.AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(E_gridfunc), bdr_attr);
  capacity.Assemble();
  // Задаем множитель перед интегралом
  GridFunction ones(&fespace_capacity);
  ones = -1./(4*3.1415);
  // Возвращаем значение ёмкости.
  return capacity(ones); 
}; 


int main(int argc, char *argv[])
{
   // 1. Инициализируем MPI и HYPRE.
   Mpi::Init();
   int num_procs = Mpi::WorldSize();
   int myid = Mpi::WorldRank();
   Hypre::Init();

   // 2. Задаем параметры запуска.
   const char *mesh_file = "./Drop.msh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   bool fa = false;
   const char *device_config = "cpu";
   bool visualization = false;
   bool algebraic_ceed = false;

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
#ifdef MFEM_USE_CEED
   args.AddOption(&algebraic_ceed, "-a", "--algebraic",
                  "-no-a", "--no-algebraic",
                  "Use algebraic Ceed solver");
#endif
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      if (myid == 0)
      {
         args.PrintUsage(cout);
      }
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(cout);
   }

   // 3. Включаем аппаратные устройства, такие как графические процессоры, и модели программирования, такие как
   // CUDA, OCCA, RAJA и OpenMP, на основе параметров командной строки.
   Device device(device_config);
   if (myid == 0) { device.Print(); }

   // 4. Считаем (последовательную) сетку из данного файла сетки. Мы
   // можем обрабатывать треугольные, четырехугольные, тетраэдрические, шестигранные, поверхностные
   // и объемные сетки с одним и тем же кодом.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 5. Refine the serial mesh on all processors to increase the resolution. In
   //    this example we do 'ref_levels' of uniform refinement. We choose
   //    'ref_levels' to be the largest number that gives a final mesh with no
   //    more than 10,000 elements.
   {
      int ref_levels =
         (int)floor(log(5000./mesh.GetNE())/log(2.)/dim);
      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
   }

   // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
   //    this mesh further in parallel to increase the resolution. Once the
   //    parallel mesh is defined, the serial mesh can be deleted.
   ParMesh pmesh(MPI_COMM_WORLD, mesh);
   mesh.Clear();
   {
      int par_ref_levels = 1;
      for (int l = 0; l < par_ref_levels; l++)
      {
         pmesh.UniformRefinement();
      }
   }

   // 7. Define a parallel finite element space on the parallel mesh. Here we
   //    use continuous Lagrange finite elements of the specified order. If
   //    order < 1, we instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   FiniteElementCollection *fecgrad = new ND_FECollection(order, dim); // создаём простравнство конечных элементов для градиента потенциала
   bool delete_fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
      delete_fec = true;
   }
   else if (pmesh.GetNodes())
   {
      fec = pmesh.GetNodes()->OwnFEC();
      delete_fec = false;

      if (myid == 0)
      {
         cout << "Using isoparametric FEs: " << fec->Name() << endl;
      }
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
      delete_fec = true;
   }
   ParFiniteElementSpace fespace(&pmesh, fec);
   ParFiniteElementSpace *fespacegrad = new ParFiniteElementSpace(&pmesh, fecgrad);
   HYPRE_BigInt size = fespace.GlobalTrueVSize();
   if (myid == 0)
   {
      cout << "Number of finite element unknowns: " << size << endl;
   }

   // 8. Determine the list of true (i.e. parallel conforming) essential
   //    boundary dofs. In this example, the boundary conditions are defined
   //    by marking all the boundary attributes from the mesh as essential
   //    (Dirichlet) and converting them to a list of true dofs.
   Array<int> ess_tdof_list;
   if (pmesh.bdr_attributes.Size())
   {
      Array<int> ess_bdr(pmesh.bdr_attributes.Max()); // определяем массив, длина которого равна количеству физических поверхностей
      ess_bdr = 1; // определяем все физические поверхности как "существенные"
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 9. Задаем параллельную линейную форму b(.), 
   //которая соответствует правой части линейной системы 
   // метода конечных элементов Ax=b.
   // В данном случае она равна (0, phi_i), 
   // где phi_i -- базисная функция.  
   ParLinearForm b(&fespace);
   ConstantCoefficient zero(0.0);
   ConstantCoefficient one(1.0);
   ConstantCoefficient antione(-1.0);
   b.AddDomainIntegrator(new DomainLFIntegrator(zero));
   b.Assemble();

   // 10. Задаём сеточную функцию и начальные условия.
   ParGridFunction x(&fespace); // потенциал поля
   /* ParGridFunction x_grad(fespacegrad); // напряженность поля
   //x = 0.0;
    {
      Array<int> ess_bdr(pmesh.bdr_attributes.Max());
      Coefficient* coeff[1]; // задаем массив значений, которые хотим присвоить как граничное условие Дирихле
      ess_bdr = 0; // делаем "несущественными" все физические поверхности
      ess_bdr[3] = 1; // делаем "существенной" поверхность первой капли
      ess_bdr[4] = 0; // делаем "существенной" поверхность второй капли
      coeff[0]=&one;
      x.ProjectBdrCoefficient(coeff, ess_bdr); // Задаем значение "1" на поверхности первой капли

      // Аналогично задаём значение потенциала на второй капле
      ess_bdr = 0;
      ess_bdr[3] = 0;
      ess_bdr[4] = 1;
      coeff[0] = &one;
      x.ProjectBdrCoefficient(coeff, ess_bdr);

      //Задаём потенциал на электродах 
      //Первый электрод
      ess_bdr = 0;
      ess_bdr[7] = 1;
      ess_bdr[8] = 0;
      coeff[0]=&one;
      x.ProjectBdrCoefficient(coeff, ess_bdr);
   
      // Второй электрод
      ess_bdr = 0;
      ess_bdr[7] = 0;
      ess_bdr[8] = 1;
      coeff[0]=&one;
      x.ProjectBdrCoefficient(coeff, ess_bdr);

      //Первый затвор
      ess_bdr = 0;
      ess_bdr[9] = 1;
      ess_bdr[10] = 0;
      coeff[0]=&one;
      x.ProjectBdrCoefficient(coeff, ess_bdr);
   
      // Второй затвор
      ess_bdr = 0;
      ess_bdr[9] = 0;
      ess_bdr[10] = 1;
      coeff[0]=&one;
      x.ProjectBdrCoefficient(coeff, ess_bdr);
   } 
 */
   // 11. Создаем параллельную билинейную форму a(.,.) в пространстве конечных элементов,
   // соответствующую оператору Лапласа, добавив DiffusionIntegrator
   ParBilinearForm a(&fespace);
   if (pa) { 
      a.SetAssemblyLevel(AssemblyLevel::PARTIAL); 
   }
   if (fa)
   {
      a.SetAssemblyLevel(AssemblyLevel::FULL);
      // Sort the matrix column indices when running on GPU or with OpenMP (i.e.
      // when Device::IsEnabled() returns true). This makes the results
      // bit-for-bit deterministic at the cost of somewhat longer run time.
      a.EnableSparseMatrixSorting(Device::IsEnabled());
   }
   // Задаём диэлектрическую проницаемость среды.
   // В воздузе eps = 1
   // В изоляторе eps = 9.4
   Vector eps(pmesh.attributes.Max());
   eps = 1.0;
   eps(9) = eps(10)*9.4;
   PWConstCoefficient eps_func(eps);
   a.AddDomainIntegrator(new DiffusionIntegrator(eps_func));

   // 12. Assemble the parallel bilinear form and the corresponding linear
   //     system, applying any necessary transformations such as: parallel
   //     assembly, eliminating boundary conditions, applying conforming
   //     constraints for non-conforming AMR, static condensation, etc.

   if (static_cond) { a.EnableStaticCondensation(); }
   a.Assemble();
   if (!pa) { a.Finalize(); }

   /* OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
   // 13. Solve the linear system A X = B.
   //     * With full assembly, use the BoomerAMG preconditioner from hypre.
   //     * With partial assembly, use Jacobi smoothing, for now.
   Solver *prec = NULL;
   if (pa)
   {
      if (UsesTensorBasis(fespace))
      {
         if (algebraic_ceed)
         {
            prec = new ceed::AlgebraicSolver(a, ess_tdof_list);
         }
         else
         {
            prec = new OperatorJacobiSmoother(a, ess_tdof_list);
         }
      }
   }
   else
   {
      prec = new HypreBoomerAMG;
   }
   CGSolver cg(MPI_COMM_WORLD);
   cg.SetRelTol(1e-12);
   cg.SetMaxIter(2000);
   cg.SetPrintLevel(1);
   if (prec) { cg.SetPreconditioner(*prec); }
   cg.SetOperator(*A);
   cg.Mult(B, X);
   delete prec;

   // 14. Recover the parallel grid function corresponding to X. This is the
   //     local finite element solution on each processor.
   a.RecoverFEMSolution(X, b, x); */

  solver(0, 0, x, pmesh, b, a, ess_tdof_list);

  // Расчет напряженности поля.

  //Задём дискретный линейный оператор градиента. 
  ParDiscreteLinearOperator grad(&fespace, fespacegrad);
  grad.AddDomainInterpolator(new GradientInterpolator);
  //grad->SetAssemblyLevel(AssemblyLevel::PARTIAL);
  grad.Assemble();
  grad.Finalize();
  //Задаём сеточную функцию, в которую запишем значения напряжености электрического поля
  ParGridFunction ugrad(fespacegrad);
  grad.Mult(x, ugrad);
  ugrad *= -1.0;


  //Рассчет емкости.
/* 
  //Создаём отделное пространство конечных элементов для расчета ёмкости.
  FiniteElementCollection *fec_capacity = new H1_FECollection(order, dim); 
  FiniteElementSpace fespace_capacity(&pmesh, fec_capacity);
  //Задаем линейную форму. Тут мы рассчитываем заряд по теореме Гаусса. Но так как потенциал на капле равен еденице, 
  //то сразу получаем ёмкость.
  LinearForm capacity(&fespace_capacity);
  // Передаём значения вектора напряженности поля из сеточной функции в векторный кожффициент,
  //который будет играть роль вектора напряженности в интеграле.
  VectorGridFunctionCoefficient E_gridfunc(&ugrad);
  //Задаём существенные поверхности. Это по верхности, по которым будет вестись интегрирование.
  // 1 -- существенная поверхность
  // 0 -- несущественная поверхность
  Array<int> bdr_attr(pmesh.bdr_attributes.Max());
  bdr_attr = 0;
  // bdr_atr[3] - это поверхность первой капли
  bdr_attr[3] = 0; 
  // bdr_atr[4] - это поверхность второй капли
  bdr_attr[4] = 0;
  //bdr_atr[7] - это поверхность первого электрода
  bdr_attr[7] = 0;
  //bdr_atr[8] - это поверхность второго электрода
  bdr_attr[8] = 1;
  //Добавляем интегратор по поверхности вида ((Е, n), v), где v - базисная функция.
  // bdr_atr - указатель на объект, по поверхности которого будет идти интеграрование.
  capacity.AddBoundaryIntegrator(new BoundaryNormalLFIntegrator(E_gridfunc), bdr_attr);
  capacity.Assemble();
  // Задаем множитель перед интегралом
  GridFunction ones(&fespace_capacity);
  ones = -1./(4*3.1415);
  // Выводим значение ёмкости в консоль.
  cout << "Ёмкость равна = " << capacity(ones) << endl;  */


  // 1 - первая капля
  // 2 - вторая капля
  // 3, 4 - электроды
  // 5, 6 - затворы
  // cap[i][j]: i - строка; j - столбец
  double cap[6][6] = {};
  for (int i = 0; i < 6; i++){
   for(int j = 0; j < 6; j++){
      if (j == i){
         cap[i][j] = Capacity(i, order, dim, pmesh, ugrad);
      }
   }
  }

  solver(1, 0, x, pmesh, b, a, ess_tdof_list);
  grad.Mult(x, ugrad);
  ugrad *= -1.0;
  cap[1][0] = Capacity(1, order, dim, pmesh, ugrad);

  // Вывод матрицы емкости
  cout << "Матрица емкости: " << endl << endl;  
   for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 6; ++j)
        {
            cout.precision(3);
            cout << cap[i][j] << ' ';
        }
        cout << endl;
    }


   // 15. Save the refined mesh and the solution in parallel. This output can
   //     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
   {
      ostringstream mesh_name, sol_name;
      mesh_name << "mesh." << setfill('0') << setw(6) << myid;
      sol_name << "sol." << setfill('0') << setw(6) << myid;

      ofstream mesh_ofs(mesh_name.str().c_str());
      mesh_ofs.precision(8);
      pmesh.Print(mesh_ofs);

      ofstream sol_ofs(sol_name.str().c_str());
      sol_ofs.precision(8);
      x.Save(sol_ofs);
   }
  {
      ostringstream gsol_name;
      gsol_name << "gsol." << setfill('0') << setw(6) << myid;

      ofstream gsol_ofs(gsol_name.str().c_str());
      gsol_ofs.precision(8);
      ugrad.Save(gsol_ofs);
   }

   // 16. Отпрвляем решение на GLVIS сервер.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "solution\n" << pmesh << x << flush;
   }
     if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock << "parallel " << num_procs << " " << myid << "\n";
      sol_sock.precision(8);
      sol_sock << "grad\n" << pmesh << ugrad << flush;
   } 

   // 17. Free the used memory.
   if (delete_fec)
   {
      delete fec;
      delete fecgrad;
   }

   return 0;
}
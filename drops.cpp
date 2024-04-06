#define M_PI 3.14159265358979323846 /* pi */
// #include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include "gmsh.h"
using namespace std;
// using namespace mfem;

int main(int argc, char *argv[])
{

    // Параметры окружающей среды
    double R_air = 200.0;   // Радиус сферы окружающей стреды
    double R_air_1 = 150.0; // Нужен для эллипса

    // Параметры изолятора
    double h_insulator = 200.0; // толщина изолятора

    // Параметры капель
    double r_drop = 15.0; // Радиус капли
    double drop_dist = 1;
    double x_drop = 0.5 * drop_dist + r_drop; // координата Х капли
    double y_drop = 0;                        // координата Y капли
    double z_drop = h_insulator;              // координата Z капли

    // Параметры электродов
    double r_electrode = r_drop * 5 / 4; // радиус электрода
    double l_electrode = 50.0;           // длина электрода
    double el_dist = 1;
    double r_gate = 2 * r_drop;
    double gate_dist = 10;

    gmsh::initialize(argc, argv);
    gmsh::model::add("Drop");

    // добавляем нижний слой
    gmsh::model::occ::addEllipse(0, 0, 0, R_air, R_air_1, 1);
    gmsh::model::occ::addCurveLoop({1}, 1);
    gmsh::model::occ::synchronize();
    // добавляем слой изолятора
    gmsh::model::occ::addEllipse(0, 0, h_insulator, R_air, R_air_1, 2);
    gmsh::model::occ::addCurveLoop({2}, 2);
    std::vector<std::pair<int, int>> out;
    gmsh::model::occ::addThruSections({2, 1}, out, 1, 1, 1);
    gmsh::model::occ::synchronize();

    // добавляем капли на изоляторе
    gmsh::model::occ::addSphere(-x_drop, y_drop, z_drop, r_drop, 2, 0, M_PI / 2);
    gmsh::model::occ::addSphere(x_drop, y_drop, z_drop, r_drop, 3, 0, M_PI / 2);
    gmsh::model::occ::synchronize();

    // добавляем сетку окружения
    gmsh::model::occ::addSphere(0, 0, h_insulator, R_air, 4, 0, M_PI / 2);
    gmsh::model::occ::dilate({{3, 4}}, 0, 0, h_insulator, 1, R_air_1/R_air, 1);
    gmsh::model::occ::synchronize();
    // Убираем объемы капель
    std::vector<std::pair<int, int>> ov;
    std::vector<std::vector<std::pair<int, int>>> ovv;
    gmsh::model::occ::cut({{3, 4}}, {{3, 2}, {3, 3}}, ov, ovv, -1, true, true);
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize();

    // добавляем первый электрод
    gmsh::model::occ::addCylinder(-x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, -l_electrode, 0, 0, r_electrode, -1, M_PI);
    gmsh::model::occ::rotate({{3, 5}}, -x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, 1, 0, 0, 3 * M_PI / 2);
    gmsh::model::occ::addSphere(-x_drop - l_electrode - r_electrode - el_dist - r_drop, 0, h_insulator, r_electrode, -1, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 6}}, -x_drop - l_electrode - r_electrode - el_dist - r_drop, 0, h_insulator, 0, 0, 1, M_PI / 2);
    gmsh::model::occ::addSphere(-x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, r_electrode, -1, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 7}}, -x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, 0, 0, 1, -M_PI / 2);
    gmsh::model::occ::fuse({{3, 5}}, {{3, 6}, {3, 7}}, ov, ovv);
    gmsh::model::occ::synchronize();

    // добавляем второй электрод
    gmsh::model::occ::copy({{3, 5}}, ov);
    gmsh::model::occ::mirror(ov, 1, 0, 0, 0);
    gmsh::model::occ::synchronize();
    
    // добавляем первый затвор
    gmsh::model::occ::addCylinder(0, 0 - r_gate - gate_dist - r_drop, h_insulator, 0, -l_electrode, 0, r_gate, -1, M_PI);
    gmsh::model::occ::rotate({{3, 7}}, 0, 0 - r_gate - gate_dist - r_drop, h_insulator, 0, 1, 0, 3 * M_PI / 2);
    gmsh::model::occ::addSphere(0, 0 - l_electrode - r_gate - gate_dist - r_drop, h_insulator, r_gate, -1, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 8}}, 0, 0 - l_electrode - r_gate - gate_dist - r_drop, h_insulator, 1, 0, 0, M_PI / 2);
    gmsh::model::occ::addSphere(0, 0 - r_gate - gate_dist - r_drop, h_insulator, r_gate, -1, 0, M_PI / 2, M_PI);
    gmsh::model::occ::fuse({{3, 7}}, {{3, 8}, {3, 9}}, ov, ovv);
    gmsh::model::occ::synchronize();

    gmsh::model::occ::copy({{3, 7}}, ov);
    gmsh::model::occ::mirror(ov, 0, 1, 0, 0);
    gmsh::model::occ::translate(ov, x_drop, 0, 0);
    gmsh::model::occ::translate({{3, 7}}, -x_drop, 0, 0);
    gmsh::model::occ::synchronize();

    gmsh::model::occ::cut({{3, 4}}, {{3, 5}, {3, 6}, {3, 7}, {3, 8}}, ov, ovv);
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize(); 

     // добавляем физические группы
    gmsh::model::addPhysicalGroup(2, {33}, 1, "Surface of insulator");
    gmsh::model::addPhysicalGroup(2, {40}, 2, "Bottom layer");
    gmsh::model::addPhysicalGroup(3, {1}, 3, "Insulator");
    gmsh::model::addPhysicalGroup(2, {31, 38}, 4, "Sphere surface_1"); 
    gmsh::model::addPhysicalGroup(2, {32, 39}, 5, "Sphere surface_2");
    gmsh::model::addPhysicalGroup(3, {4}, 6, "Surrounding space");
    gmsh::model::addPhysicalGroup(2, {29}, 7, "Surrounding space surface");
    gmsh::model::addPhysicalGroup(2, {9, 10, 12, 35}, 8, "Electrode surface_1");
    gmsh::model::addPhysicalGroup(2, {13, 14, 16, 36}, 9, "Electrode surface_2");
    gmsh::model::addPhysicalGroup(2, {25, 26, 28, 34}, 10, "Gate surface_1");
    gmsh::model::addPhysicalGroup(2, {21, 22, 24, 37}, 11, "Gate surface_1");

    // создаем поле для сетки
    gmsh::model::mesh::field::add("Distance", 1);
    // создаём поле, опираясь на то, вокруг каких объектов нужно сгущать сетку
    gmsh::model::mesh::field::setNumbers(1, "SurfacesList", {31, 38, 32, 39}); // задаём поверхности "капель"
    gmsh::model::mesh::field::setNumber(1, "Sampling", 1);                     // регулировка густоты сетки
                                                                               // создаём новое поле, которое будет менять размер сетки, опираясь на значение, возвращаемое "Полем 1"
    gmsh::model::mesh::field::add("Threshold", 2);
    gmsh::model::mesh::field::setNumber(2, "InField", 1);   // взяли значения параметров сетки из "Поля 1"
    gmsh::model::mesh::field::setNumber(2, "SizeMin", 0.1*r_drop); // минимаьный размер
    gmsh::model::mesh::field::setNumber(2, "SizeMax", 2*r_drop);  // максимальный размер
    gmsh::model::mesh::field::setNumber(2, "DistMin", 0.15);
    gmsh::model::mesh::field::setNumber(2, "DistMax", 0.5);
    gmsh::model::mesh::field::setAsBackgroundMesh(2);
    gmsh::model::mesh::generate(3);
    // cout << gmsh::model::occ::getMaxTag(3) << endl;
    
    gmsh::write("Drop.msh");
    // запуск gmsh
    std::set<std::string> args(argv, argv + argc);
    if (!args.count("-nopopup"))
        gmsh::fltk::run();

    gmsh::finalize();

    return 0;
}
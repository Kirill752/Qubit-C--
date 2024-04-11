#define M_PI 3.14159265358979323846 /* pi */
// #include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <set>
#include "gmsh.h"
using namespace std;

int main(int argc, char *argv[])
{

    // Параметры окружающей среды
    double R_air = 200.0;   // Радиус сферы окружающей стреды
    double R_air_1 = 190.0; // Нужен для эллипса

    // Параметры изолятора
    double h_insulator = 200.0; // толщина изолятора

    // Параметры капель
    double r_drop = 15.0; // Радиус капли
    double drop_dist = 1;
    double x_drop = 0; // координата Х капли
    double y_drop = 0;                        // координата Y капли
    double z_drop = h_insulator;              // координата Z капли

    // Параметры электродов
    double r_electrode = r_drop * 5 / 4; // радиус электрода
    double l_electrode_gate = 120.0;           // длина затвора
    double l_electrode = 150;
    double el_dist = 1;
    double r_gate = 2 * r_drop;
    double gate_dist = 10;

    gmsh::initialize(argc, argv);
    gmsh::model::add("Drop");

    std::vector<std::pair<int, int>> ov;
    std::vector<std::vector<std::pair<int, int>>> ovv;
    int tag_1 = gmsh::model::occ::getMaxTag(1);
    int tag_2 = gmsh::model::occ::getMaxTag(2);
    int tag_3 = gmsh::model::occ::getMaxTag(3);

    // добавляем капли на изоляторе
    gmsh::model::occ::addSphere(-x_drop, y_drop, z_drop, r_drop, 1, 0, M_PI/2);
    gmsh::model::occ::synchronize();

     // добавляем первый электрод
    gmsh::model::occ::addCylinder(x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, -l_electrode, 0, 0, r_electrode, 2, M_PI);
    gmsh::model::occ::rotate({{3, 2}}, x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, 1, 0, 0, 3 * M_PI / 2);
/*     gmsh::model::occ::addSphere(-x_drop - l_electrode - r_electrode - el_dist - r_drop, 0, h_insulator, r_electrode, 3, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 3}}, x_drop - l_electrode - r_electrode - el_dist - r_drop, 0, h_insulator, 0, 0, 1, M_PI / 2); */
    gmsh::model::occ::addSphere(x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, r_electrode, 4, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 4}}, x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, 0, 0, 1, -M_PI / 2);
    gmsh::model::occ::fuse({{3, 2}}, { {3, 4}}, ov, ovv);

    gmsh::model::occ::synchronize();

    // добавляем второй электрод
    gmsh::model::occ::copy({{3, 2}}, ov);
    gmsh::model::occ::mirror(ov, 1, 0, 0, 0);
    gmsh::model::occ::synchronize();
    
    // добавляем первый затвор
    gmsh::model::occ::addCylinder(0, 0 - r_gate - gate_dist - r_drop, h_insulator, 0, -l_electrode_gate, 0, r_gate, 4, M_PI);
    gmsh::model::occ::rotate({{3, 4}}, 0, 0 - r_gate - gate_dist - r_drop, h_insulator, 0, 1, 0, 3 * M_PI / 2);
/*     gmsh::model::occ::addSphere(0, 0 - l_electrode - r_gate - gate_dist - r_drop, h_insulator, r_gate, 5, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3,5}}, 0, 0 - l_electrode - r_gate - gate_dist - r_drop, h_insulator, 1, 0, 0, M_PI/ 2);*/
    gmsh::model::occ::addSphere(0, 0 - r_gate - gate_dist - r_drop, h_insulator, r_gate, 6, 0, M_PI / 2, M_PI);
    gmsh::model::occ::fuse({{3, 4}}, {{3, 6}}, ov, ovv);
    gmsh::model::occ::translate({{3, 4}}, -x_drop, 0, 0); 
    gmsh::model::occ::synchronize();

     // добавляем слой изолятора
    tag_1 = gmsh::model::occ::getMaxTag(1);
    tag_2 = gmsh::model::occ::getMaxTag(2);
    gmsh::model::occ::addEllipse(0, 0, 0, R_air, R_air_1, tag_1+1);
    gmsh::model::occ::addCurveLoop({tag_1+1}, tag_1+2);
    gmsh::model::occ::synchronize();

    gmsh::model::occ::addEllipse(0, 0, h_insulator, R_air, R_air_1, tag_1+3);
    gmsh::model::occ::addCurveLoop({tag_1+3}, tag_1+4);
    std::vector<std::pair<int, int>> out;
    gmsh::model::occ::addThruSections({tag_1+4, tag_1+2}, out, -1, 1, 1); 
    gmsh::model::occ::synchronize();   

    tag_1 = gmsh::model::occ::getMaxTag(1);
    tag_2 = gmsh::model::occ::getMaxTag(2);
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize(); 

     tag_3 = gmsh::model::occ::getMaxTag(3);
    gmsh::model::occ::addSphere(0, 0, h_insulator, R_air, tag_3+1, 0, M_PI / 2);
    gmsh::model::occ::dilate({{3, tag_3+1}}, 0, 0, h_insulator, 1, R_air_1/R_air, 1);
    gmsh::model::occ::synchronize(); 
    
    gmsh::model::occ::cut({{3, 6}}, {{3, 1}, {3,2}, {3, 3}, {3, 4}}, ov, ovv, -1, true, true);
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize(); 

     // добавляем физические группы     
    gmsh::model::addPhysicalGroup(2, {15, 3}, 1, "Sphere surface_1"); 
    gmsh::model::setColor({{2, 3}, {2, 15}}, 0, 255, 0);
    gmsh::model::addPhysicalGroup(2, {4, 12, 13, 14}, 3, "Electrode surface_1");
    gmsh::model::setColor({{2, 4}, {2, 12}, {2, 13}, {2,14}}, 135, 206, 235);
    gmsh::model::addPhysicalGroup(2, {5, 16, 17, 18}, 4, "Electrode surface_2");
    gmsh::model::setColor({{2, 5}, {2, 16}, {2, 17}, {2,18}}, 135, 206, 235);
    gmsh::model::addPhysicalGroup(2, {6, 9, 10, 11}, 5, "Gate surface_1");
    gmsh::model::setColor({{2, 6}, {2, 9}, {2, 10}, {2,11}}, 135, 206, 250);
    gmsh::model::addPhysicalGroup(2, {1}, 7, "Surface of insulator");
    gmsh::model::setColor({{2, 1}, {2, 2}}, 127, 127, 127);
    gmsh::model::addPhysicalGroup(2, {7}, 8, "Bottom layer");
    gmsh::model::addPhysicalGroup(3, {5}, 9, "Insulator");
    gmsh::model::addPhysicalGroup(3, {6}, 10, "Surrounding space");
    gmsh::model::setVisibility({{3, 6}}, 0);
    gmsh::model::addPhysicalGroup(2, {8}, 11, "Surrounding space surface");   
    gmsh::model::setVisibility({{2, 8}}, 0);
    gmsh::model::setVisibility({{1, 18}}, 0); 

    // создаем поле для сетки
    gmsh::model::mesh::field::add("Distance", 1);
    // создаём поле, опираясь на то, вокруг каких объектов нужно сгущать сетку
    gmsh::model::mesh::field::setNumbers(1, "SurfacesList", {3, 15}); // задаём поверхности "капель"
    gmsh::model::mesh::field::setNumber(1, "Sampling", 1);                     // регулировка густоты сетки
                                                                               // создаём новое поле, которое будет менять размер сетки, опираясь на значение, возвращаемое "Полем 1"
    gmsh::model::mesh::field::add("Threshold", 2);
    gmsh::model::mesh::field::setNumber(2, "InField", 1);   // взяли значения параметров сетки из "Поля 1"
    gmsh::model::mesh::field::setNumber(2, "SizeMin", 0.1*r_drop); // минимаьный размер
    gmsh::model::mesh::field::setNumber(2, "SizeMax", 0.8*r_drop);  // максимальный размер
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
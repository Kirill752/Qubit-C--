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
    double R_air_1 = 190.0; // Нужен для эллипса

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
    double l_electrode = 140.0;           // длина электрода
    double l_electrode_gate = 120.0;
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
    gmsh::model::occ::addSphere(x_drop, y_drop, z_drop, r_drop, 2, 0, M_PI/2);
    gmsh::model::occ::synchronize();
/* 
    // добавляем сетку окружения
    gmsh::model::occ::addSphere(0, 0, h_insulator, R_air, 4, 0, M_PI / 2);
    gmsh::model::occ::dilate({{3, 4}}, 0, 0, h_insulator, 1, R_air_1/R_air, 1);
    gmsh::model::occ::synchronize(); 
    // Убираем объемы капель
    gmsh::model::occ::cut({{3, 4}}, {{3, 2}, {3, 3}}, ov, ovv, -1, true, true);
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize(); */

     // добавляем первый электрод
    gmsh::model::occ::addCylinder(-x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, -l_electrode, 0, 0, r_electrode, 3, M_PI);
    gmsh::model::occ::rotate({{3, 3}}, -x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, 1, 0, 0, 3 * M_PI / 2);
/*     gmsh::model::occ::addSphere(-x_drop - l_electrode - r_electrode - el_dist - r_drop, 0, h_insulator, r_electrode, 4, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 4}}, -x_drop - l_electrode - r_electrode - el_dist - r_drop, 0, h_insulator, 0, 0, 1, M_PI / 2); */
    gmsh::model::occ::addSphere(-x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, r_electrode, 5, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 5}}, -x_drop - r_electrode - el_dist - r_drop, 0, h_insulator, 0, 0, 1, -M_PI / 2);
    gmsh::model::occ::fuse({{3, 3}}, { {3, 5}}, ov, ovv);

    gmsh::model::occ::synchronize();

    // добавляем второй электрод
    gmsh::model::occ::copy({{3, 3}}, ov);
    gmsh::model::occ::mirror(ov, 1, 0, 0, 0);
    gmsh::model::occ::synchronize();
    
    // добавляем первый затвор
    gmsh::model::occ::addCylinder(0, 0 - r_gate - gate_dist - r_drop, h_insulator, 0, -l_electrode_gate, 0, r_gate, 5, M_PI);
    gmsh::model::occ::rotate({{3, 5}}, 0, 0 - r_gate - gate_dist - r_drop, h_insulator, 0, 1, 0, 3 * M_PI / 2);
/*     gmsh::model::occ::addSphere(0, 0 - l_electrode - r_gate - gate_dist - r_drop, h_insulator, r_gate, 6, 0, M_PI / 2, M_PI);
    gmsh::model::occ::rotate({{3, 6}}, 0, 0 - l_electrode - r_gate - gate_dist - r_drop, h_insulator, 1, 0, 0, M_PI / 2); */
    gmsh::model::occ::addSphere(0, 0 - r_gate - gate_dist - r_drop, h_insulator, r_gate, 7, 0, M_PI / 2, M_PI);
    gmsh::model::occ::fuse({{3, 5}}, { {3, 7}}, ov, ovv);
    gmsh::model::occ::translate({{3, 5}}, -x_drop, 0, 0); 
    gmsh::model::occ::synchronize();

    gmsh::model::occ::copy({{3, 5}}, ov);
    gmsh::model::occ::mirror(ov, 0, 1, 0, 0);
    gmsh::model::occ::synchronize();
    gmsh::model::occ::translate(ov, 2*x_drop, 0, 0);
    gmsh::model::occ::synchronize();

    /*gmsh::model::occ::cut({{3, 4}}, {{3, 5}, {3, 6}, {3, 7}, {3, 8}}, ov, ovv);
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize(); */ 
    // добавляем нижний слой
    /* tag_1 = gmsh::model::occ::getMaxTag(1);
    tag_2 = gmsh::model::occ::getMaxTag(2);
    gmsh::model::occ::addEllipse(0, 0, 0, R_air, R_air_1, tag_1+1);
    tag_1 = gmsh::model::occ::getMaxTag(1);
    gmsh::model::occ::addCurveLoop({tag_1}, tag_1+1);
    gmsh::model::occ::addPlaneSurface({tag_1+1}, tag_2+1);
    gmsh::model::occ::synchronize();

    tag_1 = gmsh::model::occ::getMaxTag(1);
    tag_2 = gmsh::model::occ::getMaxTag(2);
    gmsh::model::occ::addEllipse(0, 0, h_insulator, R_air, R_air_1, tag_1+1);
    tag_1 = gmsh::model::occ::getMaxTag(1);
    gmsh::model::occ::addCurveLoop({tag_1}, tag_1+1);
    gmsh::model::occ::addCurveLoop({3}, tag_1+2);
    gmsh::model::occ::addCurveLoop({6}, tag_1+3);
    gmsh::model::occ::addCurveLoop({9, 15, 12 ,14}, tag_1+4);
    gmsh::model::occ::addCurveLoop({25, 19, 24, 22}, tag_1+5); 
    gmsh::model::occ::addCurveLoop({28, 35, 31, 34, 33}, tag_1+6);
    gmsh::model::occ::addCurveLoop({45, 38, 43, 44, 41}, tag_1+7); 
    gmsh::model::occ::addPlaneSurface({tag_1+1, tag_1+2, tag_1+3, tag_1+4, tag_1+5, tag_1+6, tag_1+7}, tag_2+1);
    gmsh::model::mesh::createGeometry();
     //gmsh::model::occ::addCylinder(300, 15 , 15, 0, -l_electrode, 0, r_gate);
    gmsh::model::occ::synchronize();*/

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
    
    gmsh::model::occ::cut({{3, 8}}, {{3, 1}, {3,2}, {3, 3}, {3, 4}, {3,5}, {3, 6}}, ov, ovv, -1, true, true);
    gmsh::model::occ::removeAllDuplicates();
    gmsh::model::occ::synchronize(); 
     // добавляем физические группы     
    gmsh::model::addPhysicalGroup(2, {17, 3}, 1, "Sphere surface_1"); 
    gmsh::model::setColor({{2, 3}, {2, 17}}, 0, 255, 0); // зеленный
    gmsh::model::addPhysicalGroup(2, {18, 4}, 2, "Sphere surface_2");
    gmsh::model::setColor({{2, 4}, {2, 18}}, 0, 255, 0); // зеленный
    gmsh::model::addPhysicalGroup(2, {5, 16, 15, 14}, 3, "Electrode surface_1");
    gmsh::model::setColor({{2, 5}, {2, 16}, {2, 14}, {2,15}}, 135, 206, 235); // голубой
    gmsh::model::addPhysicalGroup(2, {19, 20, 21, 6}, 4, "Electrode surface_2");
    gmsh::model::setColor({{2, 6}, {2, 19}, {2, 20}, {2, 21}}, 135, 206, 235); //голубой
    gmsh::model::addPhysicalGroup(2, {13, 12, 11, 7}, 5, "Gate surface_1");
    gmsh::model::setColor({{2, 7}, {2, 11}, {2, 12}, {2, 13}}, 135, 206, 250); // синий
    gmsh::model::addPhysicalGroup(2, {22, 23, 24, 8}, 6, "Gate surface_2");
    gmsh::model::setColor({{2, 8}, {2, 22}, {2, 23}, {2, 24}}, 135, 206, 250); // синий
    gmsh::model::addPhysicalGroup(2, {1}, 7, "Surface of insulator");
    gmsh::model::setColor({{2, 1}, {2, 2}}, 127, 127, 127); // серый
    gmsh::model::addPhysicalGroup(2, {9}, 8, "Bottom layer");
    gmsh::model::addPhysicalGroup(3, {7}, 9, "Insulator");
    gmsh::model::addPhysicalGroup(3, {8}, 10, "Surrounding space");
    gmsh::model::setVisibility({{3, 8}}, 0);
    gmsh::model::addPhysicalGroup(2, {10}, 11, "Surrounding space surface");  
    gmsh::model::setVisibility({{2, 10}}, 0);
    gmsh::model::setVisibility({{1, 23}}, 0);
    
    // создаем поле для сетки
    gmsh::model::mesh::field::add("Distance", 1);
    // создаём поле, опираясь на то, вокруг каких объектов нужно сгущать сетку
    gmsh::model::mesh::field::setNumbers(1, "SurfacesList", {17, 3, 18, 4}); // задаём поверхности "капель" 
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
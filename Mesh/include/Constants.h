#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants
{
    constexpr double PI = 3.14159265358979323846;
    constexpr double PI_2 = 1.57079632679489661923;
    // Параметры окружающей среды
    constexpr double R_AIR = 210.0;   // Радиус сферы окружающей среды
    constexpr double R_AIR_1 = 200.0; // Радиус для эллипса

    // Параметры изолятора
    constexpr double H_INSULATOR = 200.0; // Толщина изолятора

    // Параметры капель
    constexpr double R_DROP = 15.0;                     // Радиус капли
    constexpr double DROP_DIST = 1.0;                   // Расстояние между каплями
    constexpr double X_DROP = 0.5 * DROP_DIST + R_DROP; // Координата X капли
    constexpr double Y_DROP = 0.0;                      // Координата Y капли
    constexpr double Z_DROP = H_INSULATOR;              // Координата Z капли

    // Параметры электродов
    constexpr double R_ELECTRODE = R_DROP * 4 / 3; // Радиус электрода
    constexpr double L_ELECTRODE = 140.0;          // Длина электрода
    constexpr double L_ELECTRODE_GATE = 120.0;     // Длина затвора
    constexpr double EL_DIST = 1.0;                // Расстояние между электродами
    constexpr double R_GATE = 2 * R_DROP;          // Радиус затвора
    constexpr double GATE_DIST = 10.0;             // Расстояние до затвора
}

#endif
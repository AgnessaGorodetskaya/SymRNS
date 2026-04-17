#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "SymRNSFixed.cpp"

int main(void) {
    constexpr Positional_Int S = 10; // коэффициент масштабирования (не кратен ни одному из оснований)
    SymRnsBase base{{3, 7, 11, 13}, S}; // набор оснований ССОК
    // SymRnsBase base{{3, 7, 11, 13, 17}, S};

    Positional_Float prec = 0.5 / S;  // целевая точность 1/2S

    // тесты
    bool ok = true;

    // целочисленные проверки конвертации в позиционную ИС
    for (Positional_Int x_pos = -base.P / 2; x_pos <= base.P / 2; ++x_pos) {
        SymRnsFixed x_rns{x_pos, base};
        Positional_Float x_pos_float = static_cast<Positional_Float>(x_pos) / S;
        Positional_Float x_pos_crt = x_rns.to_positional_crt();
        Positional_Float x_pos_mrc = x_rns.to_positional_mrc();
        Positional_Float x_pos_frac = x_rns.to_positional_frac_crt();
        if ((std::abs(x_pos_crt - x_pos_float) > prec) ||
            (std::abs(x_pos_crt - x_pos_mrc) > prec) ||
            (std::abs(x_pos_crt - x_pos_frac) > prec))
        {
            std::cout << "ОШИБКА: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_crt << " MRC=" << x_pos_mrc << " MRC_F=" << x_pos_frac << std::endl;
            ok = false;
        } else {
            // std::cout << "ИНФО: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_crt << " MRC=" << x_pos_mrc << " MRC_F=" << x_pos_frac << std::endl;
        }
    }

    // // проверки fixed-point обратной конвертации в позиционную ИС
    // for (Positional_Int x_int = 0; x_int < base.P; ++x_int) {
    //     Positional_Float x = static_cast<Positional_Float>(x_int) / q;
    //     RnsNumber x_rns = RnsNumber{x, base};
    //     Positional_Int x_pos_int;
    //     Positional_Float x_pos = x_rns.to_positional(&x_pos_int);
    //     if (fabs(x - x_pos) > 0.001) {
    //         std::cout << "Error: original " << x << " != converted " << x_pos << std::endl;
    //         ok = false;
    //     } else {
    //         std::cout << "int: " << x_int << " calc: " << x << " " << x_rns << " conv: " << x_pos
    //                 << " (" << x_pos_int << ")" << " rank=" << x_rns.get_rank() << std::endl;
    //     }
    // }

    // for (Positional_Int x_int = -base.P/2; x_int <= base.P/2; ++x_int) {

    //     // проверки унарного минуса
    //     Positional_Float res = -static_cast<Positional_Float>(x_int) / q;
    //     Positional_Int res_int = -x_int;
    //     if (res_int >= -base.P/2 && res_int <= base.P/2) {
    //         RnsNumber r_rns = -RnsNumber(x_int, base, q);
    //         RnsNumber res_rns = RnsNumber(res_int, base, q);
    //         if (res_rns != r_rns) {
    //             std::cout << "ОШИБКА: x: " << static_cast<Positional_Float>(x_int) / q << "
    //             convert: "
    //                 << res_rns.to_positional() << " calc: " << r_rns.to_positional() <<
    //                 std::endl;
    //             ok = false;
    //         }
    //     }

    //     for (Positional_Int y_int = -base.P/2; y_int <= base.P/2; ++y_int) {

    //         // проверки сложения
    //         res = static_cast<Positional_Float>(x_int + y_int) / q;
    //         res_int = x_int + y_int;
    //         if (res_int >= -base.P/2 && res_int <= base.P/2) {
    //             RnsNumber r_rns = RnsNumber(x_int, base, q) + RnsNumber(y_int, base,
    //             q); RnsNumber res_rns = RnsNumber(res, base, q); if (res_rns != r_rns) {
    //                 std::cout << "ОШИБКА: " << static_cast<Positional_Float>(x_int) / q << " + "
    //                           << static_cast<Positional_Float>(y_int) / q << " convert: " <<
    //                           res_rns.to_positional()
    //                           << " calc: " << r_rns.to_positional() << std::endl;
    //                 ok = false;
    //             }
    //         }

    //         // проверки умножения
    //         res = static_cast<Positional_Float>(x_int) / q * static_cast<Positional_Float>(y_int)
    //         / q; res_int = x_int * y_int; if (res_int >= -base.P/2 && res_int <=
    //         base.P/2) {
    //             RnsNumber r_rns = RnsNumber(x_int, base, q) * RnsNumber(y_int, base,
    //             q); RnsNumber res_rns = RnsNumber(res, base, q*q); if (res_rns != r_rns) {
    //                 std::cout << "ОШИБКА: x: " << static_cast<Positional_Float>(x_int) / q << " *
    //                 y:"
    //                     << static_cast<Positional_Float>(y_int) / q << " convert: " <<
    //                     res_rns.to_positional()
    //                     << " calc: " << r_rns.to_positional() << std::endl;
    //                 ok = false;
    //             }
    //         }
    //     }
    // }
    if (ok) std::cout << "ТЕСТЫ ПРОйДЕНЫ!!!" << std::endl << std::endl;

    // примеры использования fixed point
    // умножение
    Positional_Int a_int = 36; // целочисленное представление числа A = 36/S
    Positional_Int b_int = -38; // целочисленное представление числа B = -38/S
    SymRnsFixed a{a_int, base}; // ССОК Fixed Point представление числа A
    SymRnsFixed b{b_int, base}; // ССОК Fixed Point представление числа B
    SymRnsFixed c = a * b; // Умножение в ССОК Fixed Point
    Positional_Float c_fl = c.to_positional_frac_crt();  // перевод результата C в позиционное Floating Point

    Positional_Float a_fl = static_cast<Positional_Float>(a_int) / S; // Число A как позиционное Floating Point
    Positional_Float b_fl = static_cast<Positional_Float>(b_int) / S; // Число B как позиционное Floating Point
    Positional_Float cc_fl = a_fl * b_fl; // Вычисление умножения в позиционном Floating Point
    Positional_Float eps = std::abs(c_fl - cc_fl); // Погрешность вычисления

    if (eps > prec) {  // погрешность больше целевой точности
        ok = false;
        std::cout << "ОШИБКА: ";
    } else {
        std::cout << "OK: ";
    }
    std::cout << c_fl << " ? " << cc_fl << " eps=" << eps <<
            " prec=" << prec << std::endl;

    std::cout << a.to_positional_crt() << '*' << b.to_positional_crt() << '=' << c.to_positional_crt() << ' ' << c << std::endl << std::endl;

    // деление
    a_int = 62; // целочисленное представление числа A = 62/S
    b_int = -29; // целочисленное представление числа A = -29/S
    a = SymRnsFixed{a_int, base}; // ССОК Fixed Point представление числа A
    b = SymRnsFixed{b_int, base}; // ССОК Fixed Point представление числа B
    c = a / b; // Деление в ССОК Fixed Point
    c_fl = c.to_positional_frac_crt(); // перевод результата C в позиционное Floating Point

    a_fl = static_cast<Positional_Float>(a_int) / S; // Число A как позиционное Floating Point
    b_fl = static_cast<Positional_Float>(b_int) / S; // Число B как позиционное Floating Point
    cc_fl = a_fl / b_fl; // Вычисление деления в позиционном Floating Point
    eps = std::abs(c_fl - cc_fl); // Погрешность вычисления

    if (eps > prec) { // погрешность больше целевой точности
        ok = false;
        std::cout << "ОШИБКА: ";
    } else {
        std::cout << "OK: ";
    }
    std::cout << c_fl << " ? " << cc_fl << " eps=" << eps <<
            " prec=" << prec << std::endl;

    std::cout << a.to_positional_crt() << '/' << b.to_positional_crt() << '=' << c.to_positional_crt() << ' ' << c << std::endl << std::endl;

    // сортировка через реализованную операцию сравнения
    constexpr size_t VSIZE = 10; // размер вектора
    std::vector<SymRnsFixed> v(VSIZE, SymRnsFixed{0, base}); // создаем вектор из ССОК Fixed Point
    for (size_t i = 0; i < VSIZE; ++i) {
        v[i] = SymRnsFixed(std::rand() % base.P - base.P / 2, base);  // случайные значения целоч.представления от -P/2 до P/2
    }
    std::sort(v.begin(), v.end()); // стандартный алгоритм сортировки через реализованный оператор сравнения <
    for (const auto& x_rns : v) { // вывод вектора на экран
        std::cout << x_rns << ' ' << x_rns.to_positional_frac_crt() << std::endl;
    }

    if (ok) std::cout << "ТЕСТЫ ПРОйДЕНЫ!!!" << std::endl << std::endl;
}

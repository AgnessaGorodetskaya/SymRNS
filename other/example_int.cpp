#include "RNS.cpp"
#include <iostream>

int main() {
    RnsBase base{{3, 5, 7}};
    // RnsBase base{{5, 7, 11, 13, 17, 19}};
    // RnsBase base{{5, 7, 11, 13, 17, 19}};

    // for (int x= -10; x <= 10; ++x) {
    //     std::cout << "x=" << x << " x^-1=" << RnsNumber::mod_inverse(x, 5) << std::endl;
    // }

    // примеры округления
    // RnsNumber a{{2, 3, 8}, base};
    // RnsNumber a{{2, 5, 8}, base};
    // std::cout << a.to_positional() << ' ' << a << " rank=" << a.get_rank() << std::endl;

    // RnsNumber b = a.round(10);
    // std::cout << b.to_positional() << ' ' << b <<  std::endl;


    // тесты
    // Positional_Int q = 1;  // делитель для fixed point тестов
    bool ok = true;

    // // проверки целочисленного округления
    // Positional_Int rnd = 10;
    // for (Positional_Int x_int = 110; x_int <= 115; ++x_int) {
    //     RnsNumber x_rns = RnsNumber{x_int, base};
    //     RnsNumber x_rns_round = x_rns.round(rnd);
    //     Positional_Int x_rns_round_int = x_rns_round.to_positional_int();
    //     Positional_Int x_round_int = x_int / rnd * rnd;
    //     if (x_round_int != x_rns_round_int) {
    //         std::cout << "Error: original " << x_int << x_rns << "->" << x_round_int <<  " != converted " << x_rns_round_int << std::endl;
    //         ok = false;
    //     } else {
    //         std::cout << "Info: original " << x_int << x_rns <<"->" << x_round_int <<  " == converted " << x_rns_round_int << std::endl;
    //     }
    // }

    // целочисленные проверки конвертации в позиционную ИС
    for (Positional_Int x_pos = 0; x_pos < base.P; ++x_pos) {
        RnsNumber x_rns = RnsNumber{x_pos, base};
        Positional_Int x_pos_crt = x_rns.to_positional_crt();
        Positional_Int x_pos_mrc = x_rns.to_positional_mrc();
        if (x_pos != x_pos_crt || x_pos != x_pos_mrc) {
            std::cout << "ОШИБКА: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_crt << " MRC=" << x_pos_mrc << std::endl;
            ok = false;
        } else {
            // std::cout << "ИНФО: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_crt << " MRC=" << x_pos_mrc << std::endl;
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

    // // примеры использования fixed point
    // RnsNumber a{5, base, 6};  // = 5/6
    // RnsNumber b{-8, base, 9};  // = -8/9

    // // RnsNumber a{0.82f, base, 6};   // ~= 5/6
    // // RnsNumber b{-0.86f, base, 9};  // ~= -8/9

    // Positional_Int a_int, b_int, c_int;
    // Positional_Float a_fl = a.to_positional(&a_int);
    // Positional_Float b_fl = b.to_positional(&b_int);

    // RnsNumber c = a + b;
    // Positional_Float c_fl = c.to_positional(&c_int);
    // std::cout << a << '=' << a_int << '/' << a.q_int << '=' << a_fl << " + " << b << '=' << b_int
    //           << '/' << b.q_int << '=' << b_fl << " = " << c << '=' << c_int << '/' << c.q_int
    //           << '=' << c_fl << std::endl;

    // c = a - b;
    // c_fl = c.to_positional(&c_int);
    // std::cout << a << '=' << a_int << '/' << a.q_int << '=' << a_fl << " - " << b << '=' << b_int
    //           << '/' << b.q_int << '=' << b_fl << " = " << c << '=' << c_int << '/' << c.q_int
    //           << '=' << c_fl << std::endl;

    // c = a * b;
    // c_fl = c.to_positional(&c_int);
    // std::cout << a << '=' << a_int << '/' << a.q_int << '=' << a_fl << " * " << b << '=' << b_int
    //           << '/' << b.q_int << '=' << b_fl << " = " << c << '=' << c_int << '/' << c.q_int
    //           << '=' << c_fl << std::endl;
}

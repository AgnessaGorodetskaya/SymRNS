#include "SymRNSFrac.cpp"
#include <iostream>

int main() {
    SymRnsBase base{{3, 7, 11}};

    // тесты
    bool ok = true;

    // примеры использования fixed point
    SymRnsFrac a{{5, base}, {6, base}};
    SymRnsFrac b{{-8, base}, {9, base}};
    SymRnsFrac c{{-8, base}, {6, base}};

    std::cout << "a=" << a << '=' << a.numerator.to_positional_crt() << '/' << a.denominator.to_positional_crt() << std::endl;
    std::cout << "b=" << b << '=' << b.numerator.to_positional_crt() << '/' << b.denominator.to_positional_crt() << std::endl;
    std::cout << "c=" << c << '=' << c.numerator.to_positional_crt() << '/' << c.denominator.to_positional_crt() << std::endl;

    SymRnsFrac z = a * b;
    std::cout << "a*b=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;

    z = a + b;
    std::cout << "a+b=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;

    z = a - b;
    std::cout << "a-b=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;

    z = a / b;
    std::cout << "a/b=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;

    z = a * c;
    std::cout << "a*c=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;

    z = a + c;
    std::cout << "a+c=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;

    z = a - c;
    std::cout << "a-c=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;

    z = a / c;
    std::cout << "a/c=" << z << '=' << z.numerator.to_positional_crt() << '/' << z.denominator.to_positional_crt() << std::endl;


    // Positional_Int a_int, b_int, c_int;
    // Positional_Float a_fl = a.to_positional(&a_int);
    // Positional_Float b_fl = b.to_positional(&b_int);

    // SymRnsNumber z = a + b;
    // Positional_Float c_fl = z.to_positional(&c_int);
    // std::cout << a << '=' << a_int << '/' << a.q_int << '=' << a_fl << " + " << b << '=' << b_int
    //           << '/' << b.q_int << '=' << b_fl << " = " << z << '=' << c_int << '/' << z.q_int
    //           << '=' << c_fl << std::endl;



    // for (Positional_Int x_int = -base.P/2; x_int <= base.P/2; ++x_int) {

    //     for (Positional_Int y_int = -base.P/2; y_int <= base.P/2; ++y_int) {

    //         // проверки сложения
    //         Positional_Int res_int = x_int + y_int;
    //         if (res_int >= -base.P/2 && res_int <= base.P/2) {
    //             SymRnsNumber r_rns = SymRnsNumber(x_int, base) + SymRnsNumber(y_int, base);
    //             SymRnsNumber res_rns = SymRnsNumber(res_int, base);
    //             if (res_rns != r_rns) {
    //                 std::cout << "ОШИБКА: " << x_int << " + "
    //                           << y_int << " convert: " <<
    //                           res_rns.to_positional_crt()
    //                           << " calc: " << r_rns.to_positional_crt() << std::endl;
    //                 ok = false;
    //             }
    //         }

    //         // проверки умножения
    //         res_int = x_int * y_int;
    //         if (res_int >= -base.P/2 && res_int <= base.P/2) {
    //             SymRnsNumber r_rns = SymRnsNumber(x_int, base) * SymRnsNumber(y_int, base);
    //             SymRnsNumber res_rns = SymRnsNumber(res_int, base);
    //             if (res_rns != r_rns) {
    //                 std::cout << "ОШИБКА: x: " << x_int << " * y:"
    //                     << y_int << " convert: " <<
    //                     res_rns.to_positional_crt()
    //                     << " calc: " << r_rns.to_positional_crt() << std::endl;
    //                 ok = false;
    //             }
    //         }
    //     }
    // }
    if (ok) std::cout << "ТЕСТЫ ПРОйДЕНЫ!!!" << std::endl << std::endl;


    // z = a - b;
    // c_fl = z.to_positional(&c_int);
    // std::cout << a << '=' << a_int << '/' << a.q_int << '=' << a_fl << " - " << b << '=' << b_int
    //           << '/' << b.q_int << '=' << b_fl << " = " << z << '=' << c_int << '/' << z.q_int
    //           << '=' << c_fl << std::endl;

    // z = a * b;
    // c_fl = z.to_positional(&c_int);
    // std::cout << a << '=' << a_int << '/' << a.q_int << '=' << a_fl << " * " << b << '=' << b_int
    //           << '/' << b.q_int << '=' << b_fl << " = " << z << '=' << c_int << '/' << z.q_int
    //           << '=' << c_fl << std::endl;

}

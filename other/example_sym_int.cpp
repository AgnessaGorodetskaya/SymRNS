#include "SymRNS.cpp"
#include <iostream>

int main() {
    SymRnsBase base{{3, 7, 11}};

    // тесты
    bool ok = true;

    // целочисленные проверки конвертации в позиционную ИС
    for (Positional_Int x_pos = -base.P / 2; x_pos <= base.P / 2; ++x_pos) {
        SymRnsNumber x_rns = SymRnsNumber{x_pos, base};
        Positional_Int x_pos_crt = x_rns.to_positional_crt();
        Positional_Int x_pos_mrc = x_rns.to_positional_mrc();
        if (x_pos != x_pos_crt || x_pos != x_pos_mrc) {
            std::cout << "ОШИБКА: x=" << x_pos << " ССОК=" << x_rns << " ORT=" << x_pos_crt << " MRC=" << x_pos_mrc << std::endl;
            ok = false;
        } else {
            // std::cout << "ИНФО: x=" << x_pos << " ССОК=" << x_rns << " ORT=" << x_pos_crt << " MRC=" << x_pos_mrc << std::endl;
        }
    }

    for (Positional_Int x_int = -base.P/2; x_int <= base.P/2; ++x_int) {

        for (Positional_Int y_int = -base.P/2; y_int <= base.P/2; ++y_int) {

            // проверки сложения
            Positional_Int res_int = x_int + y_int;
            if (res_int >= -base.P/2 && res_int <= base.P/2) {
                SymRnsNumber r_rns = SymRnsNumber(x_int, base) + SymRnsNumber(y_int, base);
                SymRnsNumber res_rns = SymRnsNumber(res_int, base);
                if (res_rns != r_rns) {
                    std::cout << "ОШИБКА: " << x_int << " + "
                              << y_int << " convert: " <<
                              res_rns.to_positional_crt()
                              << " calc: " << r_rns.to_positional_crt() << std::endl;
                    ok = false;
                }
            }

            // проверки умножения
            res_int = x_int * y_int;
            if (res_int >= -base.P/2 && res_int <= base.P/2) {
                SymRnsNumber r_rns = SymRnsNumber(x_int, base) * SymRnsNumber(y_int, base);
                SymRnsNumber res_rns = SymRnsNumber(res_int, base);
                if (res_rns != r_rns) {
                    std::cout << "ОШИБКА: x: " << x_int << " * y:"
                        << y_int << " convert: " <<
                        res_rns.to_positional_crt()
                        << " calc: " << r_rns.to_positional_crt() << std::endl;
                    ok = false;
                }
            }
        }
    }
    if (ok) std::cout << "ТЕСТЫ ПРОйДЕНЫ!!!" << std::endl << std::endl;
}

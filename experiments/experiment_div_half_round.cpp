#include <vector>
#include <iomanip>
#include "../SymRNSFixed.cpp"
#include "../RNSFixed.cpp"

int main(void) {
    constexpr Positional_Int S = 100; // коэффициент масштабирования (должен быть не кратен ни одному из оснований)

    SymRnsBase srns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований ССОК
    RnsBase rns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований СОК

    // ДЕЛЕНИЯ: накопление ошибки округления S/2 вверх в СОК
    constexpr Positional_Int a_int_orig = 10000;

    SymRnsFixed a_srns{a_int_orig, srns_base}; // ССОК 1/S
    RnsFixed a_rns{a_int_orig, rns_base}; // СОК 1/S

    SymRnsFixed res_srns{0, srns_base}; // результат ССОК
    RnsFixed res_rns{0, rns_base}; // результат СОК

    for (Positional_Int z_int = 2; z_int <= 10000; ++z_int) {
        Positional_Int b_int_raw = 2 * a_int_orig * S;
        if (b_int_raw % (2 * z_int + 1) == 0) {
            Positional_Int b_int = b_int_raw / (2 * z_int + 1);
            if (b_int % 2 == 0 && srns_base.has_mod_inverse_sym(b_int)) {
                std::cout << "z = " << z_int << ", b = " << b_int << std::endl;
                SymRnsFixed b_srns{b_int, srns_base};
                RnsFixed b_rns{b_int, rns_base};

                SymRnsFixed div_srns = a_srns / b_srns;
                res_srns += div_srns;
                std::cout << "ССОК Фикс.точка: " << std::setprecision(100) << div_srns.to_positional_frac_crt_unscaled() << std::endl;

                RnsFixed div_rns = a_rns / b_rns;
                res_rns += div_rns;
                std::cout << " СОК Фикс.точка: " << std::setprecision(100) << div_rns.to_positional_frac_crt_unscaled() << std::endl << std::endl;
            }
        } else {
            // std::cout << b_int_raw << " делится с остатком на " << a_int_orig << std::endl;
        }
    }

    std::cout << "СУММА ССОК Фикс.точка: " << std::setprecision(100) << res_srns.to_positional_frac_crt() << std::endl;
    std::cout << "СУММА  СОК Фикс.точка: " << std::setprecision(100) << res_rns.to_positional_frac_crt() << std::endl << std::endl;

    return 0;
}

#include <vector>
#include <iomanip>
#include "SymRNSFixed.cpp"
#include "RNSFixed.cpp"

int main(void) {
    constexpr Positional_Int S = 1000; // коэффициент масштабирования (должен быть не кратен ни одному из оснований)
    // constexpr Positional_Int a_int_orig = 1000;

    SymRnsBase srns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований ССОК
    RnsBase rns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований СОК
    // Positional_Int b_max = srns_base.P / 2 / S; // максимальное значение b для тестов
    // std::cout << "b_max = " << b_max << std::endl;

    // ССОК
    SymRnsFixed a_srns{1, srns_base}; // ССОК 1/10000
    SymRnsFixed b_srns{200, srns_base}; // ССОК 1/10000
    a_srns /= b_srns;
    a_srns *= b_srns;
    std::cout << "ССОК Фикс.точка: " << ' ' << std::setprecision(100) << a_srns.to_positional_frac_crt_unscaled() << std::endl;

    RnsFixed a_rns{1, rns_base}; // ССОК 1/10000
    RnsFixed b_rns{200, rns_base}; // ССОК 1/10000
    a_rns /= b_rns;
    a_rns *= b_rns;
    std::cout << " СОК Фикс.точка: " << ' ' << std::setprecision(100) << a_rns.to_positional_frac_crt_unscaled() << std::endl;

    // // СОК
    // RnsFixed a_rns_orig{a_int_orig, rns_base}; // СОК 1/10000
    // RnsFixed a_rns = a_rns_orig; // СОК 1/10000
    // for (Positional_Int b_int = 2; b_int <= b_max; ++b_int) {
    //     if (srns_base.has_mod_inverse_sym(b_int)) {
    //         RnsFixed b_rns{b_int, rns_base};
    //         a_rns /= b_rns;
    //         a_rns *= b_rns;
    //         if (a_rns != a_rns_orig) {
    //             std::cout << "СОК Фикс.точка: " << b_int << ' ' << std::setprecision(100) << a_rns.to_positional_frac_crt_unscaled() << std::endl;
    //             break;
    //         }
    //     }
    // }

    return 0;
}

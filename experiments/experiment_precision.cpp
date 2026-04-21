#include <vector>
#include <iomanip>
#include "../SymRNSFixed.cpp"
#include "../RNSFixed.cpp"

int main(void) {
    constexpr Positional_Int S = 10000; // коэффициент масштабирования (должен быть не кратен ни одному из оснований)
    constexpr Positional_Int a_int_orig = 1000;

    SymRnsBase srns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований ССОК
    RnsBase rns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований СОК
    Positional_Int b_max = srns_base.P / 2 / S; // максимальное значение b для тестов
    std::cout << "b_max = " << b_max << std::endl;

    // машинная плавающая точка
    Positional_Float a_fl_orig = static_cast<Positional_Float>(a_int_orig) / S; // 1/10000
    Positional_Float a_fl = a_fl_orig; // 1/10000
    for (Positional_Int b_int = 2; b_int <= b_max; ++b_int) {
        if (srns_base.has_mod_inverse_sym(b_int)) {
            Positional_Float b_fl = static_cast<Positional_Float>(b_int) / S; // 2/10000 ... N/10000
            a_fl /= b_fl;
            a_fl *= b_fl;
            if (std::abs(a_fl - a_fl_orig) > 0) {
                std::cout << "Аппаратная плав.точка: " << b_int << ' ' << std::setprecision(100) << a_fl << std::endl;
                break;
            }
        }
    }

    // ССОК
    SymRnsFixed a_srns_orig{a_int_orig, srns_base}; // ССОК 1/10000
    SymRnsFixed a_srns = a_srns_orig; // ССОК 1/10000
    // СОК
    RnsFixed a_rns_orig{a_int_orig, rns_base}; // СОК 1/10000
    RnsFixed a_rns = a_rns_orig; // СОК 1/10000
    for (Positional_Int b_int = 2; b_int <= b_max; ++b_int) {
        if (srns_base.has_mod_inverse_sym(b_int)) {
            SymRnsFixed b_srns{b_int, srns_base};
            a_srns /= b_srns;
            std::cout << "ССОК Фикс.точка: " << b_int << ' ' << std::setprecision(100) << a_srns.to_positional_frac_crt_unscaled() << std::endl;
            a_srns *= b_srns;
            // std::cout << "ССОК Фикс.точка: " << b_int << ' ' << std::setprecision(100) << a_srns.to_positional_frac_crt_unscaled() << std::endl;
            if (a_srns != a_srns_orig) {
                break;
            }
            RnsFixed b_rns{b_int, rns_base};
            a_rns /= b_rns;
            std::cout << " СОК Фикс.точка: " << b_int << ' ' << std::setprecision(100) << a_rns.to_positional_frac_crt_unscaled() << std::endl;
            a_rns *= b_rns;
            // std::cout << " СОК Фикс.точка: " << b_int << ' ' << std::setprecision(100) << a_rns.to_positional_frac_crt_unscaled() << std::endl;
            if (a_rns != a_rns_orig) {
                break;
            }
        }
    }

    return 0;
}

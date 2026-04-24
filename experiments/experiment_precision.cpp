#include <vector>
#include <iomanip>
#include "../SymRNSFixed.cpp"
#include "../RNSFixed.cpp"

int main(void) {
    constexpr Positional_Int S = 10000; // коэффициент масштабирования (должен быть не кратен ни одному из оснований)
    constexpr Positional_Int a_int_orig = 1000;
    constexpr Positional_Int PRECISION = 20; // точность для вывода результатов

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
                std::cout << "Аппаратная плав.точка: " << b_int << ' ' << std::setprecision(PRECISION) << a_fl << std::endl;
                break;
            }
        }
    }

    // ССОК
    SymRnsFixed a_srns_orig{a_int_orig, srns_base}; // ССОК A/S
    SymRnsFixed a_srns = a_srns_orig; // ССОК A/S
    // СОК
    RnsFixed a_rns_orig{a_int_orig, rns_base}; // СОК A/S
    RnsFixed a_rns = a_rns_orig; // СОК A/S
    // служебные переменные
    Positional_Int num_iter = 0; // число итераций
    Positional_Float err_rns_fp = 0, err_srns_fp = 0; // накопленная ошибка деления
    for (Positional_Int b_int = 2; b_int <= b_max; ++b_int) {
        if (srns_base.has_mod_inverse_sym(b_int) && (b_int % 2 != 0 || (a_int_orig * S) % b_int != b_int / 2)) {
            Positional_Int res_int_unscaled = static_cast<Positional_Int>(std::round(static_cast<Positional_Float>(a_int_orig) * S / b_int));
            // ССОК
            SymRnsFixed b_srns{b_int, srns_base};
            a_srns /= b_srns;
            Positional_Float a_srns_fp = a_srns.to_positional_frac_crt_unscaled();
            err_srns_fp += std::abs(a_srns_fp - res_int_unscaled);
            std::cout << "ССОК Фикс.зап: B = " << b_int << " A=" << std::fixed << std::setprecision(PRECISION) << a_srns_fp  << std::endl;
            a_srns *= b_srns;
            if (a_srns != a_srns_orig) {
                std::runtime_error("Ошибка восстановления первоначального значения ССОК!");
            }
            // СОК
            RnsFixed b_rns{b_int, rns_base};
            a_rns /= b_rns;
            Positional_Float a_rns_fp = a_rns.to_positional_frac_crt_unscaled();
            err_rns_fp += std::abs(a_rns_fp - res_int_unscaled);
            std::cout << " СОК Фикс.зап: B = " << b_int << " A=" << std::fixed << std::setprecision(PRECISION) << a_rns_fp <<std::endl;
            a_rns *= b_rns;
            if (a_rns != a_rns_orig) {
                std::runtime_error("Ошибка восстановления первоначального значения СОК!");
            }
            // std::cout << err_rns_fp - err_srns_fp << std::endl;
            ++num_iter;
        }
    }

    std::cout << std::endl;
    std::cout << "Коэффициент масштаб. : " << S << std::endl;
    std::cout << "Количество итераций  : " << num_iter << std::endl;
    std::cout << "ОШИБКА  СОК Фикс.зап.: " << std::fixed << std::setprecision(PRECISION) << err_rns_fp << std::endl;
    std::cout << "ОШИБКА ССОК Фикс.зап.: " << std::fixed << std::setprecision(PRECISION) << err_srns_fp << std::endl;
    std::cout << "Теор.погрешннсть     : " << std::fixed << std::setprecision(PRECISION) << 1.0 / (2 * S) << std::endl;
    std::cout << "Факт.погрешннсть  СОК: " << std::fixed << std::setprecision(PRECISION) << err_rns_fp / num_iter << std::endl;
    std::cout << "Факт.погрешннсть ССОК: " << std::fixed << std::setprecision(PRECISION) << err_srns_fp / num_iter << std::endl;

    return 0;
}

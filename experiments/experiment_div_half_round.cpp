#include <vector>
#include <iomanip>
#include "../SymRNSFixed.cpp"
#include "../RNSFixed.cpp"

int main(void) {
    constexpr Positional_Int S = 100; // коэффициент масштабирования (должен быть не кратен ни одному из оснований)
    constexpr Positional_Int PRECISION = 20; // точность для вывода результатов

    SymRnsBase srns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований ССОК
    RnsBase rns_base{{3, 7, 11, 13, 17, 19, 23}, S}; // набор оснований СОК

    // ДЕЛЕНИЯ: накопление ошибки округления S/2 вверх в СОК
    constexpr Positional_Int a_int = 10000;

    SymRnsFixed a_srns{a_int, srns_base}; // ССОК A/S
    RnsFixed a_rns{a_int, rns_base}; // СОК A/S

    SymRnsFixed res_srns{0, srns_base}; // накопленный результат ССОК = 0
    RnsFixed res_rns{0, rns_base}; // накопленный результат СОК = 0
    Positional_Float res_fp = 0; // накопленный результат в позиционной
    Positional_Int num_iter = 0; // число итераций

    std::vector<Positional_Int> b_int_vec;
    for (Positional_Int z_int = 2; z_int <= 10000; ++z_int) {
        Positional_Int b_int = 2 * a_int * S;
        if (b_int % (2 * z_int + 1) == 0) {
            b_int /= (2 * z_int + 1);
            if (b_int % 2 == 0 && srns_base.has_mod_inverse_sym(b_int)) {
                b_int_vec.push_back(b_int);
            }
        }
    }
    std::cout << "Найдено подходящих делителей: " << b_int_vec.size() << std::endl;

    for (int i = 0; i < 100; ++i) {  // цикл итераций
        Positional_Int b_int = b_int_vec[static_cast<size_t>(std::rand()) % b_int_vec.size()];

        std::cout << "A(" << a_int << ") * S(" << S << ") / B(" << b_int << ')' << std::endl;
        SymRnsFixed b_srns{b_int, srns_base};
        RnsFixed b_rns{b_int, rns_base};

        SymRnsFixed div_srns = a_srns / b_srns;
        res_srns += div_srns;
        std::cout << "ССОК Фикс.точка: " << std::fixed << std::setprecision(PRECISION) << div_srns.to_positional_frac_crt_unscaled() << std::endl;

        RnsFixed div_rns = a_rns / b_rns;
        res_rns += div_rns;
        std::cout << " СОК Фикс.точка: " << std::fixed << std::setprecision(PRECISION) << div_rns.to_positional_frac_crt_unscaled() << std::endl << std::endl;

        res_fp += static_cast<Positional_Float>(a_int) / b_int;
        ++num_iter;
    }

    Positional_Float res_srns_fp = res_srns.to_positional_frac_crt();
    Positional_Float res_rns_fp = res_rns.to_positional_frac_crt();

    std::cout << std::endl;
    std::cout << "Коэффициент масштаб. : " << S << std::endl;
    std::cout << "Количество итераций  : " << num_iter << std::endl;
    std::cout << "СУММА   СОК Фикс.зап.: " << std::fixed << std::setprecision(PRECISION) << res_rns_fp << std::endl;
    std::cout << "СУММА  ССОК Фикс.зап.: " << std::fixed << std::setprecision(PRECISION) << res_srns_fp << std::endl;
    std::cout << "СУММА     Позиционная: " << std::fixed << std::setprecision(PRECISION) << res_fp << std::endl;
    std::cout << "ОШИБКА  СОК Фикс.зап.: " << std::fixed << std::setprecision(PRECISION) << std::abs(res_rns_fp - res_fp) << std::endl;
    std::cout << "ОШИБКА ССОК Фикс.зап.: " << std::fixed << std::setprecision(PRECISION) << std::abs(res_srns_fp - res_fp) << std::endl;
    std::cout << "Теор.погрешннсть     : " << std::fixed << std::setprecision(PRECISION) << 1.0 / (2 * S) << std::endl;
    std::cout << "Факт.погрешннсть  СОК: " << std::fixed << std::setprecision(PRECISION) << std::abs(res_rns_fp - res_fp) / num_iter << std::endl;
    std::cout << "Факт.погрешннсть ССОК: " << std::fixed << std::setprecision(PRECISION) << std::abs(res_srns_fp - res_fp) / num_iter << std::endl;

    return 0;
}

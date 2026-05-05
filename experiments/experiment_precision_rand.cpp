#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <random>
#include "../SymRNSFixed.cpp"
#include "../RNSFixed.cpp"

int main(void) {
    constexpr Positional_Int S = 10000; // коэффициент масштабирования (должен быть не кратен ни одному из оснований)
    constexpr Positional_Int a_pos = 5555; // A начальное
    SymRnsBase srns_base{{59, 61, 67, 71, 73}, S}; // набор оснований ССОК
    RnsBase rns_base{{59, 61, 67, 71, 73}, S}; // набор оснований СОК

    Positional_Int a_max = static_cast<Positional_Int>(std::sqrt((srns_base.P - 1) / 2));
    Positional_Int b_max = std::min((srns_base.P - 1) / 2 / S, a_max); // максимальное значение b для тестов
    if (!srns_base.has_mod_inverse_sym(S)) {
        throw std::runtime_error(std::format("Коэффициент масштабирования S={} кратен основанию!", S));
    }
    if (a_pos > a_max) {
        throw std::runtime_error(std::format("Недопустимое значение a_pos={} > a_max={}!", a_pos, a_max));
    }

    // СОК
    RnsFixed a_rns_orig{a_pos, rns_base}; // СОК A/S
    RnsFixed a_rns{a_rns_orig}; // СОК A/S
    // ССОК
    SymRnsFixed a_srns_orig{a_pos, srns_base}; // ССОК A/S
    SymRnsFixed a_srns{a_srns_orig}; // ССОК A/S
    std::cout << a_srns << ' ' << a_rns << std::endl;

    // подбор подходящих делителей
    std::vector<Positional_Int> b_int_vec;
    for (Positional_Int b_int = 1; b_int <= b_max; ++b_int) {
        if (srns_base.has_mod_inverse_sym(b_int)) {
            b_int_vec.push_back(b_int);
        }
    }
    // Инициализация random_device
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, b_int_vec.size() - 1);

    // служебные переменные
    Positional_Int num_iter = 0; // число итераций
    Positional_Float err_rns = 0, err_srns = 0; // накопленная (немасштабированная) ошибка
    size_t num_iter_max = 200000;
    size_t num_iter_max_print = num_iter_max - 100;
    for (size_t i = 0; i < num_iter_max; ++i) {  // цикл итераций
        // произвольный делитель из подходящих под B/2
        Positional_Int b_int = b_int_vec[dist(rng)];
        // СОК
        RnsFixed b_rns{b_int, rns_base};
        Positional_Int a_rns_pos_pre = static_cast<Positional_Int>(std::round(a_rns.to_positional_frac_crt_unscaled()));
        if (i > num_iter_max_print) std::cout << std::format(" СОК: ({}/{} / {}/{}", a_rns_pos_pre, S, b_int, S);
        a_rns /= b_rns;
        Positional_Int a_rns_pos_div = static_cast<Positional_Int>(std::round(a_rns.to_positional_frac_crt_unscaled()));
        if (i > num_iter_max_print) std::cout << std::format(" = {}/{})", a_rns_pos_div, S);

        a_rns *= b_rns;
        Positional_Int a_rns_pos_mul = static_cast<Positional_Int>(std::round(a_rns.to_positional_frac_crt_unscaled()));
        if (i > num_iter_max_print) std::cout << std::format(" * {}/{} = {}/{}\n", b_int, S, a_rns_pos_mul, S);
        Positional_Int err_rns_iter = std::abs(a_rns_pos_mul - a_rns_pos_pre);
        Positional_Int err_theor = static_cast<Positional_Int>(std::round(static_cast<Positional_Float>(b_int) / S / 2));
        if (err_rns_iter > err_theor) {
            throw std::runtime_error(std::format("Превышение теоретической ошибки err_rns_iter={} > {}!", err_rns_iter, err_theor));
        }
        err_rns += err_rns_iter;

        // ССОК
        SymRnsFixed b_srns{b_int, srns_base};
        Positional_Int a_srns_pos_pre = static_cast<Positional_Int>(std::round(a_srns.to_positional_frac_crt_unscaled()));
        if (i > num_iter_max_print) std::cout << std::format("ССОК: ({}/{} / {}/{}", a_srns_pos_pre, S, b_int, S);
        a_srns /= b_srns;
        Positional_Int a_srns_pos_div = static_cast<Positional_Int>(std::round(a_srns.to_positional_frac_crt_unscaled()));
        if (i > num_iter_max_print) std::cout << std::format(" = {}/{})", a_srns_pos_div, S);

        a_srns *= b_srns;
        Positional_Int a_srns_pos_mul = static_cast<Positional_Int>(std::round(a_srns.to_positional_frac_crt_unscaled()));
        if (i > num_iter_max_print) std::cout << std::format(" * {}/{} = {}/{}\n", b_int, S, a_srns_pos_mul, S);
        Positional_Int err_srns_iter = std::abs(a_srns_pos_mul - a_srns_pos_pre);
        if (err_srns_iter > err_theor) {
            throw std::runtime_error(std::format("Превышение теоретической ошибки err_srns_iter={} > {}!", err_srns_iter, err_theor));
        }
        err_srns += err_srns_iter;

        char err_rns_iter_sign = (err_rns_iter < 0) ? '-' : '+';
        char err_srns_iter_sign = (err_srns_iter < 0) ? '-' : '+';
        if (i > num_iter_max_print) {
            std::cout << std::format("Ошибка СОК={}({}{})/{}, Ошибка ССОК={}({}{})/{}", err_rns, err_rns_iter_sign, std::abs(err_rns_iter), S, err_srns, err_srns_iter_sign, std::abs(err_srns_iter), S);
            if (err_rns_iter != err_srns_iter) {
                std::cout << " <-- отличие накопления ошибки СОК и ССОК!";
            }
            std::cout << std::endl << std::endl;
        }
        ++num_iter;
    }

    std::cout << std::format("Начальное значение         : {}/{}\n", a_pos, S);
    std::cout << std::format("Максимальный делитель      : {}/{}\n", b_max, S);
    std::cout << std::format("Подходящих делителей       : {}\n", b_int_vec.size());
    std::cout << std::format("Количество итераций        : {}\n", num_iter);
    std::cout << std::format("Накоп.ошибка СОК  фикс.зап.: {}/{}\n", err_rns, S);
    std::cout << std::format("Накоп.ошибка ССОК фикс.зап.: {}/{}\n", err_srns, S);
    std::cout << std::format("Теор.ошибка на операцию    : {:.10f}\n", 1.0 / (2 * S));
    std::cout << std::format("Факт.ошибка  СОК фикс.зап. : {:.10f}\n", err_rns / S / num_iter / 2);
    std::cout << std::format("Факт.ошибка ССОК фикс.зап. : {:.10f} (в {:.5} раз меньше СОК)\n", err_srns / S / num_iter / 2, err_rns / err_srns);

    return 0;
}

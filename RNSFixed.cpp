#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

typedef int32_t Positional_Int;       // позиционное целочисленное
typedef double Positional_Float;       // позиционное вещественное
typedef int16_t Module;               // целочисленный модуль разряда ССОК
typedef std::vector<Module> Modules;  // вектор модулей разрядов ССОК

class RnsNumber;

// класс оснований СОК
class RnsBase {
   public:
    Modules p;  // вектор оснований (p1, p2 .. pn)
    std::vector<Positional_Int> B;  // ортогональные базисы, соответствующие (1; 0; 0..),
                                    // (0; 1; 0...), (0; 0; 1..) для конвертации в позиционную СИ
    std::vector<Positional_Int> m;  // веса ортогональных базисов
    std::vector<std::vector<Module>> t;  // t_ij для MRC
    Positional_Int P;   // динамический диапазон P = p1 * ... * pn,
                        // кол-во чисел в представлении в данном наборе оснований
                        // ex. 3 * 5 * 7 = 105
    std::vector<Positional_Int> Pi;   // P / p[i]
    Positional_Int S;  // коэф масштабирования
    RnsBase(const Modules& p0, Positional_Int S0);  // конструктор из вектора оснований
};

// класс представления числа в СОК
class RnsNumber {
   public:
    std::reference_wrapper<const RnsBase> base;  // ссылка на класс оснований
    Modules a;                                   // вектор значений модулей (a1, a2 .. an)

    RnsNumber(Positional_Int x, const RnsBase& base0);  // из целого позиционного
    RnsNumber(Positional_Float x, const RnsBase& base0);    // из вещест. позиционного
    RnsNumber(const Modules& a0, const RnsBase& base0); // из вектора остатков
    template<typename ModType> static ModType mod(Positional_Int x, ModType p); // математический модуль числа
    template<typename ModType> static ModType mod_inverse(Positional_Int a, ModType p); // обратное по модулю число a^-1 mod p
    RnsNumber& operator+=(const RnsNumber& y);  // this += y
    RnsNumber operator+(const RnsNumber& y) const;  // z = this + y
    RnsNumber& operator-=(const RnsNumber& y);  // this -= y
    RnsNumber operator-(const RnsNumber& y) const;  // z = this - y
    RnsNumber operator-() const;  // Унарный минус z = -this
    RnsNumber& operator*=(const RnsNumber& y);  // this *= y
    RnsNumber operator*(const RnsNumber& y) const;  // z = this * y
    RnsNumber& operator/=(const RnsNumber& y);  // this /= y
    RnsNumber operator/(const RnsNumber& y) const;  // z = this / y
    bool operator==(const RnsNumber& y) const;  // (this == y)
    bool operator!=(const RnsNumber& y) const;  // (this != y)
    void div_int(const RnsNumber& y); // деление нацело
    friend std::ostream& operator<<(std::ostream& os, const RnsNumber& y);  // для вывода на cout
    Positional_Float to_positional_ort() const;  // в позиционное представление
    Positional_Int to_positional_frac_ort_raw() const;  // в позиционное представление (без масштабирования)
    Positional_Float to_positional_frac_ort() const;  // в позиционное представление
    Positional_Float to_positional_mrc() const;  // в позиционное представление
    Positional_Int get_rank() const; // получение ранга числа
    Positional_Int get_remainder(Positional_Int divisor) const; // остаток от деления
};

RnsBase::RnsBase(const Modules& p0, Positional_Int S0) :
        p{p0},
        B(p.size()),
        m(p.size()),
        t(p.size(), Modules(p.size())),
        Pi(p.size()),
        S{S0}
{
    P = 1;
    bool first = true;
    std::cout << "(";
    for (Module pi : p) {
        if (!first) std::cout << ", ";
        std::cout << pi;
        P *= pi;
        first = false;
    }
    std::cout << ") P=" << P;
    // Ms = RnsNumber::mod(P, S);
    // std::cout << " Ms=" << Ms;
    std::cout << std::endl;

    for (size_t i = 0; i < p.size(); ++i) {
        Pi[i] = P / p[i];
        m[i] = RnsNumber::mod_inverse(Pi[i], p[i]);  // вес i-го ортогонального базиса
        B[i] = m[i] * Pi[i];  // i-й ортогональный базис
        // C[i] = RnsNumber::mod(Pi, S);
        RnsNumber B_rns{B[i], *this};
        std::cout << "ORT: Ортогональный базис B" << i+1 << '=' << B[i] << ' ' << B_rns << ", вес базиса m" << i+1 << "=" << m[i]
                  << " P" << i+1 << '=' << Pi[i];
        // std::cout << " C[i]=" << C[i];
        std::cout << std::endl;
    }

    for (size_t i = 0; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            t[i][j] = RnsNumber::mod_inverse(p[i], p[j]);
            std::cout << "MRC: t[" << i+1 << "][" << j+1 <<"] = " << t[i][j] << std::endl;
        }
    }
}

Positional_Float RnsNumber::to_positional_ort() const {
    Positional_Int res_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        res_int += a[i] * base.get().B[i];
    }
    res_int = mod(res_int, base.get().P);  // по модулю максимального представления
    return (static_cast<Positional_Float>(res_int) / base.get().S);
}

Positional_Int RnsNumber::to_positional_frac_ort_raw() const {
    Positional_Float pos_float = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        pos_float += static_cast<Positional_Float>(RnsNumber::mod(a[i] * base.get().m[i], base.get().p[i])) / base.get().p[i];
        // res += static_cast<Positional_Float>(a[i] * base.get().m[i]) / base.get().p[i];
    }
    // std::cout << res;
    Positional_Float int_part;
    pos_float = std::modf(pos_float, &int_part) * base.get().P;
    Positional_Int pos_int = static_cast<Positional_Int>(std::round(pos_float));
    // std::cout << ' ' << pos_int << std::endl;
    return pos_int;
}

Positional_Float RnsNumber::to_positional_frac_ort() const {
    Positional_Float res = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        res += static_cast<Positional_Float>(RnsNumber::mod(a[i] * base.get().m[i], base.get().p[i])) / base.get().p[i];
        // res += static_cast<Positional_Float>(a[i] * base.get().m[i]) / base.get().p[i];
    }
    // std::cout << res;
    Positional_Float int_part;
    res = std::modf(res, &int_part) * base.get().P;
    // std::cout << ' ' << res << std::endl;
    Positional_Float pos_float = static_cast<Positional_Float>(to_positional_frac_ort_raw()) / base.get().S;
    return pos_float;
}


Positional_Int RnsNumber::get_rank() const {
    Positional_Float pos_float = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        pos_float += static_cast<Positional_Float>(RnsNumber::mod(a[i] * base.get().m[i], base.get().p[i])) / base.get().p[i];
        //pos_float += static_cast<Positional_Float>(a[i] * base.get().m[i]) / base.get().p[i];
    }
    // std::cout << pos_float;
    Positional_Int rank = static_cast<Positional_Int>(pos_float);
    std::cout << " Rank=" << rank << " Pos_Float=" << pos_float << std::endl;
    return rank;
}

Positional_Int RnsNumber::get_remainder(Positional_Int divisor) const {
    Positional_Int rem_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        rem_int += RnsNumber::mod(base.get().Pi[i] * RnsNumber::mod(a[i] * base.get().m[i], base.get().p[i]), divisor);
    }
    Positional_Int rank = get_rank();
    rem_int -= RnsNumber::mod(rank * base.get().P, divisor);
    rem_int = RnsNumber::mod(rem_int, divisor);
    return rem_int;
}

Positional_Float RnsNumber::to_positional_mrc() const {
    Modules x(a.size());
    x[0] = a[0];
    // std::cout << "MRC: (x1 = " << x[0] << ") +" << std::endl;
    Positional_Int Ai = x[0], Pi_1 = 1, xi;
    for (size_t i = 1; i < a.size(); ++i) {
        xi = a[i];
        // std::cout << "MRC: (x" << i+1 << " = (";
        // for (size_t j = 0; j < i; ++j) std::cout << '(';
        // std::cout << xi;
        for (size_t j = 0; j < i; ++j) {
            // std::cout << " - " << x[j] << ") * " << base.get().t[j][i];
            xi = (xi - x[j]) * base.get().t[j][i];
        }
        Pi_1 *= base.get().p[i-1];
        // std::cout << " = " << xi << ") mod " << base.get().p[i];
        x[i] = mod(xi, base.get().p[i]);
        // std::cout << " = " << x[i] << ") * " << Pi_1;
        // if (i + 1 != a.size()) std::cout << " +" << std::endl;
        Ai += Pi_1 * x[i];
    }
    // std::cout << " = " << Ai << std::endl;
    // std::cout << Ai << ';';
    // for (size_t i = 0; i < a.size(); ++i) std::cout << x[i] << ';';
    // std::cout << std::endl;
    return (static_cast<Positional_Float>(Ai) / base.get().S);
}

RnsNumber::RnsNumber(Positional_Int x, const RnsBase& base0)
    : base{base0}, a(base.get().p.size()) {
    size_t i = 0;
    for (Module pi : base.get().p) {
        a[i++] = mod(x, pi);
    }
}

RnsNumber::RnsNumber(Positional_Float x, const RnsBase& base0)
    : RnsNumber(static_cast<Positional_Int>(std::round(x * base0.S)), base0) {}

RnsNumber::RnsNumber(const Modules& a0, const RnsBase& base0)
    : base{base0}, a{a0} {}

template<typename ModType> ModType RnsNumber::mod(Positional_Int x, ModType p) {
    ModType r = x % p;
    if (r < 0) r += p;
    return r;
}

template<typename ModType> ModType RnsNumber::mod_inverse(Positional_Int x, ModType p) {
    Positional_Int p1 = p, t, q;
    Positional_Int x0 = 0, x1 = 1;
    x = mod(x, p);
    if (p == 1 || x == 0) {
        throw std::runtime_error("Division by zero in mod_inverse(" + std::to_string(x) + ", " + std::to_string(p) + ")");
    }
    // Алгоритм Евклида
    while (x > 1) {
        q = x / p1; // Частное
        t = p1;
        p1 = x % p1, x = t;
        t = x0;
        x0 = x1 - q * x0, x1 = t; // Обновляем x0 и x1
    }
    if (x1 < 0) x1 += p; // Корректируем отрицательное
    return static_cast<ModType>(x1);
}

RnsNumber& RnsNumber::operator+=(const RnsNumber& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod(a[i] + y.a[i], base.get().p[i]);
    }
    return *this;
}

RnsNumber RnsNumber::operator+(const RnsNumber& y) const {
    RnsNumber res{*this};
    res += y;
    return res;
}

RnsNumber& RnsNumber::operator-=(const RnsNumber& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod(a[i] - y.a[i], base.get().p[i]);
    }
    return *this;
}

RnsNumber RnsNumber::operator-(const RnsNumber& y) const {
    RnsNumber res{*this};
    res -= y;
    return res;
}

RnsNumber RnsNumber::operator-() const {
    RnsNumber res{*this};
    for (size_t i = 0; i < a.size(); ++i) {
        res.a[i] = -a[i];
    }
    return res;
}

void RnsNumber::div_int(const RnsNumber& y) {  // деление нацело (когда точно делится)
    for (size_t i = 0; i < a.size(); ++i) {
        Positional_Int Si = RnsNumber::mod_inverse(y.a[i], base.get().p[i]);
        std::cout << ' ' << Si;
        a[i] = mod(a[i] * Si, base.get().p[i]);
    }
    std::cout << std::endl;
}

RnsNumber& RnsNumber::operator*=(const RnsNumber& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod(a[i] * y.a[i], base.get().p[i]);
    }

    std::cout << "RAW=" << to_positional_frac_ort();
    // масштабирование на S
    Positional_Int remainder = get_remainder(base.get().S);
    std::cout << " REM=" << remainder;

    if (remainder < base.get().S / 2) {  // округление к ближайшему
        RnsNumber remainder_rns{remainder, base};
        *this -= remainder_rns;
    } else {
        RnsNumber remainder_rns{base.get().S - remainder, base};
        *this += remainder_rns;
    }
    std::cout << " WO_REM=" << to_positional_frac_ort() << ' ' << *this;

    RnsNumber S_rns{base.get().S, base};
    std::cout << " S=" << S_rns;
    div_int(S_rns);

    std::cout << std::endl;
    return *this;
}

RnsNumber& RnsNumber::operator/=(const RnsNumber& y) {
    RnsNumber S_rns{base.get().S, base};
    for (size_t i = 0; i < a.size(); ++i) { // a * S
        a[i] = mod(a[i] * S_rns.a[i], base.get().p[i]);
    }

    std::cout << "RAW=" << to_positional_frac_ort_raw();
    // масштабирование на y
    Positional_Int y_int = y.to_positional_frac_ort_raw();
    Positional_Int remainder = get_remainder(y_int);
    std::cout << " REM=" << remainder;

    if (remainder < y_int / 2) {  // округление к ближайшему
        RnsNumber remainder_rns{remainder, base};
        *this -= remainder_rns;
    } else {
        RnsNumber remainder_rns{y_int - remainder, base};
        *this += remainder_rns;
    }
    std::cout << " WO_REM=" << to_positional_frac_ort_raw() << ' ' << *this;
    div_int(y);

    std::cout << std::endl;
    return *this;
}

RnsNumber RnsNumber::operator*(const RnsNumber& y) const {
    RnsNumber res{*this};
    res *= y;
    return res;
}

RnsNumber RnsNumber::operator/(const RnsNumber& y) const {
    RnsNumber res{*this};
    res /= y;
    return res;
}

bool RnsNumber::operator==(const RnsNumber& y) const {
    if (a.size() != y.a.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != y.a[i]) return false;
    }
    return true;
}

bool RnsNumber::operator!=(const RnsNumber& y) const { return (!(*this == y)); }

std::ostream& operator<<(std::ostream& os, const RnsNumber& y) {
    size_t sz = y.a.size();
    os << '(';
    for (size_t i = 0; i < sz; ++i) {
        if (i) os << ";";
        os << y.a[i];
    }
    os << ')';
    return os;
}

int main() {
    // constexpr Positional_Int S = 100;
    // // RnsBase base{{3, 5, 7}, S};
    // RnsBase base{{3, 7, 11, 13, 17}, S};
    constexpr Positional_Int S = 10;
    RnsBase base{{3, 7, 11, 13}, S};


    // for (int x= -10; x <= 10; ++x) {
    //     std::cout << "x=" << x << " x^-1=" << RnsNumber::mod_inverse(x, 5) << std::endl;
    // }

    // тесты
    // Positional_Int q = 1;  // делитель для fixed point тестов
    bool ok = true;

    // целочисленные проверки конвертации в позиционную ИС
    for (Positional_Int x_pos = 0; x_pos < base.P; ++x_pos) {
        RnsNumber x_rns = RnsNumber{x_pos, base};
        Positional_Float x_pos_ort = x_rns.to_positional_ort();
        Positional_Float x_pos_mrc = x_rns.to_positional_mrc();
        Positional_Float x_pos_frac = x_rns.to_positional_frac_ort();
        if (x_pos != static_cast<Positional_Int>(std::round(x_pos_ort * S)) ||
            x_pos != static_cast<Positional_Int>(std::round(x_pos_mrc * S)) ||
            x_pos != static_cast<Positional_Int>(std::round(x_pos_frac * S))
        ) {
            std::cout << "ОШИБКА: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_ort << " MRC=" << x_pos_mrc << " MRC_F=" << x_pos_frac << std::endl;
            ok = false;
        } else {
        //     std::cout << "ИНФО: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_ort << " MRC=" << x_pos_mrc << " MRC_F=" << x_pos_frac << std::endl;
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

    // примеры использования fixed point
    // Positional_Int a_int = 89; // 90/S
    // Positional_Int b_int = 108; // 100/S
    Positional_Int a_int = 54; // 90/S
    Positional_Int b_int = 33; // 100/S
    RnsNumber a{a_int, base};
    RnsNumber b{b_int, base};
    RnsNumber c = a * b;
    Positional_Float c_fl = c.to_positional_frac_ort();

    Positional_Float a_fl = static_cast<Positional_Float>(a_int) / S; // 90/S
    Positional_Float b_fl = static_cast<Positional_Float>(b_int) / S; // 100/S
    Positional_Float cc_fl = a_fl * b_fl;
    Positional_Float eps = std::abs(c_fl - cc_fl);
    Positional_Float prec = 0.5 / S;

    if (eps > prec) {
        ok = false;
        std::cout << "ОШИБКА: ";
    } else {
        std::cout << "OK: ";
    }
    std::cout << c_fl << " ? " << cc_fl << " eps=" << eps <<
            " prec=" << prec << std::endl;

    std::cout << a.to_positional_ort() << '*' << b.to_positional_ort() << '=' << c.to_positional_ort() << ' ' << c << std::endl;

    if (ok) std::cout << "ТЕСТЫ ПРОйДЕНЫ!!!" << std::endl << std::endl;
}

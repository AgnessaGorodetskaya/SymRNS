#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

typedef int32_t Positional_Int;       // позиционное целочисленное
typedef float Positional_Float;       // позиционное вещественное
typedef int16_t Module;               // целочисленный модуль разряда ССОК
typedef std::vector<Module> Modules;  // вектор модулей разрядов ССОК

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
    RnsBase(const Modules& p0);  // конструктор из вектора оснований
};

// класс представления числа в СОК
class RnsNumber {
   public:
    std::reference_wrapper<const RnsBase> base;  // ссылка на класс оснований
    Modules a;                                   // вектор значений модулей (a1, a2 .. an)

    RnsNumber(Positional_Int x, const RnsBase& base0);  // из целого позиционного
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
    bool operator==(const RnsNumber& y) const;  // (this == y)
    bool operator!=(const RnsNumber& y) const;  // (this != y)
    friend std::ostream& operator<<(std::ostream& os, const RnsNumber& y);  // для вывода на cout
    Positional_Int to_positional_ort() const;  // в позиционное представление
    Positional_Int to_positional_mrc() const;  // в позиционное представление
    Positional_Int get_rank() const; // получение ранга числа
    // RnsNumber round(Positional_Int b) const;  // округление к ближайшему кратному x
};

RnsBase::RnsBase(const Modules& p0) : p{p0}, B(p.size()), m(p.size()), t(p.size(), Modules(p.size())) {
    P = 1;
    bool first = true;
    std::cout << "(";
    for (Module pi : p) {
        if (!first) std::cout << ", ";
        std::cout << pi;
        P *= pi;
        first = false;
    }
    std::cout << ") P=" << P << std::endl;

    for (size_t i = 0; i < p.size(); ++i) {
        Positional_Int Pi = P / p[i];
        m[i] = RnsNumber::mod_inverse(Pi, p[i]);  // вес i-го ортогонального базиса
        B[i] = m[i] * Pi;  // i-й ортогональный базис
        RnsNumber B_rns{B[i], *this};
        std::cout << "ORT: Ортогональный базис B" << i+1 << '=' << B[i] << ' ' << B_rns << ", вес базиса m" << i+1 << "=" << m[i] << std::endl;
    }

    for (size_t i = 0; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            t[i][j] = RnsNumber::mod_inverse(p[i], p[j]);
            std::cout << "MRC: t[" << i+1 << "][" << j+1 <<"] = " << t[i][j] << std::endl;
        }
    }
}

Positional_Int RnsNumber::to_positional_ort() const {
    Positional_Int res_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        res_int += a[i] * base.get().B[i];
    }
    return mod(res_int, base.get().P);  // по модулю максимального представления
}

Positional_Int RnsNumber::to_positional_mrc() const {
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
    std::cout << Ai << ';';
    for (size_t i = 0; i < a.size(); ++i) std::cout << x[i] << ';';
    std::cout << std::endl;
    return Ai;
}

RnsNumber::RnsNumber(Positional_Int x, const RnsBase& base0)
    : base{base0}, a(base.get().p.size()) {
    size_t i = 0;
    for (Module pi : base.get().p) {
        a[i++] = mod(x, pi);
    }
}

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
    if (p == 1 || x == 0) return 0;
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

RnsNumber& RnsNumber::operator*=(const RnsNumber& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod(a[i] * y.a[i], base.get().p[i]);
    }
    return *this;
}

RnsNumber RnsNumber::operator*(const RnsNumber& y) const {
    RnsNumber res{*this};
    res *= y;
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

Positional_Int RnsNumber::get_rank() const {
    Positional_Int pos_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        pos_int += a[i] * base.get().B[i];
    }
    Positional_Int rank = pos_int / base.get().P;
    std::cout << "rank of "<< pos_int << "=" << rank << std::endl;
    return rank;
}

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
        Positional_Int x_pos_ort = x_rns.to_positional_ort();
        Positional_Int x_pos_mrc = x_rns.to_positional_mrc();
        if (x_pos != x_pos_ort || x_pos != x_pos_mrc) {
            std::cout << "ОШИБКА: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_ort << " MRC=" << x_pos_mrc << std::endl;
            ok = false;
        } else {
            // std::cout << "ИНФО: x=" << x_pos << " СОК=" << x_rns << " ORT=" << x_pos_ort << " MRC=" << x_pos_mrc << std::endl;
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

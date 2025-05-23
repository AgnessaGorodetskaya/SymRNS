#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

typedef int32_t Positional_Int;       // позиционное целочисленное
typedef float Positional_Float;       // позиционное вещественное
typedef int16_t Module;               // целочисленный модуль разряда ССОК
typedef std::vector<Module> Modules;  // вектор модулей разрядов ССОК

// класс оснований ССОК
class SymRnsBase {
   public:
    Modules p;  // вектор оснований (p1, p2 .. pn)
    std::vector<Positional_Int> B;  // ортогональные базисы, соответствующие (1; 0; 0..),
                                    // (0; 1; 0...), (0; 0; 1..) для конвертации в позиционную СИ
    std::vector<Positional_Int> m;  // веса ортогональных базисов
    std::vector<std::vector<Module>> t;  // t_ij для MRC
    Positional_Int P;   // динамический диапазон P = p1 * ... * pn,
                        // кол-во чисел в представлении в данном наборе оснований
                        // ex. 3 * 5 * 7 = 105
    SymRnsBase(const Modules& p0);  // конструктор из вектора оснований
};

// класс представления числа в ССОК
class SymRnsNumber {
   public:
    std::reference_wrapper<const SymRnsBase> base;  // ссылка на класс оснований
    Modules a;                                      // вектор значений модулей (a1, a2 .. an)
    Positional_Int q_int;  // делитель представления с фиксированной точкой

    SymRnsNumber(Positional_Int x, const SymRnsBase& base0, Positional_Int q0 = 1);  // из целого позиционного
    SymRnsNumber(Positional_Float x, const SymRnsBase& base0, Positional_Int q0);    // из вещест. позиционного
    SymRnsNumber(const Modules& a0, const SymRnsBase& base0, Positional_Int q0 = 1); // из вектора остатков
    template<typename ModType> static ModType mod(Positional_Int x, ModType p); // математический модуль числа
    template<typename ModType> static ModType mod_sym(Positional_Int x, ModType p);  // симметричный модуль
    template<typename ModType> static ModType mod_inverse_sym(Positional_Int a, ModType p); // обратное по модулю число a^-1 mod p
    template<typename ModType> static inline ModType modsym_to_mod(ModType a, ModType p) {
        if (a < 0) a += p;
        return a;
    }
    template<typename ModType> static inline ModType mod_to_modsym(ModType a, ModType p) {
        if (a > p / 2) a -= p;
        return a;
    }
    SymRnsNumber& operator+=(const SymRnsNumber& y);  // this += y
    SymRnsNumber operator+(const SymRnsNumber& y) const;  // z = this + y
    SymRnsNumber& operator-=(const SymRnsNumber& y);  // this -= y
    SymRnsNumber operator-(const SymRnsNumber& y) const;  // z = this - y
    SymRnsNumber operator-() const;  // Унарный минус z = -this
    SymRnsNumber& operator*=(const SymRnsNumber& y);  // this *= y
    SymRnsNumber operator*(const SymRnsNumber& y) const;  // z = this * y
    bool operator==(const SymRnsNumber& y) const;  // (this == y)
    bool operator!=(const SymRnsNumber& y) const;  // (this != y)
    friend std::ostream& operator<<(std::ostream& os, const SymRnsNumber& y);  // для вывода на cout
    Positional_Int to_positional_int() const;  // в позиционное (int) представление
    Positional_Int to_positional_int_MRC() const;  // в позиционное представление
    Positional_Int to_positional_int0() const;  // в позиционное из СОК (int) представление
    Positional_Float to_positional(Positional_Int* x_int = NULL) const;  // в позиционное представление
    Positional_Int get_rank(bool* sign=NULL) const;  // получение ранга числа
    Positional_Int get_rank_x(Positional_Int x) const; // получение ранга числа * x
    SymRnsNumber round(Positional_Int b) const;  // округление к ближайшему кратному x
};

SymRnsBase::SymRnsBase(const Modules& p0) : p{p0}, B(p.size()), m(p.size()), t(p.size(), Modules(p.size())) {
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

    std::cout << "Поиск ортогональных базисов..." << std::endl;
    for (size_t i = 0; i < p.size(); ++i) {
        Positional_Int Pi = P / p[i];
        m[i] = SymRnsNumber::mod_inverse_sym(Pi, p[i]);  // вес i-го ортогонального базиса
        B[i] = m[i] * Pi;  // i-й ортогональный базис
        SymRnsNumber B_rns{B[i], *this, 1};
        std::cout << "Ортогональный базис B" << i+1 << '=' << B[i] << ' ' << B_rns << ", вес базиса m" << i+1 << "=" << m[i] << std::endl;
    }

    std::cout << "Поиск t_ij..." << std::endl;
    for (size_t i = 0; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            t[i][j] = SymRnsNumber::mod_inverse_sym(p[i], p[j]);
            std::cout << "t[" <<i<< "][" << j <<"] = " << t[i][j] << std::endl;
        }
    }
}

Positional_Int SymRnsNumber::to_positional_int() const {
    Positional_Int res_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        res_int += a[i] * base.get().B[i];
    }
    return mod_sym(res_int, base.get().P);  // по симметричному модулю динамического диапазона
}

Positional_Int SymRnsNumber::to_positional_int_MRC() const {
    Modules x(a.size());
    x[0] = a[0];
    Positional_Int Ai = x[0], Pi_1 = 1;
    for (size_t i = 1; i < a.size(); ++i) {
        x[i] = a[i];
        for (size_t j = 0; j < i; ++j) {
            x[i] = (x[i] - x[j]) * base.get().t[j][i];
        }
        x[i] = mod_sym(x[i], base.get().p[i]);
        std::cout << "x"<<i<<" = " << x[i] << std::endl;
        Pi_1 *= base.get().p[i-1];
        Ai += x[i] * Pi_1;
    }
    return Ai;
}

Positional_Int SymRnsNumber::to_positional_int0() const {  // воспринимать как несимметричное представление
    Positional_Int res_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        res_int += modsym_to_mod(a[i], base.get().p[i]) *
            modsym_to_mod(base.get().B[i], base.get().P);
    }
    return mod(res_int, base.get().P);  // по модулю динамического диапазона
}

Positional_Float SymRnsNumber::to_positional(Positional_Int* x_int) const {
    Positional_Int res_int = to_positional_int();
    if (x_int) *x_int = res_int;  // опционально целочисленный числитель возвращаем по указателю в параметре
    return static_cast<Positional_Float>(res_int) / q_int;
}

SymRnsNumber::SymRnsNumber(Positional_Int x, const SymRnsBase& base0, Positional_Int q0)
    : base{base0}, a(base.get().p.size()), q_int(q0) {
    size_t i = 0;
    for (Module pi : base.get().p) {
        a[i++] = mod_sym(x, pi);
    }
}

SymRnsNumber::SymRnsNumber(Positional_Float x, const SymRnsBase& base0, Positional_Int q0)
    : SymRnsNumber(static_cast<Positional_Int>(std::round(x * q0)), base0, q0) {}

SymRnsNumber::SymRnsNumber(const Modules& a0, const SymRnsBase& base0, Positional_Int q0)
    : base{base0}, a{a0}, q_int(q0) {}

template<typename ModType> ModType SymRnsNumber::mod(Positional_Int x, ModType p) {
    ModType r = x % p;
    if (r < 0) r += p;
    return r;
}

template<typename ModType> ModType SymRnsNumber::mod_sym(Positional_Int x, ModType p) {
    ModType p_half = p / 2;
    ModType r = x % p;
    if (r > p_half) { r -= p; }
    else if (r < -p_half) { r += p; }
    return r;
}

template<typename ModType> ModType SymRnsNumber::mod_inverse_sym(Positional_Int x, ModType p) {
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
    Positional_Int p_half = p / 2;  // вводим в симметричный диапазон
    if (x1 > p_half) { x1 -= p; }
    else if (x1 < -p_half) { x1 += p; }
    return static_cast<ModType>(x1);
}

SymRnsNumber& SymRnsNumber::operator+=(const SymRnsNumber& y) {
    Positional_Int xm, ym;
    if (q_int == y.q_int) {  // сложение с одинаковыми делителями (знаменателями)
        xm = 1;
        ym = 1;
    } else {
        Positional_Int lcm = std::lcm(q_int, y.q_int);  // делитель суммы = НОК делителей слагаемых
                                                        // 5/6 + 8/9 = (5*3 + 8*2) / НОК(6,9)=18
        xm = lcm / q_int;
        ym = lcm / y.q_int;
        q_int = lcm;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(xm * a[i] + ym * y.a[i], base.get().p[i]);
    }
    return *this;
}

SymRnsNumber SymRnsNumber::operator+(const SymRnsNumber& y) const {
    SymRnsNumber res{*this};
    res += y;
    return res;
}

SymRnsNumber& SymRnsNumber::operator-=(const SymRnsNumber& y) {
    Positional_Int xm, ym;
    if (q_int == y.q_int) {  // вычитание с одинаковыми делителями (знаменателями)
        xm = 1;
        ym = 1;
    } else {
        Positional_Int lcm =
            std::lcm(q_int, y.q_int);  // делитель разности = НОК делителей слагаемых
                                       // 5/6 - 8/9 = (5*3 - 8*2) / НОК(6,9)=18
        xm = lcm / q_int;
        ym = lcm / y.q_int;
        q_int = lcm;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(xm * a[i] - ym * y.a[i], base.get().p[i]);
    }
    return *this;
}

SymRnsNumber SymRnsNumber::operator-(const SymRnsNumber& y) const {
    SymRnsNumber res{*this};
    res -= y;
    return res;
}

SymRnsNumber SymRnsNumber::operator-() const {
    SymRnsNumber res{*this};
    for (size_t i = 0; i < a.size(); ++i) {
        res.a[i] = -a[i];
    }
    return res;
}

SymRnsNumber& SymRnsNumber::operator*=(const SymRnsNumber& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(a[i] * y.a[i], base.get().p[i]);
    }
    q_int *= y.q_int;  // при умножении fixed-делители умножаются 2/3 * 4/5 = (2*4)/(3*5)
    return *this;
}

SymRnsNumber SymRnsNumber::operator*(const SymRnsNumber& y) const {
    SymRnsNumber res{*this};
    res *= y;
    return res;
}

bool SymRnsNumber::operator==(const SymRnsNumber& y) const {
    if (a.size() != y.a.size() || q_int != y.q_int) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != y.a[i]) return false;
    }
    return true;
}

bool SymRnsNumber::operator!=(const SymRnsNumber& y) const { return (!(*this == y)); }

std::ostream& operator<<(std::ostream& os, const SymRnsNumber& y) {
    size_t sz = y.a.size();
    os << '(';
    for (size_t i = 0; i < sz; ++i) {
        if (i) os << ";";
        os << y.a[i];
    }
    if (y.q_int != 1) os << ";q=" << y.q_int;
    os << ')';
    return os;
}

Positional_Int SymRnsNumber::get_rank(bool* sign) const {
    Positional_Int pos_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        Positional_Int Bi = base.get().B[i];
        Module ai = a[i];
        Bi = modsym_to_mod(Bi, base.get().P); // выведем B[i] из симм.диапазона        
        ai = modsym_to_mod(ai, base.get().p[i]);  // выведем из симм
        pos_int += ai * Bi;
    }
    Positional_Int P = base.get().P;
    Positional_Int rank = pos_int / P;
    if (sign != NULL) {
        if (pos_int % P > P / 2) {
            *sign = true;  // negative
        } else {
            *sign = false;  // positive
        }
    }
    return rank;
}

Positional_Int SymRnsNumber::get_rank_x(Positional_Int x) const { // вычисление ранга числа * x
    Positional_Int r_a = get_rank(); // ранг текущего числа
    Positional_Int r_xa = r_a * x;
    std::cout << "r_xa = " << r_a << "*" << x;
    for (size_t i = 0; i < a.size(); ++i) {
        Module pi = base.get().p[i];
        Module ai = modsym_to_mod(a[i], pi);
        Positional_Int t = (ai * x / pi) * base.get().m[i];
        if (r_a >= 0) {
            r_xa -= t;
            std::cout << "- (" << ai << "*" << x << "/" << pi <<") * " << base.get().m[i];
        } else {
            r_xa += t;
            std::cout << "+ (" << ai << "*" << x << "/" << pi <<") * " << base.get().m[i];
        }
    }
    std::cout << " = " << r_xa << std::endl;
    return r_xa;
}

SymRnsNumber SymRnsNumber::round(Positional_Int b) const {
    bool sign;
    Positional_Int ra = get_rank(&sign);  // получение ранга и знака текущего числа
    SymRnsNumber aq{*this}; // копия текущего числа
    if (sign) aq -= SymRnsNumber{3, base}; // скорректируем само число
    std::cout << "orig=" << to_positional_int0() << ", corrected=" << aq.to_positional_int0() << std::endl;
    for (size_t i = 0; i < a.size(); ++i) {
        Module pi = base.get().p[i];
        // поразрядно делим копию числа на x
        aq.a[i] = mod_sym(aq.a[i] * mod_inverse_sym(b, pi), pi);
    }
    std::cout << "aq=" << aq.to_positional() << ' ' << aq << " ra=" << aq.get_rank() << std::endl;

    Positional_Int r_xa = aq.get_rank_x(b);  // ранг числа (aq * x)
    // k - разница в рангах числа a и 10*(a/10)
    Positional_Int k = ra - r_xa;
    std::cout << "k=" << k << " ra=" << ra << " r_xa=" << r_xa << std::endl; 
    SymRnsNumber res{*this};
    if (k) { // если есть разница в рангах
        // t - последняя цифра делимого, округляющее число
        Positional_Int t;
        t = mod(-k * base.get().P, b);
        std::cout << "t=" << t << std::endl;
        res -= SymRnsNumber{t, base}; // отнимаем округляющее число
    }
    std::cout << "rounded=" << res.to_positional_int0() << std::endl;
    // if (sign) {
    //     res += SymRnsNumber{3, base};
    // }
    return res;
}

int floor_int_to_10(int x) {
    if (x >= 0)
        return (x / 10) * 10;
    else
        return ((x - 9) / 10) * 10;  // Эквивалент floor(x / 10) * 10
}

int main() {
    SymRnsBase base{{3, 5, 7}};

    // for (int x=-10; x <=10; ++x) {
    //     std::cout << "x=" << x << " x^-1=" << SymRnsNumber::mod_inverse_sym(x, 5) << std::endl;
    // }

    // примеры округления
    //SymRnsNumber a{{2, 3, 8}, base};
    // SymRnsNumber a{101, base};
    // std::cout << a.to_positional() << ' ' << a << " rank=" << a.get_rank() << std::endl;

    // SymRnsNumber b = a.round(10);
    // std::cout << b.to_positional() << ' ' << b <<  std::endl;

    // тесты
    // Positional_Int q = 1;  // делитель для fixed point тестов
    bool ok = true;

    // // // проверки целочисленного округления
    // Positional_Int rnd = 10;
    // std::vector<Positional_Int> xx = {-11, -10, 10, 11};
    // for (Positional_Int x_int : xx) {
    // // for (Positional_Int x_int = 0; x_int <= base.P / 2; ++x_int) {
    // // for (Positional_Int x_int = -3; x_int <= 3; ++x_int) {
    // // for (Positional_Int x_int = -11; x_int <= -10; x_int += 1) {
    //     SymRnsNumber x_rns = SymRnsNumber{x_int, base};
    //     SymRnsNumber x_rns_round = x_rns.round(rnd);
    //     Positional_Int x_rns_round_int = x_rns_round.to_positional_int();
    //     Positional_Int x_round_int = floor_int_to_10(x_int);
    //     if (x_round_int != x_rns_round_int) {
    //         std::cout << "Error: original " << x_int << x_rns << "->" << x_round_int <<  " != converted " << x_rns_round_int << std::endl << std::endl;
    //         ok = false;
    //     } else {
    //         std::cout << "Info: original " << x_int << x_rns << "->" << x_round_int <<  " == converted " << x_rns_round_int << std::endl << std::endl;
    //     }
    // }

    // целочисленные проверки конвертации в позиционную ИС
    for (Positional_Int x_int = -base.P / 2; x_int <= base.P / 2; ++x_int) {
        SymRnsNumber x_rns = SymRnsNumber{x_int, base};
        Positional_Int x_pos_int = x_rns.to_positional_int_MRC();
        if (x_int != x_pos_int) {
            std::cout << "Error: x: " << x_int << " rns: " << x_rns << " != x_pos " << x_pos_int << std::endl;
            ok = false;
        } else {
            std::cout << "Info: x: " << x_int << " rns: " << x_rns << " x_pos: " << x_pos_int << std::endl;
        }
    }

    // // fixed point проверки конвертации в позиционную ИС
    // for (Positional_Int x_int = -base.P / 2; x_int <= base.P / 2; ++x_int) {
    //     Positional_Float x = static_cast<Positional_Float>(x_int) / q;
    //     SymRnsNumber x_rns = SymRnsNumber{x, base, q};
    //     Positional_Int x_pos_int;
    //     Positional_Float x_pos = x_rns.to_positional(&x_pos_int);
    //     // std::cout << "int: " << x_int << " calc: " << x << " " << x_rns << " conv: " << x_pos
    //     //           << " (" << x_pos_int << "/" << q << ")" << std::endl;

    //     if (fabs(x - x_pos) > 0.001) {
    //         std::cout << "Error: original " << x << " != converted " << x_pos << std::endl;
    //         ok = false;
    //     }
    // }

    // for (Positional_Int x_int = -base.P/2; x_int <= base.P/2; ++x_int) {

    //     // проверки унарного минуса
    //     Positional_Float res = -static_cast<Positional_Float>(x_int) / q;
    //     Positional_Int res_int = -x_int;
    //     if (res_int >= -base.P/2 && res_int <= base.P/2) {
    //         SymRnsNumber r_rns = -SymRnsNumber(x_int, base, q);
    //         SymRnsNumber res_rns = SymRnsNumber(res_int, base, q);
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
    //             SymRnsNumber r_rns = SymRnsNumber(x_int, base, q) + SymRnsNumber(y_int, base,
    //             q); SymRnsNumber res_rns = SymRnsNumber(res, base, q); if (res_rns != r_rns) {
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
    //             SymRnsNumber r_rns = SymRnsNumber(x_int, base, q) * SymRnsNumber(y_int, base,
    //             q); SymRnsNumber res_rns = SymRnsNumber(res, base, q*q); if (res_rns != r_rns) {
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
    // SymRnsNumber a{5, base, 6};  // = 5/6
    // SymRnsNumber b{-8, base, 9};  // = -8/9

    // // SymRnsNumber a{0.82f, base, 6};   // ~= 5/6
    // // SymRnsNumber b{-0.86f, base, 9};  // ~= -8/9

    // Positional_Int a_int, b_int, c_int;
    // Positional_Float a_fl = a.to_positional(&a_int);
    // Positional_Float b_fl = b.to_positional(&b_int);

    // SymRnsNumber c = a + b;
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

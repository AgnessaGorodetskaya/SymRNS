#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

typedef int32_t Positional_Int;       // позиционное целочисленное
typedef float Positional_Float;       // позиционное вещественное
typedef int16_t Module;               // целочисленный модуль разряда ССОК
typedef std::vector<Module> Modules;  // вектор модулей разрядов ССОК

// Функция возвращает обратное по модулю число a^-1 mod p
Positional_Int mod_inverse(Positional_Int a, Positional_Int p) {
    Positional_Int p0 = p, t, q;
    Positional_Int x0 = 0, x1 = 1;
    if (p == 1) return 0;
    // Алгоритм Евклида
    while (a > 1) {
        q = a / p;        // Частное
        t = p;
        p = a % p, a = t;
        t = x0;
        x0 = x1 - q * x0, x1 = t; // Обновляем x0 и x1
    }
    if (x1 < 0) x1 += p0; // Обратное может быть отрицательным
    return x1;
}

// класс оснований СОК
class RnsBase {
   public:
    Modules p;  // вектор оснований (p1, p2 .. pn)
    std::vector<Positional_Int> B;  // ортогональные базисы, соответствующие (1; 0; 0..),
                                    // (0; 1; 0...), (0; 0; 1..) для конвертации в позиционную СИ
    std::vector<Positional_Int> m;  // веса ортогональных базисов
    Positional_Int P;   // кол-во чисел в представлении в данном наборе оснований P = p1 * ... * pn,
                        // ex. 3 * 5 * 7 = 105
    RnsBase(const Modules& p0);  // конструктор из вектора оснований
};

// класс представления числа в СОК
class RnsNumber {
   public:
    std::reference_wrapper<const RnsBase> base;  // ссылка на класс оснований
    Modules a;                                   // вектор значений модулей (a1, a2 .. an)
    Positional_Int q_int;  // делитель представления с фиксированной точкой

    RnsNumber(Positional_Int x, const RnsBase& base0, Positional_Int q0 = 1);  // из целого позиционного
    RnsNumber(Positional_Float x, const RnsBase& base0, Positional_Int q0);    // из вещест. позиционного
    template <typename ModType> ModType mod(Positional_Int x, ModType p) const; // симметричный модуль
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
    Positional_Float to_positional(Positional_Int* x_int = NULL) const;  // в позиционное представление
    RnsNumber round(Positional_Int x_int);  // округление к ближайшему кратному x_int
};

RnsBase::RnsBase(const Modules& p0) : p{p0}, B(p.size()), m(p.size()) {
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
        m[i] = mod_inverse(P / p[i], p[i]);  // вес i-го ортогонального базиса
        B[i] = m[i] * P / p[i];  // i-й ортогональный базис
        RnsNumber B_rns{B[i], *this, 1};
        std::cout << "Ортогональный базис B" << i+1 << '=' << B[i] << ", вес базиса m" << i+1 << "=" << m[i] << ' ' << B_rns << std::endl;
    }

    // // поиск ортогональных базисов во всем диапазоне
    // for (Positional_Int x_int = 0; x_int <= P; ++x_int) {
    //     RnsNumber t{x_int, *this, 1};
    //     bool is_ref = false;
    //     size_t j = 0;  // позиция единичного значения
    //     for (size_t i = 0; i < t.a.size(); ++i) {
    //         if (t.a[i] == 1) {   // рассматриваемый модуль единичный (=1)
    //             if (!is_ref) {      // число еще не базис, значит единичный модуль первый
    //                 is_ref = true;  // предпологаем, что число базис
    //                 j = i;          // сохраняем позицию единичного модуля
    //             } else if (is_ref) {  // уже считали, что число базис, значит это не первый
    //                                   // единичный модуль
    //                 is_ref = false;  // и следовательно число не базис
    //                 break;           // нет смысла продолжать дальше
    //             }
    //         } else if (t.a[i] != 0) {  // наличие любого не единичного и не нулевого модуля
    //             is_ref = false;           // значит число - не базис
    //             break;                    // нет смысла продолжать дальше
    //         }
    //     }
    //     if (is_ref) {       // если в итоге число базис
    //         B[j] = x_int;   // сохраняем его в массиве ортогональных базисов по индексу,
    //                         // соответствующему позиции единичного модуля
    //         m[j] = B[j] * p[j] / P;  // вес j-го базиса
    //         std::cout << 'B' << j+1 << '=' << x_int << " m" << j+1 << "=" << m[j] << ' ' << t << std::endl;
    //     }
    // }
}

Positional_Float RnsNumber::to_positional(Positional_Int* x_int) const {
    Positional_Int res_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        res_int += a[i] * base.get().B[i];
    }
    res_int = mod(res_int, base.get().P);  // по симметричному модулю максимального представления
    if (x_int) *x_int = res_int;  // опционально целочисленный числитель возвращаем по указателю в параметре
    return static_cast<Positional_Float>(res_int) / q_int;
}

RnsNumber::RnsNumber(Positional_Int x, const RnsBase& base0, Positional_Int q0)
    : base{base0}, a(base.get().p.size()), q_int(q0) {
    size_t i = 0;
    for (Module pi : base.get().p) {
        a[i++] = mod(x, pi);
    }
}

RnsNumber::RnsNumber(Positional_Float x, const RnsBase& base0, Positional_Int q0)
    : RnsNumber(static_cast<Positional_Int>(std::round(x * q0)), base0, q0) {}

template <typename ModType>
ModType RnsNumber::mod(Positional_Int x, ModType p) const {
    ModType r = x % p;
    if (r < 0) r += p;
    return r;
}

RnsNumber& RnsNumber::operator+=(const RnsNumber& y) {
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
        a[i] = mod(xm * a[i] + ym * y.a[i], base.get().p[i]);
    }
    return *this;
}

RnsNumber RnsNumber::operator+(const RnsNumber& y) const {
    RnsNumber res{*this};
    res += y;
    return res;
}

RnsNumber& RnsNumber::operator-=(const RnsNumber& y) {
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
        a[i] = mod(xm * a[i] - ym * y.a[i], base.get().p[i]);
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
    q_int *= y.q_int;  // при умножении fixed-делители умножаются 2/3 * 4/5 = (2*4)/(3*5)
    return *this;
}

RnsNumber RnsNumber::operator*(const RnsNumber& y) const {
    RnsNumber res{*this};
    res *= y;
    return res;
}

bool RnsNumber::operator==(const RnsNumber& y) const {
    if (a.size() != y.a.size() || q_int != y.q_int) return false;
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
    if (y.q_int != 1) os << ";q=" << y.q_int;
    os << ')';
    return os;
}

RnsNumber RnsNumber::round(Positional_Int x_int) {
    if (x_int >= base.get().p[0]) {
        RnsNumber res{*this};
        Positional_Int half_x_int = x_int / 2;
        RnsNumber half_x_rns{half_x_int, base};
        res += half_x_rns;
        for (size_t i = 0; i < a.size(); ++i) {
            if (x_int % base.get().p[i] == 0) {
                x_int /= base.get().p[i];
                res.a[i] = 0;
                if (x_int == 1) break;
            }
        }
        if (x_int == 1) return res;
    }
    return RnsNumber{*this};  // не получилось округлить
}

int main() {
    RnsBase base{{3, 5, 7, 17}};

    // тесты
    Positional_Int q = 1;  // делитель для fixed point тестов
    bool ok = true;

    // проверки обратной конвертации в позиционную ИС
    for (Positional_Int x_int = 0; x_int < base.P; ++x_int) {
        Positional_Float x = static_cast<Positional_Float>(x_int) / q;
        RnsNumber x_rns = RnsNumber{x, base, q};
        Positional_Int x_pos_int;
        Positional_Float x_pos = x_rns.to_positional(&x_pos_int);
        // std::cout << "int: " << x_int << " calc: " << x << " " << x_rns << " conv: " << x_pos
        //           << " (" << x_pos_int << "/" << q << ")" << std::endl;

        if (fabs(x - x_pos) > 0.001) {
            std::cout << "Error: original " << x << " != converted " << x_pos << std::endl;
            ok = false;
        }
    }

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

    // примеры округления
    RnsNumber a{43, base};
    Positional_Int a_int;
    a.to_positional(&a_int);
    std::cout << a_int << ' ' << a << std::endl;

    RnsNumber a_round = a.round(15);
    Positional_Int a_round_int;
    a_round.to_positional(&a_round_int);
    std::cout << a_round_int << ' ' << a_round << std::endl;
}

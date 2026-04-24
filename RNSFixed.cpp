#include "RNSFixed.hpp"

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
    // Ms = RnsFixed::mod(P, S);
    // std::cout << " Ms=" << Ms;
    std::cout << std::endl;

    for (size_t i = 0; i < p.size(); ++i) {
        Pi[i] = P / p[i];
        m[i] = RnsFixed::mod_inverse(Pi[i], p[i]);  // вес i-го ортогонального базиса
        B[i] = m[i] * Pi[i];  // i-й ортогональный базис
        // C[i] = RnsFixed::mod(Pi, S);
        RnsFixed B_rns{B[i], *this};
        std::cout << "ORT: Ортогональный базис B" << i+1 << '=' << B[i] << ' ' << B_rns << ", вес базиса m" << i+1 << "=" << m[i]
                  << " P" << i+1 << '=' << Pi[i];
        // std::cout << " C[i]=" << C[i];
        std::cout << std::endl;
    }

    for (size_t i = 0; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            t[i][j] = RnsFixed::mod_inverse(p[i], p[j]);
            std::cout << "MRC: t[" << i+1 << "][" << j+1 <<"] = " << t[i][j] << std::endl;
        }
    }
}

Positional_Float RnsFixed::to_positional_crt() const {
    Positional_Int res_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        res_int += a[i] * base.get().B[i];
    }
    res_int = mod(res_int, base.get().P);  // по модулю максимального представления
    return (static_cast<Positional_Float>(res_int) / base.get().S);
}

Positional_Float RnsFixed::frac_crt_sum() const {
    Positional_Float pos_float = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        pos_float += static_cast<Positional_Float>(mod(a[i] * base.get().m[i], base.get().p[i])) / base.get().p[i];
        // res += static_cast<Positional_Float>(a[i] * base.get().m[i]) / base.get().p[i];
    }
    return pos_float;
}

Positional_Float RnsFixed::to_positional_frac_crt_unscaled() const {
    Positional_Float pos_float = frac_crt_sum(), int_part;
    pos_float = std::modf(pos_float, &int_part);
    pos_float *= base.get().P;
    // std::cout << ' ' << pos_float << std::endl;
    return pos_float;
}

Positional_Float RnsFixed::to_positional_frac_crt() const {
    Positional_Float pos_float = to_positional_frac_crt_unscaled();
    pos_float /= base.get().S;
    return pos_float;
}

Positional_Int RnsFixed::get_rank() const {
    Positional_Float pos_float = frac_crt_sum();
    Positional_Int rank = static_cast<Positional_Int>(pos_float);
    // std::cout << " Rank=" << rank << " Pos_Float=" << pos_float << std::endl;
    return rank;
}

Positional_Int RnsFixed::get_remainder(Positional_Int divisor) const {
    Positional_Int rem_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        rem_int += mod(base.get().Pi[i] * mod(a[i] * base.get().m[i], base.get().p[i]), divisor);
    }
    Positional_Int rank = get_rank();
    rem_int -= mod(rank * base.get().P, divisor);
    rem_int = mod(rem_int, divisor);
    return rem_int;
}

Positional_Float RnsFixed::to_positional_mrc() const {
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

RnsFixed::RnsFixed(Positional_Int x, const RnsBase& base0)
    : base{base0}, a(base.get().p.size()) {
    size_t i = 0;
    for (Module pi : base.get().p) {
        a[i++] = mod(x, pi);
    }
}

RnsFixed::RnsFixed(const Modules& a0, const RnsBase& base0)
    : base{base0}, a{a0} {}

template<typename ModType> ModType RnsFixed::mod(Positional_Int x, ModType p) {
    ModType r = x % p;
    if (r < 0) r += p;
    return r;
}

template<typename ModType> ModType RnsFixed::mod_inverse(Positional_Int x, ModType p) {
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

RnsFixed& RnsFixed::operator+=(const RnsFixed& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod(a[i] + y.a[i], base.get().p[i]);
    }
    return *this;
}

RnsFixed RnsFixed::operator+(const RnsFixed& y) const {
    RnsFixed res{*this};
    res += y;
    return res;
}

RnsFixed& RnsFixed::operator-=(const RnsFixed& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod(a[i] - y.a[i], base.get().p[i]);
    }
    return *this;
}

RnsFixed RnsFixed::operator-(const RnsFixed& y) const {
    RnsFixed res{*this};
    res -= y;
    return res;
}

RnsFixed RnsFixed::operator-() const {
    RnsFixed res{*this};
    for (size_t i = 0; i < a.size(); ++i) {
        res.a[i] = -a[i];
    }
    return res;
}

void RnsFixed::div_int(const RnsFixed& y) {  // деление нацело (когда точно делится)
    for (size_t i = 0; i < a.size(); ++i) {
        Positional_Int Si = RnsFixed::mod_inverse(y.a[i], base.get().p[i]);
        // std::cout << ' ' << Si;
        a[i] = mod(a[i] * Si, base.get().p[i]);
    }
    // std::cout << std::endl;
}

bool RnsFixed::operator<(const RnsFixed &y) const
{
    if (*this == y) return false;
    Positional_Float x_fl = frac_crt_sum(), y_fl = y.frac_crt_sum(), int_part;
    x_fl = std::modf(x_fl, &int_part);
    y_fl = std::modf(y_fl, &int_part);
    return (x_fl < y_fl);
}

RnsFixed& RnsFixed::operator*=(const RnsFixed& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod(a[i] * y.a[i], base.get().p[i]);
    }

    // std::cout << "RAW=" << to_positional_frac_crt_unscaled();
    // масштабирование на S
    Positional_Int remainder = get_remainder(base.get().S);
    // std::cout << " REM=" << remainder << std::endl;

    if (remainder < base.get().S / 2) {  // округление к ближайшему
        RnsFixed remainder_rns{remainder, base};
        *this -= remainder_rns;
    } else {
        if (remainder == base.get().S / 2) {
            std::cout << "!!!! СОК: Остаток == S/2" << std::endl;
        }
        RnsFixed remainder_rns{base.get().S - remainder, base};
        *this += remainder_rns;
    }
    // std::cout << " WO_REM=" << to_positional_frac_crt_unscaled();

    RnsFixed S_rns{base.get().S, base};
    // std::cout << " S=" << S_rns;
    div_int(S_rns);

    return *this;
}

RnsFixed& RnsFixed::operator/=(const RnsFixed& y) {
    RnsFixed S_rns{base.get().S, base};
    for (size_t i = 0; i < a.size(); ++i) { // a * S
        a[i] = mod(a[i] * S_rns.a[i], base.get().p[i]);
    }

    // std::cout << "RAW=" << to_positional_frac_crt_unscaled();
    // масштабирование на y
    Positional_Int y_int = static_cast<Positional_Int>(std::round(y.to_positional_frac_crt_unscaled()));
    Positional_Int remainder = get_remainder(y_int);
    // std::cout << " REM=" << remainder << std::endl;

    if (remainder < y_int / 2) {  // округление к ближайшему
        RnsFixed remainder_rns{remainder, base};
        *this -= remainder_rns;
    } else {
        if (remainder == y_int / 2) {
            std::cout << "!!!! СОК: Остаток == divider/2" << std::endl;
        }
        RnsFixed remainder_rns{y_int - remainder, base};
        *this += remainder_rns;
    }
    // std::cout << " WO_REM=" << to_positional_frac_crt_unscaled();
    div_int(y);

    return *this;
}

RnsFixed RnsFixed::operator*(const RnsFixed& y) const {
    RnsFixed res{*this};
    res *= y;
    return res;
}

RnsFixed RnsFixed::operator/(const RnsFixed& y) const {
    RnsFixed res{*this};
    res /= y;
    return res;
}

bool RnsFixed::operator==(const RnsFixed& y) const {
    if (a.size() != y.a.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != y.a[i]) return false;
    }
    return true;
}

bool RnsFixed::operator!=(const RnsFixed& y) const { return (!(*this == y)); }

std::ostream& operator<<(std::ostream& os, const RnsFixed& y) {
    size_t sz = y.a.size();
    os << '(';
    for (size_t i = 0; i < sz; ++i) {
        if (i) os << ";";
        os << y.a[i];
    }
    os << ')';
    return os;
}

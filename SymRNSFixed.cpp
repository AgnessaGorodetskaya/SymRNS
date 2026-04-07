#include "SymRNSFixed.hpp"

SymRnsBase::SymRnsBase(const Modules& p0, Positional_Int S0) :
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
    std::cout << ") P=" << P << std::endl;

    for (size_t i = 0; i < p.size(); ++i) {
        Pi[i] = P / p[i];
        m[i] = SymRnsFixed::mod_inverse_sym(Pi[i], p[i]);  // вес i-го ортогонального базиса
        B[i] = m[i] * Pi[i];  // i-й ортогональный базис
        SymRnsFixed B_rns{B[i], *this};
        std::cout << "ORT: Ортогональный базис B" << i+1 << '=' << B[i] << ' ' << B_rns << ", вес базиса m" << i+1 << "=" << m[i] << std::endl;
    }

    for (size_t i = 0; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            t[i][j] = SymRnsFixed::mod_inverse_sym(p[i], p[j]);
            std::cout << "MRC: t[" << i+1 << "][" << j+1 <<"] = " << t[i][j] << std::endl;
        }
    }
}

Positional_Float SymRnsFixed::to_positional_crt() const {
    Positional_Int sm = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sm += a[i] * base.get().B[i];
    }
    Positional_Int A = mod_sym(sm, base.get().P);
    // Module rank = static_cast<Module>((sm - A) / base.get().P);
    // std::cout << std::endl << "ORT: " << sm << " mod_sym " << base.get().P << " = " << A << ", ранг=" << rank <<std::endl;
    return (static_cast<Positional_Float>(A) / base.get().S);  // по симметричному модулю динамического диапазона
}

Positional_Float SymRnsFixed::frac_crt_sum() const {
    Positional_Float pos_float = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        pos_float += static_cast<Positional_Float>(mod_sym(a[i] * base.get().m[i], base.get().p[i])) / base.get().p[i];
        // pos_float += static_cast<Positional_Float>(a[i] * base.get().m[i]) / base.get().p[i];
    }
    return pos_float;
}

Positional_Int SymRnsFixed::to_positional_frac_crt_unscaled() const {
    Positional_Float pos_float = frac_crt_sum(), int_part;
    pos_float = std::modf(pos_float, &int_part);
    if (pos_float > 0.5) { pos_float -= 1.0; }
    else if (pos_float < -0.5) { pos_float += 1.0; }
    pos_float *= base.get().P;
    Positional_Int pos_int = static_cast<Positional_Int>(std::round(pos_float));
    return pos_int;
}

Positional_Float SymRnsFixed::to_positional_frac_crt() const {
    Positional_Float pos_float = frac_crt_sum(), int_part;
    pos_float = std::modf(pos_float, &int_part);
    if (pos_float > 0.5) { pos_float -= 1.0; }
    else if (pos_float < -0.5) { pos_float += 1.0; }
    // std::cout << pos_float << std::endl;
    pos_float *= base.get().P;
    pos_float /= base.get().S;
    return pos_float;
}

Positional_Int SymRnsFixed::get_rank() const {
    Positional_Float pos_float = frac_crt_sum(), int_part;
    std::cout << " get_rank(): " << pos_float;
    pos_float = std::modf(pos_float, &int_part);
    Positional_Int rank = static_cast<Positional_Int>(int_part);
    if (pos_float > 0.5) { pos_float -= 1.0; ++rank; }
    else if (pos_float < -0.5) { pos_float += 1.0; --rank; }
    std::cout << " Rank=" << rank << " Pos_Float=" << pos_float << std::endl;
    return rank;
}

Positional_Int SymRnsFixed::get_remainder(Positional_Int divisor) const {
    Positional_Int rem_int = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        rem_int += mod_sym(base.get().Pi[i] * mod_sym(a[i] * base.get().m[i], base.get().p[i]), divisor);
    }
    Positional_Int rank = get_rank();
    rem_int -= mod_sym(rank * base.get().P, divisor);
    rem_int = mod_sym(rem_int, divisor);
    return rem_int;
}

Positional_Float SymRnsFixed::to_positional_mrc() const {
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
        // std::cout << " = " << xi << ") mod_sym " << base.get().p[i];
        x[i] = mod_sym(xi, base.get().p[i]);
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

SymRnsFixed::SymRnsFixed(Positional_Int x, const SymRnsBase& base0)
    : base{base0}, a(base.get().p.size()) {
    size_t i = 0;
    for (Module pi : base.get().p) {
        a[i++] = mod_sym(x, pi);
    }
}

template<typename ModType> ModType SymRnsFixed::mod(Positional_Int x, ModType p) {
    ModType r = x % p;
    if (r < 0) r += p;
    return r;
}

template<typename ModType> ModType SymRnsFixed::mod_sym(Positional_Int x, ModType p) {
    ModType p_half = p / 2;
    ModType r = x % p;
    if (r > p_half) { r -= p; }
    else if (r < -p_half) { r += p; }
    return r;
}

template<typename ModType> ModType SymRnsFixed::mod_inverse_sym(Positional_Int x, ModType p) {
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

SymRnsFixed& SymRnsFixed::operator+=(const SymRnsFixed& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(a[i] + y.a[i], base.get().p[i]);
    }
    return *this;
}

SymRnsFixed SymRnsFixed::operator+(const SymRnsFixed& y) const {
    SymRnsFixed res{*this};
    res += y;
    return res;
}

SymRnsFixed& SymRnsFixed::operator-=(const SymRnsFixed& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(a[i] - y.a[i], base.get().p[i]);
    }
    return *this;
}

SymRnsFixed SymRnsFixed::operator-(const SymRnsFixed& y) const {
    SymRnsFixed res{*this};
    res -= y;
    return res;
}

SymRnsFixed SymRnsFixed::operator-() const {
    SymRnsFixed res{*this};
    for (size_t i = 0; i < a.size(); ++i) {
        res.a[i] = -a[i];
    }
    return res;
}

void SymRnsFixed::div_int(const SymRnsFixed& y) {  // деление нацело (когда точно делится)
    for (size_t i = 0; i < a.size(); ++i) {
        Positional_Int Si = mod_inverse_sym(y.a[i], base.get().p[i]);
        std::cout << ' ' << Si;
        a[i] = mod(a[i] * Si, base.get().p[i]);
    }
    std::cout << std::endl;
}

SymRnsFixed& SymRnsFixed::operator*=(const SymRnsFixed& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(a[i] * y.a[i], base.get().p[i]);
    }

    std::cout << "RAW=" << to_positional_crt();
    // масштабирование на S
    Positional_Int remainder = get_remainder(base.get().S);
    std::cout << " REM=" << remainder;

    SymRnsFixed remainder_rns{remainder, base};
    *this -= remainder_rns;

    std::cout << " WO_REM=" << to_positional_crt() << ' ' << *this;

    SymRnsFixed S_rns{base.get().S, base};
    std::cout << " S=" << S_rns;
    div_int(S_rns);

    std::cout << std::endl;

    return *this;
}

SymRnsFixed SymRnsFixed::operator*(const SymRnsFixed& y) const {
    SymRnsFixed res{*this};
    res *= y;
    return res;
}

bool SymRnsFixed::operator==(const SymRnsFixed& y) const {
    if (a.size() != y.a.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != y.a[i]) return false;
    }
    return true;
}

bool SymRnsFixed::operator!=(const SymRnsFixed& y) const { return (!(*this == y)); }

bool SymRnsFixed::operator<(const SymRnsFixed &y) const
{
    if (*this == y) return false;  // быстро

    Positional_Float x_fl = frac_crt_sum(), y_fl = y.frac_crt_sum(), int_part;

    x_fl = std::modf(x_fl, &int_part);
    if (x_fl > 0.5) { x_fl -= 1.0; }
    else if (x_fl < -0.5) { x_fl += 1.0; }

    y_fl = std::modf(y_fl, &int_part);
    if (y_fl > 0.5) { y_fl -= 1.0; }
    else if (y_fl < -0.5) { y_fl += 1.0; }

    return (x_fl < y_fl);
}

std::ostream& operator<<(std::ostream& os, const SymRnsFixed& y) {
    size_t sz = y.a.size();
    os << '(';
    for (size_t i = 0; i < sz; ++i) {
        if (i) os << ",";
        os << y.a[i];
    }
    os << ')';
    return os;
}

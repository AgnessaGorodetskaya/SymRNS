#include "SymRNS.hpp"

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

    for (size_t i = 0; i < p.size(); ++i) {
        Positional_Int Pi = P / p[i];
        m[i] = SymRnsNumber::mod_inverse_sym(Pi, p[i]);  // вес i-го ортогонального базиса
        B[i] = m[i] * Pi;  // i-й ортогональный базис
        SymRnsNumber B_rns{B[i], *this};
        std::cout << "ORT: Ортогональный базис B" << i+1 << '=' << B[i] << ' ' << B_rns << ", вес базиса m" << i+1 << "=" << m[i] << std::endl;
    }

    for (size_t i = 0; i < p.size(); ++i) {
        for (size_t j = i + 1; j < p.size(); ++j) {
            t[i][j] = SymRnsNumber::mod_inverse_sym(p[i], p[j]);
            std::cout << "MRC: t[" << i+1 << "][" << j+1 <<"] = " << t[i][j] << std::endl;
        }
    }
}

Positional_Int SymRnsNumber::to_positional_ort() const {
    Positional_Int sm = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        sm += a[i] * base.get().B[i];
    }
    Positional_Int A = mod_sym(sm, base.get().P);
    // Module rank = static_cast<Module>((sm - A) / base.get().P);
    // std::cout << std::endl << "ORT: " << sm << " mod_sym " << base.get().P << " = " << A << ", ранг=" << rank <<std::endl;
    return A;  // по симметричному модулю динамического диапазона
}

Positional_Int SymRnsNumber::to_positional_mrc() const {
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
    std::cout << Ai << ';';
    for (size_t i = 0; i < a.size(); ++i) std::cout << x[i] << ';';
    std::cout << std::endl;
    return Ai;
}

SymRnsNumber::SymRnsNumber(Positional_Int x, const SymRnsBase& base0)
    : base{base0}, a(base.get().p.size()) {
    size_t i = 0;
    for (Module pi : base.get().p) {
        a[i++] = mod_sym(x, pi);
    }
}

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
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(a[i] + y.a[i], base.get().p[i]);
    }
    return *this;
}

SymRnsNumber SymRnsNumber::operator+(const SymRnsNumber& y) const {
    SymRnsNumber res{*this};
    res += y;
    return res;
}

SymRnsNumber& SymRnsNumber::operator-=(const SymRnsNumber& y) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] = mod_sym(a[i] - y.a[i], base.get().p[i]);
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
    return *this;
}

SymRnsNumber SymRnsNumber::operator*(const SymRnsNumber& y) const {
    SymRnsNumber res{*this};
    res *= y;
    return res;
}

bool SymRnsNumber::operator==(const SymRnsNumber& y) const {
    if (a.size() != y.a.size()) return false;
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
        if (i) os << ",";
        os << y.a[i];
    }
    os << ')';
    return os;
}

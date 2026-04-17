#include "RNS.hpp"

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

Positional_Int RnsNumber::to_positional_crt() const {
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

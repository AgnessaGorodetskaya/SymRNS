#include "SymRNSFrac.hpp"

SymRnsFrac& SymRnsFrac::operator*=(const SymRnsFrac& y) {
    numerator *= y.numerator;
    denominator *= y.denominator;
    return *this;
}

SymRnsFrac SymRnsFrac::operator*(const SymRnsFrac& y) const {
    SymRnsFrac res{*this};
    res *= y;
    return res;
}

SymRnsFrac& SymRnsFrac::operator/=(const SymRnsFrac& y) {
    numerator *= y.denominator;
    denominator *= y.numerator;
    return *this;
}

SymRnsFrac SymRnsFrac::operator/(const SymRnsFrac& y) const {
    SymRnsFrac res{*this};
    res /= y;
    return res;
}

SymRnsFrac& SymRnsFrac::operator+=(const SymRnsFrac& y) {
    if (denominator == y.denominator) {
        numerator += y.numerator;
    } else {
        numerator = numerator * y.denominator + y.numerator * denominator;
        denominator *= y.denominator;
    }
    return *this;
}

SymRnsFrac SymRnsFrac::operator+(const SymRnsFrac& y) const {
    SymRnsFrac res{*this};
    res += y;
    return res;
}

SymRnsFrac& SymRnsFrac::operator-=(const SymRnsFrac& y) {
    if (denominator == y.denominator) {
        numerator -= y.numerator;
    } else {
        numerator = numerator * y.denominator - y.numerator * denominator;
        denominator *= y.denominator;
    }
    return *this;
}

SymRnsFrac SymRnsFrac::operator-(const SymRnsFrac& y) const {
    SymRnsFrac res{*this};
    res -= y;
    return res;
}

std::ostream& operator<<(std::ostream& os, const SymRnsFrac& y) {
    os << y.numerator << '/' << y.denominator;
    return os;
}

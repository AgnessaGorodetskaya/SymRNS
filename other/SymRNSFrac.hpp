#include <iostream>
#include "SymRNS.cpp"

class SymRnsFrac {
   public:
    SymRnsNumber numerator, denominator;  // числитель и знаменатель
    SymRnsFrac(SymRnsNumber num, SymRnsNumber denom) : numerator{num}, denominator{denom} {};
    SymRnsFrac& operator*=(const SymRnsFrac& y);
    SymRnsFrac operator*(const SymRnsFrac& y) const;
    SymRnsFrac& operator/=(const SymRnsFrac& y);
    SymRnsFrac operator/(const SymRnsFrac& y) const;
    SymRnsFrac& operator+=(const SymRnsFrac& y);
    SymRnsFrac operator+(const SymRnsFrac& y) const;
    SymRnsFrac& operator-=(const SymRnsFrac& y);
    SymRnsFrac operator-(const SymRnsFrac& y) const;
    friend std::ostream& operator<<(std::ostream& os, const SymRnsFrac& y);
};

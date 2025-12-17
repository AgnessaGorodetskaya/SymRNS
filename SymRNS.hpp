#include <numeric>
#include <functional>
#include <vector>
#include <iostream>

typedef int32_t Positional_Int;       // позиционное целочисленное
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

    SymRnsNumber(Positional_Int x, const SymRnsBase& base0);  // из целого позиционного
    SymRnsNumber(const Modules& a0, const SymRnsBase& base0) : base{base0}, a{a0} {}; // из вектора остатков
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
    Positional_Int to_positional_ort() const;  // в позиционное (int) представление по ORT
    Positional_Int to_positional_mrc() const;  // в позиционное представление по MRC
};

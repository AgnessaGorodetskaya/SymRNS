#include <numeric>
#include <functional>
#include <vector>
#include <iostream>

typedef int32_t Positional_Int;       // позиционное целочисленное
typedef double Positional_Float;      // позиционное вещественное
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
    std::vector<Positional_Int> Pi;   // P / p[i]
    Positional_Int S;  // коэф масштабирования
    SymRnsBase(const Modules& p0, Positional_Int S0);  // конструктор из вектора оснований
};

// класс представления числа в ССОК
class SymRnsFixed {
   public:
    std::reference_wrapper<const SymRnsBase> base;  // ссылка на класс оснований
    Modules a;                                      // вектор значений модулей (a1, a2 .. an)

    SymRnsFixed(Positional_Int x, const SymRnsBase& base0);  // из целого позиционного
    SymRnsFixed(const Modules& a0, const SymRnsBase& base0) : base{base0}, a{a0} {}; // из вектора остатков
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
    SymRnsFixed& operator+=(const SymRnsFixed& y);  // this += y
    SymRnsFixed operator+(const SymRnsFixed& y) const;  // z = this + y
    SymRnsFixed& operator-=(const SymRnsFixed& y);  // this -= y
    SymRnsFixed operator-(const SymRnsFixed& y) const;  // z = this - y
    SymRnsFixed operator-() const;  // Унарный минус z = -this
    SymRnsFixed& operator*=(const SymRnsFixed& y);  // this *= y
    SymRnsFixed operator*(const SymRnsFixed& y) const;  // z = this * y
    SymRnsFixed& operator/=(const SymRnsFixed& y);  // this /= y
    SymRnsFixed operator/(const SymRnsFixed& y) const;  // z = this / y
    bool operator==(const SymRnsFixed& y) const;  // (this == y)
    bool operator!=(const SymRnsFixed& y) const;  // (this != y)
    bool operator<(const SymRnsFixed& y) const; // сравнение
    void div_int(const SymRnsFixed& y);  // деление нацело (когда точно делится)
    friend std::ostream& operator<<(std::ostream& os, const SymRnsFixed& y);  // для вывода на cout
    Positional_Float frac_crt_sum() const; // дробная сумма дробной CRT
    Positional_Int to_positional_frac_crt_unscaled() const;  // в позиционное представление (немасштабированное целое)
    Positional_Float to_positional_frac_crt() const;  // в позиционное представление через дробную CRT
    Positional_Float to_positional_crt() const;  // в позиционное представление CRT
    Positional_Float to_positional_mrc() const;  // в позиционное представление MRC
    Positional_Int get_rank() const; // получение ранга числа
    Positional_Int get_remainder(Positional_Int divisor) const; // остаток от деления
};

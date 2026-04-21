#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

typedef int32_t Positional_Int;       // позиционное целочисленное
typedef double Positional_Float;       // позиционное вещественное
typedef int16_t Module;               // целочисленный модуль разряда ССОК
typedef std::vector<Module> Modules;  // вектор модулей разрядов ССОК

// класс оснований СОК
class RnsBase {
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
    RnsBase(const Modules& p0, Positional_Int S0);  // конструктор из вектора оснований
};

// класс представления числа в СОК
class RnsFixed {
   public:
    std::reference_wrapper<const RnsBase> base;  // ссылка на класс оснований
    Modules a;                                   // вектор значений модулей (a1, a2 .. an)

    RnsFixed(Positional_Int x, const RnsBase& base0);  // из целого позиционного
    RnsFixed(const Modules& a0, const RnsBase& base0); // из вектора остатков
    template<typename ModType> static ModType mod(Positional_Int x, ModType p); // математический модуль числа
    template<typename ModType> static ModType mod_inverse(Positional_Int a, ModType p); // обратное по модулю число a^-1 mod p
    RnsFixed& operator+=(const RnsFixed& y);  // this += y
    RnsFixed operator+(const RnsFixed& y) const;  // z = this + y
    RnsFixed& operator-=(const RnsFixed& y);  // this -= y
    RnsFixed operator-(const RnsFixed& y) const;  // z = this - y
    RnsFixed operator-() const;  // Унарный минус z = -this
    RnsFixed& operator*=(const RnsFixed& y);  // this *= y
    RnsFixed operator*(const RnsFixed& y) const;  // z = this * y
    RnsFixed& operator/=(const RnsFixed& y);  // this /= y
    RnsFixed operator/(const RnsFixed& y) const;  // z = this / y
    bool operator==(const RnsFixed& y) const;  // (this == y)
    bool operator!=(const RnsFixed& y) const;  // (this != y)
    void div_int(const RnsFixed& y); // деление нацело
    bool operator<(const RnsFixed& y) const; // сравнение
    friend std::ostream& operator<<(std::ostream& os, const RnsFixed& y);  // для вывода на cout
    Positional_Float frac_crt_sum() const;  // дробная сумма дробной CRT
    Positional_Float to_positional_frac_crt_unscaled() const;  // в позиционное представление через дробную CRT (без масштабирования)
    Positional_Float to_positional_frac_crt() const;  // в позиционное представление через дробную CRT
    Positional_Float to_positional_crt() const;  // в позиционное представление CRT
    Positional_Float to_positional_mrc() const;  // в позиционное представление MRC
    Positional_Int get_rank() const; // получение ранга числа
    Positional_Int get_remainder(Positional_Int divisor) const; // остаток от деления
};

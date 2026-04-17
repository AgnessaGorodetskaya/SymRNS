#include <functional>
#include <iostream>
#include <numeric>
#include <vector>

typedef int32_t Positional_Int;       // позиционное целочисленное
typedef float Positional_Float;       // позиционное вещественное
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
    RnsBase(const Modules& p0);  // конструктор из вектора оснований
};

// класс представления числа в СОК
class RnsNumber {
   public:
    std::reference_wrapper<const RnsBase> base;  // ссылка на класс оснований
    Modules a;                                   // вектор значений модулей (a1, a2 .. an)

    RnsNumber(Positional_Int x, const RnsBase& base0);  // из целого позиционного
    RnsNumber(const Modules& a0, const RnsBase& base0); // из вектора остатков
    template<typename ModType> static ModType mod(Positional_Int x, ModType p); // математический модуль числа
    template<typename ModType> static ModType mod_inverse(Positional_Int a, ModType p); // обратное по модулю число a^-1 mod p
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
    Positional_Int to_positional_crt() const;  // в позиционное представление
    Positional_Int to_positional_mrc() const;  // в позиционное представление
    Positional_Int get_rank() const; // получение ранга числа
    // RnsNumber round(Positional_Int b) const;  // округление к ближайшему кратному x
};

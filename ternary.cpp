// Операции над троичными числами по 64 трита (каждый трит = 2 бита: minus, plus)
// Код трита: 0 = 00, +1 = 01, -1 = 10. Код 11 запрещён

#include <iostream>

constexpr int BITS = 64;

class Ternary
{
    uint64_t _minus{};
    uint64_t _plus{};

    // single_add: складывает два набора тритов a и b (каждый задан парами битов minus/plus)
    // Возвращает результат в res_minus/res_plus и переносы (не сдвинутые) в carry_minus/carry_plus
    // Все операции — побитовые, одновременно для 64 тритов
    // Входные пары не содержат код 11, который не должен появляться в корректных данных
    void single_add(uint64_t a_minus, uint64_t a_plus, uint64_t b_minus, uint64_t b_plus, uint64_t &res_minus,
                    uint64_t &res_plus, uint64_t &carry_minus, uint64_t &carry_plus) const
    {
        // нулевые маски
        uint64_t a0 = ~(a_minus | a_plus);
        uint64_t b0 = ~(b_minus | b_plus);

        // результирующие биты
        res_minus = (a_plus & b_plus) // +1 + +1 -> -1 (перенос carry+)
                    | (a_minus & b0)  // -1 +  0 -> -1
                    | (a0 & b_minus); //  0 + -1 -> -1

        res_plus = (a_minus & b_minus) // -1 + -1 -> +1 (перенос carry-)
                   | (a0 & b_plus)     //  0 + +1 -> +1
                   | (a_plus & b0);    // +1 +  0 -> +1

        // первичные переносы, не сдвинутые
        carry_minus = a_minus & b_minus; // -1 + -1 -> carry -1 в следующий трит
        carry_plus = a_plus & b_plus;    // +1 + +1 -> carry +1 в следующий трит
    }

    // add_trits: складывает a и b (любые 64-тритовые числа) и возвращает результат
    // Все аргументы — uint64_t, где i-бит соответствует i-триту
    void add_trits(uint64_t a_minus, uint64_t a_plus, uint64_t b_minus, uint64_t b_plus, uint64_t &res_minus,
                   uint64_t &res_plus) const
    {
        // Первичная сумма и переносы
        uint64_t c_minus, c_plus;
        single_add(a_minus, a_plus, b_minus, b_plus, res_minus, res_plus, c_minus, c_plus);

        // Итеративно добавляем переносы, пока они не исчезнут
        while (c_minus | c_plus)
        {
            // На каждой итерации перенос сдвигается на 1 трит в старший разряд
            c_minus <<= 1; c_plus <<= 1;
            // Затем мы складываем его с текущим результатом (используем single_add)
            // получая новый результат и новые переносы
            single_add(res_minus, res_plus, c_minus, c_plus, res_minus, res_plus, c_minus, c_plus);
        }
    }

  public:
    Ternary() {}; // Пустой конструктор
    Ternary(int64_t decimal) // Конструктор из десятичного представления
    {
        for (int i = 0; decimal != 0 && i < BITS; ++i)
        {
            int trit = decimal % 3;
            if (trit < 0)
                trit += 3; // для отрицательных чисел корректируем остаток
            switch (trit)
            {
            case 0: // трит = 0
                decimal /= 3;
                break;
            case 1: // трит = +1
                _plus |= (1ULL << i);
                decimal = (decimal - 1) / 3;
                break;
            case 2: // трит = -1
                _minus |= (1ULL << i);
                decimal = (decimal + 1) / 3;
                break;
            default:
                break;
            }
        }
    }
    int64_t to_decimal() const // Преобразование в десятичное число
    {
        int64_t res = 0;
        uint64_t power = 1; // 3^i = 3^0
        for (int i = 0; i < BITS; ++i)
        {
            uint64_t m = (_minus >> i) & 1ULL;
            uint64_t p = (_plus >> i) & 1ULL;
            if (m)
            {
                res -= power;
            }
            else if (p)
            {
                res += power;
            }
            power *= 3;
        }
        return res;
    }
    Ternary operator+(const Ternary &y) const // res = this + y
    {
        Ternary res{};
        add_trits(_minus, _plus, y._minus, y._plus, res._minus, res._plus);
        return res;
    }
    Ternary operator-(const Ternary &y) const // res = this - y
    {
        Ternary res{};
        // сложение с -y (меняем местами y_plus и y_minus)
        add_trits(_minus, _plus, y._plus, y._minus, res._minus, res._plus);
        return res;
    }
    Ternary operator-() const  // Унарный минус res = -this
    {
        Ternary res{};
        res._minus = _plus; // меняем местами _plus и _minus
        res._plus = _minus;
        return res;
    }
    // Печать 64-тритового числа в удобочитаемой форме
    friend std::ostream &operator<<(std::ostream &os, const Ternary &t)
    {
        for (int i = BITS - 1; i >= 0; --i)
        {
            uint64_t m = (t._minus >> i) & 1ULL;
            uint64_t p = (t._plus >> i) & 1ULL;
            if (m && !p)
            {
                os << "-";
            }
            else if (!m && !p)
            {
                os << "0";
            }
            else if (!m && p)
            {
                os << "+";
            }
            else
            {
                os << "?"; // код 11 запрещён
            }
        }
        return os;
    }
};

// Пример использования
int main(void)
{
    Ternary a(57934588798797);
    Ternary b(-32432408797923);
    std::cout << "  A = " << a << ' ' << a.to_decimal() << std::endl;
    std::cout << "  B = " << b << ' ' << b.to_decimal() << std::endl;
    Ternary c = a + b;
    std::cout << "A+B = " << c << ' ' << c.to_decimal() << std::endl;
    c = a - b;
    std::cout << "A-B = " << c << ' ' << c.to_decimal() << std::endl;
    c = -b;
    std::cout << " -B = " << c << ' ' << c.to_decimal() << std::endl;
    return 0;
}

//https://github.com/timacpp/fuzzy-numbers/tree/main/src
#ifndef FUZZYNUM_H
#define FUZZYNUM_H

#include <array>
#include <compare>
#include <iostream>

/** Type to represent triangular fuzzy number parameters **/
using real_t = double;

/** Class representing a triangular fuzzy number **/
class TriFuzzyNum {
public:
    TriFuzzyNum() = delete;

    /** Constructs a triangular fuzzy number with ascending order of values */
    constexpr TriFuzzyNum(real_t first, real_t second, real_t third)
            : lower{first}, modal{second}, upper{third} {
        this->order_values();
    }

    TriFuzzyNum(const TriFuzzyNum&) = default;

    TriFuzzyNum(TriFuzzyNum&&) = default;

    ~TriFuzzyNum() = default;

    [[nodiscard]] constexpr inline real_t lower_value() const {
        return this->lower;
    }

    [[nodiscard]] constexpr inline real_t modal_value() const {
        return this->modal;
    }

    [[nodiscard]] constexpr inline real_t upper_value() const {
        return this->upper;
    }

    friend TriFuzzyNum operator+(const TriFuzzyNum&, const TriFuzzyNum&);

    friend TriFuzzyNum operator-(const TriFuzzyNum&, const TriFuzzyNum&);

    friend TriFuzzyNum operator*(const TriFuzzyNum&, const TriFuzzyNum&);

    TriFuzzyNum operator-() const;

    TriFuzzyNum operator/(real_t) const;

    TriFuzzyNum& operator=(const TriFuzzyNum&) = default;

    TriFuzzyNum& operator=(TriFuzzyNum&&) = default;

    TriFuzzyNum& operator+=(const TriFuzzyNum&);

    TriFuzzyNum& operator-=(const TriFuzzyNum&);

    TriFuzzyNum& operator*=(const TriFuzzyNum&);

    friend std::ostream& operator<<(std::ostream&, const TriFuzzyNum&);

    friend std::strong_ordering operator<=>(const TriFuzzyNum&, const TriFuzzyNum&);

    friend constexpr inline bool operator==(const TriFuzzyNum&, const TriFuzzyNum&);

    friend constexpr inline bool operator!=(const TriFuzzyNum&, const TriFuzzyNum&);

private:
    real_t lower;
    real_t modal;
    real_t upper;

    constexpr static size_t RANK_SIZE{3};

    /**
     * Computes the rank of a triangular fuzzy number.
     * It is required for the comparison of the numbers.
     * For details check the link in README.md file.
     * @return rank of a triangular number
     */
    [[nodiscard]] std::array<real_t, RANK_SIZE> rank() const;

    constexpr inline void order_values() {
        if (lower > modal)
            std::swap(lower, modal);
        if (lower > upper)
            std::swap(lower, upper);
        if (modal > upper)
            std::swap(modal, upper);
    }
};

constexpr inline bool operator==(const TriFuzzyNum& lhs, const TriFuzzyNum& rhs) {
    return lhs.lower == rhs.lower && lhs.modal == rhs.modal && lhs.upper == rhs.upper;
}

constexpr inline bool operator!=(const TriFuzzyNum& lhs, const TriFuzzyNum& rhs) {
    return !(lhs == rhs);
}

/**
 * Compile-time only default fuzzy number generator.
 * @param v single parameter
 * @return (v,v,v)
 **/
consteval inline TriFuzzyNum crisp_number(real_t v) {
    return {v, v, v};
}

/** Compile-time only addition neutral element. **/
constinit inline TriFuzzyNum crisp_zero{crisp_number(0)};

#endif //FUZZYNUM_H
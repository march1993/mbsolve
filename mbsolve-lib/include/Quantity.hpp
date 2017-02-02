#ifndef QUANTITY_H
#define QUANTITY_H

#include <string>
#include <mbsolve-lib_EXPORTS.hpp>
#include <Types.hpp>

namespace mbsolve {

class MBSOLVE_LIB_EXPORT Quantity
{
private:
    std::string m_name;
    real m_value;
public:
    Quantity(const real& value = 0.0);

    virtual const real& operator()() const;

    operator real() const {
	return m_value;
    }

    /* TODO: operator with real? */

    Quantity operator+(const Quantity& rhs) const;

    Quantity& operator+=(const Quantity& rhs);

    Quantity operator-(const Quantity& rhs) const;

    Quantity& operator-=(const Quantity& rhs);

    Quantity operator*(const Quantity& rhs) const;

    Quantity& operator*=(const Quantity& rhs);

    Quantity operator/(const Quantity& rhs) const;

    Quantity& operator/=(const Quantity& rhs);

    bool operator<(const Quantity& rhs) const;

    bool operator>(const Quantity& rhs) const;

};

const Quantity PI = 3.14159265359;
const Quantity HBAR = 1.05457266e-34;
const Quantity MU0 =  PI() * 4e-7;
const Quantity EPS0 = 8.854187817e-12;
const Quantity E0 =  1.60217733e-19;

}

#endif

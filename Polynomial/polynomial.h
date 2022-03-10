#include <vector>
#include <limits>

#ifndef POLYNOMIAL_LIBRARY_H
#define POLYNOMIAL_LIBRARY_H

namespace polynomial{
    class Polynomial{
    public:
        Polynomial(std::vector<double> nums); // constructor given coefficients in ascending degree (coefs[0] should be constant term)
        Polynomial(std::vector<double> x, std::vector<double> y, int degree); // fit least squares regression on x and y
        Polynomial operator+(const auto& num);
        Polynomial operator-(const auto& num);
        Polynomial operator*(const auto& num);
        Polynomial operator/(const auto& num);
        Polynomial operator+(Polynomial other);
        Polynomial operator-(Polynomial other);
        double operator()(const auto& x);
        Polynomial integ(); // compute the indefinite integral
        Polynomial integ(const auto& lbnd); // the resulting polynomial should be 0 at lbnd
        double integ(const auto& lb, const auto& ub); // compute the definite integral over [lb, ub]
        Polynomial deriv();
        std::vector<double> roots(const double& real_threshold = 1e-8); // find all real roots of the polynomial
        std::vector<double> roots(const auto& lb, const auto& ub, const double& real_threshold = 1e-8); // find all real roots within a given interval
        std::vector<double> budans_roots(const auto& lb = -std::numeric_limits<float>::max(), const auto& ub = std::numeric_limits<float>::max(), const int& newton_iter = 5, const double& precision = 1e-3); // faster for high degree polynomials (>20)
        std::vector<double> solve(const auto& num); // solve for all occurrences of when the polynomial equals a certain value
        double solve(const auto& num, const auto& lb, const auto& ub, const int& newton_iter = 7); //solve for single occurrence of when the polynomials equals a given value within a certain range
        void print();
        double val(const auto& x); // return f(x)
    private:
        int degree;
        std::vector<double> coefs;
    };
}

#endif //POLYNOMIAL_LIBRARY_H
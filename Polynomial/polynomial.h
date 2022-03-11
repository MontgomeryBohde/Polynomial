#include <vector>
#include <limits>


#ifndef POLYNOMIAL_LIBRARY_H
#define POLYNOMIAL_LIBRARY_H

namespace polynomial{
    class Polynomial{
    public:
        Polynomial(const std::vector<double>& nums); // constructor given coefficients in ascending degree (coefs[0] should be constant term)
        Polynomial(const std::vector<double>&, const std::vector<double>&, const int& degree); // fit least squares regression on x and y
        Polynomial operator+(const double& num);
        Polynomial operator-(const double& num);
        Polynomial operator*(const double& num);
        Polynomial operator/(const double& num);
        Polynomial operator+(const Polynomial& other);
        Polynomial operator-(const Polynomial& other);
        double operator()(const double& x);
        Polynomial integ(); // compute the indefinite integral
        Polynomial integ(const double& lbnd); // the resulting polynomial should be 0 at lbnd
        double integ(const double& lb, const double& ub); // compute the definite integral over [lb, ub]
        Polynomial deriv();
        std::vector<double> roots(const double& real_threshold = 1e-8); // find all real roots of the polynomial
        std::vector<double> roots(const double& lb, const double& ub, const double& real_threshold = 1e-8); // find all real roots within a given interval
        std::vector<double> budans_roots(const double& lb = -std::numeric_limits<float>::max(), const double& ub = std::numeric_limits<float>::max(), const int& newton_iter = 5, const double& precision = 1e-3); // faster for high degree polynomials (>20)
        std::vector<double> solve(const double& num); // solve for all occurrences of when the polynomial equals a certain value
        double solve(const double& num, const double& lb, const double& ub, const int& newton_iter = 7); //solve for single occurrence of when the polynomials equals a given value within a certain range
        void print();
        double val(const double& x); // return f(x)
    private:
        int degree;
        std::vector<double> coefs;
    };
}

#endif //POLYNOMIAL_LIBRARY_H
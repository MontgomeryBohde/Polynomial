#include "polynomial.h"
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <stack>
#include <array>
#include <complex>

using namespace polynomial;

Polynomial::Polynomial(std::vector<double> nums) {
    /* Instantiate a polynomial from a list of coefficients in increasing order. You must include all terms, including
     * a constant term at the beginning and terms with a 0 coefficient EX: x^2 = {0, 0, 1}
     */

    degree = nums.size() - 1; // subtract one because the polynomial vector also includes a constant term
    coefs = nums;
    // remove extra terms at the front if the coeficient is 0 to ensure stability with rest of the code
    while (coefs.back() == 0 and degree > 0){
        coefs.pop_back();
        degree -= 1;
    }
}

Polynomial::Polynomial(std::vector<double> x, std::vector<double> y, int order){
    /* Fit an n-th degree polynomial of y as a function of x using least squares regression. Order is the degree
     * of the resulting polynomial
     */

    //necessary assumptions in order to be able to do calculations
    assert(x.size() == y.size());
    assert(x.size() >= order + 1);

    Eigen::MatrixXd X(x.size(), order+1);
    Eigen::VectorXd Y = Eigen::VectorXd::Map(&y.front(), y.size());
    Eigen::VectorXd result;

    // populate matrix
    for (size_t i = 0; i < x.size(); i++){
        for (size_t j = 0; j < order + 1; j++){
            X(i, j) = pow(x.at(i), j);
        }
    }
    // solve matrix
    result = X.householderQr().solve(Y);
    // add coefficients to the coefs vector
    for (int i = 0; i < order+1; i++){
        coefs.push_back(result[i]);
    }
    degree = order;
}

Polynomial Polynomial::operator+(const auto&num){
    std::vector<double> new_coefs(coefs);
    new_coefs[0] += num;
    return Polynomial(new_coefs);
}

Polynomial Polynomial::operator-(const auto& num) {
    std::vector<double> new_coefs(coefs);
    new_coefs[0] -= num;
    return Polynomial(new_coefs);
}

Polynomial Polynomial::operator*(const auto& num){
    std::vector<double> new_coefs(degree+1);
    std::transform(coefs.begin(), coefs.end(), new_coefs.begin(), [&num](double elm){return elm *= num;});
    return Polynomial(new_coefs);
}

Polynomial Polynomial::operator/(const auto& num) {
    assert(num != 0);
    std::vector<double> new_coefs(degree+1);
    std::transform(coefs.begin(), coefs.end(), new_coefs.begin(), [&num](double elm){return elm /= num;});
    return Polynomial(new_coefs);
}

Polynomial Polynomial::operator+(Polynomial other){
    int new_degree = std::max(degree, other.degree);
    std::vector<double> new_coefs(new_degree+1);
    for (int i = 0; i <= degree; i++){
        new_coefs[i] += coefs[i];
    }
    for (int i = 0; i <= other.degree; i++){
        new_coefs[i] += other.coefs[i];
    }
    return Polynomial(new_coefs);
}

Polynomial Polynomial::operator-(Polynomial other){
    int new_degree = std::max(degree, other.degree);
    std::vector<double> new_coefs(new_degree+1);
    for (int i = 0; i <= degree; i++){
        new_coefs[i] += coefs[i];
    }
    for (int i = 0; i <= other.degree; i++){
        new_coefs[i] -= other.coefs[i];
    }
    return Polynomial(new_coefs);
}

double Polynomial::operator()(const auto& x){
    return val(x);
}

double Polynomial::val(const auto& x){
    /* Return the value of the polynomial at x
     */
    double res = coefs[0];
    for (int i = 1; i <= degree; i++){
        res += coefs[i] * pow(x, i);
    }
    return res;
}

std::vector<double> Polynomial::roots(const double& real_threshold){
    /* Return all real roots of the polynomial. The roots are calculated using the eigenvalues of the companion matrix
     * and will be automatically sorted.
     * @param real_threshold: The maximum imagninary threshold at which a root can still be considered real
     */
    std::vector<double> roots;

    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();

    companion(0, degree-1) = -(coefs[0] / coefs.back()); // intialize first row in companion matrix
    for (int row = 1; row < degree; row++){ // initialize rest of companion matrix
        companion(row, row-1) = 1;
        companion(row, degree-1) = -(coefs[row] / coefs.back());
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(companion, false); // false because we do not need the eigenvectors, only the eigenvalues
    for (std::complex<double> l : es.eigenvalues()){
        if (abs(l.imag()) < real_threshold){
            roots.push_back(l.real());
        }
    }
    return roots;
}

std::vector<double> Polynomial::roots(const auto& lb, const auto& ub, const double& real_threshold){
    /* Return all real roots within the specified bounds, exclusive: (lb, ub)
     */
    std::vector<double> roots;

    Eigen::MatrixXd companion(degree, degree);
    companion.setZero();

    companion(0, degree-1) = -(coefs[0] / coefs.back()); // intialize first row in companion matrix
    for (int row = 1; row < degree; row++){ // initialize rest of companion matrix
        companion(row, row-1) = 1;
        companion(row, degree-1) = -(coefs[row] / coefs.back());
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(companion, false); // false because we do not need the eigenvectors, only the values
    for (std::complex<double> l : es.eigenvalues()){
        if (abs(l.imag()) < real_threshold && l.real() > lb && l.real() < ub){
            roots.push_back(l.real());
        }
    }
    return roots;
}

std::vector<double> Polynomial::budans_roots(const auto& lb, const auto& ub, const int& newton_iter, const double& precision){
    /* This polynomial root function uses budan's theorem to aproximate a range that each root lies within, then
     * newtons method is used to find the exact value of the root. This method is faster than the other roots method when
     * the polynomial has a degree higher than 10
     * @param newton_iter is the number of iterations of newtons method, more iterations will increase precision of the roots
     * @param precision is the maximum range in the budan's approximation of the roots. Higher precision will be better but more time consuming
     */
    // IT IS VERY HIGHLY RECOMMENDED YOU PASS BOUNDS TO THIS METHOD TO ENSURE IT DOES NOT OVERFLOW ON HIGH DEGREE POLYNOMIALS
    std::vector<double> roots;
    std::vector<Polynomial> eqs;
    eqs.push_back(deriv()); // add the first derivative to the budan table
    for (int i = 1; i < degree; i++){ // populate vector with the rest of the derivatives of the polynomial
        eqs.push_back(eqs.back().deriv()); // get last polynomial, take its derivative, then add it to the end of the vector
    }

    // if degree is > 8, and bounds not passed, the doubles will overflow
    double x1 = lb;
    double x2 = ub;

    std::vector<std::array<double, 2>> one_root_regions; // final regions, each with one root, note regions are half-open: (x1, x2]
    std::stack<std::array<double, 2>> regions_stack;

    regions_stack.push({x1, x2});
    while (!regions_stack.empty()){
        std::array<double, 2> reg = regions_stack.top();
        double l = reg[0];
        double r = reg[1];
        regions_stack.pop();
        int left_sign_switches = 0;
        double last = val(l);
        for (Polynomial p : eqs){
            // increase count whenever the last and p.val(l) have opposite signs
            double val = p.val(l);
            if ((val > 0) && (last < 0)){
                last = val;
                left_sign_switches++;
            }else if((val < 0) && (last > 0)){
                last = val;
                left_sign_switches++;
            }else if(last == 0){
                last = val;
            }
        }

        int right_sign_switches = 0;
        last = val(r);
        for (Polynomial p : eqs){
            // increase count whenever the last and p.val(r) have opposite signs
            double val = p.val(r);
            if ((val > 0) && (last < 0)){
                last = val;
                right_sign_switches++;
            }else if((val < 0) && (last > 0)){
                last = val;
                right_sign_switches++;
            }else if(last == 0){
                last = val;
            }
        }

        int s = left_sign_switches - right_sign_switches;
        if (s == 0){ // no roots in region
            continue;
        }else if (s == 1){
            if ((r - l) < precision){
                one_root_regions.push_back(reg);
            }else{
                double m = (l + r)/2;
                regions_stack.push({l, m});
                regions_stack.push({m, r});
            }
        }else{ // if more then 1 roots in the region, split into left and rights sides then add both to the stack
            if ((r - l) < (precision)){
                continue;
            }else{
                double m = (l + r)/2;
                regions_stack.push({l, m});
                regions_stack.push({m, r});
            }
        }
    }

    for (std::array<double, 2> arr : one_root_regions){
        double x = (arr[0] + arr[1]) / 2; // initial geuss is the midpoint of the region
        for (int i = 0; i < newton_iter; i++){
            x = x - (val(x) / eqs[0].val(x));
        }
        if (x > arr[0] && x <= arr[1]){
            roots.push_back(x);
        }
    }
    sort(roots.begin(), roots.end());
    return roots;
}

Polynomial Polynomial::deriv(){
    /* Take the derivative of the polynomial and return it as a new polynomial
     */
    std::vector<double> res_coefs(degree);
    for (int i = 1; i <= degree; i++){
        res_coefs[i-1] = coefs[i] * i;
    }
    return Polynomial(res_coefs);
}

Polynomial Polynomial::integ(){
    /* Compute the indefinite integral of the polynomial, this does not include any integration constant
     */
    std::vector<double> new_coefs(degree+2);
    for (int i = 0; i <= degree; i++){
        new_coefs[i+1] = coefs[i] / (i+1);
    }
    return Polynomial(new_coefs);
}

Polynomial Polynomial::integ(const auto& lbnd) {
    /* Compute the indefinite integral of the polynomial, but scale it to be 0 at lbnd
     */
    std::vector<double> new_coefs(degree+2);
    for (int i = 0; i <= degree; i++){
        new_coefs[i+1] = coefs[i] / (i+1);
    }
    Polynomial temp(new_coefs);
    return (temp - temp.val(lbnd));
}

double Polynomial::integ(const auto& lb, const auto& ub) {
    /* Compute the definite integral of the polynomial from lb to ub
     */
    std::vector<double> new_coefs(degree+2);
    for (int i = 0; i <= degree; i++){
        new_coefs[i+1] = coefs[i] / (i+1);
    }
    Polynomial poly(new_coefs);
    return poly.val(ub) - poly.val(lb);
}

std::vector<double> Polynomial::solve(const auto& num){
    /* Solve for all occurrences of when the polynomial equals num, by default this uses the companion matrix method, so
     * for very high degree polynomials you should consider using budan's method
     */
    std::vector<double> new_coefs = coefs;
    new_coefs[0] -= num;
    Polynomial temp(new_coefs);
    return temp.roots();
}

double Polynomial::solve(const auto& num, const auto& lb, const auto& ub, const int& newton_iter){
    /* This version of solve is meant for when a function is strictly monotonic over the interval (lb, ub) and it is known
     * that exactly one solution to exists within the bounds, for example solving a cdf function.
     */
    std::vector<double> new_coefs = coefs;
    new_coefs[0] -= num;
    Polynomial f(new_coefs);
    Polynomial f_prime = f.deriv();
    double x = (lb+ub)/2; // initial geuss is midpoint
    for (int i = 0; i < newton_iter; i++){
        x = x - (f.val(x) / f_prime.val(x));
    }
    return x;
}

void Polynomial::print() {
    if (degree == 0){
        std::cout << coefs[0] << std::endl;
        return;
    }
    std::cout << coefs[degree] << "x^" << degree;
    for (int i = degree-1; i > 0; i--){
        if (coefs[i] > 0){
            std::cout << " + ";
        }else{
            std::cout << " - ";
        }
        std::cout << abs(coefs[i]) << "x^" << i;
    }
    if (coefs[0] > 0){
        std::cout << " + " << coefs[0] << std::endl;
    }else{
        std::cout << " - " << abs(coefs[0]) << std::endl;
    }
    return;
}
# Polynomial
C++ Polynomial Library with Least Squares Regression, Root Finding, Differentiation, Integration and Arithmetic

# Dependencies
The only dependency is Eigen, which is a C++ Linear Algebra library required for the least-squares regression and root-finding algorithims.

# Example Code: 
'''
// create a polynomial from a list of coefficients:
Polynomial poly({-4, 0, 1); // x^2 - 4

// create a 3rd degree polynomial using least-squares regression:
std::vector<double> x = {0.0192341804504395, 0.0394501686096191, ... ,0.0987141132354736,  0.119336366653442};
std::vector<double> y = {1.8, 1.86, ...  ,2.81, 2.87};
Polynomial poly2(x, y, 3);
  
// find the roots
std::vector<double> roots = poly2.roots();

// take the derivative
Polynomial derivative = poly2.deriv();
  
// compute the definite integral
double integral = poly2.integ(0, 3);
 
// see the docs for all the available functions
'''

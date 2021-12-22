#include<Eigen/Dense>

using namespace std;
using namespace Eigen;

double CHI2(MatrixXd y, MatrixXd y_measured, MatrixXd W);
MatrixXd function_y(MatrixXd para, MatrixXd x);
double rho_function(MatrixXd y, MatrixXd y_measured, MatrixXd W, MatrixXd h, MatrixXd Jacobian, double lambda, double chi2, double chi2plus);
MatrixXd jacobian_function(MatrixXd x,  MatrixXd para, MatrixXd y, MatrixXd initial_deflection);
MatrixXd levenberg_fit(MatrixXd W, MatrixXd y_measured, MatrixXd x, MatrixXd initial_deflection, MatrixXd para, double lambda);

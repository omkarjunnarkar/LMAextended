#include"myfunctions.h"
#include<iostream>
#include<Eigen/Dense>
#include<math.h>

using namespace std;
using namespace Eigen;

/*To compute the value of function y=f(x,z) ; arguments : Parameters, x,z,p */
MatrixXd function_y(MatrixXd para, MatrixXd x) {
	int number_of_data_points = x.rows();
	MatrixXd y(number_of_data_points, 1);
	for (int i = 0; i < number_of_data_points; i++) {

		//y(i, 0) = para(0, 0) * pow(x(i, 0), 2) + para(1, 0) * exp(pow(z(i, 0), 2) + 1) + para(2, 0);

		y(i, 0) = para(0, 0) * pow(x(i, 0), 1) + para(1, 0) * exp(para(2,0)*x(i, 0))  ;
	}
	return y;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*To compute the squared difference (rf. Literature EQN:2 ) */
double CHI2(MatrixXd y, MatrixXd y_measured, MatrixXd W) {

	MatrixXd val = (y_measured - y).transpose() * W * (y_measured - y);
	return val(0, 0);
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*To compute the deciding factor for parameter & damping factor upgradation (rf. Literature EQN:12 )*/
double rho_function(MatrixXd y, MatrixXd y_measured, MatrixXd W, MatrixXd h, MatrixXd Jacobian, double lambda, double chi2, double chi2plus) {
	cout << "h=" << h;
	cout << "lambda=" << lambda;
	cout << "ch2=" << chi2;
	cout << "ch2plus=" << chi2plus;
	

	MatrixXd value = h.transpose() * (lambda * h + Jacobian.transpose() * W * (y_measured - y));
	// [ chi2(p) - chi2(p+h) ] / [ hT*(lamda*h + JT*W*(measuredy-computedy) ]
	return (chi2 - chi2plus) / value(0, 0);
}

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*To compute the Jacobian Matrix by finite difference scheme*/
MatrixXd jacobian_function(MatrixXd x,  MatrixXd para, MatrixXd y, MatrixXd initial_deflection) {

	MatrixXd Jacobian_Matrix(x.rows(), para.rows());
	MatrixXd y_para = y;										//Original Value of function
	MatrixXd y_deflected(x.rows(), 1);

	for (int i = 0; i < para.rows(); i++) {

		para(i, 0) = para(i, 0) + initial_deflection(i, 0);		/*Changing the parameters one by one */

		y_deflected = function_y(para, x);				/*Computing the deflected function arrray */
		for (int j = 0; j < x.rows(); j++) {

			// [f(v, p + dp) - f(v, p) ] / [dp] 

			Jacobian_Matrix(j, i) = (y_deflected(j, 0) - y_para(j, 0)) / initial_deflection(i, 0); 
		}
		para(i, 0) = para(i, 0) - initial_deflection(i, 0);		/*Bringing back the parametes to original value*/
	}
	return Jacobian_Matrix;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*Levenberg-Marquardt Algorithm ; Arguements: Actual function value, x ,Initial deflection of parameters , Initial guess of parameters , Damping factor Lambda */
MatrixXd levenberg_fit(MatrixXd W, MatrixXd y_measured, MatrixXd x,  MatrixXd initial_deflection, MatrixXd para, double lambda) {

	MatrixXd IdentityMat = MatrixXd::Identity(para.rows(), para.rows());
	int count = 0;								//Number of Iterations
	double v = 2.0;								//damping facctor manipulator
	MatrixXd y = function_y(para, x);		//Initial value of 'y' using initial guess of parameters
	double chi2p = CHI2(y, y_measured, W);		//Initial Chi2

	do {
		count++;
		cout << "\nCount = " << count << endl;

		//cout << "Parameters=\n" << para << endl;

		y = function_y(para, x);													/*Function Array*/
		double chi2p_plus_h;															/*Defining Residual for next iteration*/
		MatrixXd J = jacobian_function(x, para,y, initial_deflection);			/*Jacobian Matrix*/
		MatrixXd intermediate = J.transpose() * W * J + lambda * IdentityMat;		
		cout << "J="<<J << endl;
		MatrixXd h = intermediate.inverse() * (J.transpose() * W * (y_measured - y));	/*Computed Change in parameters*/
		MatrixXd y_plus_h = function_y(para + h, x);
		//cout << "h =\n" << h << endl;

		chi2p = CHI2(y, y_measured, W);
		chi2p_plus_h = CHI2(y_plus_h, y_measured, W);									/*Updated Residual after changing Parameters*/

		double rho = rho_function(y, y_measured, W, h, J, lambda, chi2p, chi2p_plus_h);	/*Computing the deciding factor for Lambda change*/

		cout << "chi2p = " << chi2p << ", rho = " << rho << endl;

		//Update Criterion ( rf. Literature PG 4, 5, 7 )

		if (rho > 1e-2) {																
			para = para + h;
			v = 2.0;
			/*if ((1 - pow((2 * rho) - 1, 3)) <= 0) { lambda = lambda * (1 / 3); }
			else {
				lambda = lambda * (1 - pow((2 * rho) - 1, 3));
			}*/
			lambda = lambda / 3;
			cout << "Parameters updated, Lambda = " << lambda << endl;
		}
		else {
			lambda = lambda * v;
			v = 2 * v;
			cout << "Parameters NOT updated, Lambda = " << lambda << endl;
		}

		cout << para << endl;
	} while ((chi2p / (x.rows() - para.rows() + 1)) > 1e-2);
	

	return para;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
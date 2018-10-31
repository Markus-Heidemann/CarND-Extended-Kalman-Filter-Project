#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  	VectorXd rmse(4);
	for(int i=0; i < estimations.size(); i++)
    {
      	VectorXd residual(4);
      	residual = estimations[i] - ground_truth[i];
      	residual = residual.array() * residual.array();
		rmse += residual;
    }
  	rmse = rmse.array().sqrt();
  	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	float px = x_state[0];
  	float py = x_state[1];
  	float vx = x_state[2];
  	float vy = x_state[3];
  
  	float c1 = px*px + py*py;
  	float c2 = sqrt(c1);
  	float c3 = vx*py + vy*px;
  	float c4 = pow(c1, 1.5);
  
  	MatrixXd Hj(3,4);
	Hj << px/c2, py/c2, 0, 0,
  		  -py/c1, px/c1, 0, 0,
  		  py*c3/c4, px*c3/c4, px/c2, py/c2;
  
  	return Hj;
}

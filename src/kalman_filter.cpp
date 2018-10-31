#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, Eigen::MatrixXd &Hj_in, MatrixXd &R_in,
                        MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  Hj_ = Hj_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  VectorXd y = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  x_ += K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  VectorXd y = z - h(x_);

  // make sure phi is in a range of -pi to pi
  float pi = std::acos(-1);
  while (y(1) < -pi)
  {
    y(1) += 2 * pi;
  }
  while (y(1) > pi)
  {
    y(1) -= 2 * pi;
  }

  MatrixXd S = Hj_ * P_ * Hj_.transpose() + R_;
  MatrixXd K = P_ * Hj_.transpose() * S.inverse();

  x_ += K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * Hj_) * P_;
}

VectorXd KalmanFilter::h(const VectorXd &x)
{
  float px = x(0);
  float py = x(1);
  float vx = x(2);
  float vy = x(3);

  float rho = sqrt(px * px + py * py);
  float phi = atan2(py, px);
  float rho_dot = (px * vx + py * vy) / rho;

  VectorXd res(3);
  res << rho, phi, rho_dot;
  return res;
}
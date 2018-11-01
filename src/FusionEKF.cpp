#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() : noise_ax_(9.0), noise_ay_(9.0)
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;

  Qv_ = MatrixXd(2, 2);
  Qv_ << noise_ax_, 0,
      0, noise_ay_;

  ekf_ = KalmanFilter();

  tools = Tools();
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian
      * coordinates.
    */
    // first measurement

    previous_timestamp_ = measurement_pack.timestamp_;
    MatrixXd P_init(4, 4);
    MatrixXd F_init(4, 4);
    MatrixXd Q_init(4, 4);
    MatrixXd Hj_init(3, 4);

    F_init.setZero();
    Q_init.setZero();
    Hj_init.setZero();

    P_init << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float d = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float d_dot = measurement_pack.raw_measurements_(2);

      VectorXd cart_meas(4);
      cart_meas << d * sin(phi),
                   d * cos(phi),
                   d_dot * sin(phi),
                   d_dot * cos(phi);
      MatrixXd Hj = tools.CalculateJacobian(cart_meas);
      ekf_.Init(cart_meas, P_init, F_init, H_laser_, Hj, R_radar_, Q_init);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      /**
      Initialize state.
      */
      VectorXd meas(4);
      meas << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1],
          0,
          0;
      ekf_.Init(meas, P_init, F_init, H_laser_, Hj_init, R_laser_, Q_init);
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ << 1, 0, dt, 0,
      0, 1, 0, dt,
      0, 0, 1, 0,
      0, 0, 0, 1;

  MatrixXd G(4, 2);
  G << 0.5 * dt * dt, 0,
      0, 0.5 * dt * dt,
      dt, 0,
      0, dt;

  ekf_.Q_ = G * Qv_ * G.transpose();

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);
    if (!ekf_.Hj_.isZero(0))
    {
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    else
    {
    }
  }
  else
  {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

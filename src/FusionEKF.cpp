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
FusionEKF::FusionEKF() : noise_ax_(9.0), noise_ay_(9.0) {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

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

  Qv_ = MatrixXd(2,2);
  Qv_ << noise_ax_, 0,
       0, noise_ay_;

  VectorXd x(4);
  MatrixXd P(4,4);
  MatrixXd F(4,4);
  MatrixXd Q(4,4);

  ekf_ = KalmanFilter();
  ekf_.Init(x, P, F, H_laser_, Hj_, R_laser_, Q);

  tools = Tools();

  std::cout << "FusionEKF::FusionEKF() done!\n";
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  std::cout << "FusionEKF::ProcessMeasurement called" << "\n";

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    VectorXd init_meas(4);
    init_meas << 1, 1, 1, 1;
    ekf_.x_ = init_meas;
    ekf_.Hj_ = tools.CalculateJacobian(init_meas);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float d = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float d_dot = measurement_pack.raw_measurements_(2);

      VectorXd cart_meas(4);
      cart_meas << d*sin(phi), d*cos(phi), d_dot*sin(phi), d_dot*cos(phi);
      ekf_.x_ = cart_meas;
      ekf_.R_ = R_radar_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ = measurement_pack.raw_measurements_;
      ekf_.R_ = R_laser_;
    }

    // done initializing, no need to predict or update
    std::cout << "Initialization done!\n";
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

  int dt = measurement_pack.timestamp_ - previous_timestamp_;

  ekf_.F_ << 1, 0, dt, 0,
         0, 1, 0, dt,
         0, 0, 1, 0,
         0, 0, 0, 1;

  MatrixXd G(4,2);
  G << 0.5*dt*dt, 0,
      0, 0.5*dt*dt,
      dt, 0,
      0, dt;

  ekf_.Q_ = G*Qv_*G.transpose();

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.Hj_ = tools.CalculateJacobian(measurement_pack.raw_measurements_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

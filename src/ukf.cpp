#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.45;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = n_x_ + 2;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;

  // weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  // calculate weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if (!is_initialized_) {

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float rho = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];
            float rho_dot = meas_package.raw_measurements_[2];
            x_ << rho * cos(phi), rho * sin(phi), rho_dot, 0.0, 0.0;
        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            // Initialize the state ekf_.x_ with the first measurement.
            // set the state with the initial location and zero velocity
            float px = meas_package.raw_measurements_[0];
            float py = meas_package.raw_measurements_[1];
            x_ << px, py, 0.0, 0.0, 0.0;
        }

        P_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, 1;

        time_us_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  double time_delta = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Updates x_ and P_
  Prediction(time_delta);

  /*****************************************************************************
  *  Update
  ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else {
    // Lidar updates
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  VectorXd x_aug = VectorXd(7);
    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;

    // process noise
    P_aug(5,5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    MatrixXd sqrt_P_aug= P_aug.llt().matrixL();

    // create augmented sigma points
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
    Xsig_aug.col(0) = x_aug;
    for (int i=0; i<n_aug_; i++) {
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*sqrt_P_aug.col(i);
        Xsig_aug.col(n_aug_+i+1) = x_aug - sqrt(lambda_+n_aug_)*sqrt_P_aug.col(i);
    }

    // Predict sigma points
    for (int i=0; i<2*n_aug_+1; i++) {
      double px = Xsig_aug(0,i);
      double py = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double yaw = Xsig_aug(3,i);
      double yawRate = Xsig_aug(4,i);
      double v_acc_noise = Xsig_aug(5,i);
      double yaw_acc_noise = Xsig_aug(6,i);
      double px_pred;
      double py_pred;
      double v_pred;
      double yaw_pred;
      double yawRate_pred;

      double delta_t_sq = delta_t*delta_t;
      yaw_pred = yaw + yawRate*delta_t;
      if (fabs(yawRate) < 0.001) {
        px_pred = px + (v*cos(yaw)*delta_t);
        py_pred = py + (v*sin(yaw)*delta_t);
      } else {
        px_pred = px + ((v/yawRate) * (sin(yaw_pred)-sin(yaw)));
        py_pred = py + ((v/yawRate) * (-cos(yaw_pred)+cos(yaw)));
      }

      // add noise
      px_pred  += 0.5*delta_t_sq*cos(yaw)*v_acc_noise;
      py_pred  += 0.5*delta_t_sq*sin(yaw)*v_acc_noise;
      yaw_pred += 0.5*delta_t_sq*yaw_acc_noise;
      yawRate_pred = yawRate + delta_t*yaw_acc_noise;
      v_pred = v + delta_t*v_acc_noise;

      Xsig_pred_(0,i) = px_pred;
      Xsig_pred_(1,i) = py_pred;
      Xsig_pred_(2,i) = v_pred;
      Xsig_pred_(3,i) = yaw_pred;
      Xsig_pred_(4,i) = yawRate_pred;
    }

    // Predict state mean
    x_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {
      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    //predict state covariance matrix
    P_.fill(0.0);
    for (int i=0; i<2*n_aug_+1; i++) {

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;

      //angle normalization
      while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;
      P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/

　//set measurement dimension, lidar can measure px, py

  int n_z = 2;
  // 2 x 15
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

      // extract values for better readability
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);

      // measurement model
      Zsig(0, i) = p_x;
      Zsig(1, i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;

      //angle normalization
      while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
      while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

      S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;
  S = S + R;

  /*****************************************************************************
  *  Update
  ****************************************************************************/

  VectorXd z = VectorXd(n_z);
  // measurement for laser
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  // 5 x 2
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      //angle normalization
      while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
      while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_laser_ = (meas_package.raw_measurements_ - z_pred).transpose() * S.inverse() *
               (meas_package.raw_measurements_ - z_pred);

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/


  //set measurement dimension, lidar can measure px, py
  
  int n_z = 2;

  // 2 x 15
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

      // extract values for better readability
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);

      // measurement model
      Zsig(0, i) = p_x;
      Zsig(1, i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;

      //angle normalization
      while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
      while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

      S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0, 0, std_laspy_ * std_laspy_;
  S = S + R;

  /*****************************************************************************
  *  Update
  ****************************************************************************/

  VectorXd z = VectorXd(n_z);
  // measurement for laser
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  // 5 x 2
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      //angle normalization
      while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
      while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //angle normalization
      while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  NIS_laser_ = (meas_package.raw_measurements_ - z_pred).transpose() * S.inverse() *
               (meas_package.raw_measurements_ - z_pred);

}

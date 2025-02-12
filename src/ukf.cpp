#include "ukf.h"
#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include "measurement_package.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
   is_initialized_ = false;

   n_x_ = 5;

   n_aug_ = 7;

   lambda_ = 3 - n_aug_;

   Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

   weights_ = VectorXd(2*n_aug_+1);

   weights_(0) = lambda_ / (lambda_ + n_aug_);
   for(int i=1; i < 2*n_aug_+1; i++){
     weights_(i) = 0.5 / (lambda_ + n_aug_);
   }


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

   if (!is_initialized_) {

    x_.fill(0.0);
    
    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true){
        x_ << meas_package.raw_measurements_(0), 
        meas_package.raw_measurements_(1), 
        0.2, 
        0,
        0;
      
        // set initial covariance matrix
        P_<< 0.01, 0, 0, 0, 0,
            0, 0.01, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
 
    }

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true){
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);
      
      double p_x = rho * cos(phi);
      double p_y = rho * sin(phi);
      double v = rho_dot;

      x_ << p_x,
           p_y,
           v,
           phi,
           0;
      
       // set initial covariance matrix
       P_<< 0.01, 0, 0, 0, 0,
           0, 0.01, 0, 0, 0,
           0, 0, 0.01, 0, 0,
           0, 0, 0, 0.09, 0,
           0, 0, 0, 0, 0.09;


    }
    
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Prediction step
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  // Measurement Update 
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } 
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

  
}

void UKF::Prediction(double delta_t) {

   /* Augmented State Vector and Sigma Points */

   // create augmented mean vector
   VectorXd x_aug = VectorXd(n_aug_);

   // create augmented state covariance
   MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
   // create sigma point matrix
   MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  // create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
 
   // create square root matrix
   MatrixXd L = P_aug.llt().matrixL();
 
   // create augmented sigma points
   Xsig_aug.col(0)  = x_aug;
   for (int i = 0; i< n_aug_; ++i) {
     Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
     Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
   }


   /* Sigma Points Prediction Step */

   for (int i = 0; i< 2*n_aug_+1; ++i) {

    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double yawd     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double px_p, py_p;
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v/yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
        px_p = p_x + (v * delta_t * cos(yaw));
        py_p = p_y + (v * delta_t * sin(yaw));
    }

    
    double v_p = v;
    double yaw_p = yaw + (yawd * delta_t);
    double yawd_p = yawd;

    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + (0.5 * nu_yawdd * delta_t * delta_t);
    yawd_p = yawd_p + (nu_yawdd * delta_t);

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;

   }

   /* Predict Mean and Covariance */

  x_.fill(0.0);
  for(int i=0; i < 2*n_aug_+1; i++){
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
    
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
      // state difference
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      // angle normalization
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

      P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}


void UKF::UpdateLidar(MeasurementPackage meas_package) {

   int n_z = 2;

   // create matrix for sigma points in measurement space
   MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
   Zsig.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     Zsig(0, i) = Xsig_pred_(0, i);
     Zsig(1, i) = Xsig_pred_(1, i);
   }
   
   // mean predicted measurement
   VectorXd z_pred = VectorXd(n_z);
   z_pred.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     z_pred += weights_(i) * Zsig.col(i);
   }

   // measurement covariance matrix S
   MatrixXd S = MatrixXd(n_z, n_z);
   S.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

     S += weights_(i) * z_diff * z_diff.transpose();
   }
   
   // add measurement noise covariance matrix
   MatrixXd R = MatrixXd(n_z, n_z);
   R << std_laspx_ * std_laspx_, 0,
        0, std_laspy_ * std_laspy_;
   S += R;
   
   // create matrix for cross correlation Tc
   MatrixXd Tc = MatrixXd(n_x_, n_z);
   Tc.fill(0.0);
   for (int i = 0; i < 2 * n_aug_ + 1; i++) {
     VectorXd z_diff = Zsig.col(i) - z_pred;
     VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

     Tc += weights_(i) * x_diff * z_diff.transpose();
   }
   
   // Kalman gain K;
   MatrixXd K = Tc * S.inverse();
   // residual
   VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

   // update state mean and covariance matrix
   x_ += K * z_diff;
   P_ -= K * S * K.transpose();

  double NIS_lidar = z_diff.transpose() * S.inverse() * z_diff;
  cout << "NIS lidar: " << NIS_lidar << endl;

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

   // create matrix for sigma points in measurement space
   MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
   // mean predicted measurement
   VectorXd z_pred = VectorXd(n_z);
   
   // measurement covariance matrix S
   MatrixXd S = MatrixXd(n_z,n_z);
 
   // transform sigma points into measurement space
   Zsig.fill(0.0);
   for (int i = 0; i< 2*n_aug_+1; ++i) {
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      
      // measurement model
      double rho, rho_dot, var_phi;
      rho = sqrt(p_x*p_x + p_y * p_y);
      var_phi = atan2(p_y,p_x);
      rho_dot = (p_x * cos(yaw) * v + p_y * sin(yaw) * v)/ sqrt(p_x * p_x + p_y * p_y);
  
      Zsig(0,i) = rho;
      Zsig(1,i) = var_phi;
      Zsig(2,i) = rho_dot;
      
   }

    // mean predicted measurement
    z_pred.fill(0.0);
    for(int i=0; i < 2*n_aug_+1; i++){
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }
      
      // innovation covariance matrix S
      MatrixXd R = MatrixXd(n_z,n_z);
      R <<  std_radr_ * std_radr_,            0,                          0,
                   0,             std_radphi_* std_radphi_,               0,
                   0,                         0,                   std_radrd_ *std_radrd_;

      S.fill(0.0);
     for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
         // state difference
         VectorXd z_diff = Zsig.col(i) - z_pred;

         // angle normalization
         while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
         while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
 
         S = S + weights_(i) * z_diff * z_diff.transpose();
     }
     // add measurement noise to innovation covariance
     S = S + R;

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        // measurement difference
       VectorXd z_diff = Zsig.col(i) - z_pred;
       // angle normalization
       while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
       while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose() ;
    }

    MatrixXd K = Tc * S.inverse();

   // residual
  VectorXd z_diff =   meas_package.raw_measurements_ - z_pred;

   // angle normalization
   while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
   while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

   x_ = x_ + (K * z_diff);
   P_ = P_ - K * S * K.transpose();

   double NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
   std::cout << "NIS radar: " << NIS_radar << std::endl;

}
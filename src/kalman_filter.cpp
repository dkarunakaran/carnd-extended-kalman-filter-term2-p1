#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  * Predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  *Update the state by using Kalman Filter equations for LIDAR
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  
  UpdateMeasurement(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  *Update the state by using Extended Kalman Filter equations for RADAR
  */

  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  //Checking the value is not zero
  if(px == 0. && py == 0.)
    return;
  
  float rho = sqrt(px*px + py*py);
  float theta = atan2(py, px);
  
  //Checking the value is not zero
  if (rho < 0.0001) {
    rho = 0.0001;
  } 
  float rho_dot = (px*vx + py*vy) / rho;
  

  //Finding h(x)
  VectorXd h = VectorXd(3); // h(x_)
  h << rho, theta, rho_dot;
  VectorXd y = z-h;

  //Normalize the angle between -pi to pi
  while (y[1] < -M_PI)
    y[1] += 2 * M_PI;
  while (y[1] > M_PI)
    y[1] -= 2 * M_PI;

  UpdateMeasurement(y);

}

void KalmanFilter::UpdateMeasurement(const VectorXd &y) {
  //Mesaurement updates
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  
  MatrixXd K = PHt * Si;
  
  //New estimate
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

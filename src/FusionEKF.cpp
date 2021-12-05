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
FusionEKF::FusionEKF() {
  noise_ax = 9;
  noise_ay = 9;

  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  
  //measurement covariance matrix - LIDAR
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - RADAR
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //Measurement matrix - LIDAR - Matrix should be 2x4
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //Measurement matrix - RADAR - Matrix should be 3x4
  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1; 

  //Initial transition matrix F_ - Matrix should be 4x4
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  //State covariance matrix P - Matrix should be 4x4
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 100, 0,
             0, 0, 0, 100;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      float ro     = measurement_pack.raw_measurements_(0); //ro
      float phi    = measurement_pack.raw_measurements_(1); //phi
      float ro_dot = measurement_pack.raw_measurements_(2); //ro_dot
      //Polar coordinate to cartesian coordinate - trigonometry rules
      ekf_.x_(0) = ro     * cos(phi); //px = ro*cos(phi)
      ekf_.x_(1) = ro     * sin(phi); //py = ro*cos(phi)      
      ekf_.x_(2) = ro_dot * cos(phi); //vx = ro_dot*cos(phi)
      ekf_.x_(3) = ro_dot * sin(phi); //vy = ro_dot*cos(phi)
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_(0) = measurement_pack.raw_measurements_(0); //px
      ekf_.x_(1) = measurement_pack.raw_measurements_(1); //py
      ekf_.x_(2) = 0.0; //no velocity data for Lidar
      ekf_.x_(3) = 0.0; //no velocity data for Lidar
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // Done initializing, no need to predict or update
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

  //Compute elapsed time - classroom material
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  if( dt > 0.001 ) 
  { 
    float dt_2 = dt   * dt; //dt^2
    float dt_3 = dt_2 * dt; //dt^3
    float dt_4 = dt_3 * dt; //dt^4

    //State transition matrix F
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
   

    //Process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
              0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
              dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //Radar updates - Call for CalculateJacobian from tools.cpp
    // Use Jacobian instead of H
    Tools tools;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_; 
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  
  else {
    // Laser updates - simply use H_laser since it is linear
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
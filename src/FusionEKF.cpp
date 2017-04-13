#include <cmath>
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
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);

    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

    // acceleration variance for process covariance matrix
    noise_ax_ = 9;
    noise_ay_ = 9;
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

        // first measurement
        cout << "EKF: " << endl;

        // initialize state vector
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        // initialize covariance matrix
        ekf_.P_ = MatrixXd(4,4);
        ekf_.P_ << 1, 0, 0, 0,
                   0, 1, 0, 0,
                   0, 0, 1000, 0,
                   0, 0, 0, 1000;

        // initialize state transition matrix
        ekf_.F_ = MatrixXd(4,4);
        ekf_.F_ << 1, 0, 1, 0,
                   0, 1, 0, 1,
                   0, 0, 1, 0,
                   0, 0, 0, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            double rho = measurement_pack.raw_measurements_(0);   // radial distance
            double phi = measurement_pack.raw_measurements_(1);   // heading
            double rhod = measurement_pack.raw_measurements_(2);  // radial velocity

            ekf_.x_(0) = rho * cos(phi);    // px
            ekf_.x_(1) = rho * sin(phi);    // py
            ekf_.x_(2) = rhod * cos(phi);   // vx
            ekf_.x_(3) = rhod * sin(phi);   // vy
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            ekf_.x_(0) = measurement_pack.raw_measurements_[0];   // px
            ekf_.x_(1) = measurement_pack.raw_measurements_[1];   // py
        }

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;

    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    // Modify the F matrix so that time is integrated
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;

    // Set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4/4*noise_ax_, 0, dt_3/2*noise_ax_, 0,
               0, dt_4/4*noise_ay_, 0, dt_3/2*noise_ay_,
               dt_3/2*noise_ax_, 0, dt_2*noise_ax_, 0,
               0, dt_3/2*noise_ay_, 0, dt_2*noise_ay_;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}

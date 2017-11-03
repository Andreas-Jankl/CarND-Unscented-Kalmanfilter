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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //intialize state size
  n_x_=5;
  n_aug_=7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.16;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.16;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.95;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  //intialize lampda parameter
  lambda_ = 3 - n_aug_;

  //create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.fill(0.0);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  //intialize the weights
  weights_=VectorXd(2 * n_aug_ + 1);
  weights_(0)=lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // Lidar measurement noise covariance matrix
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << pow(std_laspx_, 2), 0,
              0,                  pow(std_laspy_, 2);

  // Radar measurement noise covariance matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << pow(std_radr_, 2), 0,                   0,
              0,                 pow(std_radphi_, 2), 0,
              0,                 0,                   pow(std_radrd_, 2);
}

UKF::~UKF() {}

/*******************************************************************************
 *  UKF overall logic
 ******************************************************************************/

/**
 * Implement the overall UKF logic
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    /*****************************************************************************
    *  Initialization
    ****************************************************************************/
    if (!is_initialized_) {

        cout << "Initialise UKF: " << endl;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

                //Initialize the states for a first measurement with a radar
                //This means converting coordinates from polar to cartesian
                float rho = meas_package.raw_measurements_[0]; // Range - radial distance from origin
                float phi = meas_package.raw_measurements_[1]; // Bearing - angle between rho and x
                float rho_dot = meas_package.raw_measurements_[2]; // Radial Velocity
                float x = rho * cos(phi);
                float y = rho * sin(phi);
                float vx = rho_dot * cos(phi);
                float vy = rho_dot * sin(phi);

                x_ << x, y, sqrt(vx * vx + vy * vy),0.0,0.0;

        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

                //Initialize the states for a first measurement with a laser
                //Just take the position but and leave the others to zero
                x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0;

        }

        //Set the previous timestamp as the initial one. and set the filter to be intialized
        previous_timestamp_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
  }

  if((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_==true) || (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_==true)){

      /*****************************************************************************
       *  Prediction
       ****************************************************************************/

      //compute the time elapsed between the current and previous measurements
      float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
      previous_timestamp_ = meas_package.timestamp_;

      cout << "Timestamp dt = " << dt << endl;


      cout << "Predict" << endl << endl;

      //Predict step in the Kalmanfilter
      Prediction(dt);

      //cout << "State variables x: " << x_ << endl;
      //cout << "State variables P: " << P_ << endl;
      //cout << "State variables Xsig_aug: " << Xsig_aug_ << endl;
      //cout << "State variables Xsig_pred: " << Xsig_pred_ << endl;


      /*****************************************************************************
       *  Update
       ****************************************************************************/


      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

            cout << "Update Radar" << endl << endl;

            // Radar update. Prepare measurement
            VectorXd z = VectorXd(3);
            z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
            UpdateRadar(z);
      } else {

            cout << "Update Lidar" << endl << endl;

            // Laser updates. Prepare measurement
            VectorXd z = VectorXd(2);
            z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
            UpdateLidar(z);
      }
  }
  //print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/*******************************************************************************
 *  Prediction Logic
 ******************************************************************************/

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  GeneratSigmaPoints();
  //cout << "Prediction Step 1" << endl;
  PredictSigmaPoints(delta_t);
  //cout << "Prediction Step 2" << endl;
  PredictStateMeanAndCovarinaceMatrix();
  //cout << "Prediction Step 3" << endl;
}

/**
 * Produces n_sig_ (2 * n_aug_ + 1) sigma points (augmented by process noise)
 * and assigns them to the Xsig_aug_ matrix.
 */
void UKF::GeneratSigmaPoints() {

    //create augmented mean vector
    VectorXd x_aug_ = VectorXd(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

    //create augmented mean state
    x_aug_.head(n_x_) = x_;
    x_aug_(n_x_) = 0;
    x_aug_(n_x_+1) = 0;

    //create augmented covariance matrix
    P_aug_.fill(0.0);
    P_aug_.topLeftCorner(n_x_,n_x_) = P_;
    P_aug_(n_x_,n_x_) = std_a_*std_a_;
    P_aug_(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd L_ = P_aug_.llt().matrixL();

    //create augmented sigma points
    Xsig_aug_.col(0)  = x_aug_;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L_.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L_.col(i);
    }
}

/**
 * Applies the state space model to each column which represents a
 * sigma point in Xsig_aug_ and assigns the results to the Xsig_pred_ matrix.
 * @param delta_t Time between k and k+1 in s
 */
void UKF::PredictSigmaPoints(double delta_t) {

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/**
 * Sets the state vector, x_, to be the weighted mean of the predicted sigma
 * points and the state covariance matrix, P_, as the weighted self similar
 * product (w * A * A^T) of the difference between each predicted sigma point
 * and new state vector x_s. This process completes the prediction step
 * of the UKF.
 */
void UKF::PredictStateMeanAndCovarinaceMatrix() {

    //create vector for predicted state
    VectorXd x_help = VectorXd(n_x_);

    //create covariance matrix for prediction
    MatrixXd P_help = MatrixXd(n_x_, n_x_);

      //predicted state mean
    x_help.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_help = x_help+ weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_help.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_help = P_help + weights_(i) * x_diff * x_diff.transpose() ;
    }
    P_=P_help;
    x_=x_help;
}

/*******************************************************************************
 *  Update Logic
 ******************************************************************************/

/**
* Translates the predicted sigma point matrix into the radar measurement space.
* Produce zpred vector which is the weighted radar measurement on the basis of the sigma points
*/
void UKF::PredictMeasurement_Radar() {
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(3, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y); //r

    // measurement model
    if(p_x==0){
        p_x=1e-10;
        Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y); //r
        cout << "Fehlerbehandlung 2 divide by zero" << endl;
    }


    if(Zsig_(0,i)==0){
        Zsig_(0,i)=1e-10;
        cout << "Fehlerbehandlung 3 divide by zero" << endl;
    }
    Zsig_(1,i) = atan2(p_y,p_x);        //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / Zsig_(0,i);   //r_dot

    }

  //mean predicted measurement
  z_pred_ = VectorXd(3);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  S_ = MatrixXd(3,3);
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  S_ = S_ + R_radar_;
}

/**
* Produce zpred vector which is the weighted radar measurement on the basis of the sigma points
*/
void UKF::PredictMeasurement_Lidar() {
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(2, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig_(0,i) = p_x;                        //p_x
    Zsig_(1,i) = p_y;                        //p_y
  }

  //mean predicted measurement
  z_pred_ = VectorXd(2);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  S_ = MatrixXd(2,2);
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  S_ = S_ + R_lidar_;

}

/**
* Updates the state and the state covariance matrix using a laser measurement
* @param z_lidar_ The measurement at k+1
*/
void UKF::UpdateLidar(VectorXd z_lidar) {

  VectorXd z_ = z_lidar;

  PredictMeasurement_Lidar();

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 2);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred_;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

  NIS_laser_=tools.CalculateNIS(z_diff, S_);
  //cout << "State variables x: " << x_ << endl;
  //cout << "State variables P: " << P_ << endl;
  //cout << "State variables Zsig_: " << Zsig_ << endl;
  //cout << "State variables S_: " << S_ << endl;
  //cout << "State variables Tc: " << Tc << endl;
  //cout << "State variables K: " << K << endl;

}

/**
* Updates the state and the state covariance matrix using a radar measurement
* @param z_radar_ The measurement at k+1
*/
void UKF::UpdateRadar(VectorXd z_radar) {

  VectorXd z_ = z_radar;

  PredictMeasurement_Radar();

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 3);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred_;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

  NIS_radar_=tools.CalculateNIS(z_diff, S_);


  //cout << "State variables Zsig_: " << Zsig_ << endl;
  //cout << "State variables S_: " << S_ << endl;
  //cout << "State variables Tc: " << Tc << endl;
  //cout << "State variables K: " << K << endl;
  //cout << "State variables x: " << x_ << endl;
  //cout << "State variables P: " << P_ << endl;
}

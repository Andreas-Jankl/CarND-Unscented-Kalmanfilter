#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"
#include <stdint.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///*previous timestamp
  int64_t previous_timestamp_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* measurement covariance matrix
  MatrixXd S_;

  ///* augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* Lidar measurement noise covariance matrix
  MatrixXd R_lidar_;

  ///* Radar measurement noise covariance matrix
  MatrixXd R_radar_;

  ///* Vector of predicted measurements
  VectorXd z_pred_;

  ///* Matrix of sigma points in measurements space
  MatrixXd Zsig_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

    ///* The Tool object encapsulates a handful of helper methods
  Tools tools;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

/*****************************************************************************
*  Overal UKF Logic
****************************************************************************/

      /**
      * Implement the overall UKF logic
      * @param {MeasurementPackage} meas_package The latest measurement data of
      * either radar or laser.
      */
      void ProcessMeasurement(MeasurementPackage meas_package);


/*****************************************************************************
*  Prediction Logic
****************************************************************************/

      /**
       * Prediction Predicts sigma points, the state, and the state covariance
       * matrix
       * @param delta_t Time between k and k+1 in s
       */
      void Prediction(double delta_t);

      /**
       * Produces n_sig_ (2 * n_aug_ + 1) sigma points (augmented by process noise)
       * and assigns them to the Xsig_aug_ matrix.
       */
      void GeneratSigmaPoints();

      /**
      * Applies the state space model to each column which represents a
      * sigma point in Xsig_aug_ and assigns the results to the Xsig_pred_ matrix.
      * @param delta_t Time between k and k+1 in s
      */
      void PredictSigmaPoints(double delta_t);

      /**
       * Sets the state vector, x_, to be the weighted mean of the predicted sigma
       * points and the state covariance matrix, P_, as the weighted self similar
       * product (w * A * A^T) of the difference between each predicted sigma point
       * and new state vector x_s. This process completes the prediction step
       * of the UKF.
       */
      void PredictStateMeanAndCovarinaceMatrix();

/*****************************************************************************
*  Update Logic
****************************************************************************/

      /**
       * Updates the state and the state covariance matrix using a laser measurement
       * @param z_lidar_ The measurement at k+1
       */
      void UpdateLidar(VectorXd z_lidar);

      /**
       * Updates the state and the state covariance matrix using a radar measurement
       * @param z_radar_ The measurement at k+1
       */
      void UpdateRadar(VectorXd z_radar);

      /**
       * Translates the predicted sigma point matrix into the radar measurement space.
       * Produce zpred vector which is the weighted radar measurement on the basis of the sigma points
       */
      void PredictMeasurement_Radar();

       /**
       * Produce zpred vector which is the weighted radar measurement on the basis of the sigma points
       */
      void PredictMeasurement_Lidar();
};

#endif /* UKF_H */

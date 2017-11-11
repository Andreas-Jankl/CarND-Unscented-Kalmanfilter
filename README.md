# Unscented Kalman Filter Project Starter Code
Self-Driving Car Engineer Nanodegree Program

This implements an Unscented Kalmanfilter that is able to do a sensor fusion of given sensor readings from radar and lidar object detections.

When its being run and compared against ground truth data its px, py, vx, and vy RMSE is less than [0.09, 0.09, 0.65, 0.65] as required for passing.

---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/sample-laser-radar-measurement-data-1.txt output.txt`

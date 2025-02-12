# Sensor Fusion of LIDAR and RADAR Data using UKF Filter

This project implements an **Unscented Kalman Filter (UKF)** to estimate the state of multiple cars on a highway using noisy **LIDAR** and **RADAR** measurements. The goal is to track and predict the motion of the vehicles in the simulation environment while dealing with sensor noise.

The simulation involves an ego car in the center, surrounded by other traffic cars. The UKF is applied to each of these vehicles, with the purpose of tracking their states (position and velocity) accurately despite noisy sensor data. The performance is evaluated using the **Root Mean Squared Error (RMSE)** metric, and the objective is to achieve RMSE values below a given threshold to pass the project.

## Project Overview

The program is built to simulate a **highway** environment with multiple cars, where:

- The **ego car** (green car) is the car you are tracking.
- **Traffic cars** (blue cars) are other vehicles on the road that will accelerate and change lanes.
- **LIDAR measurements** (red spheres) and **RADAR measurements** (purple lines) are used to track the vehicles.
- The coordinate system is centered around the ego car, and all other cars' positions are relative to it.

The **Unscented Kalman Filter (UKF)** is used to process the noisy LIDAR and RADAR data and estimate the position and velocity of each car in the simulation.

## Key Components

1. **UKF (Unscented Kalman Filter)**: Implemented in the `src/ukf.cpp` and `src/ukf.h` files. The UKF is used to filter and predict the state of each vehicle (position and velocity).
2. **Main Program (`main.cpp`)**: Handles the setup of the simulation, including generating the highway environment and initiating the UKF for each car.
3. **Visualization**: The simulation is visualized using 3D rendering. The red spheres represent LIDAR detections, while the purple lines show RADAR measurements. The cars are visualized in the XY plane.

## Build Instructions

### Prerequisites

Before building the project, ensure that you have the following dependencies installed on your system:

- **CMake** >= 3.5
- **Make** >= 4.1 (Linux/Mac) or 3.81 (Windows)
- **GCC/G++** >= 5.4
- **PCL (Point Cloud Library)** 1.2

### Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./ukf_highway`

## Expected Results

Once the program runs, you should see the simulation of the highway with multiple cars, where:

- The ego car is in the center (green).
- Other traffic cars are in different lanes and will change their velocities and lanes during the simulation.
- The LIDAR and RADAR measurements will be visualized in the scene.

### Performance Evaluation

The performance of the UKF filter is evaluated based on the RMSE (Root Mean Squared Error) metric. The goal is to achieve RMSE values lower than the provided tolerance in the project rubric. This indicates that the filter is successfully predicting the state of the vehicles with minimal error.

## Results

You can watch the simulation results in the following videos:

- **Simulation with LIDAR and RADAR Data:**  
  ![LIDAR and RADAR Data Visualization](media/output_tracking_vis.webm)

- **Simulation with Point Cloud Data (PCD):**  
  ![PCD Visualization](media/output_vis_pcd.webm)


These videos demonstrate how the LIDAR and RADAR measurements are used to track the vehicles and how the UKF filter processes the data to estimate their positions and velocities.

## Acknowledgements

This project was done as part of the **Udacity Sensor Fusion** course. Special thanks to **Udacity** for providing the resources and knowledge to implement this project. The course has been instrumental in helping me understand the concepts of sensor fusion and Kalman filters, which are crucial in autonomous vehicle systems.


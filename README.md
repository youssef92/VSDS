# VSDS

Copyright (C) 2023, by Youssef Michel.

The project provides the C++ implementation of the VSDS approach which aims to generate a force field
based on a first order DS and a desired stiffness profile in a closed loop configuration. The streamlines of the VSDS controller are illustrated below

<p float="left">
  <img src="images/vsds_1.png" width="400" />
  <img src="images/vsds_2.png" width="400" /> 
</p>

For more details, please check:

> [1] Y. Michel, M. Saveriano, and D. Lee, "A Passivity-Based Approach for Variable Stiffness Control With Dynamical Systems," IEEE Transactions on Automation Science and Engineering (TASE), 2023.

The main scripts in this folder are:

- `VSDS_base.cpp`: contains the base class implementation of the VSDS algorithm for the translation case. The `/config` folder contains the parameters used by VSDS such as the via-points definition, initial and goal orientation.

- `VSDS_qp.cpp`, `VSDS_vf.cpp`, `VSDS_org.cpp` provide the implementation of the specfic VSDS approach chosen. For more details check the above paper
 to the robot. Override the implementation based on your speficic robot.

- `VSDS_main.cpp`: This scipt implements the main control loop of the VSDS control law, handling the communication with a Kuka Robot via FRI, recieving sensor
              readings and sending force commands to the robot





To launch the nodes, use the `VSDS.launch`, for instance:
```
roslaunch VSDS VSDS.launch VSDS_name:=Trapezoid_DS VS_type:=qp Stiff_type:=v
```






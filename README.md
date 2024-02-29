# FE-GUT: Factor Graph Optimization hybrid with Extended Kalman Filter for tightly coupled GNSS/UWB Integration
We open-source the source code and simulation dataset of a bridging Factor Graph Optimization (FGO) and Extended Kalman Filter (EKF) architecture for tightly coupled GNSS/UWB integration with online temporal calibration. The main contributions are as follows:
* We introduce the **[Graphical State Space Model (GSSM)](https://github.com/shaolinbit/GraphicalStateSpaceModel)** which is a novel discretization method into the tightly coupled GNSS/UWB integration. 
* An architecture in which FGO is hybrid with EKF is designed to effectively leverage the respective merits of both methods.
## How to Use This Library
The library was built and tested in Ubuntu 20.04. The results may be different for different OS. If any problem was found by you, please propose an issue or report it to zhaoqj23@mails.tsinghua.edu.cn.
### Requirements
1) Ubuntu 20.04 with the newest compiler is recommended.
2) Eigen3
3) Ceres Solver
### Clone the repository
```git clone https://github.com/zhaoqj23/GNSS-UWB.git```
### Build the library
```
cd ~/GNSS-UWB
mkdir build && cd build
cmake ..
make -j8
```
### Run demo
```
cd ~/GNSS-UWB
./bin/gnss_uwb
```
## Acknowledgements
We would like to acknowledge Dr. Xiaoji Niu and the Integrated and Intelligent Navigation (i2Nav) group from Wuhan University for providing the OB_GINS software that was used in the library. We would also like to acknowledge Mr. Yihan Guo for his help in this project.

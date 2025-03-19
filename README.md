# FE-GUT: Factor Graph Optimization hybrid with Extended Kalman Filter for tightly coupled GNSS/UWB Integration
We open-source the source code and simulation dataset of a novel architecture in which the Factor Graph Optimization (FGO) is hybrid with the Extended Kalman Filter (EKF) for tightly coupled GNSS/UWB integration with online Temporal calibration (FE-GUT). The main contributions are as follows:
* We introduce the **[Graphical State Space Model (GSSM)](https://github.com/shaolinbit/GraphicalStateSpaceModel)** which is a novel discretization method into the tightly coupled GNSS/UWB integration. 
* An architecture in which FGO is hybrid with EKF is designed to effectively leverage the respective merits of both methods.
Furthermore, we employ a nonlinearity analysis approach to investigate the inherent reasons why GSSM can improve state estimation accuracy. We then validate its effectiveness through simulations in two case studies: radar tracking and GNSS/UWB tightly coupled integration.
## How to Use This Library
The library was built and tested in Ubuntu 20.04. The results may be different for different OS. If any problem was found by you, please propose an issue or report it to zhaoqj23@mails.tsinghua.edu.cn.
### Requirements
1) Ubuntu 20.04 with the newest compiler is recommended.
2) Eigen3
3) Ceres Solver
### Clone the repository
```git clone https://github.com/zhaoqj23/FE-GUT.git```
### Build the library
```
cd ~/FE-GUT
mkdir build && cd build
cmake ..
make -j8
```
### Run demo
```
cd ~/FE-GUT
./bin/gnss_uwb
```
## Acknowledgements
We would like to acknowledge Dr. Xiaoji Niu and the Integrated and Intelligent Navigation (i2Nav) group from Wuhan University for providing the OB_GINS software that was used in the library. We would also like to acknowledge Mr. Yihan Guo for his help in this project.

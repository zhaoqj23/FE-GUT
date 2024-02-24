#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Geometry>

using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

typedef struct PS {
    double time;
    Eigen::Matrix<double, 6, 1> pse;
} PS;

typedef struct PSRATE {
    double time;
    Eigen::Matrix<double, 6, 1> pserate;
} PSRATE;

typedef struct SATPOS {
    double time;
    Eigen::Matrix<double, 6, 3> satposition;
} SATPOS;

typedef struct SATVEL {
    double time;
    Eigen::Matrix<double, 6, 3> satvelocity;
} SATVEL;


typedef struct UWB {
    double time;
    Eigen::Matrix<double, 4, 1> range;
} UWB;

typedef struct Pose {
    Matrix3d R;
    Vector3d t;
} Pose;

typedef struct IntegrationState {
    double time;

    Eigen::Vector3d r{0, 0, 0};
    Eigen::Vector3d v{0, 0, 0};
    Eigen::Vector3d a{0, 0, 0};

    double tu;
    double fu;
    double tdk;

    Eigen::Matrix<double, 12, 12> P;

} IntegrationState;

typedef struct IntegrationStateData {
    double time;

    double xdata[12];

    Eigen::Matrix<double, 12, 12> pdata;

} IntegrationStateData;

#endif // TYPES_H

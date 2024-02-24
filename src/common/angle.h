
#include <eigen3/Eigen/Geometry>

#ifndef ANGLE_H
#define ANGLE_H

const double D2R = (M_PI / 180.0);
const double R2D = (180.0 / M_PI);

class Angle {

public:
    static double rad2deg(double rad) {
        return rad * R2D;
    }

    static double deg2rad(double deg) {
        return deg * D2R;
    }

    static float rad2deg(float rad) {
        return rad * R2D;
    }

    static float deg2rad(float deg) {
        return deg * D2R;
    }

    template <typename T, int Rows, int Cols>
    static Eigen::Matrix<T, Rows, Cols> rad2deg(const Eigen::Matrix<T, Rows, Cols> &array) {
        return array * R2D;
    }

    template <typename T, int Rows, int Cols>
    static Eigen::Matrix<T, Rows, Cols> deg2rad(const Eigen::Matrix<T, Rows, Cols> &array) {
        return array * D2R;
    }
};

#endif // ANGLE_H

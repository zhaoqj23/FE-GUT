#ifndef PREDICT_FACTOR_H
#define PREDICT_FACTOR_H

#include <Eigen/Geometry>
#include <ceres/ceres.h>

#include "src/common/rotation.h"
#include "src/common/types.h"
#include "src/integration/parameter.h"

class PredictFactor : public ceres::CostFunction {

public:
    explicit PredictFactor(IntegrationStateData statedata)
        : statedata_(std::move(statedata)) {

        *mutable_parameter_block_sizes() = std::vector<int>{11,11};
        set_num_residuals(11);

    }


    bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override {
        Eigen::Matrix<double,11,1> state0;
        Eigen::Matrix<double,11,1> state1;
        for(int i=0;i<11;i++) {
            state0(i,0) = parameters[0][i];
            state1(i,0) = parameters[1][i];
        }

        Eigen::Map<Eigen::Matrix<double, 11, 1>> error(residuals);

        error = state1 - Fg.block(0, 0, 11, 11)*state0;

        if (jacobians) {
            if (jacobians[0]) {
                Eigen::Map<Eigen::Matrix<double, 11, 11, Eigen::RowMajor>> jacobian_s1(jacobians[0]);
                jacobian_s1.setZero();
                jacobian_s1 = -Fg.block(0, 0, 11, 11);
            }
            if (jacobians[1]) {
                Eigen::Map<Eigen::Matrix<double, 11, 11, Eigen::ColMajor>> jacobian_s2(jacobians[1]);
                jacobian_s2.setIdentity();
            }
        }

        return true;
    }

private:
    IntegrationStateData statedata_;
};

#endif // PREDICT_FACTOR_H

#ifndef TC_FACTOR_H
#define TC_FACTOR_H

#include <Eigen/Geometry>
#include <ceres/ceres.h>

#include "../common/rotation.h"
#include "../common/types.h"
#include "../integration/parameter.h"

class TcFactor : public ceres::CostFunction {

public:
    explicit TcFactor(PS ps, PSRATE psrate, SATPOS satpos, UWB uwb, Eigen::Matrix<double,4,3> UWB_base_ecef)
        : ps_(std::move(ps))
        , psrate_(std::move(psrate))
        , satpos_(std::move(satpos))
        , uwb_(std::move(uwb))
        , UWB_base_ecef_(std::move(UWB_base_ecef)) {

        *mutable_parameter_block_sizes() = std::vector<int>{11,1};
        set_num_residuals(16);

    }


    bool Evaluate(const double *const *parameters, double *residuals, double **jacobians) const override {
        Vector3d r{parameters[0][0], parameters[0][1], parameters[0][2]};
        Vector3d v{parameters[0][3], parameters[0][4], parameters[0][5]};
        Vector3d a{parameters[0][6], parameters[0][7], parameters[0][8]};
        double tu = parameters[0][9];
        double fu = parameters[0][10];
        double td = parameters[1][0];


        Eigen::Map<Eigen::Matrix<double, 16, 1>> error(residuals);

        Eigen::Matrix<double,16,1> measurements;
        for(int i=0;i<6;i++) {
            measurements(i,0) = ps_.pse(i,0);
            measurements(i+6,0) = psrate_.pserate(i,0);
        }
        for(int i=0;i<4;i++) {
            measurements(i+12,0) = uwb_.range(i,0);
        }


        Eigen::Matrix<double,16,1> mea_predict;
        Eigen::Matrix<double,6,3> HGn;

        for(int i=0;i<6;i++) {
            HGn.row(i) = -(satpos_.satposition.row(i)-r.transpose()).normalized();
            mea_predict(i,0) = (satpos_.satposition.row(i) - r.transpose()).norm() + tu;
        }

        Eigen::Matrix<double, 6, 1> temp = HGn*v;
        for(int i=0;i<6;i++) {
            mea_predict(i+6,0) = temp(i,0) + fu;
        }


        Eigen::Matrix<double,4,3> HUm;
        Vector3d dr = r - v*td - 0.5*a*td*td;;
        for(int i=0;i<4;i++) {
            HUm.row(i) = -(UWB_base_ecef_.row(i) - dr.transpose()).normalized();
            mea_predict(12+i,0) = (UWB_base_ecef_.row(i) - dr.transpose()).norm();
        }

        Vector3d tem = v + a*td;
        Eigen::Matrix<double,4,1> Htdk;
        for(int i=0;i < 4;i++)
        {
            Htdk(i,0) = HUm(i,0)*tem[0] + HUm(i,1)*tem[1] + HUm(i,2)*tem[2];
            Htdk(i,0) = -Htdk(i,0);
        }
        
        error = mea_predict - measurements;
        error = R_gssm*error;

        if (jacobians) {
            if (jacobians[0]) {
                Eigen::Map<Eigen::Matrix<double, 16, 11, Eigen::RowMajor>> jacobian_state(jacobians[0]);
                Eigen::VectorXd ones = Eigen::VectorXd::Constant(6, 1);
                jacobian_state.setZero();
                jacobian_state.block(0,0,6,3) = HGn;
                jacobian_state.block(0,9,6,1) = ones;
                jacobian_state.block(6,3,6,3) = HGn;
                jacobian_state.block(6,10,6,1) = ones;
                jacobian_state.block(12,0,4,3) = HUm;
                jacobian_state.block(12,3,4,3) = -td*HUm;
                jacobian_state.block(12,6,4,3) = -0.5*td*td*HUm;
                jacobian_state = R_gssm*jacobian_state;
            }
            if (jacobians[1]) {
                Eigen::Map<Eigen::Matrix<double, 16, 1, Eigen::ColMajor>> jacobian_td(jacobians[1]);
                jacobian_td.setZero();
                for(int i=0;i<4;i++) {
                    jacobian_td(i+12,0) = Htdk(i,0);
                }
                jacobian_td = R_gssm*jacobian_td;
            }
        }

        return true;
    }

private:
    PS ps_;
    PSRATE psrate_;
    SATPOS satpos_;
    UWB uwb_;
    Eigen::Matrix<double,4,3> UWB_base_ecef_;
};

#endif // GNSS_FACTOR_H

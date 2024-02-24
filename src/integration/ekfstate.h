#ifndef EKFSTATE_H
#define EKFSTATE_H




#include <Eigen/Dense>
#include <iostream>
#include "src/integration/integration_state.h"
#include "src/integration/parameter.h"
#include "src/common/types.h"


class EkfState {

public:
    EkfState(IntegrationStateData ekfstate);
    void ekfupdate(PS ps, PSRATE psrate, UWB uwb, SATPOS satpos, SATVEL satvel, Eigen::Matrix<double,4,3> UWB_base_ecef);
    void measure_construct(PS ps, PSRATE psrate, UWB uwb);
    void measure_predict(PS ps, PSRATE psrate, UWB uwb, SATPOS satpos, SATVEL satvel, Eigen::Matrix<double,4,3> UWB_base_ecef);
    void ekfuwbupdate(UWB uwb, Eigen::Matrix<double,4,3> UWB_base_ecef);
    void state_init();
    IntegrationState state_construct();
    void dataprint();

public:
    int ynum_ = 2 * gnssnum + uwbnum;
    Eigen::MatrixXd mea_vector_, mea_predict_, H_, xdata_, pdata_, mea_uwb_, uwb_predict_, H_uwb_;
    IntegrationStateData statedatapre_, statedatacur_;
};

#endif
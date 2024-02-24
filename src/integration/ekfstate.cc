#include "ekfstate.h"

EkfState::EkfState(IntegrationStateData ekfstate)
    : statedatapre_(std::move(ekfstate)) {
}

void EkfState::ekfupdate(PS ps, PSRATE psrate, UWB uwb, SATPOS satpos, SATVEL satvel, Eigen::Matrix<double,4,3> UWB_base_ecef) {
    
    for (size_t i = 0; i < 12; i++)
    {
        xdata_(i,0) = statedatapre_.xdata[i];
    }
    pdata_ = statedatapre_.pdata;

    // Prior prediction
    xdata_ = F*xdata_;
    pdata_ = F*pdata_*F.transpose() + Q;

    // Measurement prediction
    measure_predict(ps, psrate, uwb, satpos, satvel, UWB_base_ecef);

    // Posterior update
    Eigen::MatrixXd K = pdata_*H_.transpose()*(H_*pdata_*H_.transpose() + R).inverse();
    xdata_ = xdata_ + K*(mea_vector_ - mea_predict_);
    pdata_ = (Eigen::MatrixXd::Identity(statenum,statenum) - K*H_)*pdata_;

    // Time update
    // The prestate need to be updated lately
    statedatacur_.time = statedatapre_.time + ddt_;

    // Matrix assignment
    for (size_t i = 0; i < 12; i++)
    {
        statedatacur_.xdata[i] = xdata_(i,0);
    }
    statedatacur_.pdata = pdata_;
    statedatapre_ = statedatacur_;
    // dataprint();
    return;
}
void EkfState::dataprint() {
    Eigen::IOFormat format(7);
    std::cout << "xdata_: " << xdata_.transpose().format(format) << std::endl;
    std::cout << "pdata_: " << pdata_.format(format) << std::endl;
    std::cout << "H_: " << H_.format(format) << std::endl;
    std::cout << "mea_vector_: " << mea_vector_.transpose().format(format) << std::endl;
    std::cout << "mea_predict_: " << mea_predict_.transpose().format(format) << std::endl;
    std::cout << "Q:" << Q.format(format) << std::endl;
    std::cout << "R:" << R.format(format) << std::endl;
}

void EkfState::measure_predict(PS ps, PSRATE psrate, UWB uwb, SATPOS satpos, SATVEL satvel, Eigen::Matrix<double,4,3> UWB_base_ecef) {
    // Measurement construction
    measure_construct(ps, psrate, uwb);
    
    // Measurement prediction
    IntegrationState tem = state_construct();
    Eigen::Matrix<double,6,3> HGn;
    Eigen::Matrix<double,4,3> HUm;
    Eigen::Matrix<double,4,1> Htdk;
    Vector3d dr = tem.r - tem.tdk*tem.v - 0.5*tem.a*tem.tdk*tem.tdk;

    // GNSS Pseudoranges
    for(int i=0;i<6;i++) {
        HGn.row(i) = -(satpos.satposition.row(i)-tem.r.transpose()).normalized();
        mea_predict_(i,0) = (satpos.satposition.row(i) - tem.r.transpose()).norm() + tem.tu;
    }

    // GNSS Doppler-shift
    Eigen::Matrix<double, 6, 1> vg = HGn*tem.v;
    for(int i=0;i<6;i++) {
        mea_predict_(i+6,0) = vg(i,0) + tem.fu;
    }

    // UWB Ranges
    for(int i=0;i<4;i++) {
        HUm.row(i) = -(UWB_base_ecef.row(i) - dr.transpose()).normalized();
        mea_predict_(12+i,0) = (UWB_base_ecef.row(i) - dr.transpose()).norm();
    }

    Vector3d vu = tem.v + tem.a*tem.tdk;
    for(int i=0;i < 4;i++)
    {
        Htdk(i,0) = HUm(i,0)*vu[0] + HUm(i,1)*vu[1] + HUm(i,2)*vu[2];
        Htdk(i,0) = -Htdk(i,0);
    }

    // Construct measurement matrix
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(gnssnum,1);
    H_.block(0,0,gnssnum,3) = HGn;
    H_.block(0,9,gnssnum,1) = ones;
    H_.block(gnssnum,3,gnssnum,3) = HGn;
    H_.block(gnssnum,10,gnssnum,1) = ones;
    H_.block(2*gnssnum,0,uwbnum,3) = HUm;
    H_.block(2*gnssnum,3,uwbnum,3) = -tem.tdk*HUm;
    H_.block(2*gnssnum,6,uwbnum,3) = -0.5*tem.tdk*tem.tdk*HUm;
    // H_.block(2*gnssnum,11,uwbnum,1) = Htdk;
}

void EkfState::measure_construct(PS ps, PSRATE psrate, UWB uwb) {
    // Measurement vector
    for(int i = 0; i < gnssnum; i++) {
        mea_vector_(i,0) = ps.pse(i,0);
        mea_vector_(i+gnssnum,0) = psrate.pserate(i,0);
    }
    for(int i = 0; i < uwbnum; i++) {
        mea_vector_(i+2*gnssnum,0) = uwb.range(i,0);
    }
}

void EkfState::state_init() {
    mea_vector_ = Eigen::MatrixXd::Zero(ynum_,1);
    mea_predict_ = Eigen::MatrixXd::Zero(ynum_,1);
    H_ = Eigen::MatrixXd::Zero(ynum_,statenum);
    xdata_ = Eigen::MatrixXd::Zero(statenum,1);
    pdata_ = Eigen::MatrixXd::Zero(statenum,statenum);
    mea_uwb_ = Eigen::MatrixXd::Zero(uwbnum,1);
    uwb_predict_ = Eigen::MatrixXd::Zero(uwbnum,1);
    H_uwb_ = Eigen::MatrixXd::Zero(uwbnum,statenum);
    for(int i = 0; i < statenum; i++) {
        xdata_(i,0) = statedatapre_.xdata[i];
    }
    pdata_ = statedatapre_.pdata;
}

IntegrationState EkfState::state_construct() {
    IntegrationState state;
    for (int i = 0; i < 3; i++)
    {
        state.r(i,0) = xdata_(i,0);
        state.v(i,0) = xdata_(i+3,0);
        state.a(i,0) = xdata_(i+6,0);
    }
    state.tu = xdata_(9,0);
    state.fu = xdata_(10,0);
    state.tdk = xdata_(11,0);
    return state;
}

void EkfState::ekfuwbupdate(UWB uwb, Eigen::Matrix<double,4,3> UWB_base_ecef) {
    for (size_t i = 0; i < 12; i++)
    {
        xdata_(i,0) = statedatapre_.xdata[i];
    }
    pdata_ = statedatapre_.pdata;

    // Prior prediction
    xdata_ = F*xdata_;
    pdata_ = F*pdata_*F.transpose() + Q;

    // Measurement construction
    for(int i = 0; i < uwbnum; i++) {
        mea_uwb_(i,0) = uwb.range(i,0);
    }

    // Measurement prediction
    IntegrationState tem = state_construct();
    Eigen::Matrix<double,4,3> HUm;
    Eigen::Matrix<double,4,1> Htdk;
    Vector3d dr = tem.r - tem.tdk*tem.v - 0.5*tem.a*tem.tdk*tem.tdk;

    // UWB Ranges
    for(int i=0;i<4;i++) {
        HUm.row(i) = -(UWB_base_ecef.row(i) - dr.transpose()).normalized();
        uwb_predict_(i,0) = (UWB_base_ecef.row(i) - dr.transpose()).norm();
    }

    Vector3d vu = tem.v + tem.a*tem.tdk;
    for(int i=0;i < 4;i++)
    {
        Htdk(i,0) = HUm(i,0)*vu[0] + HUm(i,1)*vu[1] + HUm(i,2)*vu[2];
        Htdk(i,0) = -Htdk(i,0);
    }

    // Construct measurement matrix
    H_uwb_.block(0,0,uwbnum,3) = HUm;
    H_uwb_.block(0,3,uwbnum,3) = -tem.tdk*HUm;
    H_uwb_.block(0,6,uwbnum,3) = -0.5*tem.tdk*tem.tdk*HUm;
    // H_uwb_.block(0,11,uwbnum,1) = Htdk;

    // Posterior update
    Eigen::MatrixXd K = pdata_*H_uwb_.transpose()*(H_uwb_*pdata_*H_uwb_.transpose() + R_uwb).inverse();
    xdata_ = xdata_ + K*(mea_uwb_ - uwb_predict_);
    pdata_ = (Eigen::MatrixXd::Identity(statenum,statenum) - K*H_uwb_)*pdata_;

    // Time update
    // The prestate need to be updated lately
    statedatacur_.time = statedatapre_.time + ddt_;

    // Matrix assignment
    for (size_t i = 0; i < 12; i++)
    {
        statedatacur_.xdata[i] = xdata_(i,0);
    }
    statedatacur_.pdata = pdata_;
    statedatapre_ = statedatacur_;
    // dataprint();
    return;
}
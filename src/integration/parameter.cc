#include "src/integration/parameter.h"


double ddt_, ddt_gnss_;
int gnssnum = 6;
int uwbnum = 4;
int statenum = 12;
int uwb_rate = 200, gnss_rate = 10;
double Sj = 0.4, St = 36, Sf = 0.01, Sdt = 0.01;
double  std_ps = 2; double std_psrate = 0.1; double std_uwb = 0.5;
double  std_psg = 2; double std_psrateg = 0.1; double std_uwbg = 0.5;
Eigen::MatrixXd F = Eigen::MatrixXd::Zero(statenum,statenum);
Eigen::MatrixXd Fg = Eigen::MatrixXd::Zero(statenum,statenum);
Eigen::MatrixXd P0 = Eigen::MatrixXd::Identity(statenum,statenum);
Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(statenum,statenum);
Eigen::MatrixXd Q_gssm = Eigen::MatrixXd::Zero(statenum,statenum);
Eigen::MatrixXd R = Eigen::MatrixXd::Zero(2*gnssnum+uwbnum,2*gnssnum+uwbnum);
Eigen::MatrixXd R_gssm = Eigen::MatrixXd::Zero(2*gnssnum+uwbnum,2*gnssnum+uwbnum);
Eigen::MatrixXd R_uwb = Eigen::MatrixXd::Zero(uwbnum,uwbnum);

void parameters_init() {
    ddt_ = 1/(double)uwb_rate;
    ddt_gnss_ = 1/(double)gnss_rate;
    F << 
    1,0,0,ddt_,0,0,0.5*ddt_*ddt_,0,0,0,0,0,
    0,1,0,0,ddt_,0,0,0.5*ddt_*ddt_,0,0,0,0,
    0,0,1,0,0,ddt_,0,0,0.5*ddt_*ddt_,0,0,0,
    0,0,0,1,0,0,ddt_,0,0,0,0,0,
    0,0,0,0,1,0,0,ddt_,0,0,0,0,
    0,0,0,0,0,1,0,0,ddt_,0,0,0,
    0,0,0,0,0,0,1,0,0,0,0,0,
    0,0,0,0,0,0,0,1,0,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,0,
    0,0,0,0,0,0,0,0,0,1,ddt_,0,
    0,0,0,0,0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,0,0,0,0,0,1;
    Fg <<
    1,0,0,ddt_gnss_,0,0,0.5*ddt_gnss_*ddt_gnss_,0,0,0,0,0,
    0,1,0,0,ddt_gnss_,0,0,0.5*ddt_gnss_*ddt_gnss_,0,0,0,0,
    0,0,1,0,0,ddt_gnss_,0,0,0.5*ddt_gnss_*ddt_gnss_,0,0,0,
    0,0,0,1,0,0,ddt_gnss_,0,0,0,0,0,
    0,0,0,0,1,0,0,ddt_gnss_,0,0,0,0,
    0,0,0,0,0,1,0,0,ddt_gnss_,0,0,0,
    0,0,0,0,0,0,1,0,0,0,0,0,
    0,0,0,0,0,0,0,1,0,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,0,
    0,0,0,0,0,0,0,0,0,1,ddt_gnss_,0,
    0,0,0,0,0,0,0,0,0,0,1,0,
    0,0,0,0,0,0,0,0,0,0,0,1;
    Q.block(0,0,3,3) = 1/20*Sj*pow(ddt_,5)*Eigen::Matrix3d::Identity();
    Q.block(3,3,3,3) = 1/3*Sj*pow(ddt_,3)*Eigen::Matrix3d::Identity();
    Q.block(6,6,3,3) = Sj*ddt_*Eigen::Matrix3d::Identity();
    Q.block(3,0,3,3) = 1/8*Sj*pow(ddt_,4)*Eigen::Matrix3d::Identity(); Q.block(0,3,3,3) = 1/8*Sj*pow(ddt_,4)*Eigen::Matrix3d::Identity();
    Q.block(6,3,3,3) = 1/2*Sj*pow(ddt_,2)*Eigen::Matrix3d::Identity(); Q.block(3,6,3,3) = 1/2*Sj*pow(ddt_,2)*Eigen::Matrix3d::Identity();
    Q.block(6,0,3,3) = 1/6*Sj*pow(ddt_,3)*Eigen::Matrix3d::Identity(); Q.block(0,6,3,3) = 1/6*Sj*pow(ddt_,3)*Eigen::Matrix3d::Identity();
    Q(9,9) = St*ddt_+1/3*Sf*pow(ddt_,3);
    Q(9,10) = 1/2*Sf*pow(ddt_,2); Q(10,9) = 1/2*Sf*pow(ddt_,2);
    Q(10,10) = Sf*ddt_; Q(11,11) = Sdt*ddt_;
    Q_gssm.block(0,0,3,3) = 1/20*Sj*pow(ddt_gnss_,5)*Eigen::Matrix3d::Identity();
    Q_gssm.block(3,3,3,3) = 1/3*Sj*pow(ddt_gnss_,3)*Eigen::Matrix3d::Identity();
    Q_gssm.block(6,6,3,3) = Sj*ddt_*Eigen::Matrix3d::Identity();
    Q_gssm.block(3,0,3,3) = 1/8*Sj*pow(ddt_gnss_,4)*Eigen::Matrix3d::Identity(); Q_gssm.block(0,3,3,3) = 1/8*Sj*pow(ddt_gnss_,4)*Eigen::Matrix3d::Identity();
    Q_gssm.block(6,3,3,3) = 1/2*Sj*pow(ddt_gnss_,2)*Eigen::Matrix3d::Identity(); Q_gssm.block(3,6,3,3) = 1/2*Sj*pow(ddt_gnss_,2)*Eigen::Matrix3d::Identity();
    Q_gssm.block(6,0,3,3) = 1/6*Sj*pow(ddt_gnss_,3)*Eigen::Matrix3d::Identity(); Q_gssm.block(0,6,3,3) = 1/6*Sj*pow(ddt_gnss_,3)*Eigen::Matrix3d::Identity();
    Q_gssm(9,9) = St*ddt_gnss_+1/3*Sf*pow(ddt_gnss_,3);
    Q_gssm(9,10) = 1/2*Sf*pow(ddt_gnss_,2); Q_gssm(10,9) = 1/2*Sf*pow(ddt_gnss_,2);
    Q_gssm(10,10) = Sf*ddt_gnss_; Q_gssm(11,11) = Sdt*ddt_gnss_;
    R.block(0,0,gnssnum,gnssnum) = std_ps*std_ps*Eigen::MatrixXd::Identity(gnssnum,gnssnum);
    R.block(gnssnum,gnssnum,gnssnum,gnssnum) = std_psrate*std_psrate*Eigen::MatrixXd::Identity(gnssnum,gnssnum);
    R.block(2*gnssnum,2*gnssnum,uwbnum,uwbnum) = std_uwb*std_uwb*Eigen::MatrixXd::Identity(uwbnum,uwbnum);
    R_gssm.block(0,0,gnssnum,gnssnum) = 1/std_psg*Eigen::MatrixXd::Identity(gnssnum,gnssnum);
    R_gssm.block(gnssnum,gnssnum,gnssnum,gnssnum) = 1/std_psrateg*Eigen::MatrixXd::Identity(gnssnum,gnssnum);
    R_gssm.block(2*gnssnum,2*gnssnum,uwbnum,uwbnum) = 1/std_uwbg*Eigen::MatrixXd::Identity(uwbnum,uwbnum);
    R_uwb = std_uwb*std_uwb*Eigen::MatrixXd::Identity(uwbnum,uwbnum);
}

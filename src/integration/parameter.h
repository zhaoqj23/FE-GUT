#ifndef PARAMETER_H
#define PARAMETER_H



#include <Eigen/Dense>

extern int gnssnum, uwbnum, statenum;
extern double ddt_, Sj, St, Sf, Sdt, ddt_gnss_;
extern Eigen::MatrixXd F, P0, Q, R, R_gssm, Fg, R_uwb, Q_gssm;
extern double ddt_, ddt_gnss_;
extern int gnssnum;
extern int uwbnum;
extern int statenum;
extern int uwb_rate, gnss_rate;







#endif // PARAMETER_H
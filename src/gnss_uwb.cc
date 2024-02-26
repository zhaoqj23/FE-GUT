#include "common/earth.h"
#include "common/angle.h"

#include "fileio/filesaver.h"
#include "fileio/psfileloader.h"
#include "fileio/psratefileloader.h"
#include "fileio/uwbfileloader.h"
#include "fileio/satposfileloader.h"
#include "fileio/satvelfileloader.h"


#include "integration/integration_state.h"
#include "integration/ekfstate.h"
#include "integration/parameter.h"

#include "factors/tc_factor.h"
#include "factors/marginalization_factor.h"
#include "factors/predict_factor.h"

#include <ceres/ceres.h>
#include <absl/strings/str_format.h>
#include <absl/time/clock.h>
#include <iostream>
#include <deque>
#include <iomanip>
#include "glog/logging.h"

#define tole 0.001


void writeResult(double time, IntegrationState state, FileSaver &navfile, FileSaver &tdfile, double td, Vector3d station_origin);
void parameters_init();

int main(int argc, char *argv[]) {

    google::InitGoogleLogging(argv[0]);                // glog init
    google::ParseCommandLineFlags(&argc, &argv, false);

    // running time recorder
    auto ts = absl::Now();

    parameters_init();

    std::shared_ptr<MarginalizationInfo> last_marginalization_info;
    std::vector<double *> last_marginalization_parameter_blocks;

    // TODO: replace it with the real path
    std::string psepath = "/home/zhaoqj23/Documents/GNSS-UWB/dataset/psdata.txt";
    std::string psratepath = "/home/zhaoqj23/Documents/GNSS-UWB/dataset/psratedata.txt";
    std::string satpospath = "/home/zhaoqj23/Documents/GNSS-UWB/dataset/satposdata.txt";
    std::string satvelpath = "/home/zhaoqj23/Documents/GNSS-UWB/dataset/satveldata.txt";
    std::string uwbpath = "/home/zhaoqj23/Documents/GNSS-UWB/dataset/uwbdata.txt";
    
    // file loader and saver
    PSFileLoader psefile(psepath);
    PSRATEFileLoader psratefile(psratepath);
    SATPOSFileLoader satposfile(satpospath);
    SATVELFileLoader satvelfile(satvelpath);
    UWBFileLoader uwbfile(uwbpath);
    FileSaver navfile("/home/zhaoqj23/Documents/GNSS-UWB/dataset/navdata.txt",7,FileSaver::TEXT);
    FileSaver tdfile("/home/zhaoqj23/Documents/GNSS-UWB/dataset/tddata.txt",2,FileSaver::TEXT);

    if (!psefile.isOpen() || !psratefile.isOpen() || !satposfile.isOpen() || !satvelfile.isOpen() || !uwbfile.isOpen() || !navfile.isOpen() || !tdfile.isOpen()) {
        std::cout << "Failed to open data file" << std::endl;
        return -1;
    }

    int datalength = 12000;
    int window = 30;
    int datanow = 0;
    std::vector<IntegrationState> statelist(window);
    std::vector<IntegrationStateData> statedatalist(window);
    std::vector<double> timelist(window);
    std::vector<PS> pslist(window);
    std::vector<PSRATE> psratelist(window);
    std::vector<SATPOS> satposlist(window);
    std::vector<SATVEL> satvellist(window);
    std::vector<UWB> uwblist(window);

    double td = 0;
    double endtime = 1200;
    Eigen::Vector3d recPos0{39.904987,116.405289,60.0352};
    recPos0[0] *= D2R;
    recPos0[1] *= D2R;
    Vector3d station_origin = recPos0;

    PS ps; ps = psefile.next(); pslist[datanow] = ps;
    PSRATE psrate; psrate = psratefile.next(); psratelist[datanow] = psrate;
    SATPOS satpos; satpos = satposfile.next(); satposlist[datanow] = satpos;
    SATVEL satvel; satvel = satvelfile.next(); satvellist[datanow] = satvel;
    UWB uwb; uwb = uwbfile.next(); uwblist[datanow] = uwb;
    Vector3d vini = {5,0,0}; 

    Eigen::Matrix<double,4,3> UWB_base_ned;
    UWB_base_ned << 0,150,-5,50,100,-5,0,50,-5,-50,100,-5;
    
    Eigen::Matrix<double,4,3> UWB_base_ecef;
    for(int i=0;i<4;i++){
        Vector3d blh = Earth::local2global(station_origin,UWB_base_ned.row(i));
        UWB_base_ecef.row(i) = Earth::blh2ecef(blh);
    }

    IntegrationState state_curr = {
        .time = ps.time,
        .r = Earth::blh2ecef(recPos0),
        .v = Earth::cne(recPos0)*vini,//cne is the transformation matrix from n frame to e frame!
        .a = {0,0,0},
        .tu = 0,
        .fu = 0,
        .P  = P0
    };

    IntegrationStateData statedata_curr = stateToData(state_curr);
    writeResult(state_curr.time,state_curr,navfile,tdfile,td,station_origin);

    // Initialization
    EkfState ekfx(statedata_curr);
    ekfx.state_init();
    statedatalist[datanow] = ekfx.statedatapre_;
    statelist[datanow] = stateFromData(ekfx.statedatapre_);
    datanow++;
    double sow = 0;

    for(int i = 0 ; i < datalength; i++) {

        while(1) {
            uwb = uwbfile.next();
            ekfx.ekfuwbupdate(uwb, UWB_base_ecef);
            sow += ddt_;
            if((fabs(uwb.time - ps.time - ddt_gnss_ + ddt_) < tole) || uwbfile.isEof()) {
                uwb = uwbfile.next();
                uwblist[datanow] = uwb;
                break;
            }
        }

        if((sow >= endtime)  || uwbfile.isEof()) {
            break;
        }

        ps = psefile.next(); pslist[datanow] = ps;
        psrate = psratefile.next(); psratelist[datanow] = psrate;
        satpos = satposfile.next(); satposlist[datanow] = satpos;
        satvel = satvelfile.next(); satvellist[datanow] = satvel;
        timelist[datanow] = ps.time;

        

        // Initialized by EKF
        ekfx.ekfupdate(ps,psrate,uwb,satpos,satvel,UWB_base_ecef);
        statedatalist[datanow] = ekfx.statedatacur_;
        
        datanow++;

        // Construct optimiation problem
        ceres::Problem::Options problem_options;
        problem_options.enable_fast_removal = true;

        ceres::Problem problem(problem_options);
        ceres::Solver solver;
        ceres::Solver::Summary summary;
        ceres::Solver::Options options;
        options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
        options.linear_solver_type         = ceres::SPARSE_NORMAL_CHOLESKY;
        options.num_threads                = 8;
        options.use_nonmonotonic_steps     = true;
        options.dynamic_sparsity           = true;

        // Add the first 11 variables in statedatalist[k].xdata to the problem as parameters
        for (int i = 0; i < datanow; i++)
        {
            problem.AddParameterBlock(statedatalist[i].xdata,11);
        }

        problem.AddParameterBlock(&td,1);

        // measurements factor: gnss/uwb
        for(int k = 0; k < datanow; k++) {
            auto factor = new TcFactor(pslist[k], psratelist[k], satposlist[k], uwblist[k], UWB_base_ecef);
            problem.AddResidualBlock(factor,nullptr,statedatalist[k].xdata,&td);
        }

        // predict factor
        // ceres::LossFunction *loss_function = new ceres::HuberLoss(1.0);
        for(int k = 0; k < datanow-1; k++) {
            auto factor = new PredictFactor(statedatalist[k]);
            problem.AddResidualBlock(factor, nullptr, statedatalist[k].xdata, statedatalist[k+1].xdata);
        }

        // prior factor
        if (last_marginalization_info && last_marginalization_info->isValid()) {
            auto factor = new MarginalizationFactor(last_marginalization_info);
            problem.AddResidualBlock(factor, nullptr, last_marginalization_parameter_blocks);
        }

        options.max_num_iterations = 200;
        if(i > 50) solver.Solve(options, &problem, &summary);
        // std::cout << summary.BriefReport() << "\n";


        int percent = ((int)ps.time)*100/endtime;
        static int lastpercent = 0;
        if (abs(percent - lastpercent) >= 1) {
                lastpercent = percent;
                std::cout << "Percentage: " << std::setw(3) << percent << "%\r";
                flush(std::cout);
        }

        state_curr = stateFromData(ekfx.statedatacur_);
        statelist[datanow-1] = state_curr;
        // ekfx.statedatapre_ = statedatalist[datanow-1];
        ekfx.statedatapre_.xdata[11] = td;
        writeResult(state_curr.time,state_curr,navfile,tdfile,td,station_origin);

        if(datanow >= window) {
            {
                // marginalization
                std::shared_ptr<MarginalizationInfo> marginalization_info = std::make_shared<MarginalizationInfo>();
                if (last_marginalization_info && last_marginalization_info->isValid()) {
                    std::vector<int> marginilized_index;
                    for (size_t k = 0; k < last_marginalization_parameter_blocks.size(); k++) {
                        if (last_marginalization_parameter_blocks[k] == statedatalist[0].xdata) {
                            marginilized_index.push_back(static_cast<int>(k));
                        }
                    }

                    auto factor   = std::make_shared<MarginalizationFactor>(last_marginalization_info);
                    auto residual = std::make_shared<ResidualBlockInfo>(
                        factor, nullptr, last_marginalization_parameter_blocks, marginilized_index);
                    marginalization_info->addResidualBlockInfo(residual);
                }

                {
                    auto factor   = std::make_shared<TcFactor>(pslist[0], psratelist[0], satposlist[0], uwblist[0], UWB_base_ecef);
                    auto residual = std::make_shared<ResidualBlockInfo>(
                        factor, nullptr, 
                        std::vector<double *>{statedatalist[0].xdata,&td}, std::vector<int>{0});
                    marginalization_info->addResidualBlockInfo(residual);
                }

                {
                    auto factor = std::make_shared<PredictFactor>(statedatalist[0]);
                    auto residual = std::make_shared<ResidualBlockInfo>(
                        factor, nullptr, 
                        std::vector<double *>{statedatalist[0].xdata, statedatalist[1].xdata}, std::vector<int>{0});
                    marginalization_info->addResidualBlockInfo(residual);
                }

                // do marginalization
                marginalization_info->marginalization();

                // get new pointers
                std::unordered_map<long, double *> address;
                for (int k = 1; k < window; k++) {
                    address[reinterpret_cast<long>(statedatalist[k].xdata)] = statedatalist[k - 1].xdata;
                }
                address[reinterpret_cast<long>(&td)] = &td;
                last_marginalization_parameter_blocks = marginalization_info->getParamterBlocks(address);
                last_marginalization_info             = std::move(marginalization_info);
            }

            // sliding window
            {
                for (int k = 0; k < window - 1; k++) {
                    statedatalist[k] = statedatalist[k + 1];
                    statelist[k]     = statelist[k + 1];
                    timelist[k] = timelist[k + 1];
                    pslist[k] = pslist[k + 1];;
                    psratelist[k] = psratelist[k + 1];
                    satposlist[k] = satposlist[k + 1];
                    satvellist[k] = satvellist[k + 1];
                    uwblist[k] = uwblist[k + 1];
                }
                datanow--;
            }
        }

    }
    
    // std::cout << "Endtime is : " << timelist[datalength] << " s" << std::endl;

    navfile.close();
    psefile.close();
    psratefile.close();
    satposfile.close();
    satvelfile.close();
    uwbfile.close();
    tdfile.close();

    auto te = absl::Now();
    std::cout << std::endl << std::endl << "Cost " << absl::ToDoubleSeconds(te - ts) << " s in total" << std::endl;
    
    return 0;

}


void writeResult(double time,IntegrationState state,FileSaver &navfile,FileSaver &tdfile, double td, Vector3d station_origin) {
    state.r = Earth::ecef2blh(state.r);
    state.r = Earth::global2local(station_origin,state.r);
    state.v = Earth::cne(station_origin).transpose()*state.v;
    vector<double> result;
    {
        result.clear();
        result.push_back(time);
        result.push_back(state.r[0]);
        result.push_back(state.r[1]);
        result.push_back(state.r[2]);
        result.push_back(state.v[0]);
        result.push_back(state.v[1]);
        result.push_back(state.v[2]);
        navfile.dump(result);
    }
    {
        result.clear();
        result.push_back(time);
        result.push_back(td*1000);
        tdfile.dump(result);
    }
}



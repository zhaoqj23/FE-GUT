#ifndef LOGGING_H
#define LOGGING_H

#include <Eigen/Geometry>
#include <absl/strings/str_format.h>
#include <glog/logging.h>
#include <glog/stl_logging.h>
#include <iostream>

using std::string;

#define LOGI (LOG(INFO))
#define LOGW (LOG(WARNING))
#define LOGE (LOG(ERROR))
#define LOGF (LOG(FATAL))

#if !DCHECK_IS_ON()
#define DLOGI (static_cast<void>(0), true ? (void) 0 : google::LogMessageVoidify() & LOG(INFO))
#define DLOGW (static_cast<void>(0), true ? (void) 0 : google::LogMessageVoidify() & LOG(WARNING))
#define DLOGE (static_cast<void>(0), true ? (void) 0 : google::LogMessageVoidify() & LOG(ERROR))
#define DLOGF (static_cast<void>(0), true ? (void) 0 : google::LogMessageVoidify() & LOG(FATAL))
#else
#define DLOGI LOGI
#define DLOGW LOGW
#define DLOGE LOGE
#define DLOGF LOGF
#endif

class Logging {

public:
    static void initialization(char **argv, bool logtostderr = true, bool logtofile = false) {
        if (logtostderr & logtofile) {
            FLAGS_alsologtostderr = true;
        } else if (logtostderr) {
            FLAGS_logtostderr = true;
        }

        if (logtostderr) {
            // glog init
            FLAGS_colorlogtostderr = true;
        }
        
        google::InitGoogleLogging(argv[0]);
    }

    template <typename T, int Rows, int Cols>
    static void printMatrix(const Eigen::Matrix<T, Rows, Cols> &matrix, const string &prefix = "Matrix: ") {
        std::cout << prefix << matrix.rows() << "x" << matrix.cols() << std::endl;
        if (matrix.cols() == 1) {
            std::cout << matrix.transpose() << std::endl;
        } else {
            std::cout << matrix << std::endl;
        }
    }

    static string doubleData(double data) {
        return absl::StrFormat("%0.6lf", data);
    }
};

#endif // LOGGING_H

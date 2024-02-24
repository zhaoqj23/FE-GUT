#ifndef SATPOSFILELOADER_H
#define SATPOSFILELOADER_H

#include "fileloader.h"
#include "../common/types.h"

class SATPOSFileLoader : public FileLoader {

public:
    SATPOSFileLoader() = delete;
    explicit SATPOSFileLoader(const string &filename, int columns = 19) {
        open(filename, columns, FileLoader::TEXT);
    }

    const SATPOS &next() {
        data_ = load();
        satpos_.time = data_[0];
        for(int i = 0;i < 6;i++) {
            for(int j=0;j<3;j++) {
                satpos_.satposition(i,j) = data_[1+3*i+j];
            }
        }
        return satpos_;
    }

private:
    SATPOS satpos_;
    vector<double> data_;
};

#endif // SATPOSFILELOADER_H

#ifndef SATVELFILELOADER_H
#define SATVELFILELOADER_H

#include "fileloader.h"
#include "../common/types.h"

class SATVELFileLoader : public FileLoader {

public:
    SATVELFileLoader() = delete;
    explicit SATVELFileLoader(const string &filename, int columns = 19) {
        open(filename, columns, FileLoader::TEXT);
    }

    const SATVEL &next() {
        data_ = load();
        satvel_.time = data_[0];
        for(int i = 0;i < 6;i++) {
            for(int j=0;j<3;j++) {
                satvel_.satvelocity(i,j) = data_[1+3*i+j];
            }
        }
        return satvel_;
    }

private:
    SATVEL satvel_;
    vector<double> data_;
};

#endif // SATVELFILELOADER_H

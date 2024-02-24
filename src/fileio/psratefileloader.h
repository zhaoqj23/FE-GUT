#ifndef PSRATEFILELOADER_H
#define PSRATEFILELOADER_H

#include "fileloader.h"
#include "../common/types.h"

class PSRATEFileLoader : public FileLoader {

public:
    PSRATEFileLoader() = delete;
    explicit PSRATEFileLoader(const string &filename, int columns = 7) {
        open(filename, columns, FileLoader::TEXT);
    }

    const PSRATE &next() {
        data_ = load();
        psrate_.time = data_[0];
        for(int i=0;i<6;i++) {
            psrate_.pserate(i,0) = data_[1+i];
        }
        return psrate_;
    }

private:
    PSRATE psrate_;
    vector<double> data_;
};

#endif // PSRATEFILELOADER_H

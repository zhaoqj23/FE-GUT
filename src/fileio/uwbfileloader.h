#ifndef UWBFILELOADER_H
#define UWBFILELOADER_H

#include "fileloader.h"
#include "../common/types.h"

class UWBFileLoader : public FileLoader {

public:
    UWBFileLoader() = delete;
    explicit UWBFileLoader(const string &filename, int columns = 5) {
        open(filename, columns, FileLoader::TEXT);
    }

    const UWB &next() {
        data_ = load();
        uwb_.time = data_[0];
        for(int i = 0;i<4;i++) {
            uwb_.range(i,0) = data_[1+i];
        }
        return uwb_;
    }

private:
    UWB uwb_;
    vector<double> data_;
};

#endif // PSFILELOADER_H

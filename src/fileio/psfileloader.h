#ifndef PSFILELOADER_H
#define PSFILELOADER_H

#include "fileloader.h"
#include "../common/types.h"

class PSFileLoader : public FileLoader {

public:
    PSFileLoader() = delete;
    explicit PSFileLoader(const string &filename, int columns = 7) {
        open(filename, columns, FileLoader::TEXT);
    }

    const PS &next() {
        data_ = load();
        ps_.time = data_[0];
        for(int i = 0;i<6;i++) {
            ps_.pse(i,0) = data_[1+i];
        }
        return ps_;
    }

private:
    PS ps_;
    vector<double> data_;
};

#endif // PSFILELOADER_H

#ifndef INTEGRATION_DEFINE_H
#define INTEGRATION_DEFINE_H

#include <Eigen/Geometry>
#include <vector>

#include "src/common/types.h"


static IntegrationStateData stateToData(const IntegrationState &state) {
    IntegrationStateData data;
    data.time = state.time;

    for (size_t i = 0; i < 3; i++)
    {
        data.xdata[i] = state.r(i);
        data.xdata[i + 3] = state.v(i);
        data.xdata[i + 6] = state.a(i);
    }
    data.xdata[9] = state.tu;
    data.xdata[10] = state.fu;
    data.xdata[11] = state.tdk;
    data.pdata = state.P;
    return data;
}

static IntegrationState stateFromData(const IntegrationStateData &data) {
    IntegrationState state;
    state.time = data.time;
    for (size_t i = 0; i < 3; i++)
    {
        state.r(i) = data.xdata[i];
        state.v(i) = data.xdata[i + 3];
        state.a(i) = data.xdata[i + 6];
    }
    state.tu = data.xdata[9];
    state.fu = data.xdata[10];
    state.tdk = data.xdata[11];
    state.P = data.pdata;
    return state;
}




#endif // INTEGRATION_DEFINE_H

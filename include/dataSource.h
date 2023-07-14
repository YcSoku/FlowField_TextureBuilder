//
// Created by Yucheng Soku on 2023/1/28.
//

#ifndef FLOWFIELD_BUILDER_DATASOURCE_H
#define FLOWFIELD_BUILDER_DATASOURCE_H

#endif //FLOWFIELD_BUILDER_DATASOURCE_H

#include <iostream>
#include <ogr_core.h>
#include "gdal_utils.h"
#include "ogrsf_frmts.h"
#include "gdal_priv.h"
#include "ogr_geometry.h"

class GeoDataSource
{
private:
    const char* Path;
    GDALDataset* Content;
public:
    GeoDataSource(const char* path);
    ~GeoDataSource();

    OGRGeometry*                        GetGeometry(int iLayer=0, int iFeature=0);
    OGRSpatialReference*                GetSpatialReference(int iLayer=0);
    char*                               GetSpatialReferenceWKT(int iLayer=0);
    OGREnvelope                         GetEnvelope(int iLayer=0, int iFeature=0);
    void                                BuildMask(const char* path, int width, int height, int iLayer=0, int iField=0);
};
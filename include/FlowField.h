//
// Created by Yucheng Soku on 2023/7/14.
//

#ifndef FLOWFIELD_BUILDER_FLOWFIELD_H
#define FLOWFIELD_BUILDER_FLOWFIELD_H

#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>

#include "nlohmann/json.hpp"
using Json = nlohmann::json;
#include <ogr_core.h>
#include "gdal_utils.h"
#include "ogrsf_frmts.h"
#include "gdal_priv.h"
#include "ogr_geometry.h"
#include "gridxtyp.h"
#include "scatdata.h"
#include "surfgrid.h"
#include "xpand.h"
#include "contour.h"

enum Platform_Types {
    Mapbox = 0,
    CESIUM,
    EPSG,

    PLATFORM_COUNT
};

struct FlowFiledFile
{
    std::string in; // Input File Name
    std::string sn; // Serial Number
};

union Texel {
    float f32;
    uint8_t rgba8[4];
};

class FlowField
{
private:
    std::string m_descriptionPath;
    std::string m_inUrl;
    std::string m_outUrl;
    std::string m_xName;
    std::string m_yName;
    std::string m_uName;
    std::string m_vName;
    std::string m_vMask;
    int m_attributeIndex[4];
    int m_headingNum;
    int m_columnNum;
    double m_outlier;
    double m_resolution;
    double m_boundary[4];
    double m_flowVelocityRange[4];

    std::string m_sourceSpace;
    std::vector<std::string> m_targetSpaces;

    GDALDataset* m_ds;
    OGRSpatialReference* m_sRef;
    OGRGeometry* m_vectorMask;
    char** m_maskOptions;
    std::vector<char> m_bmpMask;

    std::vector<FlowFiledFile> m_files;
    Json& m_result;

public:
    std::string urlPath;
    std::string resultPath;

private:
    void updateFlowBoundary(double u, double v);
    bool toTextureGRID_uv(ScatData& sc_u, ScatData& sc_v, SurfaceGrid& grd_u, SurfaceGrid& grd_v, const char* textureFile, const char* maskFile, const char* poolFile);
    void MercatorCoordsFromLonLat(OGRCoordinateTransformation* poCT, double& lon, double& lat);
    void EPSGFromLonLat(OGRCoordinateTransformation* poCT, double& lon, double& lat);
    void CesiumWgs84ToCartesian(OGRCoordinateTransformation* poCT, double& lon, double& lat, double& altitude);
    OGRGeometry* createVectorMask();
    void buildBmpMask(int width, int height);
    void readBmpMask();
    void closeDataSet();

public:
    explicit FlowField(const char* descriptionPath, Json& resultJson);
    ~FlowField();
    static FlowField* Create(const char* descriptionPath, Json& resultJson);
    void preprocess();
    void buildProjectionTexture(int width, int height, size_t index);
    void buildTextures();
    size_t getTargetSpaceNum() {return m_targetSpaces.size();}
};

#endif //FLOWFIELD_BUILDER_FLOWFIELD_H

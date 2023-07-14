#include "FlowField.h"

int main(int argc, char *argv[]) {

//    auto* test = new Cartesian3();
//    test = Cartesian3::FromDegrees(120.980697, 31.684162, 0.0);
//    auto m = float(test->x);

    // DebugPath: /Users/soku/Desktop/FlowField_Builder/resource/description.json
    Json resultJson;

    auto ff = FlowField::Create(argv[1], resultJson);
    ff->preprocess();
    for(auto i = 0; i < ff->getTargetSpaceNum(); ++i)
    {
        ff->buildProjectionTexture(1024, 2048, i);
    }
    ff->buildTextures();

    std::ofstream oStream(ff->resultPath.c_str());
    oStream << std::setw(4) << resultJson << std::endl;

    return 0;
}

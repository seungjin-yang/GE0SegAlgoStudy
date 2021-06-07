#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MuonTriggering/GE0Segment/plugins/GE0DatasetWriter.h"
#include "MuonTriggering/GE0Segment/plugins/GE0GeometryDumper.h"

DEFINE_FWK_MODULE(GE0DatasetWriter);
DEFINE_FWK_MODULE(GE0GeometryDumper);

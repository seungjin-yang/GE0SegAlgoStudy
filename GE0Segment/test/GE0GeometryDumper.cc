#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"

#include <TTree.h>

class GE0GeometryDumper : public edm::EDAnalyzer {
 public:
  explicit GE0GeometryDumper(const edm::ParameterSet&);
  ~GE0GeometryDumper();
  static void fillDescriptions(edm::ConfigurationDescriptions &);

 private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> kGEMToken_;
  const TString kTreeName_;

  TFileService* file_service_;
  TTree* tree_;
  long b_region_;
  long b_station_;
  long b_chamber_;
  long b_layer_;
  long b_ieta_;
  long b_strip_;
  float b_global_x_;
  float b_global_y_;
  float b_global_z_;
  float b_global_r_;
  float b_global_eta_;
  float b_global_phi_;
  float b_superchamber_x_;
  float b_superchamber_y_;
  float b_superchamber_z_;
  float b_eta_partition_x_;
  float b_eta_partition_y_;
  float b_eta_partition_z_;
  float b_eta_partition_phi_;

  const std::string kLogCategory_ = "GE0GeometryDumper";
};

GE0GeometryDumper::GE0GeometryDumper(const edm::ParameterSet& pset)
     : kGEMToken_(esConsumes<GEMGeometry, MuonGeometryRecord>()),
       kTreeName_(static_cast<TString>(pset.getUntrackedParameter<std::string>("treename"))) {

  tree_ = file_service_->make<TTree>(kTreeName_, kTreeName_);

  #define BRANCH_(name, suffix) tree_->Branch(#name, &b_##name##_, #name "/" #suffix);
  #define BRANCH_F(name) BRANCH_(name, F);
  #define BRANCH_L(name) BRANCH_(name, L);

  BRANCH_L(region)
  BRANCH_L(station)
  BRANCH_L(chamber)
  BRANCH_L(layer)
  BRANCH_L(ieta)
  BRANCH_L(strip)
  BRANCH_F(global_x)
  BRANCH_F(global_y)
  BRANCH_F(global_z)
  BRANCH_F(global_r)
  BRANCH_F(global_eta)
  BRANCH_F(global_phi)
  BRANCH_F(superchamber_x)
  BRANCH_F(superchamber_y)
  BRANCH_F(superchamber_z)
  BRANCH_F(eta_partition_x)
  BRANCH_F(eta_partition_y)
  BRANCH_F(eta_partition_z)
  BRANCH_F(eta_partition_phi)
}

void GE0GeometryDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.addUntracked<std::string>("treename", "Geometry");
  descriptions.add("GE0GeometryDumper", desc);
}

GE0GeometryDumper::~GE0GeometryDumper() {
  std::cout << "dtor begin" << std::endl;
  std::cout << "dtor end" << std::endl;
}


void GE0GeometryDumper::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  edm::ESHandle<GEMGeometry> gem = setup.getHandle(kGEMToken_);
  if (not gem.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMGeometry";
  }

  for (const GEMStation* station : gem->stations()) {
    if (station->station() != 0) {
      continue;
    }

    b_region_ = station->region();
    b_station_ = station->station();

    for (const GEMSuperChamber* superchamber : station->superChambers()) {
      for (const GEMChamber* chamber : superchamber->chambers()) {
        for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
          const GEMDetId&& gem_id = eta_partition->id();

          b_chamber_ = gem_id.chamber();
          b_layer_ = gem_id.layer();
          b_ieta_ = gem_id.ieta();

          for (const int strip : {0, eta_partition->nstrips()}) {
            const LocalPoint&& ieta_pos = eta_partition->centreOfStrip(strip);
            const GlobalPoint&& global_pos = eta_partition->toGlobal(ieta_pos);
            const LocalPoint&& superchamber_pos = superchamber->toLocal(global_pos);

            const StripTopology& topology = eta_partition->specificTopology();

            b_strip_ = strip;
            // global point
            b_global_x_ = global_pos.x();
            b_global_y_ = global_pos.y();
            b_global_z_ = global_pos.z();
            b_global_r_ = global_pos.perp();
            b_global_eta_ = global_pos.eta();
            b_global_phi_ = global_pos.barePhi();
            // on the superchamber
            b_superchamber_x_ = superchamber_pos.x();
            b_superchamber_y_ = superchamber_pos.y();
            b_superchamber_z_ = superchamber_pos.z();
            // on the eta partition
            b_eta_partition_x_ = ieta_pos.x();
            b_eta_partition_y_ = ieta_pos.y();
            b_eta_partition_z_ = ieta_pos.z();
            b_eta_partition_phi_ = topology.stripAngle(strip);

            tree_->Fill();

            std::cout << gem_id << " strip " << strip << " : " << superchamber_pos << std::endl;

          } // strip
        } // eta_partition
      } // chamber
    } // superchamber

  } // station
}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GE0GeometryDumper);

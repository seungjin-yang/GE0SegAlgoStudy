#ifndef MuonTriggering_MuonGEMDigis_plugins_GE0SegmentAnalyser_h
#define MuonTriggering_MuonGEMDigis_plugins_GE0SegmentAnalyser_h

#include "MuonTriggering/MuonGEMDigis/plugins/GE0SimSegment.cc"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

class TTree;
class TH1F;
class TH2F;

using namespace std;

class GE0SegmentAnalyser : public edm::EDAnalyzer {
 public:
  explicit GE0SegmentAnalyser(const edm::ParameterSet&);
  ~GE0SegmentAnalyser();
  static void fillDescriptions(edm::ConfigurationDescriptions &);

 private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  void beginFileService(); // with histograms ..
  void bookHistogram();
  void resetBranch();

  long get3DImageIndexWindow(long layer, long roll, long strip);
  
  bool isSimTrackGood(edm::SimTrackContainer::const_iterator);

  bool isSimHitGood(edm::PSimHitContainer::const_iterator);

  bool matchWithRecHit(const int, const GEMRecHit&);

  tuple<bool, uint32_t, bool> isSimSegmentGood(
      const vector<edm::PSimHitContainer::const_iterator> &);

  GE0SimSegmentCollection reconstructSimSegment(
      const edm::Handle<edm::SimTrackContainer>&,
      const edm::Handle<edm::PSimHitContainer>&,
      const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >&,
      const edm::ESHandle<GEMGeometry>&);

  bool associateRecSegToSimSeg(
      const GEMSegmentCollection::const_iterator&,
      const GE0SimSegment*,
      const edm::ESHandle<GEMGeometry>&);


  // ----------member data ---------------------------
  //////////////////////////////////////////////////////////////////////////////
  // NOTE Parameter
  //////////////////////////////////////////////////////////////////////////////
  edm::EDGetTokenT<edm::SimTrackContainer>             sim_track_token_;
  edm::EDGetTokenT<edm::PSimHitContainer>              sim_hit_token_;
  edm::EDGetTokenT<GEMDigiCollection>                  gem_digi_token_;
  edm::EDGetTokenT<edm::DetSetVector<GEMDigiSimLink> > gem_link_token_;
  edm::EDGetTokenT<GEMRecHitCollection>                gem_rechit_token_;
  edm::EDGetTokenT<GEMSegmentCollection>               gem_segment_token_;
  // edm::EDGetTokenT<reco::GenParticleCollection>        gen_particle_token_;

  double min_pt_;
  double min_quality_;
  unsigned int min_num_layers_;
  unsigned int min_digis_;  
  unsigned int max_muons_;  

  // NOTE FileService
  edm::Service<TFileService> file_service_;

  TTree* tree_;
  TTree* tree_win_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Branches
  // DL-aware branch
  // features -> float
  // index or label -> long
  //////////////////////////////////////////////////////////////////////////////
  long b_digi_size_;
  std::vector<long> b_digi_layer_; // for 
  std::vector<long> b_digi_ieta_;
  std::vector<long> b_digi_strip_;
  std::vector<long> b_digi_label_;
  std::vector<long> b_digi_particle_type_;
  std::vector<long> b_digi_track_id_;

  long b_muon_size_;
  std::vector<float> b_muon_pt_;
  std::vector<float> b_muon_eta_;
  std::vector<float> b_muon_phi_;

  long b_rechit_size_;
  std::vector<long> b_rechit_layer_;
  std::vector<long> b_rechit_ieta_;
  std::vector<long> b_rechit_strip_; // from LocalPoint
  std::vector<long> b_rechit_first_strip_;
  std::vector<long> b_rechit_cls_;
  std::vector<long> b_rechit_bx_;
  std::vector<float> b_rechit_x_;
  std::vector<float> b_rechit_y_;
  std::vector<float> b_rechit_z_;
  std::vector<long> b_rechit_label_;

  // GEMSegments built by RU (Road Usage) algorithm
  // https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_1/RecoLocalMuon/GEMSegment/python/GEMSegmentsRU_cfi.py
  long b_ru_size_;
  std::vector<long> b_ru_muon_idx_;
  std::vector<float> b_ru_norm_chi2_; // normalized chi2
  std::vector<long> b_ru_rechit_size_;
  std::vector<std::vector<long> > b_ru_rechit_layer_;
  std::vector<std::vector<long> > b_ru_rechit_ieta_;
  std::vector<std::vector<long> > b_ru_rechit_strip_;
  std::vector<std::vector<long> > b_ru_rechit_first_strip_;
  std::vector<std::vector<long> > b_ru_rechit_cls_;
  std::vector<std::vector<long> > b_ru_rechit_bx_;

  // additional
  long b_region_;
  long b_chamber_;

  // NOTE window
  bool  b_win_digi_image_[180]; // [layer][ieta 3][strip 10] --> 6 * 3 * 10
  bool  b_win_digi_image_label_[180]; // actually SimTrack with abs(PID) == 6
  long   b_win_strip_;
  long   b_win_ieta_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE histograms for summary & monitoring
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Constants
  //////////////////////////////////////////////////////////////////////////////
  const std::string kLogCategory_ = "GE0SegmentAnalyser";
  const int num_digi_ = 6 * 8 * 384;
};

#endif // MuonTriggering_MuonGEMDigis_plugins_GE0SegmentAnalyser_h

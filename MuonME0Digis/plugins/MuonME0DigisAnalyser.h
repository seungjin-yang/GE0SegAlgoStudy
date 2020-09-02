#ifndef MuonTriggering_MuonME0Digis_plugins_MuonME0DigisAnalyser_h
#define MuonTriggering_MuonME0Digis_plugins_MuonME0DigisAnalyser_h

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
#include "SimDataFormats/GEMDigiSimLink/interface/ME0DigiSimLink.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"
#include "DataFormats/GEMDigi/interface/ME0DigiCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0RecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticle.h"
// #include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"


class ME0MuonData;

class TTree;
class TH1F;
class TH2F;

using namespace std;

class MuonME0DigisAnalyser : public edm::EDAnalyzer {
 public:
  explicit MuonME0DigisAnalyser(const edm::ParameterSet&);
  ~MuonME0DigisAnalyser();

 private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  void setBranch(); // with histograms ..
  void bookHistogram();
  void resetBranch();

  int get3DImageIndex(int layer, int roll, int strip);
  int get3DImageIndexWindow(int layer, int roll, int strip);
  
  bool isSimTrackGood(edm::SimTrackContainer::const_iterator);

  bool isSimHitGood(edm::PSimHitContainer::const_iterator);

  tuple<bool, uint32_t, bool> isSimSegmentGood(
      const vector<edm::PSimHitContainer::const_iterator> &);

  map<pair<unsigned int, uint32_t>, ME0MuonData> buildDB(
      const edm::Handle<edm::SimTrackContainer>&,
      const edm::Handle<edm::PSimHitContainer>&,
      const edm::Handle<ME0DigiCollection> &,
      const edm::Handle<edm::DetSetVector<ME0DigiSimLink> >&,
      const edm::Handle<ME0SegmentCollection>&,
      const ME0Geometry*);


  // ----------member data ---------------------------
  //////////////////////////////////////////////////////////////////////////////
  // NOTE ParameterSet
  //////////////////////////////////////////////////////////////////////////////
  edm::EDGetTokenT<edm::SimTrackContainer>             sim_track_token_;
  edm::EDGetTokenT<edm::PSimHitContainer>              sim_hit_token_;
  edm::EDGetTokenT<ME0DigiCollection>                  me0_digi_token_;
  edm::EDGetTokenT<edm::DetSetVector<ME0DigiSimLink> > me0_link_token_;
  edm::EDGetTokenT<ME0RecHitCollection>                me0_rechit_token_;
  edm::EDGetTokenT<ME0SegmentCollection>               me0_segment_token_;
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
  //////////////////////////////////////////////////////////////////////////////
  int b_digi_size_;
  std::vector<int> b_digi_layer_;
  std::vector<int> b_digi_ieta_;
  std::vector<int> b_digi_strip_;
  std::vector<int> b_digi_label_;
  std::vector<int> b_digi_particle_type_;
  std::vector<int> b_digi_track_id_;
  bool b_digi_image_[18432]; // [layer][ieta][strip] --> 6 * 8 * 384
  bool b_digi_image_label_[18432]; // [layer][ieta][strip] --> 6 * 8 * 384

  int b_muon_size_;
  std::vector<float> b_muon_pt_;
  std::vector<float> b_muon_eta_;
  std::vector<float> b_muon_phi_;
  std::vector<int> b_muon_charge_;

  int b_rechit_size_;
  std::vector<int> b_rechit_layer_;
  std::vector<int> b_rechit_ieta_;
  std::vector<int> b_rechit_strip_;
  std::vector<float> b_rechit_x_;
  std::vector<float> b_rechit_y_;
  std::vector<float> b_rechit_z_;
  std::vector<float> b_rechit_tof_;
  std::vector<int> b_rechit_label_;

  // ME0Segments built by RU (Road Usage) algorithm
  // https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_1/RecoLocalMuon/GEMSegment/python/ME0SegmentsRU_cfi.py
  int b_ru_size_;
  std::vector<int> b_ru_muon_idx_;
  std::vector<float> b_ru_reduced_chi2_;
  std::vector<int> b_ru_rechit_size_;
  std::vector<std::vector<int> > b_ru_rechit_layer_;
  std::vector<std::vector<int> > b_ru_rechit_ieta_;
  std::vector<std::vector<int> > b_ru_rechit_strip_;

  // additional
  int b_region_;
  int b_chamber_;

  // NOTE window
  bool  b_win_digi_image_[180]; // [layer][ieta 3][strip 10] --> 6 * 3 * 10
  bool  b_win_digi_image_label_[180]; // actually SimTrack with abs(PID) == 6
  int   b_win_strip_;
  int   b_win_ieta_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE histograms for summary & monitoring
  //////////////////////////////////////////////////////////////////////////////
  TH1F* h_stats_;
  TH1F* h_sim_seg_nhits_;
  TH1F* h_sim_track_pt_;
  TH1F* h_sim_track_eta_;
  TH1F* h_sim_track_phi_;


  //////////////////////////////////////////////////////////////////////////////
  // NOTE Constants
  //////////////////////////////////////////////////////////////////////////////

  const std::string kLogCategory_ = "MuonME0DigisAnalyser";
};

#endif // MuonTriggering_MuonME0Digis_plugins_MuonME0DigisAnalyser_h

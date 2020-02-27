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

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartition.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"


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
  void resetBranch();


  int getIndex(int layer, int roll, int strip);
  int getIndexWindow(int layer, int roll, int strip);
  
  int getUniqueId(int region, int chamber,
                  int layer, int roll, int strip);

  int getUniqueId(const ME0DetId & det_id, int strip);

  bool isSimTrackGood(edm::SimTrackContainer::const_iterator);
  bool isSimHitGood(edm::PSimHitContainer::const_iterator);
  bool isSimSegmentGood(const vector<edm::PSimHitContainer::const_iterator> &);

  // ----------member data ---------------------------

  // NOTE ParameterSet
  edm::EDGetTokenT<edm::SimTrackContainer>             sim_track_token_;
  edm::EDGetTokenT<edm::PSimHitContainer>              sim_hit_token_;
  edm::EDGetTokenT<ME0DigiCollection>                  me0_digi_token_;
  edm::EDGetTokenT<edm::DetSetVector<ME0DigiSimLink> > me0_link_token_;
  edm::EDGetTokenT<ME0SegmentCollection>               me0_segment_token_;

  double pt_min_;

  // NOTE FileService
  edm::Service<TFileService> file_service_;

  TTree* tree_ch_;
  TTree* tree_win_;
  TTree* tree_multi_;

  // NOTE tree_ch_
  int b_num_digi_;
  bool b_digi_[18432]; // [layer][ieta][strip] --> 6 * 8 * 384
  std::vector<int> b_digi_layer_;
  std::vector<int> b_digi_ieta_;
  std::vector<int> b_digi_strip_;
  std::vector<int> b_digi_is_muon_;

  // debug
  std::vector<int> b_digi_particle_type_;
  std::vector<int> b_digi_track_id_;

  // Muon SimTrack & it's digi
  bool b_has_muon_;
  int b_num_muon_digi_;
  bool b_muon_digi_[18432];
  std::vector<int> b_muon_digi_layer_;
  std::vector<int> b_muon_digi_ieta_;
  std::vector<int> b_muon_digi_strip_;
  float b_muon_pt_;
  float b_muon_eta_;
  float b_muon_phi_;

  // ME0Segments built by RU (Road Usage) algorithm
  // https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_1/RecoLocalMuon/GEMSegment/python/ME0SegmentsRU_cfi.py
  bool b_has_ru_;
  int b_num_ru_;

  // has an associated rec segment
  bool b_has_ru_asso_;
  int b_ru_asso_nhits_;
  std::vector<int> b_ru_asso_rechit_layer_;
  std::vector<int> b_ru_asso_rechit_ieta_;
  std::vector<int> b_ru_asso_rechit_strip_;

  // fake means ME0Segments that are not associated with a muon.
  int b_num_ru_fake_;
  std::vector<int> b_ru_fake_nhits_;

  // additional
  int b_region_;
  int b_chamber_;

  // NOTE window
  bool  b_win_digi_[180]; // [layer][ieta 3][strip 10] --> 6 * 3 * 10
  bool  b_win_muon_digi_[180]; // actually SimTrack with abs(PID) == 6
  int   b_win_strip_;
  int   b_win_ieta_;

  // NOTE multiple muons
  int b_multi_num_muon_;

  float b_multi_muon_pt_[5];
  float b_multi_muon_eta_[5];
  float b_multi_muon_phi_[5];
  int b_multi_muon_num_digi_[5];
  int b_multi_muon_digi_idx_[5][20];

  int b_multi_num_ru_asso_;
  std::vector<int> b_multi_ru_asso_nhits_;
  std::vector<int> b_multi_ru_asso_muon_idx_;

  int b_multi_num_ru_fake_;
  std::vector<int> b_multi_ru_fake_nhits_;

  //////////////////////////////////////////////////////////////////////////////
  // histograms for summary & monitoring
  //////////////////////////////////////////////////////////////////////////////
  TH1F* h_sim_track_pt_;
  TH1F* h_sim_track_eta_;
  TH1F* h_sim_track_phi_;
  TH1F* h_matched_sim_track_pt_;
  TH1F* h_matched_sim_track_eta_;
  TH1F* h_matched_sim_track_phi_;

  // TODO
  // TH1F* h_segment_num_rechits_;
  // TH1F* h_matched_segment_num_rechits_;

  TH1F* h_num_simhit_;
  TH1F* h_num_muon_;
  TH1F* h_stats_;
  TH1F* h_num_rec_seg_;
  TH2F* h_num_sim_rec_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Constants
  //////////////////////////////////////////////////////////////////////////////
  const int kNumChambers_ = 18; // per region
  const int kNumLayers_ = 6; // per chamber
  const int kNumEtaPartitions_ = 8; // per layer
  const int kNumStrips_ = 384; // per eta partition
  const int kMuonPDGId_ = 13;

  const int kMaxNumMuons_ = 5;
  const unsigned int kMaxNumDigisPerMuon_ = 20;

  const std::string kLogCategory_ = "MuonME0DigisAnalyser";
};

#endif // MuonTriggering_MuonME0Digis_plugins_MuonME0DigisAnalyser_h

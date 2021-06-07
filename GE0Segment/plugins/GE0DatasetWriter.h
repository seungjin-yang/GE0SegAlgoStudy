#ifndef MuonTriggering_GE0Segment_GE0DatasetWriter_h
#define MuonTriggering_GE0Segment_GE0DatasetWriter_h

#include "MuonTriggering/GE0Segment/interface/GE0SimSegment.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMPadDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

class GE0DatasetWriter : public edm::EDAnalyzer {
 public:
  explicit GE0DatasetWriter(const edm::ParameterSet&);
  ~GE0DatasetWriter();
  static void fillDescriptions(edm::ConfigurationDescriptions &);

  template <typename T>
  inline edm::EDGetTokenT<T> getToken(const edm::ParameterSet&, const std::string&);
  
 private:
  using Digi2Index = std::map<std::tuple<int, int, int>, unsigned int>;

  // method
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  void beginFileService(); // with histograms ..
  void resetBranch();

  long get3DImageIndexWindow(long layer, long roll, long strip);
  
  bool isSimTrackGood(edm::SimTrackContainer::const_iterator);

  std::vector<GE0SimSegment::LinkData> getGE0LinksFromSimTrack(
      const GE0SimSegment::SimTrackData &,
      const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >&);

  std::tuple<bool, uint32_t, bool> checkSimSegment(
      const std::vector<GE0SimSegment::LinkData>&);

  std::vector<GE0SimSegment> buildGE0SimSegments(
      const edm::Handle<edm::SimTrackContainer> &,
      const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >&,
      const edm::Handle<edm::PSimHitContainer> &,
      const edm::ESHandle<GEMGeometry>&);

  bool matchWithRecHit(const int, const int, const GEMRecHit&); 
  bool matchWithRecHit(const GE0SimSegment::LinkData&, const GEMRecHit&);

  std::pair<float, unsigned int> computeEfficiency(
      const GEMSegmentCollection::const_iterator&,
      const GE0SimSegment*);

  float computeFakeHitRate(
      const GEMSegmentCollection::const_iterator&,
      const GE0SimSegment*);

  std::vector<const GE0SimSegment*> getSimSegmentsInSuperChamber(
      const std::vector<GE0SimSegment>&,
      const GEMDetId &);

  LocalPoint getSuperChamberPosition(
      const LocalPoint&,
      const GEMSuperChamber*,
      const GEMEtaPartition*);

  void analyzeSuperChamber(const GEMSuperChamber*);

  void analyzeMuon(const std::vector<const GE0SimSegment*>&,
                   const edm::ESHandle<GEMGeometry>&);

  bool analyzeDigi(
      const edm::Handle<GEMDigiCollection>&,
      const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >&,
      const std::vector<const GE0SimSegment*>&,
      const GEMSuperChamber*);

  void analyzeWindow(const Digi2Index&); 

  bool analyzePad(
      const edm::Handle<GEMPadDigiCollection>&,
      const std::vector<const GE0SimSegment*>&,
      const GEMSuperChamber*,
      const edm::ESHandle<GEMGeometry>&);

  bool analyzeRecHit(
      const edm::Handle<GEMRecHitCollection>&,
      const std::vector<const GE0SimSegment*>&,
      const GEMSuperChamber*,
      const edm::ESHandle<GEMGeometry>&);

  bool analyzeSegmentRU(
      const GEMSegmentCollection::range&,
      const std::vector<const GE0SimSegment*>&,
      const edm::ESHandle<GEMGeometry>&);

  template <typename T>
  inline unsigned int countUnique(const std::vector<T>&);

  // ----------member data ---------------------------
  // Parameter
  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord>     kGEMToken_;
  const edm::EDGetTokenT<edm::SimTrackContainer>             kSimTrackToken_;
  const edm::EDGetTokenT<edm::PSimHitContainer>              kSimHitToken_;
  const edm::EDGetTokenT<GEMDigiCollection>                  kDigiToken_;
  const edm::EDGetTokenT<GEMPadDigiCollection>               kPadToken_;
  const edm::EDGetTokenT<edm::DetSetVector<GEMDigiSimLink> > kLinkToken_;
  const edm::EDGetTokenT<GEMRecHitCollection>                kRecHitToken_;
  const edm::EDGetTokenT<GEMSegmentCollection>               kSegmentToken_;
  const double kMinMuonPT_;
  const unsigned int kMinNumLayers_;
  const unsigned int kMaxNumMuons_;  
  // FileService
  edm::Service<TFileService> file_service_;
  TTree* tree_;

  TH1F* h_num_ge0_links_;
  TH1F* h_num_simseg_;
  TH1F* h_simseg_num_layers_;

  TH1F* h_muon_bending_x_;
  TH1F* h_muon_bending_y_;
  TH1F* h_muon_bending_ieta_;
  TH1F* h_muon_bending_strip_;

  TH1F* h_hit_dx_;
  TH1F* h_hit_dy_;
  TH1F* h_hit_sx_;
  TH1F* h_hit_sy_;

  // CUDA-aware branch, features -> float, index or label -> long
  // GEMSuperChamber
  long b_region_;
  long b_station_;
  long b_chamber_;
  // Muon (GE0SimSegment)
  long b_muon_size_;
  std::vector<float> b_muon_pt_;
  std::vector<float> b_muon_eta_;
  std::vector<float> b_muon_phi_;
  std::vector<long> b_muon_charge_; 
  std::vector<long> b_muon_hit_size_;
  std::vector<std::vector<float> > b_muon_hit_x_; 
  std::vector<std::vector<float> > b_muon_hit_y_; 
  std::vector<std::vector<float> > b_muon_hit_z_; 
  std::vector<std::vector<long> > b_muon_hit_layer_; 
  std::vector<std::vector<long> > b_muon_hit_ieta_; 
  std::vector<std::vector<long> > b_muon_hit_strip_;
  std::vector<float> b_muon_bending_x_;
  std::vector<float> b_muon_bending_y_;
  std::vector<long> b_muon_bending_ieta_;
  std::vector<long> b_muon_bending_strip_;
  // GEMDigi
  long b_digi_size_;
  long b_digi_layer_count_;
  std::vector<long> b_digi_layer_;
  std::vector<long> b_digi_ieta_;
  std::vector<long> b_digi_strip_;
  std::vector<long> b_digi_label_;
  // GEMPadDigi
  long b_pad_size_;
  long b_pad_layer_count_;
  std::vector<long> b_pad_layer_;
  std::vector<long> b_pad_ieta_;
  std::vector<long> b_pad_pad_;
  std::vector<long> b_pad_label_;
  // GEMRecHit
  long b_rechit_size_;
  long b_rechit_layer_count_;
  std::vector<long> b_rechit_layer_;
  std::vector<long> b_rechit_ieta_;
  std::vector<long> b_rechit_strip_; // from LocalPoint
  std::vector<long> b_rechit_first_strip_;
  std::vector<long> b_rechit_cls_;
  std::vector<long> b_rechit_bx_;
  std::vector<float> b_rechit_x_;
  std::vector<float> b_rechit_y_;
  std::vector<float> b_rechit_z_;
  std::vector<float> b_rechit_err_xx_;
  std::vector<float> b_rechit_err_xy_;
  std::vector<float> b_rechit_err_yy_;
  std::vector<long> b_rechit_label_;
  // GEMSegments built by RU (Road Usage) algorithm
  long b_ru_size_;
  std::vector<long> b_ru_muon_idx_;
  std::vector<long> b_ru_rechit_size_;
  std::vector<long> b_ru_num_matched_;
  std::vector<float> b_ru_eff_;
  std::vector<float> b_ru_fake_hit_rate_;
  std::vector<float> b_ru_norm_chi2_; // normalized chi2
  std::vector<std::vector<long> > b_ru_rechit_layer_;
  std::vector<std::vector<long> > b_ru_rechit_ieta_;
  std::vector<std::vector<long> > b_ru_rechit_strip_;
  std::vector<std::vector<long> > b_ru_rechit_first_strip_;
  std::vector<std::vector<long> > b_ru_rechit_cls_;
  std::vector<std::vector<long> > b_ru_rechit_bx_;
  // Window
  bool b_has_window_;
  bool b_window_image_[180]; // [layer][ieta 3][strip 10] --> 6 * 3 * 10
  bool b_window_label_[180]; // actually SimTrack with abs(PID) == 6
  long b_window_ieta_;
  long b_window_strip_;

  const std::string kLogCategory_ = "GE0DatasetWriter";
  const int kMuonPID_ = 13;
};

template <typename T>
inline edm::EDGetTokenT<T> GE0DatasetWriter::getToken(const edm::ParameterSet& pset, const std::string& name) {
  return consumes<T>(pset.getParameter<edm::InputTag>(name));
}

template <typename T>
inline unsigned int GE0DatasetWriter::countUnique(const std::vector<T>& v) {
  const std::unordered_set us(v.begin(), v.end());
  return us.size();
}

#endif // MuonTriggering_GE0Segment_GE0DatasetWriter_h

#include "MuonTriggering/GE0Segment/plugins/GE0DatasetWriter.h"
#include "MuonTriggering/GE0Segment/interface/MinMaxScaler.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"

#include <numeric> // iota
#include <algorithm> // minmax_element

using namespace std;


GE0DatasetWriter::GE0DatasetWriter(const edm::ParameterSet& pset)
    : kGEMToken_(esConsumes<GEMGeometry, MuonGeometryRecord>()),
      kSimTrackToken_(getToken<edm::SimTrackContainer>(pset, "simTrackTag")),
      kSimHitToken_(getToken<edm::PSimHitContainer>(pset, "gemSimHitTag")),
      kDigiToken_(getToken<GEMDigiCollection>(pset, "gemDigiTag")),
      kPadToken_(getToken<GEMPadDigiCollection>(pset, "gemPadTag")),
      kLinkToken_(getToken<edm::DetSetVector<GEMDigiSimLink> >(pset, "gemLinkTag")),
      kRecHitToken_(getToken<GEMRecHitCollection>(pset, "gemRecHitTag")),
      kSegmentToken_(getToken<GEMSegmentCollection>(pset, "gemSegmentTag")),
      kMinMuonPT_(pset.getParameter<double>("minMuonPt")),
      kMinNumLayers_(pset.getParameter<unsigned int>("minNumLayers")),
      kMaxNumMuons_(pset.getParameter<unsigned int>("maxNumMuons")) {

  cout << "ctor begin" << endl;
  beginFileService();
  cout << "ctor end" << endl;
}

void GE0DatasetWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("simTrackTag", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("gemSimHitTag", edm::InputTag("g4SimHits", "MuonGEMHits"));
  desc.add<edm::InputTag>("gemDigiTag", edm::InputTag("simMuonGEMDigis"));
  desc.add<edm::InputTag>("gemLinkTag", edm::InputTag("simMuonGEMDigis", "GEM"));
  desc.add<edm::InputTag>("gemPadTag", edm::InputTag("simMuonGEMPadDigis")); // FIXME
  desc.add<edm::InputTag>("gemRecHitTag", edm::InputTag("gemRecHits"));
  desc.add<edm::InputTag>("gemSegmentTag", edm::InputTag("gemSegments"));
  desc.add<double>("minMuonPt", 5.0); // GeV
  desc.add<unsigned int>("minNumLayers", 4);
  desc.add<unsigned int>("maxNumMuons", 3);
  descriptions.add("GE0", desc);
}

GE0DatasetWriter::~GE0DatasetWriter() {
  cout << "dtor begin" << endl;
  cout << "dtor end" << endl;
}

void GE0DatasetWriter::beginFileService() {
  cout << "beginFileService begin" << endl;

  // NOTE TTree
  // FIXME
  gInterpreter->GenerateDictionary("vector<vector<long> >", "vector");
  tree_ = file_service_->make<TTree>("chamber", "chamber");

  // for the single numerical variables and the nested vector 
  #define BRANCH(name) tree_->Branch(#name, &b_##name##_);
  #define BRANCH_(name, suffix) tree_->Branch(#name, &b_##name##_, #name "/" #suffix);
  #define BRANCH_F(name) BRANCH_(name, F);
  #define BRANCH_L(name) BRANCH_(name, L);
  // the array of numerical variables
  #define BRANCH_A_(name, size, suffix) tree_->Branch(#name, &b_##name##_, #name"["#size"]/"#suffix);
  #define BRANCH_AO(name, size)  BRANCH_A_(name, size, O);
  //
  #define BRANCH_V_(name, type) tree_->Branch(#name, "vector<"#type">", &b_##name##_);
  #define BRANCH_VF(name) BRANCH_V_(name, Float_t);
  #define BRANCH_VL(name) BRANCH_V_(name, Long_t);

  // GEMSuperChamber
  BRANCH_L(region)
  BRANCH_L(station)
  BRANCH_L(chamber)
  // Muon SimTrack
  BRANCH_L(muon_size)
  BRANCH_VF(muon_pt)
  BRANCH_VF(muon_eta)
  BRANCH_VF(muon_phi)
  BRANCH_VL(muon_charge)
  BRANCH_VL(muon_hit_size)
  BRANCH(muon_hit_x)
  BRANCH(muon_hit_y)
  BRANCH(muon_hit_z)
  BRANCH(muon_hit_layer)
  BRANCH(muon_hit_ieta)
  BRANCH(muon_hit_strip)
  BRANCH(muon_hit_entry_strip)
  BRANCH(muon_hit_exit_strip)
  BRANCH(muon_hit_first_strip)
  BRANCH(muon_hit_cls)
  BRANCH(muon_hit_bx)
  BRANCH_VF(muon_delta_x);
  BRANCH_VF(muon_delta_y);
  BRANCH_VL(muon_delta_ieta);
  BRANCH_VL(muon_delta_strip);
  // GEMDigi
  BRANCH_L(digi_size)
  BRANCH_L(digi_layer_count)
  BRANCH_VL(digi_layer)
  BRANCH_VL(digi_ieta)
  BRANCH_VL(digi_strip)
  BRANCH_VL(digi_label)
  // GEMPadDigi
  BRANCH_L(pad_size)
  BRANCH_L(pad_layer_count)
  BRANCH_VL(pad_layer)
  BRANCH_VL(pad_ieta)
  BRANCH_VL(pad_pad)
  BRANCH_VL(pad_label)
  // GEMRecHit
  BRANCH_L(rechit_size)
  BRANCH_L(rechit_layer_count)
  BRANCH_VL(rechit_layer)
  BRANCH_VL(rechit_ieta)
  BRANCH_VL(rechit_strip)
  BRANCH_VL(rechit_first_strip)
  BRANCH_VL(rechit_cls)
  BRANCH_VL(rechit_bx)
  BRANCH_VF(rechit_x)
  BRANCH_VF(rechit_y)
  BRANCH_VF(rechit_z)
  BRANCH_VF(rechit_x_minmax)
  BRANCH_VF(rechit_y_minmax)
  BRANCH_VF(rechit_z_minmax)
  BRANCH_VF(rechit_err_xx)
  BRANCH_VF(rechit_err_xy)
  BRANCH_VF(rechit_err_yy)
  BRANCH_VL(rechit_label)
  BRANCH_VL(rechit_label_first_ieta)
  BRANCH_VL(rechit_label_first_strip)
  BRANCH_VL(rechit_label_last_ieta)
  BRANCH_VL(rechit_label_last_strip)
  // GEMSegment by Road Usage
  BRANCH_L(ru_size)
  BRANCH_VL(ru_rechit_size)
  BRANCH_VL(ru_muon_idx)
  BRANCH_VL(ru_num_matched)
  BRANCH_VF(ru_eff)
  BRANCH_VF(ru_fake_hit_rate)
  BRANCH_VF(ru_norm_chi2)
  BRANCH(ru_rechit_layer)
  BRANCH(ru_rechit_ieta)
  BRANCH(ru_rechit_strip)
  BRANCH(ru_rechit_first_strip)
  BRANCH(ru_rechit_cls)
  BRANCH(ru_rechit_bx)
  // Window
  BRANCH(has_window)
  BRANCH_AO(window_image, 180)
  BRANCH_AO(window_label, 180)
  BRANCH_L(window_ieta)
  BRANCH_L(window_strip)

  // Histograms
  h_num_ge0_links_ = file_service_->make<TH1F>("h_num_ge0_links", "", 11, -0.5, 10.5);
  h_num_simseg_ = file_service_->make<TH1F>("h_num_simseg", "", 11, -0.5, 10.5);
  h_simseg_num_layers_ = file_service_->make<TH1F>("h_simseg_num_layers", "", 8, -0.5, 7.5);

  h_muon_delta_x_ = file_service_->make<TH1F>("h_muon_delta_x", "", 20, -2, 2);
  h_muon_delta_y_ = file_service_->make<TH1F>("h_muon_delta_y", "", 20, 0, 10);
  h_muon_delta_ieta_ = file_service_->make<TH1F>("h_muon_delta_ieta", "", 5, -2.5, 2.5);
  h_muon_delta_strip_ = file_service_->make<TH1F>("h_muon_delta_strip", "", 21, -10.5, 10.5);

  h_hit_dx_ = file_service_->make<TH1F>("h_hit_dx", "", 20, -2, 2);
  h_hit_dy_ = file_service_->make<TH1F>("h_hit_dy", "", 20, -10, 10);
  h_hit_sx_ = file_service_->make<TH1F>("h_hit_sx", "", 30, -3, 3);
  h_hit_sy_ = file_service_->make<TH1F>("h_hit_sy", "", 30, -3, 3);

  cout << "beginFileService end" << endl;
}

void GE0DatasetWriter::resetBranch() {
  // GEMSuperChamber
  b_region_ = -100L;
  b_station_ = -100L;
  b_chamber_ = -100L;
  // SimTrack, muon 
  b_muon_size_ = -1L;
  b_muon_pt_.clear();
  b_muon_eta_.clear();
  b_muon_phi_.clear();
  b_muon_charge_.clear();
  b_muon_hit_size_.clear();
  b_muon_hit_x_.clear();
  b_muon_hit_y_.clear();
  b_muon_hit_z_.clear();
  b_muon_hit_layer_.clear();
  b_muon_hit_ieta_.clear();
  b_muon_hit_strip_.clear();
  b_muon_delta_x_.clear();
  b_muon_delta_y_.clear();
  b_muon_delta_ieta_.clear();
  b_muon_delta_strip_.clear();
  // GEMDigi
  b_digi_size_ = -1L;
  b_digi_layer_count_ = -1L;
  b_digi_layer_.clear();
  b_digi_ieta_.clear();
  b_digi_strip_.clear();
  b_digi_label_.clear();
  // GEMPadDigi
  b_pad_size_ = -1L;
  b_pad_layer_count_ = -1L;
  b_pad_layer_.clear();
  b_pad_ieta_.clear();
  b_pad_label_.clear();
  // GEMRecHit
  b_rechit_size_ = -1L;
  b_rechit_layer_count_ = -1L;
  b_rechit_layer_.clear();
  b_rechit_ieta_.clear();
  b_rechit_strip_.clear();
  b_rechit_first_strip_.clear();
  b_rechit_cls_.clear();
  b_rechit_bx_.clear();
  b_rechit_x_.clear();
  b_rechit_y_.clear();
  b_rechit_z_.clear();
  b_rechit_x_minmax_.clear();
  b_rechit_y_minmax_.clear();
  b_rechit_z_minmax_.clear();
  b_rechit_err_xx_.clear();
  b_rechit_err_xy_.clear();
  b_rechit_err_yy_.clear();
  b_rechit_label_.clear();
  b_rechit_label_first_ieta_.clear();
  b_rechit_label_first_strip_.clear();
  b_rechit_label_last_ieta_.clear();
  b_rechit_label_last_strip_.clear();
  // GEMSegment by Road Usage
  b_ru_size_ = -1L;
  b_ru_muon_idx_.clear();
  b_ru_rechit_size_.clear();
  b_ru_num_matched_.clear();
  b_ru_eff_.clear();
  b_ru_fake_hit_rate_.clear();
  b_ru_norm_chi2_.clear();
  b_ru_rechit_layer_.clear();
  b_ru_rechit_ieta_.clear();
  b_ru_rechit_strip_.clear();
  b_ru_rechit_first_strip_.clear();
  b_ru_rechit_cls_.clear();
  b_ru_rechit_bx_.clear();
  // Window 
  b_has_window_ = false;
  std::fill_n(b_window_image_, 180, false);
  std::fill_n(b_window_label_, 180, false);
  b_window_strip_ = -1L;
  b_window_ieta_ = -1L;
}

long GE0DatasetWriter::get3DImageIndexWindow(long layer, long ieta, long strip) {
  return 30L * layer + 10L * ieta + strip - 41L;
}

bool GE0DatasetWriter::isSimTrackGood(edm::SimTrackContainer::const_iterator sim_track) {
  const EncodedEventId & event_id = sim_track->eventId();
  if (event_id.event() != 0) return false;
  if (event_id.bunchCrossing() != 0) return false;

  if (abs(sim_track->type()) != kMuonPID_) return false;
  if (sim_track->momentum().pt() < kMinMuonPT_) return false;

  return true;
}

std::vector<ge0::LinkData> GE0DatasetWriter::getGE0LinksFromSimTrack(
    const ge0::SimTrackData& sim_track ,
    const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >& link_collection) {

  const EncodedEventId&& event_id = sim_track->eventId();
  const unsigned int track_id = sim_track->trackId();

  std::vector<ge0::LinkData> output;
  for (auto link_set = link_collection->begin(); link_set != link_collection->end(); link_set++) {
    if (link_set->empty()) continue;

    for (auto link = link_set->begin(); link != link_set->end(); link++) {
      if (track_id != link->getTrackId()) continue;
      if (event_id != link->getEventId()) continue;

      const GEMDetId gem_id{link->getDetUnitId()};
      // TODO if (gem->idToDet(gem_id) == nullptr) {
      if (gem_id.station() != 0) continue;

      output.push_back(link);
    }
  }

  return output;
}


std::tuple<bool, uint32_t, bool> GE0DatasetWriter::checkSimSegment(
    const std::vector<ge0::LinkData>& link_vector) {

  bool is_good = false;
  uint32_t primary_superchamber{0};
  bool need_to_prune = false;
 
  if (link_vector.size() < kMinNumLayers_) {
    return std::make_tuple(is_good, primary_superchamber, need_to_prune);
  }

  // <superchamber_rawid, layer>
  map<uint32_t, std::set<int> > layers_per_superchamber;

  for (const auto & link : link_vector) {
    const GEMDetId gem_id{link->getDetUnitId()};
    const uint32_t superchamber_id = gem_id.superChamberId().rawId();

    if (layers_per_superchamber.find(superchamber_id) == layers_per_superchamber.end()) {
      layers_per_superchamber.emplace(superchamber_id, std::set<int>{});
    }

    layers_per_superchamber.at(superchamber_id).insert(gem_id.layer());
  }

  set<int> layers;
  if (layers_per_superchamber.size() == 1) {
    primary_superchamber = layers_per_superchamber.begin()->first;
    layers = layers_per_superchamber.begin()->second;

  } else {
    auto tmp_max_elem = std::max_element(
        layers_per_superchamber.begin(),
        layers_per_superchamber.end(),
        [](const pair<uint32_t, std::set<int> >& lhs, const pair<uint32_t, std::set<int> >&  rhs) {
            return lhs.second.size() < rhs.second.size();});

    primary_superchamber = tmp_max_elem->first;
    layers = tmp_max_elem->second;

    need_to_prune = true;
  }

  h_simseg_num_layers_->Fill(static_cast<int>(layers.size()));

  is_good = layers.size() >= kMinNumLayers_;

  return std::make_tuple(is_good, primary_superchamber, need_to_prune);
}

std::vector<ge0::GE0SimSegment> GE0DatasetWriter::buildGE0SimSegments(
    const edm::Handle<edm::SimTrackContainer> & sim_track_container,
    const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >& link_collection,
    const edm::Handle<edm::PSimHitContainer> & sim_hit_container,
    const edm::ESHandle<GEMGeometry>& gem) {

  // key == {track_id, superchamber raw id}
  std::vector<ge0::GE0SimSegment> sim_segment_collection;
  for (auto sim_track = sim_track_container->begin(); sim_track != sim_track_container->end(); sim_track++) {
    if (not isSimTrackGood(sim_track)) {
      continue;
    }

    // TODO rename
    auto&& ge0_links = getGE0LinksFromSimTrack(sim_track, link_collection);
    // FIXME
    h_num_ge0_links_->Fill(std::min(static_cast<int>(ge0_links.size()), 10));

    auto [is_sim_seg_good, primary_superchamber, need_to_prune] = checkSimSegment(ge0_links);
    if (not is_sim_seg_good) {
      continue;
    }

    if (need_to_prune) {
      // Remove hits if they aren't on the primary superchamber.
      ge0_links.erase(
          std::remove_if(
              ge0_links.begin(),
              ge0_links.end(),
              [&primary_superchamber](const ge0::LinkData& link) {
                  return GEMDetId{link->getDetUnitId()}.superChamberId().rawId() != primary_superchamber;}),
          ge0_links.end());
    }

    const GEMDetId superchamber_id{primary_superchamber};

    // only digitized hits are stored.
    // https://github.com/cms-sw/cmssw/blob/CMSSW_12_0_0_pre1/SimMuon/GEMDigitizer/src/GEMSignalModel.cc#L130-L162
    std::vector<ge0::GE0SimHit> hits;
    for (auto sim_hit = sim_hit_container->begin(); sim_hit != sim_hit_container->end(); sim_hit++) {
      if (sim_track->trackId() != sim_hit->trackId()) continue;
      if (sim_track->eventId() != sim_hit->eventId()) continue;

      const GEMDetId hit_gem_id{sim_hit->detUnitId()};
      if (hit_gem_id.station() != 0) continue;

      const LocalPoint hit_entry_point = sim_hit->entryPoint();
      const float hit_tof = sim_hit->timeOfFlight();

      std::vector<ge0::LinkData> cluster;
      for (const ge0::LinkData& link : ge0_links) {
        const GEMDetId link_gem_id{link->getDetUnitId()};
        if (hit_gem_id != link_gem_id) continue;
        if (not (hit_entry_point == link->getEntryPoint())) continue;
        if (hit_tof != link->getTimeOfFlight()) continue;
        cluster.push_back(link);
      } // links

      if (not cluster.empty()) {
        const GEMEtaPartition* eta_partition = gem->etaPartition(hit_gem_id);
        const int strip = static_cast<int>(eta_partition->strip(sim_hit->localPosition()));

        // sort links by strip
        std::sort(cluster.begin(), cluster.end());

        // FIXME sanity check? bx, adjacency..
        hits.emplace_back(sim_hit, hit_gem_id.rawId(), strip,
                          cluster.front()->getBx(), cluster);
      }
    } // simhit

    sim_segment_collection.emplace_back(sim_track, superchamber_id, hits);
  } // SimTrackContainer

  return sim_segment_collection;
}

bool GE0DatasetWriter::matchWithRecHit(const int bx, const int strip, const GEMRecHit& rechit) {
  // temporal matching
  if (bx != rechit.BunchX()) {
    return false;
  }

  // spatial matching
  const int rechit_first_strip = rechit.firstClusterStrip();
  const int rechit_last_strip = rechit_first_strip + rechit.clusterSize() - 1;
  return (strip >= rechit_first_strip) and (strip <= rechit_last_strip);
}

bool GE0DatasetWriter::matchWithRecHit(
    const ge0::GE0SimHit& simhit, const GEMRecHit& rechit) {
  return matchWithRecHit(simhit.getBx(), simhit.getStrip(), rechit);
}


std::pair<float, unsigned int> GE0DatasetWriter::computeEfficiency(
    const GEMSegmentCollection::const_iterator& rec_segment,
    const ge0::GE0SimSegment* sim_segment) {

  const std::vector<GEMRecHit>& rechit_collection = rec_segment->specificRecHits();
  const std::vector<ge0::GE0SimHit>& simhit_collection = sim_segment->hits();

  // Compute the segment-wise efficiency.
  // Efficiency = (# of matched simhits) / (# of simhits)
  unsigned int num_matched = 0;

  for (const auto & simhit : simhit_collection) {
    unsigned int simhit_det_id = simhit.getDetUnitId();

    for (const GEMRecHit& rechit : rechit_collection) {
      if (simhit_det_id != rechit.gemId().rawId()) {
        continue;
      }

      if (matchWithRecHit(simhit, rechit)) {
        num_matched++;
        break;
      }
    } // rechit
  } // digi

  const float efficiency = static_cast<float>(num_matched) / simhit_collection.size();
  return std::make_pair(efficiency, num_matched);
}

float GE0DatasetWriter::computeFakeHitRate(
    const GEMSegmentCollection::const_iterator& rec_segment,
    const ge0::GE0SimSegment* sim_segment) {

  const std::vector<GEMRecHit>& rechit_collection = rec_segment->specificRecHits();
  const std::vector<ge0::GE0SimHit>&& simhit_collection = sim_segment->hits();

  const int num_total = rec_segment->nRecHits();
  int num_fake = 0;

  for (const auto & rechit : rechit_collection) {
    const uint32_t rechit_id = rechit.gemId().rawId();

    bool is_not_matched = true;
    for (const auto & simhit : simhit_collection) {
      if (rechit_id != simhit.getDetUnitId()) continue;
      if (is_not_matched) {
        if (matchWithRecHit(simhit, rechit)) {
          is_not_matched = false;
          break;
        }
      }
    } // link

    if (is_not_matched) {
      num_fake++;
    }
  } // rechit

  // return fake rate
  return static_cast<float>(num_fake) / num_total;
}

std::vector<const ge0::GE0SimSegment*> GE0DatasetWriter::getSimSegmentsInSuperChamber(
    const std::vector<ge0::GE0SimSegment>& sim_segment_collection,
    const GEMDetId & superchamber_id) {

  std::vector<const ge0::GE0SimSegment*> gemini_simsegs;

  for (const auto& sim_segment : sim_segment_collection) {
    if (sim_segment.superChamberId() == superchamber_id) {
      gemini_simsegs.push_back(&sim_segment);
    }
  }

  // in the decreasing pt order
  if (gemini_simsegs.size() > 1) {
    std::sort(gemini_simsegs.begin(), gemini_simsegs.end(),
              [](const ge0::GE0SimSegment* lhs, const ge0::GE0SimSegment* rhs) {
                  return lhs->pt() > rhs->pt();});
  }

  return gemini_simsegs;
}

LocalPoint GE0DatasetWriter::getSuperChamberPosition(
    const LocalPoint& ieta_pos,
    const GEMSuperChamber* superchamber,
    const GEMEtaPartition* eta_partition) {
  return superchamber->toLocal(eta_partition->toGlobal(ieta_pos));
}

void GE0DatasetWriter::analyzeSuperChamber(
    const GEMSuperChamber* superchamber) {
  const GEMDetId & superchamber_id = superchamber->id();
  b_region_ = static_cast<long>(superchamber_id.region());
  b_station_ = static_cast<long>(superchamber_id.station());
  b_chamber_ = static_cast<long>(superchamber_id.chamber());
}

void GE0DatasetWriter::analyzeMuon(
    const std::vector<const ge0::GE0SimSegment*>& sim_segment_collection,
    const edm::ESHandle<GEMGeometry>& gem) {

  b_muon_size_ = static_cast<long>(sim_segment_collection.size());
  if (b_muon_size_ >= 1) {
    b_muon_hit_x_.resize(b_muon_size_);
    b_muon_hit_y_.resize(b_muon_size_);
    b_muon_hit_z_.resize(b_muon_size_);
    b_muon_hit_layer_.resize(b_muon_size_);
    b_muon_hit_ieta_.resize(b_muon_size_);
    b_muon_hit_strip_.resize(b_muon_size_);
    b_muon_hit_entry_strip_.resize(b_muon_size_);
    b_muon_hit_exit_strip_.resize(b_muon_size_);
    b_muon_hit_first_strip_.resize(b_muon_size_);
    b_muon_hit_cls_.resize(b_muon_size_);
    b_muon_hit_bx_.resize(b_muon_size_);
  }

  for (unsigned int muon_idx = 0; muon_idx < b_muon_size_; muon_idx++) {
    const ge0::GE0SimSegment* sim_segment = sim_segment_collection.at(muon_idx);
    const GEMSuperChamber* superchamber = gem->superChamber(sim_segment->superChamberId());
    const auto&& muon_ge0_hits = sim_segment->hits();
    const long nhits = static_cast<long>(muon_ge0_hits.size());

    b_muon_pt_.push_back(sim_segment->pt());
    b_muon_eta_.push_back(sim_segment->eta());
    b_muon_phi_.push_back(sim_segment->phi());
    b_muon_charge_.push_back(static_cast<long>(sim_segment->charge()));
    b_muon_hit_size_.push_back(nhits);

    b_muon_hit_x_.at(muon_idx).reserve(nhits);
    b_muon_hit_y_.at(muon_idx).reserve(nhits);
    b_muon_hit_z_.at(muon_idx).reserve(nhits);
    b_muon_hit_layer_.at(muon_idx).reserve(nhits);
    b_muon_hit_ieta_.at(muon_idx).reserve(nhits);
    b_muon_hit_strip_.at(muon_idx).reserve(nhits);
    b_muon_hit_entry_strip_.at(muon_idx).reserve(nhits);
    b_muon_hit_exit_strip_.at(muon_idx).reserve(nhits);
    b_muon_hit_first_strip_.at(muon_idx).reserve(nhits);
    b_muon_hit_cls_.at(muon_idx).reserve(nhits);
    b_muon_hit_bx_.at(muon_idx).reserve(nhits);

    for (const auto simhit : muon_ge0_hits) {
      const GEMDetId gem_id{simhit.getDetUnitId()};

      const GEMEtaPartition* eta_partition = gem->etaPartition(gem_id);
      const LocalPoint&& hit_pos = getSuperChamberPosition(
          simhit.getLocalPosition(), superchamber, eta_partition);
      const long strip = static_cast<long>(simhit.getStrip());
      const long entry_strip = static_cast<long>(eta_partition->strip(simhit.data()->entryPoint()));
      const long exit_strip = static_cast<long>(eta_partition->strip(simhit.data()->exitPoint()));
      const long first_strip = static_cast<long>(simhit.getCluster().front()->getStrip());
      const long cls = static_cast<long>(simhit.getClusterSize());

      b_muon_hit_x_.at(muon_idx).push_back(hit_pos.x());
      b_muon_hit_y_.at(muon_idx).push_back(hit_pos.y());
      b_muon_hit_z_.at(muon_idx).push_back(hit_pos.z());
      b_muon_hit_layer_.at(muon_idx).push_back(gem_id.layer());
      b_muon_hit_ieta_.at(muon_idx).push_back(gem_id.ieta());
      b_muon_hit_strip_.at(muon_idx).push_back(strip);
      b_muon_hit_entry_strip_.at(muon_idx).push_back(entry_strip);
      b_muon_hit_exit_strip_.at(muon_idx).push_back(exit_strip);
      b_muon_hit_first_strip_.at(muon_idx).push_back(first_strip);
      b_muon_hit_cls_.at(muon_idx).push_back(cls);
      b_muon_hit_bx_.at(muon_idx).push_back(cls);
    }

    std::vector<unsigned int> indices(b_muon_hit_x_.at(muon_idx).size());
    std::iota(indices.begin(), indices.end(), 0u);

    const std::vector<long>& layers = b_muon_hit_layer_.at(muon_idx);
    const auto layer_argminmax = std::minmax_element(
        indices.begin(), indices.end(),
        [&layers](const long lhs, const long rhs) -> bool {
            return layers.at(lhs) < layers.at(rhs);});
    const auto layer_argmin = *(layer_argminmax.first);
    const auto layer_argmax = *(layer_argminmax.second);

    const float delta_x = b_muon_hit_x_.at(muon_idx).at(layer_argmax) - b_muon_hit_x_.at(muon_idx).at(layer_argmin);
    const float delta_y = b_muon_hit_y_.at(muon_idx).at(layer_argmax) - b_muon_hit_y_.at(muon_idx).at(layer_argmin);
    const long delta_ieta = b_muon_hit_ieta_.at(muon_idx).at(layer_argmax) - b_muon_hit_ieta_.at(muon_idx).at(layer_argmin);
    const long delta_strip = b_muon_hit_strip_.at(muon_idx).at(layer_argmax) - b_muon_hit_strip_.at(muon_idx).at(layer_argmin);

    b_muon_delta_x_.push_back(delta_x);
    b_muon_delta_y_.push_back(delta_y);
    b_muon_delta_ieta_.push_back(delta_ieta);
    b_muon_delta_strip_.push_back(delta_strip);

    h_muon_delta_x_->Fill(delta_x);
    h_muon_delta_y_->Fill(delta_y);
    h_muon_delta_ieta_->Fill(std::clamp(delta_ieta, -2L, 2L));
    h_muon_delta_strip_->Fill(std::clamp(delta_strip, -10L, 10L));

  } // muon
}

bool GE0DatasetWriter::analyzeDigi(
    const edm::Handle<GEMDigiCollection>& digi_collection,
    const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >& link_collection,
    const std::vector<const ge0::GE0SimSegment*>& sim_segment_collection,
    const GEMSuperChamber* superchamber) {

  Digi2Index digi2idx;
  for (const GEMChamber* chamber : superchamber->chambers()) {
    for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
      const GEMDetId & gem_id = eta_partition->id();
      const long layer = static_cast<long>(gem_id.layer());
      const long ieta = static_cast<long>(gem_id.ieta());

      auto digi_range = digi_collection->get(gem_id);
      for (auto digi = digi_range.first; digi != digi_range.second; ++digi) {
        if (not digi->isValid()) {
          edm::LogError(kLogCategory_) << "got an invalid digi";
          return false;
        }

        const long strip = static_cast<long>(digi->strip());
        // NOTE ignore bx
        // const long bx = static_cast<long>(digi->bx());

        const Digi2Index::key_type key{layer, ieta, strip};
        if (digi2idx.find(key) == digi2idx.end()) {
          digi2idx.emplace(key, b_digi_layer_.size());

          b_digi_layer_.push_back(layer);
          b_digi_ieta_.push_back(ieta);
          b_digi_strip_.push_back(strip);

        } else {
          // ignore bx
          continue;

        }
      } // digi
    } // eta partition
  } // layer

  b_digi_layer_count_ = countUnique(b_digi_layer_);
  b_digi_size_ = b_digi_layer_.size();

  if (b_digi_layer_count_ < kMinNumLayers_) {
    edm::LogInfo(kLogCategory_) << "Too few digis. Skip this superchamber";
    return false;
  }

  // Labeling
  b_digi_label_.resize(b_digi_size_, 0L);
  for (unsigned int muon_idx = 0; muon_idx < sim_segment_collection.size(); muon_idx++) {
    auto sim_segment = sim_segment_collection.at(muon_idx);
    // label 0 means strips fired by background hits or noise
    const long label = static_cast<long>(muon_idx) + 1L;

    for (const auto& simhit : sim_segment->hits()) {
      const GEMDetId simhit_gem_id{simhit.getDetUnitId()};
      const Digi2Index::key_type digi_key{simhit_gem_id.layer(), simhit_gem_id.ieta(),
                                          simhit.getStrip()};

      if (digi2idx.find(digi_key) == digi2idx.end()) {
        edm::LogError(kLogCategory_) << "can't find a muon GEMDigiSimLink in GEMDigis";
        return false;
      }

      const unsigned int digi_idx = digi2idx.at(digi_key);
      b_digi_label_.at(digi_idx) = label;
    }
  }

  analyzeWindow(digi2idx);

  return true;
}

bool GE0DatasetWriter::analyzePad(
    const edm::Handle<GEMPadDigiCollection>& pad_collection,
    const std::vector<const ge0::GE0SimSegment*>& sim_segment_collection,
    const GEMSuperChamber* superchamber,
    const edm::ESHandle<GEMGeometry>& gem) {

  // for labeling
  // map<tuple<layer, ieta, pad>, index>
  Digi2Index pad2idx;

  for (const GEMChamber* chamber : superchamber->chambers()) {
    for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
      const GEMDetId & gem_id = eta_partition->id();
      const long layer = static_cast<long>(gem_id.layer());
      const long ieta = static_cast<long>(gem_id.ieta());

      auto pad_range = pad_collection->get(gem_id);
      for (auto pad_digi = pad_range.first; pad_digi != pad_range.second; ++pad_digi) {
        if (not pad_digi->isValid()) {
          edm::LogError(kLogCategory_) << "got an invalid pad digi";
          return false;
        }

        const long pad = static_cast<long>(pad_digi->pad());
        // NOTE ignore bx 
        // const long bx = static_cast<long>(pad_digi->bx());

        const Digi2Index::key_type key{layer, ieta, pad};
        if (pad2idx.find(key) == pad2idx.end()) {
          pad2idx.emplace(key, b_pad_layer_.size());

          b_pad_layer_.push_back(layer);
          b_pad_ieta_.push_back(ieta);
          b_pad_pad_.push_back(pad);

        } else {
          // NOTE ignore bx
          continue;

        }
      } // pad
    } // eta partition
  } // layer

  b_pad_layer_count_ = countUnique(b_pad_layer_);
  b_pad_size_ = static_cast<long>(b_pad_layer_.size());

  // Labeling
  b_pad_label_.resize(b_pad_size_, 0L);

  for (unsigned int muon_idx = 0; muon_idx < sim_segment_collection.size(); muon_idx++) {
    auto sim_segment = sim_segment_collection.at(muon_idx);

    // label 0 means strips fired by background hits or noise
    const long label = static_cast<long>(muon_idx) + 1L;

    for (const auto& simhit : sim_segment->hits()) {
      const GEMDetId simhit_gem_id{simhit.getDetUnitId()};
      const GEMEtaPartition* eta_partition = gem->etaPartition(simhit_gem_id);
      const int pad = static_cast<int>(eta_partition->padOfStrip(simhit.getStrip()));

      const Digi2Index::key_type pad_key{simhit_gem_id.layer(), simhit_gem_id.ieta(), pad};

      if (pad2idx.find(pad_key) == pad2idx.end()) {
        edm::LogError(kLogCategory_) << "no GEMPad associated to a muon digi";
        return false;
      }

      const unsigned int pad_idx = pad2idx.at(pad_key);
      if (b_pad_label_.at(pad_idx) == 0) {
        b_pad_label_.at(pad_idx) = label;

      } else if (b_pad_label_.at(pad_idx) == label) {
        continue;

      } else {
        edm::LogInfo(kLogCategory_) << "found multiple muons passing the same pad. skip this superchamber";
        return false;

      }

    } // muon digis
  } // GE0SimSegment

  return true;
}

bool GE0DatasetWriter::analyzeRecHit(
    const edm::Handle<GEMRecHitCollection>& rechit_collection,
    const std::vector<const ge0::GE0SimSegment*>& sim_segment_collection,
    const GEMSuperChamber* superchamber,
    const edm::ESHandle<GEMGeometry>& gem) {

   
  GE0Scaler scaler{superchamber};

  for (const GEMChamber* chamber : superchamber->chambers()) {
    for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
      const GEMDetId & gem_id = eta_partition->id();
      const long layer = static_cast<long>(gem_id.layer());
      const long ieta = static_cast<long>(gem_id.ieta());

      auto rechit_range = rechit_collection->get(gem_id);
      for (auto rechit = rechit_range.first; rechit != rechit_range.second; rechit++) {
        const LocalPoint&& local_pos = rechit->localPosition();
        const LocalPoint&& superchamber_pos = getSuperChamberPosition(
            local_pos, superchamber, eta_partition);
        // it should be fine to use LocalError on EtaPartition but I do transform
        // const LocalError&& local_error = rechit->localPositionError();
        const GlobalError&& global_err = ErrorFrameTransformer::transform(
            rechit->localPositionError(), eta_partition->surface());
        const LocalError&& superchamber_err = ErrorFrameTransformer::transform(
            global_err, superchamber->surface());
        const long strip = static_cast<long>(eta_partition->strip(local_pos));

        const float rechit_x = superchamber_pos.x();
        const float rechit_y = superchamber_pos.y();
        const float rechit_z = superchamber_pos.z();
        const float rechit_x_minmax = scaler.transformX(rechit_x);
        const float rechit_y_minmax = scaler.transformY(rechit_y);
        const float rechit_z_minmax = scaler.transformZ(rechit_z);

        b_rechit_layer_.push_back(layer);
        b_rechit_ieta_.push_back(ieta);
        b_rechit_strip_.push_back(strip);
        b_rechit_bx_.push_back(static_cast<long>(rechit->BunchX()));
        b_rechit_first_strip_.push_back(static_cast<long>(rechit->firstClusterStrip()));
        b_rechit_cls_.push_back(static_cast<long>(rechit->clusterSize()));
        b_rechit_x_.push_back(rechit_x);
        b_rechit_y_.push_back(rechit_y);
        b_rechit_z_.push_back(rechit_z);
        b_rechit_x_minmax_.push_back(rechit_x_minmax);
        b_rechit_y_minmax_.push_back(rechit_y_minmax);
        b_rechit_z_minmax_.push_back(rechit_z_minmax);
        b_rechit_err_xx_.push_back(superchamber_err.xx());
        b_rechit_err_xy_.push_back(superchamber_err.xy());
        b_rechit_err_yy_.push_back(superchamber_err.yy());

        // label 0 means strips fired by background hits or noise
        std::set<long> label_candidate;
        float dx = 1e4f, dy = 1e4f, sx = 1e4f, sy = 1e4f;
        for (unsigned int muon_idx = 0; muon_idx < sim_segment_collection.size(); muon_idx++) {
          auto muon = sim_segment_collection.at(muon_idx);
          const long label = static_cast<long>(muon_idx) + 1L; 

          for (const auto& simhit : muon->hits()) {
            const GEMDetId simhit_gem_id{simhit.getDetUnitId()};
            if (gem_id != simhit_gem_id) continue;

            if (matchWithRecHit(simhit, *rechit)) {
              label_candidate.insert(label);

              const LocalPoint& simhit_pos = getSuperChamberPosition(
                  simhit.getLocalPosition(), superchamber, eta_partition);

              dx = superchamber_pos.x() - simhit_pos.x();
              dy = superchamber_pos.y() - simhit_pos.y();
              sx = dx / std::sqrt(superchamber_err.xx());
              sy = dy / std::sqrt(superchamber_err.yy());

              break;
            }
          } // muon digis
        } // muon

        if (label_candidate.empty()) {
          b_rechit_label_.push_back(0L);
          b_rechit_label_first_ieta_.push_back(0L);
          b_rechit_label_first_strip_.push_back(0L);
          b_rechit_label_last_ieta_.push_back(0L);
          b_rechit_label_last_strip_.push_back(0L);

        } else if (label_candidate.size() == 1) {
          const long label = *label_candidate.begin();
          const long muon_idx = label - 1;

          const std::vector<long>& layers = b_muon_hit_layer_.at(muon_idx);
          const std::vector<long>& ietas = b_muon_hit_ieta_.at(muon_idx);
          const std::vector<long>& strips = b_muon_hit_strip_.at(muon_idx);

          std::vector<unsigned int> indices(layers.size());
          std::iota(indices.begin(), indices.end(), 0u);
          const auto layer_argminmax = std::minmax_element(
              indices.begin(), indices.end(),
              [&layers](const long lhs, const long rhs) -> bool {
                  return layers.at(lhs) < layers.at(rhs);});
          const auto layer_argmin = *(layer_argminmax.first);
          const auto layer_argmax = *(layer_argminmax.second);

          b_rechit_label_.push_back(label);
          b_rechit_label_first_ieta_.push_back(ietas.at(layer_argmin));
          b_rechit_label_first_strip_.push_back(strips.at(layer_argmin));
          b_rechit_label_last_ieta_.push_back(ietas.at(layer_argmax));
          b_rechit_label_last_strip_.push_back(strips.at(layer_argmax));

          h_hit_dx_->Fill(dx);
          h_hit_dy_->Fill(dy);
          h_hit_sx_->Fill(sx);
          h_hit_sy_->Fill(sy);

        } else {
          auto msg = Form("a rechit is matched with %lu muons",
                          label_candidate.size());
          edm::LogInfo(kLogCategory_) << msg;
          return false;

        }
      } // rechit
    } // eta partition
  } // layer

  b_rechit_size_ = static_cast<long>(b_rechit_layer_.size());
  b_rechit_layer_count_ = countUnique(b_rechit_layer_);

  return true;
}

bool GE0DatasetWriter::analyzeSegmentRU(
    const GEMSegmentCollection::range& rec_seg_range,
    const std::vector<const ge0::GE0SimSegment*>& sim_segment_collection,
    const edm::ESHandle<GEMGeometry>& gem) {

  b_ru_size_ = std::distance(rec_seg_range.first, rec_seg_range.second);

  b_ru_muon_idx_.reserve(b_ru_size_);
  b_ru_norm_chi2_.reserve(b_ru_size_);
  b_ru_rechit_size_.reserve(b_ru_size_);
  b_ru_num_matched_.reserve(b_ru_size_);
  b_ru_eff_.reserve(b_ru_size_);
  b_ru_fake_hit_rate_.reserve(b_ru_size_);
  b_ru_rechit_layer_.reserve(b_ru_size_);
  b_ru_rechit_ieta_.reserve(b_ru_size_);
  b_ru_rechit_strip_.reserve(b_ru_size_);
  b_ru_rechit_first_strip_.reserve(b_ru_size_);
  b_ru_rechit_cls_.reserve(b_ru_size_);
  b_ru_rechit_bx_.reserve(b_ru_size_);

  for (auto rec_segment = rec_seg_range.first; rec_segment != rec_seg_range.second; rec_segment++) {
    if (not rec_segment->isValid()) {
      edm::LogError(kLogCategory_) << "found a invalid GEMSegment. skip this superchamber";
      return false;
    }

    const long ru_rechit_size = static_cast<long>(rec_segment->nRecHits());
    const float norm_chi2 = static_cast<float>(rec_segment->chi2()) / rec_segment->degreesOfFreedom();

    b_ru_norm_chi2_.push_back(norm_chi2);
    b_ru_rechit_size_.push_back(ru_rechit_size);

    // the index of an associated muon (GE0SimSegment)
    long asso_muon_idx = -1;
    long num_matched = -1;
    float efficiency = 0.0f;
    float fake_hit_rate = 1.0f;

    for (unsigned int muon_idx = 0; muon_idx < sim_segment_collection.size(); muon_idx++) {
      const ge0::GE0SimSegment* sim_segment = sim_segment_collection.at(muon_idx);
      auto eff_result = computeEfficiency(rec_segment, sim_segment);
      if (eff_result.second >= kMinNumLayers_) {
        asso_muon_idx = static_cast<long>(muon_idx);
        std::tie(efficiency, num_matched) = eff_result;
        fake_hit_rate = computeFakeHitRate(rec_segment, sim_segment);
        break;
      }
    }

    b_ru_muon_idx_.push_back(asso_muon_idx);
    b_ru_num_matched_.push_back(static_cast<long>(num_matched));
    b_ru_eff_.push_back(efficiency);
    b_ru_fake_hit_rate_.push_back(fake_hit_rate);

    std::vector<long> rechit_layer;
    std::vector<long> rechit_ieta;
    std::vector<long> rechit_strip;
    std::vector<long> rechit_first_strip;
    std::vector<long> rechit_cls;
    std::vector<long> rechit_bx;

    rechit_layer.reserve(ru_rechit_size);
    rechit_ieta.reserve(ru_rechit_size);
    rechit_strip.reserve(ru_rechit_size);
    rechit_first_strip.reserve(ru_rechit_size);
    rechit_cls.reserve(ru_rechit_size);
    rechit_bx.reserve(ru_rechit_size);

    for (const auto & rechit : rec_segment->specificRecHits()) {
      const auto&& gem_id = rechit.gemId();
      const GEMEtaPartition* eta_partition = gem->etaPartition(gem_id);
      const long strip = static_cast<long>(eta_partition->strip(rechit.localPosition()));

      rechit_layer.push_back(gem_id.layer());
      rechit_ieta.push_back(gem_id.ieta());
      rechit_strip.push_back(strip);
      rechit_first_strip.push_back(rechit.firstClusterStrip());
      rechit_cls.push_back(rechit.clusterSize());
      rechit_bx.push_back(rechit.BunchX());
    }

    b_ru_rechit_layer_.push_back(rechit_layer);
    b_ru_rechit_ieta_.push_back(rechit_ieta);
    b_ru_rechit_strip_.push_back(rechit_strip);
    b_ru_rechit_first_strip_.push_back(rechit_first_strip);
    b_ru_rechit_cls_.push_back(rechit_cls);
    b_ru_rechit_bx_.push_back(rechit_bx);
  } // rec segment

  return true;
}

void GE0DatasetWriter::analyzeWindow(const Digi2Index& digi2idx) {
  // scaning for windows
  // finding the one with most hits in 3 strip window
  int max_ieta = 0;
  int max_nstrip = 0;
  int max_hits = 0;
  for (int ieta = 1; ieta <= 8; ++ieta) {
    for (int nstrip = 2; nstrip <= 383; ++nstrip) {

      int current_nhits=0;
      for (int win_nstrip = -1; win_nstrip < 2; ++win_nstrip) {
        for (int nlayer = 1; nlayer <= 6; ++nlayer) {            
          const Digi2Index::key_type key{nlayer, ieta, nstrip + win_nstrip};
          if (digi2idx.find(key) != digi2idx.end()) {
            current_nhits++;
          }
        }
      }

      if (current_nhits > max_hits){
        max_hits = current_nhits;
        max_nstrip = nstrip;
        max_ieta = ieta;
      }
    } // strip
  } // eta partition

  b_has_window_ = max_hits >= 1;
  if (b_has_window_) {
    // found max strip window center
    // saving window
    for (int win_ieta = 1; win_ieta < 4; ++win_ieta) {
      for (int win_nstrip = 1; win_nstrip < 11; ++win_nstrip) {
        for (int win_nlayer = 1; win_nlayer < 7; ++win_nlayer) {
          const int test_ieta = max_ieta + win_ieta -1;
          const int test_nstrip = max_nstrip + win_nstrip - 5;
          const int index_win = get3DImageIndexWindow(win_nlayer, win_ieta, win_nstrip);
          bool has_hit = false;
          bool has_muon_hit = false;
          // for padding
          if ((test_ieta > 0 and test_ieta < 9) and (test_nstrip > 0 and test_nstrip < 385) ) {
            const Digi2Index::key_type key{win_nlayer, test_ieta, test_nstrip};
            has_hit = digi2idx.find(key) != digi2idx.end();
            has_muon_hit = has_hit ? b_digi_label_.at(digi2idx.at(key)) > 0 : false;
          }

          b_window_image_[index_win] = has_hit;
          b_window_label_[index_win] = has_muon_hit;

        } //layer
      } // strip
    } // ieta

    b_window_strip_ = max_nstrip;
    b_window_ieta_ = max_ieta;
  }
}

void GE0DatasetWriter::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  edm::Handle<edm::SimTrackContainer> sim_track_container;
  event.getByToken(kSimTrackToken_, sim_track_container);
  if (not sim_track_container.isValid()) {
    edm::LogError(kLogCategory_) << "invalid SimTrackContainer";
    return;
  }

  edm::Handle<edm::PSimHitContainer> sim_hit_container;
  event.getByToken(kSimHitToken_, sim_hit_container);
  if (not sim_hit_container.isValid()) {
    edm::LogError(kLogCategory_) << "invalid PSimHitContainer" << endl;
    return;
  }

  edm::Handle<GEMDigiCollection> digi_collection;
  event.getByToken(kDigiToken_, digi_collection);
  if (not digi_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMDigiCollection";
    return;
  }

  edm::Handle<edm::DetSetVector<GEMDigiSimLink>> link_collection;
  event.getByToken(kLinkToken_, link_collection);
  if (not link_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMDigiSimLink";
    return;
  }

  edm::Handle<GEMPadDigiCollection> pad_collection;
  event.getByToken(kPadToken_, pad_collection);
  if (not pad_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMPadDigiCollection";
    return;
  }

  edm::Handle<GEMRecHitCollection> rechit_collection;
  event.getByToken(kRecHitToken_, rechit_collection);
  if (not rechit_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMRecHitCollection";
    return;
  }

  edm::Handle<GEMSegmentCollection> rec_segment_collection;
  event.getByToken(kSegmentToken_, rec_segment_collection);
  if (not rec_segment_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMSegmentCollection";
    return;
  }

  edm::ESHandle<GEMGeometry> gem = setup.getHandle(kGEMToken_);
  if (not gem.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMGeometry";
  }

  const auto&& sim_segment_collection = buildGE0SimSegments(
      sim_track_container, link_collection, sim_hit_container, gem);
  h_num_simseg_->Fill(std::min(static_cast<int>(sim_segment_collection.size()), 10));

  std::map<int, MinMaxScaler> scalers;

  for (const GEMStation* station : gem->stations()) {
    if (station->station() != 0) {
      continue;
    }

    for (const GEMSuperChamber* superchamber : station->superChambers()) {
      const auto&& gemini_simsegs = getSimSegmentsInSuperChamber(
          sim_segment_collection, superchamber->id());

      if (gemini_simsegs.size() > kMaxNumMuons_) {
        auto msg = Form("found %zu muons (capacity = %u). skip this superchamber.",
                        gemini_simsegs.size(), kMaxNumMuons_);
        edm::LogInfo(kLogCategory_) << msg;
        continue;
      }

      resetBranch();

      analyzeSuperChamber(superchamber);

      analyzeMuon(gemini_simsegs, gem);

      const bool digi_okay = analyzeDigi(digi_collection, link_collection,
                                         gemini_simsegs, superchamber);
      if (not digi_okay) continue;

      const bool pad_okay = analyzePad(pad_collection, gemini_simsegs,
                                       superchamber, gem);
      if (not pad_okay) continue;

      const bool rechit_okay = analyzeRecHit(rechit_collection, gemini_simsegs,
                                             superchamber, gem);
      if (not rechit_okay) continue;

      const bool ru_okay = analyzeSegmentRU(
          rec_segment_collection->get(superchamber->id()), gemini_simsegs, gem);
      if (not ru_okay) continue;

      tree_->Fill();

    } // chamber loop
  } // station
}

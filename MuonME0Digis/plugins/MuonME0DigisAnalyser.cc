#include "MuonTriggering/MuonME0Digis/plugins/MuonME0DigisAnalyser.h"
#include "MuonTriggering/MuonME0Digis/plugins/ME0MuonData.cc"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;


MuonME0DigisAnalyser::MuonME0DigisAnalyser(const edm::ParameterSet& pset) { 
  cout << "ctor begin" << endl;

  auto sim_track_tag = pset.getParameter<edm::InputTag>("simTrackTag");
  sim_track_token_ = consumes<edm::SimTrackContainer>(sim_track_tag);

  auto sim_hit_tag = pset.getParameter<edm::InputTag>("simHitTag");
  sim_hit_token_ = consumes<edm::PSimHitContainer>(sim_hit_tag);

  auto digi_tag = pset.getParameter<edm::InputTag>("me0DigiTag");
  me0_digi_token_ = consumes<ME0DigiCollection>(digi_tag);

  auto link_tag = pset.getParameter<edm::InputTag>("me0DigiSimLinkTag");
  me0_link_token_ = consumes<edm::DetSetVector<ME0DigiSimLink>>(link_tag);

  auto rechit_tag = pset.getParameter<edm::InputTag>("me0RecHitTag");
  me0_rechit_token_ = consumes<ME0RecHitCollection>(rechit_tag);

  auto segment_tag = pset.getParameter<edm::InputTag>("me0SegmentTag");
  me0_segment_token_ = consumes<ME0SegmentCollection>(segment_tag);

  // auto gen_tag = pset.getParameter<edm::InputTag>("genParticleTag");
  // gen_particle_token_ = consumes<reco::GenParticleCollection>(gen_tag);

  min_pt_ = pset.getParameter<double>("min_pt");
  min_quality_ = pset.getParameter<double>("min_quality");
  min_num_layers_ = pset.getParameter<unsigned int>("min_num_layers");
  min_digis_ = pset.getParameter<unsigned int>("min_digis");
  max_muons_ = pset.getParameter<unsigned int>("max_muons");

  setBranch();
  
  cout << "ctor end" << endl;
}


MuonME0DigisAnalyser::~MuonME0DigisAnalyser() {
  cout << "dtor begin" << endl;
  cout << "dtor end" << endl;
}


void MuonME0DigisAnalyser::setBranch() {
  cout << "setBranch end" << endl;

  // NOTE
  tree_ = file_service_->make<TTree>("chamber", "chamber");

  // digi
  tree_->Branch("digi_size", &b_digi_size_, "digi_size/I");
  tree_->Branch("digi_layer", "vector<int>", &b_digi_layer_);
  tree_->Branch("digi_ieta", "vector<int>", &b_digi_ieta_);
  tree_->Branch("digi_strip", "vector<int>", &b_digi_strip_);
  tree_->Branch("digi_label", "vector<int>", &b_digi_label_);
  tree_->Branch("digi_particle_type", "vector<int>", &b_digi_particle_type_);
  tree_->Branch("digi_track_id", "vector<int>", &b_digi_track_id_);
  tree_->Branch("digi_image", &b_digi_image_, "digi_image[18432]/O");
  tree_->Branch("digi_image_label", &b_digi_image_label_, "digi_image_label[18432]/I");

  // muon
  tree_->Branch("muon_size", &b_muon_size_, "muon_size/I"); // for compaitibility with multi_muon tree
  tree_->Branch("muon_pt", "vector<float>", &b_muon_pt_);
  tree_->Branch("muon_eta", "vector<float>", &b_muon_eta_);
  tree_->Branch("muon_phi", "vector<float>", &b_muon_phi_);
  tree_->Branch("muon_charge", "vector<int>", &b_muon_charge_);

  // rechit
  tree_->Branch("rechit_size", &b_rechit_size_, "rechit_size/I");
  tree_->Branch("rechit_layer", "vector<int>", &b_rechit_layer_);
  tree_->Branch("rechit_ieta", "vector<int>", &b_rechit_ieta_);
  tree_->Branch("rechit_strip", "vector<int>", &b_rechit_strip_);
  tree_->Branch("rechit_x", "vector<float>", &b_rechit_x_);
  tree_->Branch("rechit_y", "vector<float>", &b_rechit_y_);
  tree_->Branch("rechit_z", "vector<float>", &b_rechit_z_);
  tree_->Branch("rechit_tof", "vector<float>", &b_rechit_tof_);
  tree_->Branch("rechit_label", "vector<int>", &b_rechit_label_);

  // RoadUsage Reconstructed ME0 Segment
  tree_->Branch("ru_size", &b_ru_size_, "ru_size/I");
  tree_->Branch("ru_muon_idx", "vector<int>", &b_ru_muon_idx_);
  tree_->Branch("ru_reduced_chi2", "vector<float>", &b_ru_reduced_chi2_);
  tree_->Branch("ru_rechit_size", "vector<int>", &b_ru_rechit_size_);
  tree_->Branch("ru_rechit_layer", &b_ru_rechit_layer_);
  tree_->Branch("ru_rechit_ieta", &b_ru_rechit_ieta_);
  tree_->Branch("ru_rechit_strip", &b_ru_rechit_strip_);

  // additional position information
  tree_->Branch("region", &b_region_, "region/I");
  tree_->Branch("chamber", &b_chamber_, "chamber/I");

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Window
  //////////////////////////////////////////////////////////////////////////////
  tree_win_ = file_service_->make<TTree>("window", "window");

  tree_win_->Branch("digi_layer", "vector<int>", &b_digi_layer_);
  tree_win_->Branch("digi_ieta", "vector<int>", &b_digi_ieta_);
  tree_win_->Branch("digi_strip", "vector<int>", &b_digi_strip_);
  tree_win_->Branch("digi_label", "vector<int>", &b_digi_label_);

  tree_win_->Branch("muon_pt", &b_muon_pt_, "muon_pt/F");
  tree_win_->Branch("muon_eta", &b_muon_eta_, "muon_eta/F");
  tree_win_->Branch("muon_phi", &b_muon_phi_, "muon_phi/F");

  tree_win_->Branch("region", &b_region_, "region/I");
  tree_win_->Branch("chamber", &b_chamber_, "chamber/I");

  tree_win_->Branch("digi_image", &b_win_digi_image_, "digi_image[180]/O");
  tree_win_->Branch("digi_image_label", &b_win_digi_image_label_, "digi_image_label[180]/O");
  tree_win_->Branch("strip", &b_win_strip_, "strip/I");
  tree_win_->Branch("ieta", &b_win_ieta_, "ieta/I");

  cout << "setBranch end" << endl;
}

void MuonME0DigisAnalyser::bookHistogram() {
}


void MuonME0DigisAnalyser::resetBranch() {
  //
  b_digi_size_ = -1;
  b_digi_layer_.clear();
  b_digi_ieta_.clear();
  b_digi_strip_.clear();
  b_digi_label_.clear();
  b_digi_particle_type_.clear();
  b_digi_track_id_.clear();
  fill_n(b_digi_image_, 18432, false);
  fill_n(b_digi_image_label_, 18432, 0);

  // Muon
  b_muon_size_ = - 1;
  b_muon_pt_.clear();
  b_muon_eta_.clear();
  b_muon_phi_.clear();
  b_muon_charge_.clear();

  // RecHit
  b_rechit_size_ = -1;
  b_rechit_layer_.clear();
  b_rechit_ieta_.clear();
  b_rechit_strip_.clear();
  b_rechit_x_.clear();
  b_rechit_y_.clear();
  b_rechit_z_.clear();
  b_rechit_tof_.clear();
  b_rechit_label_.clear();

  // Road Usage
  b_ru_size_ = -1;
  b_ru_muon_idx_.clear();
  b_ru_reduced_chi2_.clear();
  b_ru_rechit_size_.clear();
  b_ru_rechit_layer_.clear();
  b_ru_rechit_ieta_.clear();
  b_ru_rechit_strip_.clear();

  b_region_ = -100;
  b_chamber_ = -100;

  // NOTE tree_window_
  // NOTE these branches doesn't need to be reset
  // fill_n(b_win_digi_, 180, false);
  // fill_n(b_win_muon_digi_, 180, false);
  // b_win_strip_ = -100;
  // b_win_ieta_ = -100;

}


int MuonME0DigisAnalyser::get3DImageIndex(int layer, int roll, int strip) {
  // return (8 * 384) * (layer - 1) + 384 * (roll - 1) + (strip - 1);
  // L: layer, R: roll, S: strip
  // index = (8 * 384)*(L - 1) + 384*(R - 1) + (S - 1)
  //       = 3072*L + 384*R + S - (3072 + 384 +1)
  //       = 3072*L + 384*R + S - 3457
  return 3072 * layer + 384 * roll + strip - 3457;
}


int MuonME0DigisAnalyser::get3DImageIndexWindow(int layer, int roll, int strip) {
  return 30 * layer + 10 * roll + strip - 41;
}


bool MuonME0DigisAnalyser::isSimTrackGood(
    edm::SimTrackContainer::const_iterator sim_track) {
  // TODO
  // if ((*t).noVertex() && !isMuonGun_)
  //   return false;
  // if ((*t).noGenpart() && !isMuonGun_)
  //   return false;

  if (abs(sim_track->type()) != 13) return false;
  if (sim_track->momentum().pt() < min_pt_) return false;

  return true;
}


bool MuonME0DigisAnalyser::isSimHitGood(edm::PSimHitContainer::const_iterator sim_hit) {
  if (abs(sim_hit->particleType()) != 13) return false;

  const EncodedEventId & event_id = sim_hit->eventId();
  if (event_id.event() != 0) return false;
  if (event_id.bunchCrossing() != 0) return false;

  if (sim_hit->processType() != 0) return false;

  return true;
}


std::tuple<bool, uint32_t, bool> MuonME0DigisAnalyser::isSimSegmentGood(
    const vector<edm::PSimHitContainer::const_iterator> & sim_segment) {
  bool is_good = false;
  uint32_t primary_chamber = 0; // FIXME
  bool need_to_prune = false;
  
  if (sim_segment.size() < min_num_layers_) {
    return std::make_tuple(is_good, primary_chamber, need_to_prune);
  }

  // 
  map<uint32_t, std::set<int> > layers_per_chamber;

  for (const auto & hit : sim_segment) {
    ME0DetId me0_id{hit->detUnitId()};
    layers_per_chamber[me0_id.chamberId().rawId()].insert(me0_id.layer());
  }

  set<int> layers;
  if (layers_per_chamber.size() == 1) {
    primary_chamber = layers_per_chamber.begin()->first;
    layers = layers_per_chamber.begin()->second;

  } else {
    auto tmp_max_elem = max_element(
        layers_per_chamber.begin(),
        layers_per_chamber.end(),
        [](const pair<uint32_t, std::set<int> >& p1, const pair<uint32_t, std::set<int> >&  p2) {
            return p1.second.size() < p2.second.size();});

    primary_chamber = tmp_max_elem->first;
    layers = tmp_max_elem->second;

    need_to_prune = true;
  }

  is_good = layers.size() >= min_num_layers_;

  return std::make_tuple(is_good, primary_chamber, need_to_prune);
}


map<pair<unsigned int, uint32_t>, ME0MuonData>
MuonME0DigisAnalyser::buildDB(
    const edm::Handle<edm::SimTrackContainer> & sim_track_container,
    const edm::Handle<edm::PSimHitContainer> & sim_hit_container,
    const edm::Handle<ME0DigiCollection> & me0_digi_collection,
    const edm::Handle<edm::DetSetVector<ME0DigiSimLink> >& link_set_vector,
    const edm::Handle<ME0SegmentCollection> & me0_segment_collection,
    const ME0Geometry* me0) {

  map<pair<unsigned int, uint32_t>, ME0MuonData> me0_muon_db;

  // NOTE
  for (auto sim_track = sim_track_container->begin(); sim_track != sim_track_container->end(); sim_track++) {
    if (not isSimTrackGood(sim_track)) continue;

    vector<edm::PSimHitContainer::const_iterator> sim_segment;

    for (auto sim_hit = sim_hit_container->begin(); sim_hit != sim_hit_container->end(); sim_hit++) {
      if (sim_track->trackId() != sim_hit->trackId()) continue;
      if (not isSimHitGood(sim_hit)) continue;

      ME0DetId me0_id{sim_hit->detUnitId()};
      auto eta_partition = me0->etaPartition(me0_id);
      int sim_hit_strip = static_cast<int>(ceil(eta_partition->strip(sim_hit->localPosition())));

      bool has_link = false;
      auto link_set = link_set_vector->find(me0_id);
      for (auto link = link_set->begin(); link != link_set->end(); link++) {
        if (sim_hit_strip != static_cast<int>(link->getStrip())) continue;
        if (sim_hit->particleType() != link->getParticleType()) continue;
        if (sim_hit->trackId() != link->getTrackId()) continue;

        has_link = true;
        break;
      }

      if (not has_link) continue;

      sim_segment.push_back(sim_hit);
    } // PSimHitContainer

    auto [is_sim_seg_good, primary_chamber, need_to_prune] = isSimSegmentGood(sim_segment);
    if (not is_sim_seg_good) continue;

    // TODO pruneSimSegment
    if (need_to_prune) {
      // Remove hits if they aren't on the primary chamber.
      // I don't think a single SimTrack will not have more good segments than one.
      sim_segment.erase(
          remove_if(
              sim_segment.begin(),
              sim_segment.end(),
              [&primary_chamber](const edm::PSimHitContainer::const_iterator& sim_hit) {
                  return ME0DetId(sim_hit->detUnitId()).chamberId().rawId() != primary_chamber;}),

          sim_segment.end());
    }

    auto chamber_id = ME0DetId(sim_segment[0]->detUnitId()).chamberId();

    ////////////////////////////////////////////////////////////////////////////
    // NOTE find reconstructed segment
    ////////////////////////////////////////////////////////////////////////////

    bool is_reconstructed = false;
    auto matched_rec_segment = me0_segment_collection->end();

    int num_matched = 0;
    set<int> found_layer;

    auto rec_segment_range = me0_segment_collection->get(chamber_id);
    for (auto rec_segment = rec_segment_range.first; rec_segment != rec_segment_range.second; rec_segment++) {

      for (const auto & sim_hit : sim_segment) {
        ME0DetId me0_id{sim_hit->detUnitId()};
        auto eta_partition = me0->etaPartition(me0_id);
        int strip = ceil(eta_partition->strip(sim_hit->localPosition()));

        for (const auto &  rechit : rec_segment->specificRecHits()) {
          // hit-wise matching condition
          if (me0_id != rechit.me0Id()) continue;
          int rechit_strip = ceil(eta_partition->strip(rechit.localPosition()));
          // matched ME0RecHit found
          if (strip == rechit_strip) {
            num_matched++;
            found_layer.insert(me0_id.layer());
            // exit specificRecHits loop
            break;
          }
        } // rechit
      } // sim_hit

      // NOTE quality of RecSegment
      double quality = static_cast<double>(num_matched) / sim_segment.size();
      if ((quality > min_quality_) and (found_layer.size() >= min_num_layers_)) {
        is_reconstructed = true;
        matched_rec_segment = rec_segment;
        // exit ME0SegmentCollection loop
        break;
      }
    } // ME0SegmentCollection

    unsigned int track_id = sim_track->trackId();

    ME0MuonData data{sim_track, sim_segment, is_reconstructed, matched_rec_segment};
    me0_muon_db.insert({{track_id, chamber_id.rawId()}, data});
  } // SimTrackContainer

  return me0_muon_db;
}


void MuonME0DigisAnalyser::analyze(const edm::Event& event, const edm::EventSetup& event_setup) {
  edm::Handle<edm::SimTrackContainer> sim_track_container;
  event.getByToken(sim_track_token_, sim_track_container);
  if (not sim_track_container.isValid()) {
    edm::LogError(kLogCategory_) << "invalid SimTrackContainer" << endl;
    return;
  }

  edm::Handle<edm::PSimHitContainer> sim_hit_container;
  event.getByToken(sim_hit_token_, sim_hit_container);
  if (not sim_hit_container.isValid()) {
    edm::LogError(kLogCategory_) << "invalid PSimHitContainer" << endl;
    return;
  }

  edm::Handle<ME0DigiCollection> me0_digi_collection;
  event.getByToken(me0_digi_token_, me0_digi_collection);
  if (not me0_digi_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0DigiCollection" << endl;
    return;
  }

  edm::Handle<edm::DetSetVector<ME0DigiSimLink>> link_set_vector;
  event.getByToken(me0_link_token_, link_set_vector);
  if (not link_set_vector.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0DigiSimLink" << endl;
    return;
  }

  edm::Handle<ME0RecHitCollection> me0_rechit_collection;
  event.getByToken(me0_rechit_token_, me0_rechit_collection);
  if (not me0_rechit_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0RecHitCollection" << endl;
    return;
  }

  edm::Handle<ME0SegmentCollection> me0_segment_collection;
  event.getByToken(me0_segment_token_, me0_segment_collection);
  if (not me0_segment_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0SegmentCollection" << endl;
    return;
  }

  edm::ESHandle<ME0Geometry> me0_handle;
  event_setup.get<MuonGeometryRecord>().get(me0_handle);
  if (not me0_handle.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0Geometry" << endl;
  }
  const ME0Geometry* me0 = &*me0_handle;

  auto me0_muon_db = buildDB(sim_track_container, sim_hit_container,
                            me0_digi_collection, link_set_vector,
                            me0_segment_collection, me0);

  for (const auto & chamber : me0->chambers()) {
    resetBranch();

    const ME0DetId & chamber_id = chamber->id();
    b_region_ = chamber_id.region();
    b_chamber_ = chamber_id.chamber();

    set<pair<unsigned int, uint32_t> > db_key_set;

    // for rechit labeling
    // map<tuple<layer, ieta, strip>, rechit_idx>
    map<tuple<int, int, int>, unsigned int> rechit_det2idx;

    for (const auto & me0_layer : chamber->layers()) {
      for (const auto & eta_partition : me0_layer->etaPartitions()) {
        const ME0DetId & me0_id = eta_partition->id();
        int layer = me0_id.layer();
        int ieta = me0_id.roll();

        // NOTE Fill ME0Digi
        auto digi_range = me0_digi_collection->get(me0_id);
        for (auto digi = digi_range.first; digi != digi_range.second; ++digi) {
          int strip = digi->strip();
          int index = get3DImageIndex(layer, ieta, strip);
          b_digi_image_[index] = true;
          b_digi_layer_.push_back(layer);
          b_digi_ieta_.push_back(ieta);
          b_digi_strip_.push_back(strip);

          // TODO findLink(digi, link)
          bool has_link = false;
          auto link_set = link_set_vector->find(me0_id);
          edm::DetSet<ME0DigiSimLink>::const_iterator link;
          for (link = link_set->begin(); link != link_set->end(); link++) {
            if (strip != static_cast<int>(link->getStrip())) continue;
            if (abs(link->getParticleType()) != 13) continue;

            has_link = true;
            break;
          }

          if (has_link) {
            unsigned int track_id = link->getTrackId();

            // debug
            b_digi_particle_type_.push_back(link->getParticleType());
            b_digi_track_id_.emplace_back(track_id);

            // check if digi 
            pair<unsigned int, uint32_t> db_key{track_id, chamber_id.rawId()};

            if (me0_muon_db.find(db_key) != me0_muon_db.end()) {
              db_key_set.insert(db_key); 
              me0_muon_db.at(db_key).appendDigiIndex(b_digi_layer_.size() - 1);
            }

          } else {
            // intrinsic noise or simulated bkg contribution
            b_digi_particle_type_.push_back(0);
            b_digi_track_id_.push_back(-1);
          }
        } // digi

        // NOTE Fill ME0RecHit
        auto rechit_range = me0_rechit_collection->get(me0_id);
        for (auto rechit = rechit_range.first; rechit != rechit_range.second; ++rechit) {
          if (me0->idToDet(rechit->me0Id()) == nullptr) {
            edm::LogError(kLogCategory_) << "ME0RecHit didn't matched with ME0Geometry." << endl;
            continue;
          }

          const LocalPoint&& local_pos = rechit->localPosition();
          const GlobalPoint&& global_pos = eta_partition->toGlobal(local_pos);
          const LocalPoint&& chamber_pos = chamber->toLocal(global_pos);
          int strip = static_cast<int>(ceil(eta_partition->strip(local_pos)));

          b_rechit_x_.push_back(chamber_pos.x());
          b_rechit_y_.push_back(chamber_pos.y());
          b_rechit_z_.push_back(chamber_pos.z());
          b_rechit_tof_.push_back(rechit->tof());

          b_rechit_layer_.push_back(layer);
          b_rechit_ieta_.push_back(ieta);
          b_rechit_strip_.push_back(strip);

          rechit_det2idx.insert({{layer, ieta, strip}, b_rechit_layer_.size()});

        }
      } // eta partition
    } // layer

    // NOTE
    if (b_digi_layer_.size() < min_digis_) continue;
    if (db_key_set.size() > max_muons_) continue;

    // NOTE
    if (b_digi_layer_.size() < min_digis_) continue;
    b_digi_size_ = static_cast<int>(b_digi_layer_.size());
    b_rechit_size_ = static_cast<int>(b_rechit_layer_.size());

    b_digi_label_.resize(b_digi_size_, 0);
    b_rechit_label_.resize(b_rechit_size_, 0);

    // NOTE muon
    vector<ME0MuonData> muons;
    for (const auto & key : db_key_set) {
      auto muon = me0_muon_db.at(key);
      muons.push_back(muon);
    }

    // in decreasing pT order
    if (muons.size() > 1) {
      sort(muons.begin(), muons.end(), [](ME0MuonData mu0, ME0MuonData mu1) {
        return mu0.Pt() > mu1.Pt();
      });
    }

    for (unsigned int muon_idx = 0; muon_idx < muons.size(); muon_idx++) {
      auto mu = muons[muon_idx];
      b_muon_pt_.push_back(mu.Pt());
      b_muon_eta_.push_back(mu.Eta());
      b_muon_phi_.push_back(mu.Phi());

      int muon_label = (muon_idx + 1);

      // digi labeling
      for (unsigned int digi_idx : mu.digi_indices()) {
        b_digi_label_[digi_idx] = muon_label;

        int layer = b_digi_layer_[digi_idx];
        int ieta = b_digi_ieta_[digi_idx];
        int strip = b_digi_strip_[digi_idx];

        int img_idx = get3DImageIndex(layer, ieta, strip);
        b_digi_image_label_[img_idx] = muon_label;

        auto rechit_det2idx_key = make_tuple(layer, ieta, strip);
        if (rechit_det2idx.find(rechit_det2idx_key) != rechit_det2idx.end()) {
          unsigned int rechit_idx = rechit_det2idx[rechit_det2idx_key];
          b_rechit_label_[rechit_idx] = muon_label;
        }

      }
    } // muons

    // NOTE RecSegment
    auto rec_seg_range = me0_segment_collection->get(chamber_id);
    b_ru_size_ = std::distance(rec_seg_range.first, rec_seg_range.second);

    b_ru_muon_idx_.reserve(b_ru_size_);
    b_ru_reduced_chi2_.reserve(b_ru_size_);
    b_ru_muon_idx_.reserve(b_ru_size_);
    b_ru_rechit_layer_.reserve(b_ru_size_);
    b_ru_rechit_ieta_.reserve(b_ru_size_);
    b_ru_rechit_strip_.reserve(b_ru_size_);

    for (auto rec_segment = rec_seg_range.first; rec_segment != rec_seg_range.second; rec_segment++) {
      int ru_rechit_size = rec_segment->nRecHits();
      float reduced_chi2 = static_cast<float>(rec_segment->chi2()) / rec_segment->degreesOfFreedom();

      b_ru_reduced_chi2_.push_back(reduced_chi2);
      b_ru_rechit_size_.push_back(ru_rechit_size);

      int asso_muon_idx = -1;
      for (unsigned int muon_idx = 0; muon_idx < muons.size(); muon_idx++) {
        auto mu = muons[muon_idx];
        if (mu.is_reconstructed() and ((rec_segment == mu.rec_segment()))) {
          asso_muon_idx = static_cast<int>(muon_idx);
        }
      }
      b_ru_muon_idx_.push_back(asso_muon_idx);

      //
      std::vector<int> rechit_layer;
      std::vector<int> rechit_ieta;
      std::vector<int> rechit_strip;
      rechit_layer.reserve(ru_rechit_size);
      rechit_ieta.reserve(ru_rechit_size);
      rechit_strip.reserve(ru_rechit_size);

      for (const auto &  rechit : rec_segment->specificRecHits()) {
        const auto&& me0_id = rechit.me0Id();
        auto eta_partition = me0->etaPartition(me0_id);
        int strip = static_cast<int>(ceil(eta_partition->strip(rechit.localPosition())));

        rechit_layer.push_back(me0_id.layer());
        rechit_ieta.push_back(me0_id.layer());
        rechit_strip.push_back(strip);
      }

      b_ru_rechit_layer_.push_back(rechit_layer);
      b_ru_rechit_ieta_.push_back(rechit_ieta);
      b_ru_rechit_strip_.push_back(rechit_strip);
    } // rec segment

    tree_->Fill();

    if (muons.size() <= 1) {
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
              int index = get3DImageIndex(nlayer, ieta, nstrip+win_nstrip);
              if (b_digi_image_[index]) current_nhits++;
            }
          }
          if (current_nhits > max_hits){
            max_hits = current_nhits;
            max_nstrip = nstrip;
            max_ieta = ieta;
          }
        } // strip
      } // eta partition

      if (max_hits > 0) {
        // found max strip window center
        // saving window
        for (int win_ieta = 1; win_ieta < 4; ++win_ieta) {
          for (int win_nstrip = 1; win_nstrip < 11; ++win_nstrip) {
            for (int win_nlayer = 1; win_nlayer < 7; ++win_nlayer) {
              int test_ieta = max_ieta + win_ieta -1;
              int test_nstrip = max_nstrip + win_nstrip - 5;
              int index_win = get3DImageIndexWindow(win_nlayer, win_ieta, win_nstrip);
              bool has_hit = false;
              bool has_muon_hit = false;
              // for padding
              if ((test_ieta > 0 and test_ieta < 9) and
                  (test_nstrip > 0 and test_nstrip < 385) ) {
                int index = get3DImageIndex(win_nlayer, test_ieta, test_nstrip);
                has_hit = b_digi_image_[index];
                has_muon_hit = b_digi_image_label_[index];
              }
              b_win_digi_image_[index_win] = has_hit;
              b_win_digi_image_label_[index_win] = has_muon_hit;
            }
          }
        }

        b_win_strip_ = max_nstrip;
        b_win_ieta_ = max_ieta;
        tree_win_->Fill();
      } // if (max_hits > 0)
    } // b_num <= 1
  } // chamber loop
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonME0DigisAnalyser);

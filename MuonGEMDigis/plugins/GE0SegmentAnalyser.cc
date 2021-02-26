#include "MuonTriggering/MuonGEMDigis/plugins/GE0SegmentAnalyser.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

GE0SegmentAnalyser::GE0SegmentAnalyser(const edm::ParameterSet& pset) { 
  cout << "ctor begin" << endl;

  auto sim_track_tag = pset.getParameter<edm::InputTag>("simTrackTag");
  sim_track_token_ = consumes<edm::SimTrackContainer>(sim_track_tag);

  auto sim_hit_tag = pset.getParameter<edm::InputTag>("simHitTag");
  sim_hit_token_ = consumes<edm::PSimHitContainer>(sim_hit_tag);

  auto digi_tag = pset.getParameter<edm::InputTag>("gemDigiTag");
  gem_digi_token_ = consumes<GEMDigiCollection>(digi_tag);

  auto link_tag = pset.getParameter<edm::InputTag>("gemDigiSimLinkTag");
  gem_link_token_ = consumes<edm::DetSetVector<GEMDigiSimLink>>(link_tag);

  auto gem_rechit_tag = pset.getParameter<edm::InputTag>("gemRecHitTag");
  gem_rechit_token_ = consumes<GEMRecHitCollection>(gem_rechit_tag);

  auto segment_tag = pset.getParameter<edm::InputTag>("gemSegmentTag");
  gem_segment_token_ = consumes<GEMSegmentCollection>(segment_tag);

  min_pt_ = pset.getParameter<double>("minPt");
  min_quality_ = pset.getParameter<double>("minQuality");
  min_num_layers_ = pset.getParameter<unsigned int>("minNumLayers");
  min_digis_ = pset.getParameter<unsigned int>("minDigis");
  max_muons_ = pset.getParameter<unsigned int>("maxMuons");

  beginFileService();
  
  cout << "ctor end" << endl;
}


void GE0SegmentAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("simTrackTag", edm::InputTag("g4SimHits"));
  desc.add<edm::InputTag>("simHitTag", edm::InputTag("g4SimHits", "MuonGEMHits"));
  desc.add<edm::InputTag>("gemDigiTag", edm::InputTag("simMuonGEMDigis"));
  desc.add<edm::InputTag>("gemDigiSimLinkTag", edm::InputTag("simMuonGEMDigis", "GEM"));
  desc.add<edm::InputTag>("gemRecHitTag", edm::InputTag("gemRecHits"));
  desc.add<edm::InputTag>("gemSegmentTag", edm::InputTag("gemSegments"));

  desc.add<double>("minPt", 5.0); // GeV
  desc.add<double>("minQuality", 0.6);
  desc.add<unsigned int>("minNumLayers", 4);
  desc.add<unsigned int>("minDigis", 3);
  desc.add<unsigned int>("maxMuons", 3);

  descriptions.add("GE0SegmentAnalyserDefault", desc);
}


GE0SegmentAnalyser::~GE0SegmentAnalyser() {
  cout << "dtor begin" << endl;
  cout << "dtor end" << endl;
}


void GE0SegmentAnalyser::beginFileService() {
  cout << "beginFileService end" << endl;

  // NOTE
  tree_ = file_service_->make<TTree>("chamber", "chamber");

  // digi
  tree_->Branch("digi_size", &b_digi_size_, "digi_size/L");
  tree_->Branch("digi_layer", "vector<long>", &b_digi_layer_);
  tree_->Branch("digi_ieta", "vector<long>", &b_digi_ieta_);
  tree_->Branch("digi_strip", "vector<long>", &b_digi_strip_);
  tree_->Branch("digi_label", "vector<long>", &b_digi_label_);
  tree_->Branch("digi_particle_type", "vector<long>", &b_digi_particle_type_);
  tree_->Branch("digi_track_id", "vector<long>", &b_digi_track_id_);

  // muon
  tree_->Branch("muon_size", &b_muon_size_, "muon_size/L"); // for compaitibility with multi_muon tree
  tree_->Branch("muon_pt", "vector<float>", &b_muon_pt_);
  tree_->Branch("muon_eta", "vector<float>", &b_muon_eta_);
  tree_->Branch("muon_phi", "vector<float>", &b_muon_phi_);

  // rechit
  tree_->Branch("rechit_size", &b_rechit_size_, "rechit_size/L");
  tree_->Branch("rechit_layer", "vector<long>", &b_rechit_layer_);
  tree_->Branch("rechit_ieta", "vector<long>", &b_rechit_ieta_);
  tree_->Branch("rechit_strip", "vector<long>", &b_rechit_strip_);
  tree_->Branch("rechit_first_strip", "vector<long>", &b_rechit_first_strip_);
  tree_->Branch("rechit_cls", "vector<long>", &b_rechit_cls_);
  tree_->Branch("rechit_x", "vector<float>", &b_rechit_x_);
  tree_->Branch("rechit_y", "vector<float>", &b_rechit_y_);
  tree_->Branch("rechit_z", "vector<float>", &b_rechit_z_);
  tree_->Branch("rechit_label", "vector<long>", &b_rechit_label_);

  // Road Usage Segment
  tree_->Branch("ru_size", &b_ru_size_, "ru_size/L");
  tree_->Branch("ru_muon_idx", "vector<long>", &b_ru_muon_idx_);
  tree_->Branch("ru_norm_chi2", "vector<float>", &b_ru_norm_chi2_);
  tree_->Branch("ru_rechit_size", "vector<long>", &b_ru_rechit_size_);
  tree_->Branch("ru_rechit_layer", &b_ru_rechit_layer_);
  tree_->Branch("ru_rechit_ieta", &b_ru_rechit_ieta_);
  tree_->Branch("ru_rechit_strip", &b_ru_rechit_strip_);
  tree_->Branch("ru_rechit_first_strip", &b_ru_rechit_first_strip_);
  tree_->Branch("ru_rechit_cls", &b_ru_rechit_cls_);

  // additional position information
  tree_->Branch("region", &b_region_, "region/L");
  tree_->Branch("station", &b_station_, "station/L");
  tree_->Branch("chamber", &b_chamber_, "chamber/L");

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Window
  //////////////////////////////////////////////////////////////////////////////
  tree_win_ = file_service_->make<TTree>("window", "window");

  tree_win_->Branch("digi_layer", "vector<long>", &b_digi_layer_);
  tree_win_->Branch("digi_ieta", "vector<long>", &b_digi_ieta_);
  tree_win_->Branch("digi_strip", "vector<long>", &b_digi_strip_);
  tree_win_->Branch("digi_label", "vector<long>", &b_digi_label_);

  tree_win_->Branch("muon_pt", &b_muon_pt_, "muon_pt/F");
  tree_win_->Branch("muon_eta", &b_muon_eta_, "muon_eta/F");
  tree_win_->Branch("muon_phi", &b_muon_phi_, "muon_phi/F");

  tree_win_->Branch("region", &b_region_, "region/L");
  tree_win_->Branch("station", &b_station_, "station/L");
  tree_win_->Branch("chamber", &b_chamber_, "chamber/L");

  tree_win_->Branch("digi_image", &b_win_digi_image_, "digi_image[180]/O");
  tree_win_->Branch("digi_image_label", &b_win_digi_image_label_, "digi_image_label[180]/O");
  tree_win_->Branch("strip", &b_win_strip_, "strip/L");
  tree_win_->Branch("ieta", &b_win_ieta_, "ieta/L");

}

void GE0SegmentAnalyser::bookHistogram() {
}

void GE0SegmentAnalyser::resetBranch() {
  //
  b_digi_size_ = -1L;
  b_digi_layer_.clear();
  b_digi_ieta_.clear();
  b_digi_strip_.clear();
  b_digi_label_.clear();
  b_digi_particle_type_.clear();
  b_digi_track_id_.clear();

  // Muon
  b_muon_size_ = -1L;
  b_muon_pt_.clear();
  b_muon_eta_.clear();
  b_muon_phi_.clear();

  // RecHit
  b_rechit_size_ = -1L;
  b_rechit_layer_.clear();
  b_rechit_ieta_.clear();
  b_rechit_strip_.clear();
  b_rechit_first_strip_.clear();
  b_rechit_cls_.clear();
  b_rechit_x_.clear();
  b_rechit_y_.clear();
  b_rechit_z_.clear();
  b_rechit_label_.clear();

  // Road Usage
  b_ru_size_ = -1L;
  b_ru_muon_idx_.clear();
  b_ru_norm_chi2_.clear();
  b_ru_rechit_size_.clear();
  b_ru_rechit_layer_.clear();
  b_ru_rechit_ieta_.clear();
  b_ru_rechit_strip_.clear();
  b_ru_rechit_first_strip_.clear();
  b_ru_rechit_cls_.clear();

  b_region_ = -100L;
  b_station_ = -100L;
  b_chamber_ = -100L;

  // NOTE tree_window_
  // NOTE these branches doesn't need to be reset
  // fill_n(b_win_digi_, 180, false);
  // fill_n(b_win_muon_digi_, 180, false);
  // b_win_strip_ = -100;
  // b_win_ieta_ = -100;

}

long GE0SegmentAnalyser::get3DImageIndexWindow(long layer, long roll, long strip) {
  return 30L * layer + 10L * roll + strip - 41L;
}

bool GE0SegmentAnalyser::isSimTrackGood(
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

bool GE0SegmentAnalyser::isSimHitGood(edm::PSimHitContainer::const_iterator sim_hit) {
  if (abs(sim_hit->particleType()) != 13) return false;
  if (sim_hit->processType() != 0) return false;

  const EncodedEventId & event_id = sim_hit->eventId();
  if (event_id.event() != 0) return false;
  if (event_id.bunchCrossing() != 0) return false;

  const GEMDetId gem_id{sim_hit->detUnitId()};
  if (gem_id.station() != 0) return false;

  return true;
}

// FIXME rename
std::tuple<bool, uint32_t, bool> GE0SegmentAnalyser::areSimSegmentHitsGood(
    const vector<edm::PSimHitContainer::const_iterator> & sim_segment_hits) {
  bool is_good = false;
  uint32_t primary_superchamber = 0; // FIXME
  bool need_to_prune = false;
 
  if (sim_segment_hits.size() < min_num_layers_) {
    return std::make_tuple(is_good, primary_superchamber, need_to_prune);
  }

  // <superchamber_rawid, layer>
  map<uint32_t, std::set<int> > layers_per_superchamber;

  for (const auto & hit : sim_segment_hits) {
    const GEMDetId gem_id{hit->detUnitId()};
    layers_per_superchamber[gem_id.superChamberId().rawId()].insert(gem_id.layer());
  }

  set<int> layers;
  if (layers_per_superchamber.size() == 1) {
    primary_superchamber = layers_per_superchamber.begin()->first;
    layers = layers_per_superchamber.begin()->second;

  } else {
    auto tmp_max_elem = max_element(
        layers_per_superchamber.begin(),
        layers_per_superchamber.end(),
        [](const pair<uint32_t, std::set<int> >& lhs, const pair<uint32_t, std::set<int> >&  rhs) {
            return lhs.second.size() < rhs.second.size();});

    primary_superchamber = tmp_max_elem->first;
    layers = tmp_max_elem->second;

    need_to_prune = true;
  }

  is_good = layers.size() >= min_num_layers_;

  return std::make_tuple(is_good, primary_superchamber, need_to_prune);
}

bool GE0SegmentAnalyser::matchWithRecHit(const int strip, const GEMRecHit& rechit) {
  const int first = rechit.firstClusterStrip();
  const int last = first + rechit.clusterSize() - 1;
  return (strip >= first) and (strip <= last);
}

GE0SimSegmentCollection GE0SegmentAnalyser::reconstructSimSegment(
    const edm::Handle<edm::SimTrackContainer> & sim_track_container,
    const edm::Handle<edm::PSimHitContainer> & sim_hit_container,
    const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >& link_set_vector,
    const edm::ESHandle<GEMGeometry>& gem) {

  // key == {track_id, superchamber raw id}
  GE0SimSegmentCollection sim_segment_collection;
  for (auto sim_track = sim_track_container->begin(); sim_track != sim_track_container->end(); sim_track++) {
    if (not isSimTrackGood(sim_track)) {
      continue;
    }

    std::vector<GE0SimSegment::SimHitData> sim_segment_hits;
    std::vector<GE0SimSegment::DigiData> sim_segment_digis; // = digitized simhits

    for (auto sim_hit = sim_hit_container->begin(); sim_hit != sim_hit_container->end(); sim_hit++) {
      if (sim_track->trackId() != sim_hit->trackId()) continue;
      if (not isSimHitGood(sim_hit)) continue;

      sim_segment_hits.push_back(sim_hit);
    } // PSimHitContainer

    auto [is_sim_seg_good, primary_superchamber, need_to_prune] = areSimSegmentHitsGood(sim_segment_hits);
    if (not is_sim_seg_good) {
      continue;
    }

    // move it to isSimSegmentGood
    // TODO pruneSimSegment
    if (need_to_prune) {
      // Remove hits if they aren't on the primary superchamber.
      // I don't think a single SimTrack will not have more good segments than one.
      // but need to check it
      sim_segment_hits.erase(
          remove_if(
              sim_segment_hits.begin(),
              sim_segment_hits.end(),
              [&primary_superchamber](const edm::PSimHitContainer::const_iterator& sim_hit) {
                  return GEMDetId(sim_hit->detUnitId()).superChamberId().rawId() != primary_superchamber;}),

          sim_segment_hits.end());
    }

    const auto superchamber_id = GEMDetId(sim_segment_hits[0]->detUnitId()).superChamberId();

    // NOTE 
    std::set<int> digi_layer_count;
    for (const auto sim_hit : sim_segment_hits) {
      const GEMDetId gem_id{sim_hit->detUnitId()};
      const auto eta_partition = gem->etaPartition(gem_id);
      const int sim_hit_strip = static_cast<int>(eta_partition->strip(sim_hit->localPosition()));

      const auto link_set = link_set_vector->find(gem_id);
      for (auto link = link_set->begin(); link != link_set->end(); link++) {
        if (sim_hit_strip != static_cast<int>(link->getStrip())) continue;
        if (sim_hit->particleType() != link->getParticleType()) continue;
        if (sim_hit->trackId() != link->getTrackId()) continue;

        digi_layer_count.insert(gem_id.layer());

        const GEMDigi digi(link->getStrip(), link->getBx());
        sim_segment_digis.emplace_back(gem_id, digi);

        break;
      }
    }

    if (digi_layer_count.size() < min_num_layers_) {
      // skip
      continue;
    }

    sim_segment_collection.emplace_back(sim_track, superchamber_id, sim_segment_hits, sim_segment_digis);
  } // SimTrackContainer

  return sim_segment_collection;
}

bool GE0SegmentAnalyser::associateRecSegToSimSeg(
    const GEMSegmentCollection::const_iterator& rec_segment,
    const GE0SimSegment* sim_segment,
    const edm::ESHandle<GEMGeometry>& gem) { // FIXME GEMGeometry to Handle

  int num_matched = 0; 
  std::set<int> found_layer;

  for (const auto & rechit : rec_segment->specificRecHits()) {
    const GEMDetId rechit_id = rechit.gemId();
    found_layer.insert(rechit_id.layer());
  }

  for (const auto & [digi_id, digi] : sim_segment->digis()) {
    for (const auto & rechit : rec_segment->specificRecHits()) {
      const GEMDetId rechit_id = rechit.gemId();

      if (rechit_id != digi_id) {
        continue;
      }

      if (matchWithRecHit(digi.strip(), rechit)) {
        num_matched++;
        break;
      }
    }
  }

  // NOTE quality of RecSegment
  const double quality = static_cast<double>(num_matched) / sim_segment->simHits().size();
  return (found_layer.size() >= min_num_layers_) and (quality >= min_quality_);
}

void GE0SegmentAnalyser::analyze(const edm::Event& event, const edm::EventSetup& event_setup) {
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

  edm::Handle<GEMDigiCollection> digi_collection;
  event.getByToken(gem_digi_token_, digi_collection);
  if (not digi_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMDigiCollection" << endl;
    return;
  }

  edm::Handle<edm::DetSetVector<GEMDigiSimLink>> link_set_vector;
  event.getByToken(gem_link_token_, link_set_vector);
  if (not link_set_vector.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMDigiSimLink" << endl;
    return;
  }

  edm::Handle<GEMRecHitCollection> rechit_collection;
  event.getByToken(gem_rechit_token_, rechit_collection);
  if (not rechit_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMRecHitCollection" << endl;
    return;
  }

  edm::Handle<GEMSegmentCollection> rec_segment_collection;
  event.getByToken(gem_segment_token_, rec_segment_collection);
  if (not rec_segment_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMSegmentCollection" << endl;
    return;
  }

  edm::ESHandle<GEMGeometry> gem;
  event_setup.get<MuonGeometryRecord>().get(gem);
  if (not gem.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMGeometry" << endl;
  }

  //
  const GE0SimSegmentCollection&& sim_segment_collection = reconstructSimSegment(
      sim_track_container, sim_hit_container, link_set_vector, gem);

  for (const GEMStation* station : gem->stations()) {
    if (station->station() != 0) {
      continue;
    }
    for (const GEMSuperChamber* superchamber : station->superChambers()) {
      const GEMDetId & superchamber_id = superchamber->id();

      // NOTE SimSegment
      std::vector<const GE0SimSegment*> gemini_sim_segment_collection;
      for (const auto& sim_segment : sim_segment_collection) {
        if (sim_segment.detId().rawId() == superchamber_id.rawId()) {
          gemini_sim_segment_collection.push_back(&sim_segment);
        }
      }
      if (gemini_sim_segment_collection.size() > max_muons_) {
        // TODO logging
        continue;
      }
      // in the decreasing pt order
      if (gemini_sim_segment_collection.size() > 1) {
        sort(gemini_sim_segment_collection.begin(),
             gemini_sim_segment_collection.end(),
             [](const GE0SimSegment* lhs, const GE0SimSegment* rhs) {
                return lhs->pt() > rhs->pt();
             });
      }

      resetBranch();

      b_region_ = static_cast<long>(superchamber_id.region());
      b_station_ = static_cast<long>(superchamber_id.station());
      b_chamber_ = static_cast<long>(superchamber_id.chamber());
      b_muon_size_ = static_cast<long>(gemini_sim_segment_collection.size());

      // for labeling
      // map<tuple<layer, ieta, strip>, digi_idx>
      map<tuple<int, int, int>, unsigned int> digi2idx;

      // Fill digis and rechits
      for (const GEMChamber* chamber : superchamber->chambers()) {
        for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
          const GEMDetId & gem_id = eta_partition->id();
          const long layer = static_cast<long>(gem_id.layer());
          const long ieta = static_cast<long>(gem_id.roll());

          // NOTE GEMDigi
          auto digi_range = digi_collection->get(gem_id);
          for (auto digi = digi_range.first; digi != digi_range.second; ++digi) {
            if (not digi->isValid()) {
              edm::LogError(kLogCategory_) << "got an invalid digi" << std::endl;
              continue;
            }
            const long strip = static_cast<long>(digi->strip());

            tuple<int, int, int> key(layer, ieta, strip);
            if (digi2idx.find(key) == digi2idx.end()) {
              digi2idx.insert({key, b_digi_layer_.size()});

              b_digi_layer_.push_back(layer);
              b_digi_ieta_.push_back(ieta);
              b_digi_strip_.push_back(strip);
            }
          } // digi

          // NOTE GEMRecHit
          auto rechit_range = rechit_collection->get(gem_id);
          for (auto rechit = rechit_range.first; rechit != rechit_range.second; rechit++) {
            if (gem->idToDet(rechit->gemId()) == nullptr) {
              edm::LogError(kLogCategory_) << "GEMRecHit didn't matched with GEMGeometry." << endl;
              continue;
            }

            const LocalPoint&& local_pos = rechit->localPosition();
            const GlobalPoint&& global_pos = eta_partition->toGlobal(local_pos);
            const LocalPoint&& superchamber_pos = superchamber->toLocal(global_pos);
            const long strip = static_cast<long>(eta_partition->strip(local_pos));

            b_rechit_layer_.push_back(layer);
            b_rechit_ieta_.push_back(ieta);
            b_rechit_strip_.push_back(strip);
            b_rechit_first_strip_.push_back(static_cast<long>(rechit->firstClusterStrip()));
            b_rechit_cls_.push_back(static_cast<long>(rechit->clusterSize()));

            b_rechit_x_.push_back(superchamber_pos.x());
            b_rechit_y_.push_back(superchamber_pos.y());
            b_rechit_z_.push_back(superchamber_pos.z());
          }
        } // eta partition
      } // layer


      //
      if (b_digi_layer_.size() < min_digis_) continue;
      b_digi_size_ = static_cast<long>(b_digi_layer_.size());
      b_rechit_size_ = static_cast<long>(b_rechit_layer_.size());

      //////////////////////////////////////////////////////////////////////////
      // Fill the particle type and the track id of GEMDigi using
      // GEMDigiSimLink
      //////////////////////////////////////////////////////////////////////////
      b_digi_particle_type_.resize(b_digi_size_, 0L);
      b_digi_track_id_.resize(b_digi_size_, 0L);

      for (const GEMChamber* chamber : superchamber->chambers()) {
        for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
          const auto link_set = link_set_vector->find(eta_partition->id());
          if (not link_set->empty()) {
            continue;
          }

          for (auto link = link_set->begin(); link != link_set->end(); link++) {
            const GEMDetId id{link->getDetUnitId()};
            const std::tuple<int, int, int> det{id.layer(), id.roll(), link->getStrip()};
            const unsigned int index = digi2idx[det];
            b_digi_particle_type_[index] = link->getParticleType();
            b_digi_track_id_[index] = link->getTrackId();
          }
        }
      }

      //////////////////////////////////////////////////////////////////////////
      // Put labels on digis and fill muon informations using GE0SimSegment
      //////////////////////////////////////////////////////////////////////////
      b_digi_label_.resize(b_digi_size_, 0L);

      for (unsigned int idx = 0; idx < gemini_sim_segment_collection.size(); idx++) {
        auto sim_segment = gemini_sim_segment_collection[idx];
        const long particle_type = static_cast<long>(sim_segment->type());
        const long track_id = static_cast<long>(sim_segment->trackId());

        b_muon_pt_.push_back(sim_segment->pt());
        b_muon_eta_.push_back(sim_segment->eta());
        b_muon_phi_.push_back(sim_segment->phi());

        // TODO documentation
        const long label = static_cast<long>(idx) + 1L;
        for (const auto& [id, digi] : sim_segment->digis()) {
          const std::tuple<int, int, int> det{id.layer(), id.roll(), digi.strip()};
          // FIXME check if digi2idx has det as a key.
          const unsigned int index = digi2idx[det];
          b_digi_label_[index] = label;
          b_digi_particle_type_[index] = particle_type;
          b_digi_track_id_[index] = track_id;
        }
      } // GE0SimSegment

      //////////////////////////////////////////////////////////////////////////
      // Put labels on rechits
      //////////////////////////////////////////////////////////////////////////
      b_rechit_label_.resize(b_rechit_size_, 0L);

      for (unsigned int rechit_index = 0; rechit_index < b_rechit_layer_.size(); rechit_index++) {
        const int layer = static_cast<int>(b_rechit_layer_[rechit_index]);
        const int ieta = static_cast<int>(b_rechit_ieta_[rechit_index]);
        const int first_strip = static_cast<int>(b_rechit_first_strip_[rechit_index]);
        const int last_strip = first_strip + static_cast<int>(b_rechit_cls_[rechit_index]);

        std::set<int> label_set;
        for (int strip = first_strip; strip <= last_strip; strip++) {
          const std::tuple<int, int, int> key(layer, ieta, strip);
          const unsigned int digi_index = digi2idx[key];
          const long digi_label = b_digi_label_[digi_index];
          if (digi_label != 0) {
            label_set.insert(digi_label);
          }
        }

        if (label_set.size() == 1) {
          b_rechit_label_[rechit_index] = (*label_set.begin());

        } else {
          // In the case of ambiguity, put 0 label on the rechit.
          b_rechit_label_[rechit_index] = 0;

        }
      } // rechit


      //////////////////////////////////////////////////////////////////////////
      // NOTE RU
      //////////////////////////////////////////////////////////////////////////
      auto rec_seg_range = rec_segment_collection->get(superchamber_id);
      b_ru_size_ = std::distance(rec_seg_range.first, rec_seg_range.second);

      b_ru_muon_idx_.reserve(b_ru_size_);
      b_ru_norm_chi2_.reserve(b_ru_size_);
      b_ru_muon_idx_.reserve(b_ru_size_);
      b_ru_rechit_layer_.reserve(b_ru_size_);
      b_ru_rechit_ieta_.reserve(b_ru_size_);
      b_ru_rechit_strip_.reserve(b_ru_size_);
      b_ru_rechit_first_strip_.reserve(b_ru_size_);
      b_ru_rechit_cls_.reserve(b_ru_size_);

      for (auto rec_segment = rec_seg_range.first; rec_segment != rec_seg_range.second; rec_segment++) {
        const long ru_rechit_size = static_cast<long>(rec_segment->nRecHits());
        const float norm_chi2 = static_cast<float>(rec_segment->chi2()) / rec_segment->degreesOfFreedom();

        b_ru_norm_chi2_.push_back(norm_chi2);
        b_ru_rechit_size_.push_back(ru_rechit_size);

        // the index of an associated muon (GE0SimSegment)
        long asso_muon_idx = -1;
        for (unsigned int idx = 0; idx < gemini_sim_segment_collection.size(); idx++) {
          const GE0SimSegment* sim_segment = gemini_sim_segment_collection[idx];
          if (associateRecSegToSimSeg(rec_segment, sim_segment, gem)) {
            asso_muon_idx = static_cast<long>(idx);
            break;
          }
        }
        b_ru_muon_idx_.push_back(asso_muon_idx);

        std::vector<int> rechit_layer;
        std::vector<int> rechit_ieta;
        std::vector<int> rechit_strip;
        std::vector<int> rechit_first_strip;
        std::vector<int> rechit_cls;
        rechit_layer.reserve(ru_rechit_size);
        rechit_ieta.reserve(ru_rechit_size);
        rechit_strip.reserve(ru_rechit_size);
        rechit_first_strip.reserve(ru_rechit_size);
        rechit_cls.reserve(ru_rechit_size);

        for (const auto &  rechit : rec_segment->specificRecHits()) {
          const auto&& gem_id = rechit.gemId();
          auto eta_partition = gem->etaPartition(gem_id);
          const int strip = static_cast<int>(eta_partition->strip(rechit.localPosition()));

          rechit_layer.push_back(gem_id.layer());
          rechit_ieta.push_back(gem_id.layer());
          rechit_strip.push_back(strip);
          rechit_first_strip.push_back(rechit.firstClusterStrip());
          rechit_cls.push_back(rechit.clusterSize());
        }

        b_ru_rechit_layer_.push_back(rechit_layer);
        b_ru_rechit_ieta_.push_back(rechit_ieta);
        b_ru_rechit_strip_.push_back(rechit_strip);
        b_ru_rechit_first_strip_.push_back(rechit_first_strip);
        b_ru_rechit_cls_.push_back(rechit_cls);
      } // rec segment

      tree_->Fill();

      //////////////////////////////////////////////////////////////////////////
      // window algorithm
      //////////////////////////////////////////////////////////////////////////
      if (b_muon_size_ <= 1) {
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
                const std::tuple<int, int, int> key(nlayer, ieta, nstrip + win_nstrip);
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

        if (max_hits < 1) {
          continue;
        }

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
                const std::tuple<int, int, int> key(win_nlayer, test_ieta, test_nstrip);
                has_hit = digi2idx.find(key) != digi2idx.end();
                has_muon_hit = has_hit ? b_digi_label_[digi2idx[key]] > 0 : false;
              }

              b_win_digi_image_[index_win] = has_hit;
              b_win_digi_image_label_[index_win] = has_muon_hit;

            } //layer
          } // strip
        } // ieta

        b_win_strip_ = max_nstrip;
        b_win_ieta_ = max_ieta;
        tree_win_->Fill();
      } // b_num <= 1, window algorithm


    } // chamber loop
  } // station
}

//define this as a plug-in
DEFINE_FWK_MODULE(GE0SegmentAnalyser);

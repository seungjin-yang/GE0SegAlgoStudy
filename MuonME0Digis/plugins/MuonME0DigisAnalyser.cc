#include "MuonTriggering/MuonME0Digis/plugins/MuonME0DigisAnalyser.h"
#include "MuonTriggering/MuonME0Digis/plugins/ME0MuonData.h"

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

  auto segment_tag = pset.getParameter<edm::InputTag>("me0SegmentTag");
  me0_segment_token_ = consumes<ME0SegmentCollection>(segment_tag);

  pt_min_ = pset.getParameter<double>("pt_min");

  setBranch();
  
  cout << "ctor end" << endl;
}


MuonME0DigisAnalyser::~MuonME0DigisAnalyser() {
  cout << "dtor begin" << endl;
  cout << "dtor end" << endl;
}


void MuonME0DigisAnalyser::setBranch() {
  cout << "setBranch end" << endl;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE
  //////////////////////////////////////////////////////////////////////////////
  tree_ch_ = file_service_->make<TTree>("chamber", "chamber");

  // digi
  tree_ch_->Branch("num_digi", &b_num_digi_, "num_digi/I");
  tree_ch_->Branch("digi", &b_digi_, "digi[18432]/O");
  tree_ch_->Branch("digi_layer", "vector<int>", &b_digi_layer_);
  tree_ch_->Branch("digi_ieta", "vector<int>", &b_digi_ieta_);
  tree_ch_->Branch("digi_strip", "vector<int>", &b_digi_strip_);
  tree_ch_->Branch("digi_is_muon", "vector<int>", &b_digi_is_muon_);
  tree_ch_->Branch("digi_particle_type", "vector<int>", &b_digi_particle_type_);
  tree_ch_->Branch("digi_track_id", "vector<int>", &b_digi_track_id_);

  // muon digi
  tree_ch_->Branch("has_muon", &b_has_muon_, "has_muon/O");
  tree_ch_->Branch("num_muon_digi", &b_num_muon_digi_, "num_muon_digi/I");
  tree_ch_->Branch("muon_digi", &b_muon_digi_, "muon_digi[18432]/O");
  tree_ch_->Branch("muon_digi_layer", "vector<int>", &b_muon_digi_layer_);
  tree_ch_->Branch("muon_digi_ieta", "vector<int>", &b_muon_digi_ieta_);
  tree_ch_->Branch("muon_digi_strip", "vector<int>", &b_muon_digi_strip_);
  tree_ch_->Branch("muon_pt", &b_muon_pt_, "muon_pt/F");
  tree_ch_->Branch("muon_eta", &b_muon_eta_, "muon_eta/F");
  tree_ch_->Branch("muon_phi", &b_muon_phi_, "muon_phi/F");

  // Reconstructed ME0 Segment
  tree_ch_->Branch("has_ru", &b_has_ru_, "has_ru/O");
  tree_ch_->Branch("num_ru", &b_num_ru_, "num_ru/I");
  // 
  tree_ch_->Branch("has_ru_asso", &b_has_ru_asso_, "has_ru_asso/O");
  tree_ch_->Branch("ru_asso_nhits", &b_ru_asso_nhits_, "ru_asso_nhits/I");
  tree_ch_->Branch("ru_asso_rechit_layer", "vector<int>", &b_ru_asso_rechit_layer_);
  tree_ch_->Branch("ru_asso_rechit_ieta", "vector<int>", &b_ru_asso_rechit_ieta_);
  tree_ch_->Branch("ru_asso_rechit_strip", "vector<int>", &b_ru_asso_rechit_strip_);
  // fake
  tree_ch_->Branch("num_ru_fake", &b_num_ru_fake_, "num_ru_fake/I");
  tree_ch_->Branch("ru_fake_nhits", "vector<int>", &b_ru_fake_nhits_);

  // additional position information
  tree_ch_->Branch("region", &b_region_, "region/I");
  tree_ch_->Branch("chamber", &b_chamber_, "chamber/I");

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Window
  //////////////////////////////////////////////////////////////////////////////
  tree_win_ = file_service_->make<TTree>("window", "window");

  tree_win_->Branch("digi_layer", "vector<int>", &b_digi_layer_);
  tree_win_->Branch("digi_ieta", "vector<int>", &b_digi_ieta_);
  tree_win_->Branch("digi_strip", "vector<int>", &b_digi_strip_);

  tree_win_->Branch("muon_digi_layer", "vector<int>", &b_muon_digi_layer_);
  tree_win_->Branch("muon_digi_ieta", "vector<int>", &b_muon_digi_ieta_);
  tree_win_->Branch("muon_digi_strip", "vector<int>", &b_muon_digi_strip_);

  tree_win_->Branch("muon_pt", &b_muon_pt_, "muon_pt/F");
  tree_win_->Branch("muon_eta", &b_muon_eta_, "muon_eta/F");
  tree_win_->Branch("muon_phi", &b_muon_phi_, "muon_phi/F");

  tree_win_->Branch("region", &b_region_, "region/I");
  tree_win_->Branch("chamber", &b_chamber_, "chamber/I");

  tree_win_->Branch("digi", &b_win_digi_, "digi[180]/O");
  tree_win_->Branch("muon_digi", &b_win_muon_digi_, "muon_digi[180]/O");
  tree_win_->Branch("strip", &b_win_strip_, "strip/I");
  tree_win_->Branch("ieta", &b_win_ieta_, "ieta/I");


  //////////////////////////////////////////////////////////////////////////////
  // NOTE multi-muon
  //////////////////////////////////////////////////////////////////////////////
  tree_multi_ = file_service_->make<TTree>("multi_muon", "multi_muon");

  tree_multi_->Branch("digi", &b_digi_, "digi[18432]/O");
  tree_multi_->Branch("digi_layer", "vector<int>", &b_digi_layer_);
  tree_multi_->Branch("digi_ieta", "vector<int>", &b_digi_ieta_);
  tree_multi_->Branch("digi_strip", "vector<int>", &b_digi_strip_);
  tree_multi_->Branch("num_digi", &b_num_digi_, "num_digi/I");

  tree_multi_->Branch("num_muon", &b_multi_num_muon_, "num_muon/I");
  tree_multi_->Branch("muon_pt", &b_multi_muon_pt_, "muon_pt[5]/F");
  tree_multi_->Branch("muon_eta", &b_multi_muon_eta_, "muon_eta[5]/F");
  tree_multi_->Branch("muon_phi", &b_multi_muon_phi_, "muon_phi[5]/F");
  tree_multi_->Branch("muon_num_digi", &b_multi_muon_num_digi_,
                      "muon_num_digi[5]/I");
  tree_multi_->Branch("muon_digi_idx", &b_multi_muon_digi_idx_,
                      "muon_digi_idx[5][20]/I");

  tree_multi_->Branch("num_ru_asso", &b_multi_num_ru_asso_, "num_ru_asso/I");
  tree_multi_->Branch("ru_asso_nhits", "vector<int>", &b_multi_ru_asso_nhits_);
  tree_multi_->Branch("ru_asso_muon_idx", "vector<int>",
                      &b_multi_ru_asso_muon_idx_);


  tree_multi_->Branch("num_ru_fake", &b_multi_num_ru_fake_, "num_ru_fake/I");
  tree_multi_->Branch("ru_fake_nhits", "vector<int>", &b_multi_ru_fake_nhits_);


  //////////////////////////////////////////////////////////////////////////////
  // NOTE Histograms
  //////////////////////////////////////////////////////////////////////////////
  h_sim_track_pt_ = file_service_->make<TH1F>(
    "h_sim_track_pt", "", 20, 0.0, 100.0);

  h_sim_track_eta_ = file_service_->make<TH1F>(
    "h_sim_track_eta", "", 8, 2.0, 2.8);

  // ME0 Chamber1 Phi -10 rad ~ + 10 rad
  h_sim_track_phi_ = file_service_->make<TH1F>(
    "h_sim_track_phi", "", 36, -M_PI, M_PI);

  h_matched_sim_track_pt_ = file_service_->make<TH1F>(
    "h_matched_sim_track_pt", "", 20, 0.0, 100.0);

  h_matched_sim_track_eta_ = file_service_->make<TH1F>(
    "h_matched_sim_track_eta", "", 8, 2.0, 2.8);

  h_matched_sim_track_phi_ = file_service_->make<TH1F>(
    "h_matched_sim_track_phi", "", 36, -M_PI, M_PI);

  h_num_simhit_ = file_service_->make<TH1F>(
    "h_num_simhit", "", 17, -0.5, 16.5);

  h_num_muon_ = file_service_->make<TH1F>(
    "h_num_muon", "", 11, -0.5, 11.5);

  h_num_rec_seg_ = file_service_->make<TH1F>(
    "h_num_rec_seg", "", 21, -0.5, 21.5); 

  h_num_sim_rec_ = file_service_->make<TH2F>(
    "h_num_sim_seg", ";sim;rec",
    11, -0.5, 10.5,
    11, -0.5, 10.5);


  // TODO
  // vs rechit
  // min hits
  //   - Average : 3
  //   - RU: 4

  // 0: fine
  // 1: too many muons (> 5)
  // 2: too many digis per muon (> 20)
  //
  h_stats_ = file_service_->make<TH1F>(
      "h_stats", "", 3, -0.5, 2.5);
  h_stats_->GetXaxis()->SetBinLabel(1, "Good");
  h_stats_->GetXaxis()->SetBinLabel(2, "N_{#mu} > 5");
  h_stats_->GetXaxis()->SetBinLabel(3, "N_{DIGI} > 20");

  cout << "setBranch end" << endl;
}


void MuonME0DigisAnalyser::resetBranch() {
  // NOTE tree_ch_
  fill_n(b_digi_, 18432, false);
  b_digi_layer_.clear();
  b_digi_ieta_.clear();
  b_digi_strip_.clear();
  b_num_digi_ = 0;

  b_digi_particle_type_.clear();
  b_digi_track_id_.clear();
  b_digi_is_muon_.clear();

  b_has_muon_ = false;
  fill_n(b_muon_digi_, 18432, false);
  b_muon_digi_layer_.clear();
  b_muon_digi_ieta_.clear();
  b_muon_digi_strip_.clear();
  b_num_muon_digi_ = 0;
  b_muon_pt_ = -100.0f;
  b_muon_eta_ = -100.0f;
  b_muon_phi_ = -100.0f;

  b_has_ru_ = false;
  b_num_ru_ = 0;

  b_has_ru_asso_ = false;
  b_ru_asso_nhits_ = 0;
  b_ru_asso_rechit_layer_.clear();
  b_ru_asso_rechit_ieta_.clear();
  b_ru_asso_rechit_strip_.clear();

  b_num_ru_fake_ = 0;
  b_ru_fake_nhits_.clear();

  b_region_ = -100;
  b_chamber_ = -100;

  // NOTE tree_window_
  // NOTE these branches doesn't need to be reset
  // fill_n(b_win_digi_, 180, false);
  // fill_n(b_win_muon_digi_, 180, false);
  // b_win_strip_ = -100;
  // b_win_ieta_ = -100;


  // NOTE tree_multi_
  b_multi_num_muon_ = 0;

  fill_n(b_multi_muon_pt_, 5, -100.0f);
  fill_n(b_multi_muon_eta_, 5, -100.0f);
  fill_n(b_multi_muon_phi_, 5, -100.0f);
  fill_n(b_multi_muon_num_digi_, 5, -1);

  fill_n(*b_multi_muon_digi_idx_, 100, -1);

  b_multi_num_ru_asso_ = 0;
  b_multi_ru_asso_nhits_.clear();
  b_multi_ru_asso_muon_idx_.clear();

  b_multi_num_ru_fake_ = 0;
  b_multi_ru_fake_nhits_.clear();

}


int MuonME0DigisAnalyser::getIndex(int layer, int roll, int strip) {
  // return (8 * 384) * (layer - 1) + 384 * (roll - 1) + (strip - 1);
  // L: layer, R: roll, S: strip
  // index = (8 * 384)*(L - 1) + 384*(R - 1) + (S - 1)
  //       = 3072*L + 384*R + S - (3072 + 384 +1)
  //       = 3072*L + 384*R + S - 3457
  return 3072 * layer + 384 * roll + strip - 3457;
}


int MuonME0DigisAnalyser::getIndexWindow(int layer, int roll, int strip) {
  return 30 * layer + 10 * roll + strip - 41;
}


int MuonME0DigisAnalyser::getUniqueId(int region, int chamber, int layer,
                                      int roll, int strip) {
  // TODO make it efficient
  int id = 0;
  id += (18 * 6 * 8 * 384) * ((region + 1) / 2);
  id += (6 * 8 * 384) * (chamber - 1);
  id += (8 * 384) * (layer - 1);
  id += 384 * (roll - 1);
  id += (strip - 1);

  return id;
}


int MuonME0DigisAnalyser::getUniqueId(const ME0DetId & det_id, int strip) {
  return getUniqueId(det_id.region(), det_id.chamber(), det_id.layer(),
                     det_id.roll(), strip);
}


bool MuonME0DigisAnalyser::isSimTrackGood(
    edm::SimTrackContainer::const_iterator sim_track) {
  // TODO
  // if ((*t).noVertex() && !isMuonGun_)
  //   return false;
  // if ((*t).noGenpart() && !isMuonGun_)
  //   return false;

  if (abs(sim_track->type()) != kMuonPDGId_) return false;
  // default min pt is 5 GeV
  if (sim_track->momentum().pt() < pt_min_) return false;

  return true;
}


bool MuonME0DigisAnalyser::isSimHitGood(
    edm::PSimHitContainer::const_iterator sim_hit) {

  if (abs(sim_hit->particleType()) != kMuonPDGId_) return false;

  const EncodedEventId & event_id = sim_hit->eventId();
  if (event_id.event() != 0) return false;
  if (event_id.bunchCrossing() != 0) return false;

  if (sim_hit->processType() != 0) return false;

  return true;
}


bool MuonME0DigisAnalyser::isSimSegmentGood(
    const vector<edm::PSimHitContainer::const_iterator> & sim_segment) {

  if (sim_segment.size() < 4) return false;

  set<int> chambers;
  set<int> layers;
  for (const auto & hit : sim_segment) {
    ME0DetId me0_id{hit->detUnitId()};

    chambers.insert(me0_id.chamber());
    layers.insert(me0_id.layer());
  }

  if (chambers.size() != 1) return false;
  if (layers.size() < 4) return false;

  return true;
}


void MuonME0DigisAnalyser::analyze(const edm::Event& event,
                                   const edm::EventSetup& event_setup) {
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


  //////////////////////////////////////////////////////////////////////////////
  // NOTE SimTrack
  //////////////////////////////////////////////////////////////////////////////
  map<unsigned int, ME0MuonData> me0_muon_db;

  // NOTE
  for (auto sim_track = sim_track_container->begin();
            sim_track != sim_track_container->end();
            sim_track++) {
    if (not isSimTrackGood(sim_track)) continue;

    vector<edm::PSimHitContainer::const_iterator> sim_segment;

    for (auto sim_hit = sim_hit_container->begin();
              sim_hit != sim_hit_container->end();
              sim_hit++) {
      if (sim_track->trackId() != sim_hit->trackId()) continue;
      if (not isSimHitGood(sim_hit)) continue;

      sim_segment.push_back(sim_hit);
    } // PSimHitContainer

    if (not isSimSegmentGood(sim_segment)) continue;

    h_num_simhit_->Fill(sim_segment.size());
    auto chamber_id = ME0DetId(sim_segment[0]->detUnitId()).chamberId();

    ////////////////////////////////////////////////////////////////////////////
    // NOTE find reconstructed segment
    ////////////////////////////////////////////////////////////////////////////

    bool found = false;
    auto matched_rec_segment = me0_segment_collection->end();
    int num_matched = 0;
    set<int> found_layer;

    auto range = me0_segment_collection->get(chamber_id);
    for (auto segment = range.first; segment != range.second; segment++) {

      for (const auto & sim_hit : sim_segment) {
        ME0DetId me0_id{sim_hit->detUnitId()};
        auto eta_partition = me0->etaPartition(me0_id);
        int strip = ceil(eta_partition->strip(sim_hit->localPosition()));

        for (const auto &  rechit : segment->specificRecHits()) {

          // NOTE hit-wise matching condition
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
      } // simhit

      // NOTE quality of RecSegment
      float quality = static_cast<float>(num_matched) / sim_segment.size();
      if ((quality > 0.6f) and (found_layer.size() >= 4)) {
        found = true;
        matched_rec_segment = segment;
        // exit ME0SegmentCollection loop
        break;
      }
    } // ME0SegmentCollection

    // NOTE
    unsigned int track_id = sim_track->trackId();
    // do we neet to check return of insert?
    me0_muon_db.insert(
        {track_id, {sim_track, sim_segment, found, matched_rec_segment, {}}});

    // NOTE
    const auto & momentum = sim_track->momentum();
    h_sim_track_pt_->Fill(momentum.Pt());
    h_sim_track_eta_->Fill(fabs(momentum.Eta()));
    h_sim_track_phi_->Fill(momentum.Phi());
    if (found) {
      h_matched_sim_track_pt_->Fill(momentum.Pt());
      h_matched_sim_track_eta_->Fill(fabs(momentum.Eta()));
      h_matched_sim_track_phi_->Fill(momentum.Phi());
    }

  } // SimTrackContainer


  //////////////////////////////////////////////////////////////////////////////
  // NOTE ME0DigiSimLink
  // It seems that DetSetVector<ME0DigiSimLink> has duplicate ME0DigiSimLink
  //////////////////////////////////////////////////////////////////////////////

  map<int, edm::DetSet<ME0DigiSimLink>::const_iterator> link_map;
  for (auto link_set = link_set_vector->begin();
            link_set != link_set_vector->end();
            link_set++) {
    ME0DetId me0_id{link_set->detId()};

    // edm::DetSet<ME0DigiSimLink>::const_iterator
    for (auto link = link_set->data.begin();
              link != link_set->data.end();
              link++) {
      int strip = ceil(link->getStrip());
      int unique_id = getUniqueId(me0_id, strip);
      link_map[unique_id] = link;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // NOTE ME0DigiCollection
  // Fill a tree for each chamber with digis
  //////////////////////////////////////////////////////////////////////////////
  for (const auto & chamber : me0->chambers()) {
    resetBranch();

    set<unsigned int> track_id_set;
    const ME0DetId & chamber_id = chamber->id();
    
    // NOTE Fill digi and muon_digi
    for (const auto & layer : chamber->layers()) {
      for (const auto & eta_partition : layer->etaPartitions()) {
        const ME0DetId & id = eta_partition->id();

        b_region_ = id.region();
        b_chamber_ = id.chamber();
        int layer_id = id.layer();
        int roll_id = id.roll();

        ME0DigiCollection::Range range = me0_digi_collection->get(id);

        for (auto digi = range.first; digi != range.second; ++digi) {
          int strip = ceil(digi->strip());
          int index = getIndex(layer_id, roll_id, strip);
          b_digi_[index] = true;
          b_digi_layer_.push_back(layer_id);
          b_digi_ieta_.push_back(roll_id);
          b_digi_strip_.push_back(strip);
          b_num_digi_++;

          int unique_id = getUniqueId(b_region_, b_chamber_, layer_id, roll_id,
                                      strip);

          if (link_map.find(unique_id) != link_map.end()) {
            auto link = link_map[unique_id];

            b_digi_particle_type_.push_back(link->getParticleType());
            b_digi_track_id_.emplace_back(link->getTrackId());

            unsigned int track_id = link->getTrackId();
            if (me0_muon_db.find(track_id) != me0_muon_db.end()) {
              track_id_set.insert(track_id); 

              b_muon_digi_[index] = true;
              b_muon_digi_layer_.push_back(layer_id);
              b_muon_digi_ieta_.push_back(roll_id);
              b_muon_digi_strip_.push_back(strip);
              b_num_muon_digi_++;

              me0_muon_db[track_id].digi_idx.push_back(b_digi_layer_.size() - 1);

              b_digi_is_muon_.push_back(1);
            } else {
              b_digi_is_muon_.push_back(0);
           }

          } else {
            // intrinsic noise or simulated bkg contribution
            b_digi_particle_type_.push_back(0);
            b_digi_track_id_.push_back(-1);
            b_digi_is_muon_.push_back(0);
          }

        } // digi
      } // eta partition
    } // layer

    if (b_num_digi_ < 1) continue;

    // NOTE
    int num_muon = track_id_set.size();
    h_num_muon_->Fill(num_muon);
    b_has_muon_ = num_muon > 0;

    auto seg_range = me0_segment_collection->get(chamber_id);
    b_num_ru_ = std::distance(seg_range.first, seg_range.second);

    ////////////////////////////////////////////////////////////////////////////
    // NOTE
    ////////////////////////////////////////////////////////////////////////////
    if (num_muon <= 1) {
      if (b_has_muon_) {
        unsigned int track_id = *(track_id_set.begin());
        const ME0MuonData & muon = me0_muon_db[track_id];
        const math::XYZTLorentzVectorD & momentum = muon.sim_track->momentum();
        b_muon_pt_ = momentum.Pt();
        b_muon_eta_ = fabs(momentum.Eta());
        b_muon_phi_ = momentum.Phi();

        auto associated_segment = me0_segment_collection->end();
        if (muon.is_reconstructed) {
          b_has_ru_ = true;
          b_ru_asso_nhits_ = muon.rec_segment->nRecHits();
          associated_segment = muon.rec_segment;

          for (const auto & rechit : muon.rec_segment->specificRecHits()) {
            const ME0DetId & me0_id = rechit.me0Id();
            auto eta_partition = me0->etaPartition(me0_id);
            int strip = eta_partition->strip(rechit.localPosition());

            b_ru_asso_rechit_layer_.push_back(me0_id.layer());
            b_ru_asso_rechit_ieta_.push_back(me0_id.roll());
            b_ru_asso_rechit_strip_.push_back(strip);
          } // specificRecHits
        } // if muon.is_reconstructed

        for (auto segment = seg_range.first; segment != seg_range.second; segment++) {
          if (segment == associated_segment) continue; 
          b_ru_fake_nhits_.push_back(segment->nRecHits());
        }
        b_num_ru_fake_ = b_ru_fake_nhits_.size();

      } // if has muon

      tree_ch_->Fill();

      //////////////////////////////////////////////////////////////////////////
      // NOTE scaning for windows
      // finding the one with most hits in 3 strip window
      //////////////////////////////////////////////////////////////////////////
      int max_ieta = 0;
      int max_nstrip = 0;
      int max_hits = 0;
      for (int ieta = 1; ieta <= 8; ++ieta) {
        for (int nstrip = 2; nstrip <= 383; ++nstrip) {

          int current_nhits=0;
          for (int win_nstrip = -1; win_nstrip < 2; ++win_nstrip) {
            for (int nlayer = 1; nlayer <= 6; ++nlayer) {            
              int index = getIndex(nlayer, ieta, nstrip+win_nstrip);
              if (b_digi_[index]) current_nhits++;
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
              int index_win = getIndexWindow(win_nlayer, win_ieta, win_nstrip);
              bool has_hit = false;
              bool has_muon_hit = false;
              // for padding
              if ((test_ieta > 0 and test_ieta < 9) and
                  (test_nstrip > 0 and test_nstrip < 385) ) {
                int index = getIndex(win_nlayer, test_ieta, test_nstrip);
                has_hit = b_digi_[index];
                has_muon_hit = b_muon_digi_[index];
              }
              b_win_digi_[index_win] = has_hit;
              b_win_muon_digi_[index_win] = has_muon_hit;
            }
          }
        }

        b_win_strip_ = max_nstrip;
        b_win_ieta_ = max_ieta;
        tree_win_->Fill();
      } // if (max_hits > 0)

      h_stats_->Fill(0);

    
    } else if (num_muon <= kMaxNumMuons_) {
      bool too_many_digi_found = false;

      vector<ME0MuonData> muons;
      for (unsigned int track_id : track_id_set) {
        auto mu = me0_muon_db[track_id];
        if (mu.digi_idx.size() > kMaxNumDigisPerMuon_) {
          too_many_digi_found = true;
          break;
        }
        muons.push_back(mu);
      }

      if (too_many_digi_found) {
        h_stats_->Fill(1);
        continue;
      }

      b_multi_num_muon_ = num_muon;

      // in decreasing pT order
      sort(muons.begin(), muons.end(), [](ME0MuonData mu0, ME0MuonData mu1) {
        return mu0.sim_track->momentum().Pt() > mu1.sim_track->momentum().Pt();
      });

      std::map<ME0SegmentCollection::const_iterator, int> seg_to_muon_idx;
      for (auto seg = seg_range.first; seg != seg_range.second; seg++) {
        seg_to_muon_idx.emplace(seg, -1);
      }

      for (unsigned int muon_idx = 0; muon_idx < muons.size(); muon_idx++) {
        auto mu = muons[muon_idx];
        auto momentum = mu.sim_track->momentum();
        b_multi_muon_pt_[muon_idx] = momentum.Pt();
        b_multi_muon_eta_[muon_idx] = momentum.Eta();
        b_multi_muon_phi_[muon_idx] = momentum.Phi();
        b_multi_muon_num_digi_[muon_idx] = mu.digi_idx.size();

        for (unsigned int idx = 0; idx < mu.digi_idx.size(); idx++) {
          b_multi_muon_digi_idx_[muon_idx][idx] = mu.digi_idx[idx];
        }

        if (mu.is_reconstructed) {
          seg_to_muon_idx[mu.rec_segment] = muon_idx;
        }
      } // muons

      for (auto [segment, muon_idx] : seg_to_muon_idx) {
        int nhits = segment->nRecHits();
        if (muon_idx >= 0) {
          b_multi_ru_asso_nhits_.push_back(nhits);
          b_multi_ru_asso_muon_idx_.push_back(muon_idx);
        } else {
          b_multi_ru_fake_nhits_.push_back(nhits);
        }
      }

      b_multi_num_ru_asso_ = b_multi_ru_asso_nhits_.size();
      b_multi_num_ru_fake_ = b_multi_ru_fake_nhits_.size();

      h_stats_->Fill(0);
      tree_multi_->Fill();

    } else {
      // (num_muon <= kMaxNumMuons_)
      h_stats_->Fill(2);
    }

  } // chamber loop
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonME0DigisAnalyser);

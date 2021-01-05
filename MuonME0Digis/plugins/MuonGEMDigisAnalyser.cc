#include "MuonTriggering/MuonME0Digis/plugins/MuonGEMDigisAnalyser.h"
#include "MuonTriggering/MuonME0Digis/plugins/GEMMuonData.cc"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;


MuonGEMDigisAnalyser::MuonGEMDigisAnalyser(const edm::ParameterSet& pset) { 
  cout << "ctor begin" << endl;

  auto sim_track_tag = pset.getParameter<edm::InputTag>("simTrackTag");
  sim_track_token_ = consumes<edm::SimTrackContainer>(sim_track_tag);

  auto sim_hit_tag = pset.getParameter<edm::InputTag>("simHitTag");
  sim_hit_token_ = consumes<edm::PSimHitContainer>(sim_hit_tag);

  auto digi_tag = pset.getParameter<edm::InputTag>("gemDigiTag");
  gem_digi_token_ = consumes<GEMDigiCollection>(digi_tag);

  auto link_tag = pset.getParameter<edm::InputTag>("gemDigiSimLinkTag");
  gem_link_token_ = consumes<edm::DetSetVector<GEMDigiSimLink>>(link_tag);

  auto segment_tag = pset.getParameter<edm::InputTag>("gemSegmentTag");
  gem_segment_token_ = consumes<GEMSegmentCollection>(segment_tag);

  // auto gen_tag = pset.getParameter<edm::InputTag>("genParticleTag");
  // gen_particle_token_ = consumes<reco::GenParticleCollection>(gen_tag);

  min_pt_ = pset.getParameter<double>("min_pt");
  min_quality_ = pset.getParameter<double>("min_quality");
  min_num_layers_ = pset.getParameter<unsigned int>("min_num_layers");

  setBranch();
  bookHistogram();
  
  cout << "ctor end" << endl;
}


MuonGEMDigisAnalyser::~MuonGEMDigisAnalyser() {
  cout << "dtor begin" << endl;
  cout << "dtor end" << endl;
}


void MuonGEMDigisAnalyser::setBranch() {
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

  // Reconstructed GEM Segment
  tree_ch_->Branch("has_ru", &b_has_ru_, "has_ru/O");
  tree_ch_->Branch("num_ru", &b_num_ru_, "num_ru/I");
  // 
  tree_ch_->Branch("has_ru_asso", &b_has_ru_asso_, "has_ru_asso/O");
  tree_ch_->Branch("ru_asso_nhits", &b_ru_asso_nhits_, "ru_asso_nhits/I");
  tree_ch_->Branch("ru_asso_rechit_layer", "vector<int>", &b_ru_asso_rechit_layer_);
  tree_ch_->Branch("ru_asso_rechit_ieta", "vector<int>", &b_ru_asso_rechit_ieta_);
  tree_ch_->Branch("ru_asso_rechit_strip", "vector<int>", &b_ru_asso_rechit_strip_);

  tree_ch_->Branch("ru_asso_chi2", &b_ru_asso_chi2_, "ru_asso_chi2/F");
  tree_ch_->Branch("ru_asso_reduced_chi2", &b_ru_asso_reduced_chi2_, "ru_asso_reduced_chi2/F");

  // fake
  tree_ch_->Branch("num_ru_fake", &b_num_ru_fake_, "num_ru_fake/I");
  tree_ch_->Branch("ru_fake_nhits", "vector<int>", &b_ru_fake_nhits_);

  tree_ch_->Branch("ru_fake_chi2", "vector<float>", &b_ru_fake_chi2_);
  tree_ch_->Branch("ru_fake_reduced_chi2", "vector<float>", &b_ru_fake_reduced_chi2_);

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
  tree_multi_->Branch("digi_muon_idx", "vector<int>", &b_multi_digi_muon_idx_);
  tree_multi_->Branch("muon_pt", &b_multi_muon_pt_, "muon_pt[5]/F");
  tree_multi_->Branch("muon_eta", &b_multi_muon_eta_, "muon_eta[5]/F");
  tree_multi_->Branch("muon_phi", &b_multi_muon_phi_, "muon_phi[5]/F");
  tree_multi_->Branch("muon_num_digi", &b_multi_muon_num_digi_, "muon_num_digi[5]/I");
  tree_multi_->Branch("muon_digi_idx", &b_multi_muon_digi_idx_, "muon_digi_idx[5][20]/I");

  tree_multi_->Branch("digi_label", &b_multi_digi_label_, "digi_label[18432]/I");

  tree_multi_->Branch("num_ru_asso", &b_multi_num_ru_asso_, "num_ru_asso/I");
  tree_multi_->Branch("ru_asso_nhits", "vector<int>", &b_multi_ru_asso_nhits_);
  tree_multi_->Branch("ru_asso_muon_idx", "vector<int>",
                      &b_multi_ru_asso_muon_idx_);

  tree_multi_->Branch("num_ru_fake", &b_multi_num_ru_fake_, "num_ru_fake/I");
  tree_multi_->Branch("ru_fake_nhits", "vector<int>", &b_multi_ru_fake_nhits_);

  cout << "setBranch end" << endl;
}

void MuonGEMDigisAnalyser::bookHistogram() {
  //////////////////////////////////////////////////////////////////////////////
  // NOTE Histograms
  //////////////////////////////////////////////////////////////////////////////
  h_sim_seg_nhits_ = file_service_->make<TH1F>("h_sim_seg_nhits", "", 11, 2.5, 13.5);
  // FIXME
  h_sim_seg_nhits_->GetXaxis()->SetBinLabel(1, "#leq 3");
  for (int bin = 2; bin <= 10; bin++) {
    h_sim_seg_nhits_->GetXaxis()->SetBinLabel(bin, Form("%d", bin + 2));
  }
  h_sim_seg_nhits_->GetXaxis()->SetBinLabel(11, "#geq 13");

  h_stats_ = file_service_->make<TH1F>("h_stats", "", 3, -0.5, 2.5);
  h_stats_->GetXaxis()->SetBinLabel(1, "Good");
  h_stats_->GetXaxis()->SetBinLabel(2, "N_{#mu} > 5");
  h_stats_->GetXaxis()->SetBinLabel(3, "N_{DIGI} > 20");

  h_sim_track_pt_ = file_service_->make<TH1F>("h_sim_track_pt", ";p_{T};", 44, 0.0, 220.0);
  h_sim_track_eta_ = file_service_->make<TH1F>("h_sim_track_eta", ";|#eta|;", 8, 2.0, 2.8);
  // GEM Chamber1 Phi -10 rad ~ + 10 rad
  h_sim_track_phi_ = file_service_->make<TH1F>("h_sim_track_phi", ";#phi;", 72, -M_PI, M_PI);

}


void MuonGEMDigisAnalyser::resetBranch() {
  // NOTE tree_ch_
  b_num_digi_ = 0;
  fill_n(b_digi_, 18432, false);
  b_digi_layer_.clear();
  b_digi_ieta_.clear();
  b_digi_strip_.clear();
  b_digi_is_muon_.clear();
  b_digi_particle_type_.clear();
  b_digi_track_id_.clear();

  b_has_muon_ = false;
  b_num_muon_digi_ = 0;
  fill_n(b_muon_digi_, 18432, false);
  b_muon_digi_layer_.clear();
  b_muon_digi_ieta_.clear();
  b_muon_digi_strip_.clear();
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
  b_ru_asso_chi2_ = -1.0f;
  b_ru_asso_reduced_chi2_ = -1.0f;

  b_num_ru_fake_ = 0;
  b_ru_fake_nhits_.clear();

  b_ru_fake_chi2_.clear();
  b_ru_fake_reduced_chi2_.clear();

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
  b_multi_digi_muon_idx_.clear();

  fill_n(b_multi_muon_pt_, 5, -100.0f);
  fill_n(b_multi_muon_eta_, 5, -100.0f);
  fill_n(b_multi_muon_phi_, 5, -100.0f);
  fill_n(b_multi_muon_num_digi_, 5, -1);

  fill_n(*b_multi_muon_digi_idx_, 100, -1);

  fill_n(b_multi_digi_label_, 18432, 0);

  b_multi_num_ru_asso_ = 0;
  b_multi_ru_asso_nhits_.clear();
  b_multi_ru_asso_muon_idx_.clear();

  b_multi_num_ru_fake_ = 0;
  b_multi_ru_fake_nhits_.clear();

}


int MuonGEMDigisAnalyser::get3DImageIndex(int layer, int roll, int strip) {
  // return (8 * 384) * (layer - 1) + 384 * (roll - 1) + (strip - 1);
  // L: layer, R: roll, S: strip
  // index = (8 * 384)*(L - 1) + 384*(R - 1) + (S - 1)
  //       = 3072*L + 384*R + S - (3072 + 384 +1)
  //       = 3072*L + 384*R + S - 3457
  return 3072 * layer + 384 * roll + strip - 3457;
}


int MuonGEMDigisAnalyser::get3DImageIndexWindow(int layer, int roll, int strip) {
  return 30 * layer + 10 * roll + strip - 41;
}


bool MuonGEMDigisAnalyser::isSimTrackGood(
    edm::SimTrackContainer::const_iterator sim_track) {
  // TODO
  // if ((*t).noVertex() && !isMuonGun_)
  //   return false;
  // if ((*t).noGenpart() && !isMuonGun_)
  //   return false;

  if (abs(sim_track->type()) != kMuonPDGId_) return false;
  if (sim_track->momentum().pt() < min_pt_) return false;

  return true;
}


bool MuonGEMDigisAnalyser::isSimHitGood(
    edm::PSimHitContainer::const_iterator sim_hit) {

  if (abs(sim_hit->particleType()) != kMuonPDGId_) return false;

  const EncodedEventId & event_id = sim_hit->eventId();
  if (event_id.event() != 0) return false;
  if (event_id.bunchCrossing() != 0) return false;

  if (sim_hit->processType() != 0) return false;

  return true;
}


std::tuple<bool, uint32_t, bool> MuonGEMDigisAnalyser::isSimSegmentGood(
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
    GEMDetId gem_id{hit->detUnitId()};
    layers_per_chamber[gem_id.chamberId().rawId()].insert(gem_id.layer());
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


map<pair<unsigned int, uint32_t>, GEMMuonData>
MuonGEMDigisAnalyser::buildDB(
    const edm::Handle<edm::SimTrackContainer> & sim_track_container,
    const edm::Handle<edm::PSimHitContainer> & sim_hit_container,
    const edm::Handle<GEMDigiCollection> & gem_digi_collection,
    const edm::Handle<edm::DetSetVector<GEMDigiSimLink> >& link_set_vector,
    const edm::Handle<GEMSegmentCollection> & gem_segment_collection,
    const GEMGeometry* gem) {

  map<pair<unsigned int, uint32_t>, GEMMuonData> gem_muon_db;

  // NOTE
  for (auto sim_track = sim_track_container->begin(); sim_track != sim_track_container->end(); sim_track++) {
    if (not isSimTrackGood(sim_track)) continue;

    vector<edm::PSimHitContainer::const_iterator> sim_segment;

    for (auto sim_hit = sim_hit_container->begin(); sim_hit != sim_hit_container->end(); sim_hit++) {
      if (sim_track->trackId() != sim_hit->trackId()) continue;
      if (not isSimHitGood(sim_hit)) continue;

      GEMDetId gem_id{sim_hit->detUnitId()};
      auto eta_partition = gem->etaPartition(gem_id);
      int sim_hit_strip = ceil(eta_partition->strip(sim_hit->localPosition()));


      // checi if sim_hit is digitized
      // Step1: check if there is a strip that matches a given simhit.
      /*

      bool is_digitized = false;
      auto digi_range = gem_digi_collection->get(gem_id);
      for (auto digi = digi_range.first; digi != digi_range.second; ++digi) {
        if (sim_hit_strip == digi->strip()) {
          is_digitized = true;
          break;
        }
      }

      if (not is_digitized) continue;
      */

      // Step2: Although there is a strip that geometrically matches the sim_hit,
      // it may be due to the backgroud nise. So we check link.
      bool has_link = false;
      auto link_set = link_set_vector->find(gem_id);
      for (auto link = link_set->begin(); link != link_set->end(); link++) {
        // TODO areSimHitAndLinkMatched

        if (sim_hit_strip != static_cast<int>(link->getStrip())) continue;
        if (sim_hit->particleType() != link->getParticleType()) continue;

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
                  return GEMDetId(sim_hit->detUnitId()).chamberId().rawId() != primary_chamber;}),

          sim_segment.end());
    }

    auto chamber_id = GEMDetId(sim_segment[0]->detUnitId()).chamberId();
    ////////////////////////////////////////////////////////////////////////////
    // NOTE find reconstructed segment
    ////////////////////////////////////////////////////////////////////////////

    bool is_reconstructed = false;
    auto matched_rec_segment = gem_segment_collection->end();

    int num_matched = 0;
    set<int> found_layer;

    auto rec_segment_range = gem_segment_collection->get(chamber_id);
    for (auto rec_segment = rec_segment_range.first; rec_segment != rec_segment_range.second; rec_segment++) {

      for (const auto & sim_hit : sim_segment) {
        GEMDetId gem_id{sim_hit->detUnitId()};
        auto eta_partition = gem->etaPartition(gem_id);
        int strip = ceil(eta_partition->strip(sim_hit->localPosition()));

        for (const auto &  rechit : rec_segment->specificRecHits()) {
          // hit-wise matching condition
          if (gem_id != rechit.gemId()) continue;
          int rechit_strip = ceil(eta_partition->strip(rechit.localPosition()));
          // matched GEMRecHit found
          if (strip == rechit_strip) {
            num_matched++;
            found_layer.insert(gem_id.layer());
            // exit specificRecHits loop
            break;
          }
        } // rechit
      } // sim_hit

      // quality of RecSegment
      double quality = static_cast<double>(num_matched) / sim_segment.size();
      if ((quality > min_quality_) and (found_layer.size() >= min_num_layers_)) {
        is_reconstructed = true;
        matched_rec_segment = rec_segment;
        // exit GEMSegmentCollection loop
        break;
      }
    } // GEMSegmentCollection

    unsigned int track_id = sim_track->trackId();

    GEMMuonData data{sim_track, sim_segment,
                           is_reconstructed, matched_rec_segment};

    gem_muon_db.insert({{track_id, chamber_id.rawId()}, data});


  } // SimTrackContainer

  return gem_muon_db;
}


void MuonGEMDigisAnalyser::analyze(const edm::Event& event,
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

  edm::Handle<GEMDigiCollection> gem_digi_collection;
  event.getByToken(gem_digi_token_, gem_digi_collection);
  if (not gem_digi_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMDigiCollection" << endl;
    return;
  }

  edm::Handle<edm::DetSetVector<GEMDigiSimLink>> link_set_vector;
  event.getByToken(gem_link_token_, link_set_vector);
  if (not link_set_vector.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMDigiSimLink" << endl;
    return;
  }

  edm::Handle<GEMSegmentCollection> gem_segment_collection;
  event.getByToken(gem_segment_token_, gem_segment_collection);
  if (not gem_segment_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMSegmentCollection" << endl;
    return;
  }

  edm::ESHandle<GEMGeometry> gem_handle;
  event_setup.get<MuonGeometryRecord>().get(gem_handle);
  if (not gem_handle.isValid()) {
    edm::LogError(kLogCategory_) << "invalid GEMGeometry" << endl;
  }
  const GEMGeometry* gem = &*gem_handle;

  auto gem_muon_db = buildDB(sim_track_container, sim_hit_container,
                            gem_digi_collection, link_set_vector,
                            gem_segment_collection, gem);


  for (auto [key, muon] : gem_muon_db) {
    h_sim_track_pt_->Fill(muon.Pt());
    h_sim_track_eta_->Fill(muon.Eta());
    h_sim_track_phi_->Fill(muon.Phi());
  }


  //////////////////////////////////////////////////////////////////////////////
  // Fill a tree for each chamber with digis
  //////////////////////////////////////////////////////////////////////////////
  for (const auto & chamber : gem->superChambers()) {
    resetBranch();

    const GEMDetId & chamber_id = chamber->id();
    b_region_ = chamber_id.region();
    b_chamber_ = chamber_id.chamber();

    set<pair<unsigned int, uint32_t> > db_key_set;

    // FIXME
    // for (const auto & eta_partition : chamber->etaPartitions()) {
    // chamber->etaPartitions() do not return anything

    for (const auto & gem_layer : chamber->chambers()) {
      for (const auto & eta_partition : gem_layer->etaPartitions()) {
        const GEMDetId & gem_id = eta_partition->id();
          
        int layer = gem_id.layer();
        int roll = gem_id.roll();

        auto digi_range = gem_digi_collection->get(gem_id);
        for (auto digi = digi_range.first; digi != digi_range.second; ++digi) {
          int strip = digi->strip();
          int index = get3DImageIndex(layer, roll, strip);
          b_digi_[index] = true;
          b_digi_layer_.push_back(layer);
          b_digi_ieta_.push_back(roll);
          b_digi_strip_.push_back(strip);
          b_num_digi_++;

          // TODO findLink(digi, link)
          bool has_link = false;
          auto link_set = link_set_vector->find(gem_id);
          edm::DetSet<GEMDigiSimLink>::const_iterator link;
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

            if (gem_muon_db.find(db_key) != gem_muon_db.end()) {
              db_key_set.insert(db_key); 

              b_muon_digi_[index] = true;
              b_muon_digi_layer_.push_back(layer);
              b_muon_digi_ieta_.push_back(roll);
              b_muon_digi_strip_.push_back(strip);
              b_num_muon_digi_++;

              gem_muon_db.at(db_key).appendDigiIndex(b_digi_layer_.size() - 1);

              b_digi_is_muon_.push_back(1);

            } else {
              // not good segment
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

    // NOTE
    if (b_num_digi_ < 1) continue;

    // NOTE
    int num_muon = db_key_set.size();
    b_has_muon_ = num_muon > 0;

    auto rec_seg_range = gem_segment_collection->get(chamber_id);
    b_num_ru_ = std::distance(rec_seg_range.first, rec_seg_range.second);
    b_has_ru_ = b_num_ru_ >= 1;

    ////////////////////////////////////////////////////////////////////////////
    // NOTE
    ////////////////////////////////////////////////////////////////////////////
    if (num_muon <= 1) {
      if (b_has_muon_) {
        const GEMMuonData & muon = gem_muon_db.at(*(db_key_set.begin()));
        b_muon_pt_ = muon.Pt();
        b_muon_eta_ = muon.Eta();
        b_muon_phi_ = muon.Phi();

        auto asso_rec_seg = gem_segment_collection->end();
        if (muon.is_reconstructed()) {
          b_has_ru_asso_ = true;
          asso_rec_seg = muon.rec_segment();
          b_ru_asso_nhits_ = asso_rec_seg->nRecHits();

          b_ru_asso_chi2_ = static_cast<float>(asso_rec_seg->chi2());
          b_ru_asso_reduced_chi2_ = b_ru_asso_chi2_ / asso_rec_seg->degreesOfFreedom();

          for (const auto & rechit : asso_rec_seg->specificRecHits()) {
            const GEMDetId & gem_id = rechit.gemId();
            auto eta_partition = gem->etaPartition(gem_id);

            // FIXME
            int strip = ceil(eta_partition->strip(rechit.localPosition()));

            b_ru_asso_rechit_layer_.push_back(gem_id.layer());
            b_ru_asso_rechit_ieta_.push_back(gem_id.roll());
            b_ru_asso_rechit_strip_.push_back(strip);
          } // specificRecHits
        } // if muon.is_reconstructed

        for (auto rec_segment = rec_seg_range.first; rec_segment != rec_seg_range.second; rec_segment++) {
          if (rec_segment == asso_rec_seg) continue; 
          b_ru_fake_nhits_.push_back(rec_segment->nRecHits());

          float chi2 = static_cast<float>(rec_segment->chi2());
          float reduced_chi2 = chi2 / rec_segment->degreesOfFreedom();

          b_ru_fake_chi2_.push_back(chi2); 
          b_ru_fake_reduced_chi2_.push_back(reduced_chi2);
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
              int index = get3DImageIndex(nlayer, ieta, nstrip+win_nstrip);
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
              int index_win = get3DImageIndexWindow(win_nlayer, win_ieta, win_nstrip);
              bool has_hit = false;
              bool has_muon_hit = false;
              // for padding
              if ((test_ieta > 0 and test_ieta < 9) and
                  (test_nstrip > 0 and test_nstrip < 385) ) {
                int index = get3DImageIndex(win_nlayer, test_ieta, test_nstrip);
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

      vector<GEMMuonData> muons;
      for (const auto & key : db_key_set) {
        auto muon = gem_muon_db.at(key);
        if (muon.digi_indices().size() > kMaxNumDigisPerMuon_) {
          too_many_digi_found = true;
          break;
        }
        muons.push_back(muon);
      }

      if (too_many_digi_found) {
        h_stats_->Fill(1);
        continue;
      }

      b_multi_num_muon_ = num_muon;

      // in decreasing pT order
      sort(muons.begin(), muons.end(), [](GEMMuonData mu0, GEMMuonData mu1) {
        return mu0.Pt() > mu1.Pt();
      });

      std::map<GEMSegmentCollection::const_iterator, int> rec_seg_to_muon_idx;
      for (auto rec_seg = rec_seg_range.first; rec_seg != rec_seg_range.second; rec_seg++) {
        rec_seg_to_muon_idx.emplace(rec_seg, -1);
      }

      b_multi_digi_muon_idx_.resize(b_digi_layer_.size(), -1);
      for (unsigned int muon_idx = 0; muon_idx < muons.size(); muon_idx++) {
        auto mu = muons[muon_idx];
        b_multi_muon_pt_[muon_idx] = mu.Pt();
        b_multi_muon_eta_[muon_idx] = mu.Eta();
        b_multi_muon_phi_[muon_idx] = mu.Phi();
        b_multi_muon_num_digi_[muon_idx] = mu.digi_indices().size();

        for (unsigned int idx = 0; idx < mu.digi_indices().size(); idx++) {
          unsigned int digi_idx = mu.getDigiIndex(idx);
          b_multi_muon_digi_idx_[muon_idx][idx] = digi_idx;
          b_multi_digi_muon_idx_[digi_idx] = muon_idx;

          int img_idx = get3DImageIndex(b_digi_layer_[digi_idx],
                                        b_digi_ieta_[digi_idx],
                                        b_digi_strip_[digi_idx]);

          b_multi_digi_label_[img_idx] = (muon_idx + 1);


        }

        if (mu.is_reconstructed()) {
          rec_seg_to_muon_idx[mu.rec_segment()] = muon_idx;
        }
      } // muons

      // NOTE
      for (auto [rec_segment, muon_idx] : rec_seg_to_muon_idx) {
        int nhits = rec_segment->nRecHits();
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
DEFINE_FWK_MODULE(MuonGEMDigisAnalyser);

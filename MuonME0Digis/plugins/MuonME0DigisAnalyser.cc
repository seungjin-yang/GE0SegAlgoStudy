#include "MuonTriggering/MuonME0Digis/plugins/MuonME0DigisAnalyser.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TTree.h"
using namespace std;

MuonME0DigisAnalyser::MuonME0DigisAnalyser(const edm::ParameterSet& pset) { 
  cout << "ctor begin" << endl;

  auto digi_tag = pset.getParameter<edm::InputTag>("me0DigiToken");
  me0_digi_token_ = consumes<ME0DigiCollection>(digi_tag);

  auto link_tag = pset.getParameter<edm::InputTag>("me0DigiSimLinkToken");
  me0_digi_sim_link_token_ = consumes<edm::DetSetVector<ME0DigiSimLink>>(link_tag);

  auto sim_track_tag = pset.getParameter<edm::InputTag>("simTrackCollection");
  sim_track_token_ = consumes<edm::SimTrackContainer>(sim_track_tag);

  tree_ch_ = file_service_->make<TTree>("chamber", "chamber");
 
  tree_ch_->Branch("digi", &b_digi_, "digi[18432]/O");
  tree_ch_->Branch("muon_digi", &b_muon_digi_, "muon_digi[18432]/O");
  
  tree_ch_->Branch("num_fired_digi", &b_num_fired_digi_, "num_fired_digi/I");
  tree_ch_->Branch("num_muon_digi", &b_num_muon_digi_, "num_muon_digi/I");

  tree_ch_->Branch("digi_layer", "std::vector<Int_t>", &b_digi_layer_);
  tree_ch_->Branch("digi_ieta", "std::vector<Int_t>", &b_digi_ieta_);
  tree_ch_->Branch("digi_strip", "std::vector<Int_t>", &b_digi_strip_);
  
  tree_ch_->Branch("muon_digi_layer", "std::vector<Int_t>", &b_muon_digi_layer_);
  tree_ch_->Branch("muon_digi_ieta", "std::vector<Int_t>", &b_muon_digi_ieta_);
  tree_ch_->Branch("muon_digi_strip", "std::vector<Int_t>", &b_muon_digi_strip_);

  tree_ch_->Branch("muon_pt", &b_muon_pt_, "muon_pt/F");
  tree_ch_->Branch("muon_eta", &b_muon_eta_, "muon_eta/F");
  tree_ch_->Branch("muon_phi", &b_muon_phi_, "muon_phi/F");

  tree_ch_->Branch("simhit_pt", &b_simhit_pt_, "simhit_pt[6]/F");
  tree_ch_->Branch("simhit_eta", &b_simhit_eta_, "simhit_eta[6]/F");
  tree_ch_->Branch("simhit_phi", &b_simhit_phi_, "simhit_phi[6]/F");

  tree_ch_->Branch("region", &b_region_, "region/I");
  tree_ch_->Branch("chamber", &b_chamber_, "chamber/I");

  tree_win_ = file_service_->make<TTree>("window", "window");
  tree_win_->Branch("digi", &b_win_digi_, "digi[180]/O");
  tree_win_->Branch("muon_digi", &b_win_muon_digi_, "muon_digi[180]/O");
  tree_win_->Branch("strip", &b_win_strip_, "strip/I");
  tree_win_->Branch("ieta", &b_win_ieta_, "ieta/I");
  
  tree_win_->Branch("digi_layer", "std::vector<Int_t>", &b_digi_layer_);
  tree_win_->Branch("digi_ieta", "std::vector<Int_t>", &b_digi_ieta_);
  tree_win_->Branch("digi_strip", "std::vector<Int_t>", &b_digi_strip_);
  
  tree_win_->Branch("muon_digi_layer", "std::vector<Int_t>", &b_muon_digi_layer_);
  tree_win_->Branch("muon_digi_ieta", "std::vector<Int_t>", &b_muon_digi_ieta_);
  tree_win_->Branch("muon_digi_strip", "std::vector<Int_t>", &b_muon_digi_strip_);


  tree_win_->Branch("muon_pt", &b_muon_pt_, "muon_pt/F");
  tree_win_->Branch("muon_eta", &b_muon_eta_, "muon_eta/F");
  tree_win_->Branch("muon_phi", &b_muon_phi_, "muon_phi/F");

  tree_win_->Branch("simhit_pt", &b_simhit_pt_, "simhit_pt[6]/F");
  tree_win_->Branch("simhit_eta", &b_simhit_eta_, "simhit_eta[6]/F");
  tree_win_->Branch("simhit_phi", &b_simhit_phi_, "simhit_phi[6]/F");

  tree_win_->Branch("region", &b_region_, "region/I");
  tree_win_->Branch("chamber", &b_chamber_, "chamber/I");
  
  cout << "ctor end" << endl;
}

MuonME0DigisAnalyser::~MuonME0DigisAnalyser() {
  cout << "dtor begin" << endl;
  cout << "dtor end" << endl;
}

void MuonME0DigisAnalyser::resetBranch() {
  fill_n(b_digi_, 18432, false);
  fill_n(b_muon_digi_, 18432, false);

  b_digi_layer_.clear();
  b_digi_ieta_.clear();
  b_digi_strip_.clear();

  b_muon_digi_layer_.clear();
  b_muon_digi_ieta_.clear();
  b_muon_digi_strip_.clear();

  b_num_fired_digi_ = 0;
  b_num_muon_digi_ = 0;

  b_muon_pt_ = -100;
  b_muon_eta_ = -100;
  b_muon_phi_ = -100;

  fill_n(b_simhit_pt_, 6, -100);
  fill_n(b_simhit_eta_, 6, -100);
  fill_n(b_simhit_phi_, 6, -100);

  b_region_ = -100;
  b_chamber_ = -100;

  // FIXME these branches doesn't need to be reset
  // fill_n(b_win_digi_, 180, false);
  // fill_n(b_win_muon_digi_, 180, false);
  // b_win_strip_ = -100;
  // b_win_ieta_ = -100;
}

Int_t MuonME0DigisAnalyser::getUniqueId(Int_t region, Int_t chamber,
                                        Int_t layer, Int_t roll, Int_t strip) {

  Int_t id = 0;
  id += (18 * 6 * 8 * 384) * ((region + 1) / 2);
  id += (6 * 8 * 384) * (chamber - 1);
  id += (8 * 384) * (layer - 1);
  id += 384 * (roll - 1);
  id += (strip - 1);
  return id;
}



void MuonME0DigisAnalyser::analyze(const edm::Event& event,
                                   const edm::EventSetup& event_setup) {

  edm::Handle<ME0DigiCollection> me0_digi_collection;
  event.getByToken(me0_digi_token_, me0_digi_collection);
  if (not me0_digi_collection.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0DigiCollection" << endl;
    return;
  }

  edm::Handle<edm::DetSetVector<ME0DigiSimLink>> link_set_vector;
  event.getByToken(me0_digi_sim_link_token_, link_set_vector);
  if (not link_set_vector.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0DigiSimLink" << endl;
    return;
  }

  edm::Handle<edm::SimTrackContainer> sim_track_container_handle;
  event.getByToken(sim_track_token_, sim_track_container_handle);
  if (not sim_track_container_handle.isValid()) {
    edm::LogError(kLogCategory_) << "invalid SimTrackContainer" << endl;
    return;
  }
  // const edm::SimTrackContainer &
  auto sim_track_container = *sim_track_container_handle.product();

  edm::ESHandle<ME0Geometry> me0;
  event_setup.get<MuonGeometryRecord>().get(me0);
  if (not me0.isValid()) {
    edm::LogError(kLogCategory_) << "invalid ME0Geometry" << endl;
  }

  //////////////////////////////////////////////////////////////////////////////
  // NOTE SimTrack
  //////////////////////////////////////////////////////////////////////////////
  map<UInt_t, SimTrack> sim_track_map;
  for (const auto & sim_track : sim_track_container) {
    if (abs(sim_track.type()) != kMuonPDGId_) continue;
    UInt_t track_id = sim_track.trackId();
    sim_track_map[track_id] = sim_track;
  }

  //////////////////////////////////////////////////////////////////////////////
  // NOTE ME0DigiSimLink
  // It seems that DetSetVector<ME0DigiSimLink> has duplicate ME0DigiSimLink
  // TODO explain the reason. yechan said PSimHit correspond to the Geant4 step.
  //////////////////////////////////////////////////////////////////////////////
  map<Int_t, edm::DetSet<ME0DigiSimLink>::const_iterator> link_map;
  for (auto link_set = link_set_vector->begin(); link_set != link_set_vector->end(); link_set++) {
    ME0DetId me0_det_id{link_set->detId()};
    Int_t region_id = me0_det_id.region();
    Int_t chamber_id = me0_det_id.chamber();
    Int_t layer_id = me0_det_id.layer();
    Int_t roll_id = me0_det_id.roll();

    // edm::DetSet<ME0DigiSimLink>::const_iterator
    for (auto link = link_set->data.begin(); link != link_set->data.end(); ++link) {
      // NOTE https://github.com/cms-sw/cmssw/blob/88fc5373215fdad5d4ba0540fe6f92b718e94afe/Geometry/GEMGeometry/interface/ME0EtaPartition.h#L47-L50
      Int_t strip = ceil(link->getStrip());
      Int_t unique_id = getUniqueId(region_id, chamber_id, layer_id, roll_id, strip);
      link_map[unique_id] = link;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // NOTE ME0DigiCollection
  // Fill a tree for each chamber with digis
  //////////////////////////////////////////////////////////////////////////////
  for (const auto & chamber : me0->chambers()) {
    resetBranch();

    set<Int_t> track_id_set;
    Bool_t has_digi = false;

    for (const auto & layer : chamber->layers()) {
      for (const auto & eta_partition : layer->etaPartitions()) {
        const ME0DetId & id = eta_partition->id();

        b_region_ = id.region();
        b_chamber_ = id.chamber();
        Int_t layer_id = id.layer();
        Int_t roll_id = id.roll();

        ME0DigiCollection::Range range = me0_digi_collection->get(id);
        if ((not has_digi) and (range.first != range.second)) has_digi = true;

        for (auto digi = range.first; digi != range.second; ++digi) {
          Int_t strip = ceil(digi->strip());
          Int_t index = getIndex(layer_id, roll_id, strip);
          b_digi_[index] = true;
          b_digi_layer_.push_back(layer_id);
          b_digi_ieta_.push_back(roll_id);
          b_digi_strip_.push_back(strip);
          b_num_fired_digi_++;

          Int_t unique_id = getUniqueId(b_region_, b_chamber_, layer_id, roll_id, strip);
          auto link = link_map[unique_id];

          if (abs(link->getParticleType()) == kMuonPDGId_) {
            b_muon_digi_[index] = true;
            b_muon_digi_layer_.push_back(layer_id);
            b_muon_digi_ieta_.push_back(roll_id);
            b_muon_digi_strip_.push_back(strip);
            b_num_muon_digi_++;

            // they are same..
            LocalVector simhit_momentum = link->getMomentumAtEntry();
            // GlobalVector simhit_global_momentum = eta_partition->toGlobal(simhit_momentum);

            Int_t layer_index = layer_id - 1;
            b_simhit_pt_[layer_index] = simhit_momentum.perp();
            b_simhit_eta_[layer_index] = simhit_momentum.eta();
            b_simhit_phi_[layer_index] = simhit_momentum.phi();

            Int_t track_id = link->getTrackId();
            track_id_set.insert(track_id);            
          }

        } // digi
      } // eta partition
    } // layer

    // FIXME
    // consider only chamber having one SimTrack.
    if (has_digi and (track_id_set.size() == 1)) {
      Int_t track_id = *(track_id_set.begin());
      if (sim_track_map.find(track_id) != sim_track_map.end()) {
        auto sim_track = sim_track_map[track_id];
        const math::XYZTLorentzVectorD& momentum = sim_track.momentum();

        // https://root.cern.ch/root/html/ROOT__Math__PxPyPzE4D_double_.html#ROOT__Math__PxPyPzE4D_double_
        b_muon_pt_ = momentum.Pt();
        b_muon_eta_ = momentum.Eta();
        b_muon_phi_ = momentum.Phi();

        tree_ch_->Fill();
      }
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
              Int_t index = getIndex(nlayer, ieta, nstrip+win_nstrip);
              if (b_digi_[index]) current_nhits++;
            }
          }
          if (current_nhits > max_hits){
            max_hits = current_nhits;
            max_nstrip = nstrip;
            max_ieta = ieta;
          }
        }
      }
      if (max_hits) {
        // found max strip window center
        
        // saving window
        for (int win_ieta = 1; win_ieta < 4; ++win_ieta) {
          for (int win_nstrip = 1; win_nstrip < 11; ++win_nstrip) {
            for (int win_nlayer = 1; win_nlayer < 7; ++win_nlayer) {
              int test_ieta = max_ieta + win_ieta -1;
              int test_nstrip = max_nstrip + win_nstrip - 5;
              Int_t index_win = getIndexWindow(win_nlayer, win_ieta, win_nstrip);
              bool has_hit = false;
              bool has_muon_hit = false;
              // for padding
              if ((test_ieta > 0 and test_ieta < 9) and
                  (test_nstrip > 0 and test_nstrip < 385) ) {
                Int_t index = getIndex(win_nlayer, test_ieta, test_nstrip);
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
        // cout << "window ieta " << b_win_ieta_
        //      << " st " << b_win_strip_
        //      << " nhits " << max_hits
        //      << endl;
        
      }
    }
    
  } // chamber

}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonME0DigisAnalyser);

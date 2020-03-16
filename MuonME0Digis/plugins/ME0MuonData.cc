#ifndef MuonTriggering_MuonME0Digis_plugins_ME0MuonData_cc
#define MuonTriggering_MuonME0Digis_plugins_ME0MuonData_cc

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"


class ME0MuonData {
 public:
  typedef edm::SimTrackContainer::const_iterator SimTrackData;
  typedef std::vector<edm::PSimHitContainer::const_iterator> SimSegmentData;
  typedef ME0SegmentCollection::const_iterator RecSegmentData;

  ME0MuonData(SimTrackData sim_track,
                    SimSegmentData sim_segment,
                    bool is_reconstructed,
                    RecSegmentData rec_segment) :
    sim_track_(sim_track),
    sim_segment_(sim_segment),
    is_reconstructed_(is_reconstructed),
    rec_segment_(rec_segment) { }

  ~ME0MuonData() {}

  Float_t Pt(void) const {return sim_track_->momentum().Pt();}
  Float_t Eta(void) const {return sim_track_->momentum().Eta();}
  Float_t Phi(void) const {return sim_track_->momentum().Phi();}

  SimTrackData sim_track(void) const {return sim_track_;}
  SimSegmentData sim_segment(void) const {return sim_segment_;}
  bool is_reconstructed(void) const {return is_reconstructed_;}
  RecSegmentData rec_segment(void) const {return rec_segment_;}
  std::vector<unsigned int> digi_indices(void) const {return digi_indices_;}

  void appendDigiIndex(unsigned int idx) {digi_indices_.push_back(idx);}
  unsigned int getDigiIndex(unsigned int idx) {return digi_indices_[idx];}
  unsigned int getNumDigi(void) {return digi_indices_.size();}

 private:
  SimTrackData sim_track_;
  SimSegmentData sim_segment_;
  bool is_reconstructed_;
  RecSegmentData rec_segment_;
  std::vector<unsigned int> digi_indices_;


};


#endif // MuonTriggering_MuonME0Digis_plugins_ME0MuonData_h


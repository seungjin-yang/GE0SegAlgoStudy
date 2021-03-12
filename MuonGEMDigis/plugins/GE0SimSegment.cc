#ifndef MuonTriggering_MuonGEMDigis_plugins_GE0SimSegment_cc
#define MuonTriggering_MuonGEMDigis_plugins_GE0SimSegment_cc

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"

#include <map>
#include <utility>

class GE0SimSegment {
 public:
  using SimTrackData = edm::SimTrackContainer::const_iterator;
  using SimHitData = edm::PSimHitContainer::const_iterator;
  using DigiData = std::pair<GEMDetId, GEMDigi>;

  GE0SimSegment(SimTrackData sim_track,
                GEMDetId det_id,
                std::vector<SimHitData> simhits,
                std::vector<DigiData> digis)
      : sim_track_(sim_track),
        det_id_(det_id),
        simhits_(simhits), // TODO rename
        digis_(digis) {}

  ~GE0SimSegment() {}

  SimTrackData simTrack(void) const {return sim_track_;}
  GEMDetId detId() const {return det_id_;}
  std::vector<SimHitData> simHits(void) const {return simhits_;}
  std::vector<DigiData> digis(void) const {return digis_;}

  // helper
  float pt(void) const {return sim_track_->momentum().Pt();}
  float eta(void) const {return sim_track_->momentum().Eta();}
  float phi(void) const {return sim_track_->momentum().Phi();}
  float charge(void) const {return sim_track_->charge();} // charge is float!!!
  unsigned int trackId() const {return sim_track_->trackId();}
  int type() const {return sim_track_->type();} // pid

 private:
  SimTrackData sim_track_;
  GEMDetId det_id_;
  std::vector<SimHitData> simhits_;
  std::vector<DigiData> digis_;
};

using GE0SimSegmentCollection = std::vector<GE0SimSegment>;
// FIXME is RnageMap overkill?


#endif // MuonTriggering_MuonGEMDigis_plugins_GE0SimSegment_h

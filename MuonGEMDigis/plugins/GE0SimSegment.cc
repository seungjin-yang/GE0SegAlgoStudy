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

  Float_t pt(void) const {return sim_track_->momentum().Pt();}
  Float_t eta(void) const {return sim_track_->momentum().Eta();}
  Float_t phi(void) const {return sim_track_->momentum().Phi();}

  SimTrackData simTrack(void) const {return sim_track_;}
  GEMDetId detId() const {return det_id_;}
  std::vector<SimHitData> simHits(void) const {return simhits_;}
  std::vector<DigiData> digis(void) const {return digis_;}

 private:
  SimTrackData sim_track_;
  GEMDetId det_id_;
  std::vector<SimHitData> simhits_;
  std::vector<DigiData> digis_;
};

using GE0SimSegmentCollection = std::vector<GE0SimSegment>;
// FIXME is RnageMap overkill?


#endif // MuonTriggering_MuonGEMDigis_plugins_GE0SimSegment_h

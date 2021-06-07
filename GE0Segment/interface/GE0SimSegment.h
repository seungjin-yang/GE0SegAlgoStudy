#ifndef MuonTriggering_GE0Segment_plugins_GE0SimSegment_cc
#define MuonTriggering_GE0Segment_plugins_GE0SimSegment_cc

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"

class GE0SimSegment {
 public:
  using SimTrackData = edm::SimTrackContainer::const_iterator;
  using LinkData = edm::DetSet<GEMDigiSimLink>::const_iterator;
  using SimHitData = edm::PSimHitContainer::const_iterator;

  GE0SimSegment(SimTrackData sim_track, GEMDetId det_id,
                std::vector<LinkData> links, std::vector<SimHitData> hits,
                std::map<LinkData, SimHitData> link2hit)
      : sim_track_(sim_track), det_id_(det_id), links_(links), hits_(hits),
        link2hit_(link2hit) {}
  ~GE0SimSegment() {}

  SimTrackData simTrack(void) const {return sim_track_;}
  GEMDetId superChamberId() const {return det_id_;}
  std::vector<LinkData> links() const {return links_;}
  std::vector<SimHitData> hits() const {return hits_;}

  // helper
  float pt(void) const {return sim_track_->momentum().Pt();}
  float eta(void) const {return sim_track_->momentum().Eta();}
  float phi(void) const {return sim_track_->momentum().Phi();}
  float charge(void) const {return sim_track_->charge();} // charge is float!!!
  unsigned int trackId() const {return sim_track_->trackId();}
  int type() const {return sim_track_->type();} // pid
  SimHitData getSimHit(const LinkData& link) const {return link2hit_.at(link);}

 private:
  SimTrackData sim_track_;
  GEMDetId det_id_;
  std::vector<LinkData> links_;
  std::vector<SimHitData> hits_;
  std::map<LinkData, SimHitData> link2hit_;
};

#endif // MuonTriggering_GE0Segment_plugins_GE0SimSegment_h

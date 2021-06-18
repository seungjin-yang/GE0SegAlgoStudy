#ifndef MuonTriggering_GE0Segment_plugins_GE0SimSegment_cc
#define MuonTriggering_GE0Segment_plugins_GE0SimSegment_cc

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/GEMDigiSimLink/interface/GEMDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"

#include <vector>

namespace ge0 {
using SimTrackData = edm::SimTrackContainer::const_iterator;
using LinkData = edm::DetSet<GEMDigiSimLink>::const_iterator;
using SimHitData = edm::PSimHitContainer::const_iterator;


class GE0SimHit {
 public:
  GE0SimHit(SimHitData hit, unsigned int det_unit_id, int strip, int bx,
            std::vector<LinkData> cluster) 
      : hit_(hit), det_unit_id_(det_unit_id), strip_(strip), bx_(bx),
        cluster_(cluster) {}
  ~GE0SimHit() {}
  SimHitData data() const {return hit_;}
  unsigned int getDetUnitId() const {return det_unit_id_;}
  std::vector<LinkData> getCluster() const {return cluster_;}

  int getStrip() const {return strip_;}
  int getBx() const {return bx_;};
  int getClusterSize() const {return static_cast<int>(cluster_.size());}
  LocalPoint getLocalPosition() const {return hit_->localPosition();}

 private:
  SimHitData hit_;
  unsigned int det_unit_id_;
  int strip_;
  int bx_;
  std::vector<LinkData> cluster_;
};


class GE0SimSegment {
 public:


  GE0SimSegment(SimTrackData sim_track, GEMDetId det_id,
                std::vector<GE0SimHit> hits)
      : sim_track_(sim_track), det_id_(det_id), hits_(hits) {}
  ~GE0SimSegment() {}

  SimTrackData simTrack(void) const {return sim_track_;}
  GEMDetId superChamberId() const {return det_id_;}
  std::vector<GE0SimHit> hits() const {return hits_;}

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
  std::vector<GE0SimHit> hits_;
};

};

#endif // MuonTriggering_GE0Segment_plugins_GE0SimSegment_h

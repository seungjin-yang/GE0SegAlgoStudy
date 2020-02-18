#ifndef MuonTriggering_MuonME0Digis_plugins_ME0SegmentData_h
#define MuonTriggering_MuonME0Digis_plugins_ME0SegmentData_h

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/GEMRecHit/interface/ME0SegmentCollection.h"


struct ME0MuonData {
  edm::SimTrackContainer::const_iterator sim_track;
  std::vector<edm::PSimHitContainer::const_iterator> sim_hits;

  bool is_reconstructed;
  ME0SegmentCollection::const_iterator rec_segment;

  std::vector<unsigned int> digi_idx;
};


#endif // MuonTriggering_MuonME0Digis_plugins_ME0SegmentData_h

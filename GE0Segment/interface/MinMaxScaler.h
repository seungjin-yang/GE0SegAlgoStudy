#ifndef MuonTriggering_GE0Segment_MinMaxScaler_h
#define MuonTriggering_GE0Segment_MinMaxScaler_h

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"
#include "Geometry/GEMGeometry/interface/GEMChamber.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"

#include <limits> // numeric_limits
#include <cassert> // assert

class MinMaxScaler {
 public:
  MinMaxScaler(float x_in_min,
               float x_in_max,
               float x_out_min=-1.0f,
               float x_out_max=1.0f)
    : x_in_min_(x_in_min),
      x_in_max_(x_in_max),
      x_in_scale_(x_in_max - x_in_min),
      x_out_min_(x_out_min),
      x_out_max_(x_out_max),
      x_out_scale_(x_out_max - x_out_min) {}
  ~MinMaxScaler() {}

  float transform(const float x) {
    const float x_std = (x - x_in_min_) / x_in_scale_;
    const float x_scaled = x_std * x_out_scale_ + x_out_min_;
    return x_scaled;
  }

 private:
   const float x_in_min_;
   const float x_in_max_;
   const float x_in_scale_;
   const float x_out_min_;
   const float x_out_max_;
   const float x_out_scale_;
};

class GE0Scaler {
 public:
  GE0Scaler(const GEMSuperChamber* superchamber) {
    float x_min = std::numeric_limits<float>::max(); 
    float x_max = std::numeric_limits<float>::min(); 

    float y_min = std::numeric_limits<float>::max(); 
    float y_max = std::numeric_limits<float>::min(); 

    float z_min = std::numeric_limits<float>::max(); 
    float z_max = std::numeric_limits<float>::min(); 

    for (const GEMChamber* chamber : superchamber->chambers()) {
      for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
        for (const int strip : {0, eta_partition->nstrips() - 1}) {

          const LocalPoint&& position = superchamber->toLocal(
              eta_partition->toGlobal(eta_partition->centreOfStrip(strip)));

          x_min = std::min(x_min, position.x());
          x_max = std::max(x_max, position.x());

          y_min = std::min(y_min, position.y());
          y_max = std::max(y_max, position.y());

          z_min = std::min(z_min, position.z());
          z_max = std::max(z_max, position.z());

        }
      }
    }

    assert(x_min < x_max);
    assert(y_min < y_max);
    assert(z_min < z_max);

    x_scaler_ = std::make_unique<MinMaxScaler>(x_min, x_max, -1.0f, 1.0f);
    y_scaler_ = std::make_unique<MinMaxScaler>(y_min, y_max, -1.0f, 1.0f);
    z_scaler_ = std::make_unique<MinMaxScaler>(z_min, z_max, -1.0f, 1.0f);
  }

  ~GE0Scaler() {}

  float transformX(float x) { return x_scaler_->transform(x); }
  float transformY(float y) { return y_scaler_->transform(y); }
  float transformZ(float z) { return z_scaler_->transform(z); }

 private:
  std::unique_ptr<MinMaxScaler> x_scaler_;
  std::unique_ptr<MinMaxScaler> y_scaler_;
  std::unique_ptr<MinMaxScaler> z_scaler_;

};


#endif // MuonTriggering_GE0Segment_MinMaxScaler_h

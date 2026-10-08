#pragma once
// Stage 2: fixed timestep with an accumulator ("Fix Your Timestep", Glenn Fiedler).
//
// Rendering runs at whatever rate the display gives us, but physics must step by the SAME dt every
// time: variable steps make results non-reproducible and can destabilise stiff systems. So we bank
// real elapsed time and spend it in dt-sized chunks.
#include "math.hpp"

namespace phys {

class FixedTimestep {
public:
    // max_frame caps the time banked per frame. Without it a long stall (debugger, window drag)
    // would demand hundreds of catch-up steps, each slower than real time: the "spiral of death".
    explicit FixedTimestep(Real dt, Real max_frame = static_cast<Real>(0.25))
        : dt_(dt), max_frame_(max_frame) {}

    // Banks `frame_time` and calls step_fn(dt) as many times as it now owes. Returns the step count.
    template <class StepFn>
    int advance(Real frame_time, StepFn&& step_fn) {
        acc_ += frame_time < max_frame_ ? frame_time : max_frame_;
        int n = 0;
        while (acc_ >= dt_) {
            step_fn(dt_);
            acc_ -= dt_;
            ++n;
        }
        return n;
    }

    Real dt() const { return dt_; }
    // Leftover time as a fraction of dt in [0,1): used to interpolate rendering between two states.
    Real alpha() const { return acc_ / dt_; }

private:
    Real dt_;
    Real max_frame_;
    Real acc_ = 0;
};

}  // namespace phys

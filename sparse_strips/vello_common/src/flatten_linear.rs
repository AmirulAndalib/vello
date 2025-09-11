use crate::flatten::{TOL, TOL_2};
use crate::flatten_simd::Callback;
use crate::kurbo::{CubicBez, PathEl};

pub fn flatten_cubic(curve: &CubicBez, callback: &mut impl Callback) {
    if cubic_is_a_point(&curve) {
        return;
    }

    let mut rem = *curve;
    let mut from = rem.p0;

    let mut split = 0.5;
    loop {
        // Only check flatness of the entire remaining chunk if
        // we are not in the process of doing very fine subdivision.
        // This gives a 9% improvement.
        if split >= 0.25 && cubic_is_flat(&rem) {
            callback.callback(PathEl::LineTo(rem.p3));
            return;
        }

        loop {
            let sub = before_split(&rem, split);
            if cubic_is_flat(&sub) {
                callback.callback(PathEl::LineTo(sub.p3));
                from = sub.p3;
                rem = after_split(&rem, split);
                let next_split = split * 2.0;
                if next_split < 1.0 {
                    split = next_split;
                }
                break;
            }
            split *= 0.5;
        }
    }
}


/// Returns true if the curve can be approximated with a single line segment, given
/// a tolerance threshold.
// Note: The inline annotation here makes a huge difference for `linear`.
#[inline]
pub fn cubic_is_flat(curve: &CubicBez) -> bool {
    // Similar to Line::square_distance_to_point, except we keep
    // the sign of c1 and c2 to compute tighter upper bounds as we
    // do in fat_line_min_max.
    let baseline = curve.p3 - curve.p0;
    let v1 = curve.p1 - curve.p0;
    let v2 = curve.p2 - curve.p0;
    let v3 = curve.p2 - curve.p3;

    let c1 = baseline.cross(v1);
    let c2 = baseline.cross(v2);
    // TODO: it is faster to multiply the threshold with baseline_len2
    // instead of dividing d1 and d2, as done in lyon, but it changes
    // the behavior when the baseline length is zero in ways that breaks
    // some of the cubic intersection tests.
    let baseline_len2 = baseline.hypot2();
    let d1 = c1 * c1;
    let d2 = c2 * c2;

    let factor = if (c1 * c2) > 0.0 {
        3.0 / 4.0
    } else {
        4.0 / 9.0
    };

    let f2 = factor * factor;
    let threshold = baseline_len2 * TOL_2;

    (d1 * f2 <= threshold)
        & (d2 * f2 <= threshold)
        // TODO: This check is missing from CubicBezierSegment::is_linear, which
        // misses flat-ish curves with control points that aren't between the
        // endpoints.
        && (baseline.dot(v1) > -TOL)
        && (baseline.dot(v3) < TOL)
}

#[inline]
pub(crate) fn cubic_is_a_point(&curve: &CubicBez) -> bool {
    // Use <= so that tolerance can be zero.
    (curve.p0 - curve.p1).hypot2() <= TOL_2
        && (curve.p0 - curve.p1).hypot2() <= TOL_2
        && (curve.p3 - curve.p2).hypot2() <= TOL_2
}

fn before_split(curve: &CubicBez, t: f64) -> CubicBez {
    let p1a = curve.p0 + (curve.p1 - curve.p0) * t;
    let p2a = curve.p1 + (curve.p2 - curve.p1) * t;
    let p1aa = p1a + (p2a - p1a) * t;
    let ctrl3a = curve.p2 + (curve.p3 - curve.p2) * t;
    let p2aa = p2a + (ctrl3a - p2a) * t;
    let p1aaa = p1aa + (p2aa - p1aa) * t;

    CubicBez {
        p0: curve.p0,
        p1: p1a,
        p2: p1aa,
        p3: p1aaa,
    }
}

fn after_split(curve: &CubicBez, t: f64) -> CubicBez {
    let p1a = curve.p0 + (curve.p1 - curve.p0) * t;
    let p2a = curve.p1 + (curve.p2 - curve.p1) * t;
    let p1aa = p1a + (p2a - p1a) * t;
    let ctrl3a = curve.p2 + (curve.p3 - curve.p2) * t;
    let p2aa = p2a + (ctrl3a - p2a) * t;

    CubicBez {
        p0: p1aa + (p2aa - p1aa) * t,
        p1: p2a + (ctrl3a - p2a) * t,
        p2: ctrl3a,
        p3: curve.p3,
    }
}
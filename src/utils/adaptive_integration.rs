use bio::stats::probs::LogProb;

#[derive(Debug, Clone, Copy)]
struct GridPoint<T> 
where
    T: std::cmp::Ord
{
    point: T,
    prob: LogProb,
}

impl Ord for GridPoint {
    fn cmp(&self, other: &Self) -> Ordering {
        self.point.cmp(&other.point)
    }
}

/// Integrate over an interval of type T with a given density function while trying to minimize
/// the number of grid points evaluated and still hit the maximum likelihood point.
/// This is achieved via a binary search over the grid points.
/// The assumption is that the density is unimodal. If that is not the case,
/// the binary search will not find the maximum and the integral can miss features.
pub(crate) fn ln_integrate_exp<T, F>(density: F, min_point: T, max_point: T, max_resolution: T)
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Div<Output = T> + Mul<Output = T> + Float,
    F: FnMut(T) -> LogProb,
    f64: From<T>,
{
    let grid_point = |point| {
        GridPoint {
            point: point,
            prob: density(point),
        }
    };
    let middle_grid_point = |left, right| {
        grid_point((right.point + left.point) / NotNan::new_unchecked(2.0))
    };
    // METHOD:
    // Step 1: perform binary search for maximum likelihood point
    // Remember all points.
    let mut left = grid_point(min_point);
    let mut right = grid_point(max_point);

    let mut grid_points = vec![left, right];

    while (*(right.point - left.point) >= max_resolution) || left < right {
        let middle = middle_grid_point(&left, &right);

        if left.prob > right.prob {
            // investigate left window more closely
            right = middle;
        } else {
            // investiage right window more closely
            left = middle;
        }

        grid_points.push(middle);
    }
    grid_points.sort();

    // METHOD:
    // Step 2: integrate over grid points visited during the binary search.
    LogProb::ln_trapezoidal_integrate_grid_exp(|i, g| grid_points[i].prob, &grid_points)
}
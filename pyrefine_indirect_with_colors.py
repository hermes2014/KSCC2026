"""
pyrefine_indirect.py

A Python implementation of a refineR-like indirect reference interval estimator.

This module is intentionally not a line-by-line port of the R refineR package.
It reimplements the same main ideas in Python:
  * one-parameter Box-Cox model for the non-pathological distribution
  * inverse/original-scale histogram fitting
  * multi-level grid search over lambda, mu, sigma, and NP fraction P
  * asymmetric confidence-band bin selection
  * optional bootstrap confidence intervals
  * optional modified Box-Cox-like shift search

Dependencies:
    numpy, scipy, matplotlib
Optional:
    pandas, joblib

Clinical note:
    This is a research/validation aid. Do not deploy derived reference intervals
    clinically without local review, pre-analytical/analytical validation,
    and comparison with manufacturer/literature/direct reference data when possible.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, Sequence, Union, Dict, Any
import math
import warnings

import numpy as np
from scipy.stats import norm
from scipy.special import gammaln
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

try:
    import pandas as pd  # type: ignore
except Exception:  # pragma: no cover
    pd = None

try:
    from joblib import Parallel, delayed  # type: ignore
except Exception:  # pragma: no cover
    Parallel = None
    delayed = None


ModelName = Literal["BoxCox", "modBoxCoxFast", "modBoxCox"]
PointEstimate = Literal["full", "median_bootstrap"]


_EPS = 1e-12


def boxcox(x: np.ndarray, lam: float) -> np.ndarray:
    """One-parameter Box-Cox transformation."""
    x = np.asarray(x, dtype=float)
    x = np.maximum(x, _EPS)
    if abs(lam) < 1e-10:
        return np.log(x)
    return (np.power(x, lam) - 1.0) / lam


def inv_boxcox(y: np.ndarray, lam: float) -> np.ndarray:
    """Inverse one-parameter Box-Cox transformation."""
    y = np.asarray(y, dtype=float)
    if abs(lam) < 1e-10:
        return np.exp(y)
    inside = lam * y + 1.0
    out = np.full_like(y, np.nan, dtype=float)
    ok = inside > 0
    out[ok] = np.power(inside[ok], 1.0 / lam)
    return out


def _clean_data(data: Sequence[float], allow_zero: bool = True) -> np.ndarray:
    x = np.asarray(data, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x >= 0]
    if not allow_zero:
        x = x[x > 0]
    if x.size <= 10:
        raise ValueError("At least 11 finite, non-negative observations are required.")
    return x


def _estimate_rounding_base(x: np.ndarray) -> Optional[float]:
    """Roughly estimate decimal rounding base; returns None if data appear continuous."""
    xu = np.unique(x[np.isfinite(x)])
    if xu.size < 5:
        return None
    if xu.size > 3000:
        return None
    diffs = np.diff(np.sort(xu))
    diffs = diffs[diffs > _EPS]
    if diffs.size == 0:
        return None
    q10 = np.quantile(diffs, 0.1)
    if q10 <= 0:
        return None
    # A conservative rounded-data flag: many adjacent gaps are near one small base.
    ratio_near = np.mean(np.abs(diffs / q10 - np.round(diffs / q10)) < 1e-3)
    if ratio_near > 0.7:
        return float(q10)
    return None


@dataclass
class HistData:
    counts: np.ndarray
    break_l: np.ndarray
    break_r: np.ndarray
    mids: np.ndarray
    n_data: int
    n_bins: int
    overlap_factor: int
    ab_original: tuple[float, float]
    full_counts: np.ndarray
    full_breaks: np.ndarray


def _make_overlapping_histogram(x: np.ndarray, ab: tuple[float, float], rounding_base: Optional[float] = None) -> HistData:
    """Generate an overlapping histogram similar in spirit to refineR::generateHistData."""
    lo, hi = float(ab[0]), float(ab[1])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = np.quantile(x, [0.005, 0.995])
        if hi <= lo:
            lo, hi = float(np.min(x)), float(np.max(x) + _EPS)

    n_inside = int(np.sum((x >= lo) & (x <= hi)))
    n_bins = int(round(min(40, max(9, 11 * (max(n_inside, 1) / 5000.0) ** 0.2))))
    overlap_factor = max(1, int(round(256 / n_bins)))
    n_bins_total = max(8, n_bins * overlap_factor)

    if rounding_base is not None and rounding_base > 0:
        step = math.ceil((hi - lo) / max(256, 1) / rounding_base) * rounding_base
        step = max(step, rounding_base)
        overlap_factor = max(1, int(round((hi - lo) / step / max(n_bins, 1))))
        n_bins_total = max(8, n_bins * overlap_factor)
        lo = max(0.5 * rounding_base, round(lo / rounding_base) * rounding_base - 0.5 * rounding_base)
        hi = lo + (n_bins_total + overlap_factor - 1) * step
        breaks = lo + np.arange(n_bins_total + overlap_factor) * step
    else:
        breaks = np.linspace(lo, hi, n_bins_total + overlap_factor)

    if breaks.size < 3 or breaks[-1] <= breaks[0]:
        breaks = np.linspace(np.min(x), np.max(x) + _EPS, 256)
        n_bins_total = 256 - overlap_factor

    mask = (x > breaks[0]) & (x <= breaks[-1])
    full_counts, full_breaks = np.histogram(x[mask], bins=breaks)

    counts = []
    break_l = []
    break_r = []
    for i in range(0, n_bins_total):
        counts.append(np.sum(full_counts[i : i + overlap_factor]))
        break_l.append(full_breaks[i])
        break_r.append(full_breaks[i + overlap_factor])

    # Boundary bins catch observations outside the central histogram region.
    if np.min(x) <= lo:
        counts.append(np.sum((x > -_EPS) & (x <= lo)))
        break_l.append(max(float(np.min(x)), _EPS))
        break_r.append(lo)
        full_counts = np.r_[np.sum((x > -_EPS) & (x <= lo)), full_counts]
        full_breaks = np.r_[max(float(np.min(x)), _EPS), full_breaks]
    else:
        counts.append(0)
        break_l.append(lo)
        break_r.append(lo)
        full_counts = np.r_[0, full_counts]
        full_breaks = np.r_[lo, full_breaks]

    if np.max(x) > hi:
        counts.append(np.sum((x > hi) & (x <= np.max(x))))
        break_l.append(hi)
        break_r.append(float(np.max(x)))
        full_counts = np.r_[full_counts, np.sum((x > hi) & (x <= np.max(x)))]
        full_breaks = np.r_[full_breaks, float(np.max(x))]
    else:
        counts.append(0)
        break_l.append(hi)
        break_r.append(hi)
        full_counts = np.r_[full_counts, 0]
        full_breaks = np.r_[full_breaks, hi]

    counts = np.asarray(counts, dtype=float)
    break_l = np.asarray(break_l, dtype=float)
    break_r = np.asarray(break_r, dtype=float)
    mids = 0.5 * (break_l + break_r)
    order = np.argsort(mids)

    return HistData(
        counts=counts[order],
        break_l=break_l[order],
        break_r=break_r[order],
        mids=mids[order],
        n_data=int(x.size),
        n_bins=n_bins,
        overlap_factor=overlap_factor,
        ab_original=(lo, hi),
        full_counts=full_counts.astype(float),
        full_breaks=full_breaks.astype(float),
    )


def _smooth_density_grid(y: np.ndarray, n_grid: int = 512, smooth_sigma: float = 4.0) -> tuple[np.ndarray, np.ndarray]:
    """Fast smoothed histogram density on a grid."""
    y = y[np.isfinite(y)]
    if y.size < 20:
        raise ValueError("Too few finite values for density estimation.")
    lo, hi = np.quantile(y, [0.005, 0.995])
    if hi <= lo:
        lo, hi = np.min(y), np.max(y) + _EPS
    counts, edges = np.histogram(y, bins=n_grid, range=(lo, hi), density=False)
    dx = edges[1] - edges[0]
    dens = gaussian_filter1d(counts.astype(float), sigma=smooth_sigma, mode="nearest")
    area = np.sum(dens) * dx
    if area > 0:
        dens = dens / area
    grid = 0.5 * (edges[:-1] + edges[1:])
    return grid, dens


def _find_main_peak_index(grid: np.ndarray, dens: np.ndarray) -> int:
    peaks, _ = find_peaks(dens)
    if peaks.size == 0:
        return int(np.argmax(dens))
    # Choose peak with largest area between surrounding valleys, with a small height tie-breaker.
    valleys, _ = find_peaks(-dens)
    scores = []
    for p in peaks:
        left_candidates = valleys[valleys < p]
        right_candidates = valleys[valleys > p]
        left = int(left_candidates[-1]) if left_candidates.size else 0
        right = int(right_candidates[0]) if right_candidates.size else len(dens) - 1
        area = np.trapz(dens[left : right + 1], grid[left : right + 1])
        scores.append(area + 0.05 * dens[p] / max(np.max(dens), _EPS))
    return int(peaks[int(np.argmax(scores))])


def _estimate_mu_sigma_region(x_shifted: np.ndarray, lam: float) -> tuple[tuple[float, float], tuple[float, float], float, float]:
    """Estimate starting search ranges for mu and sigma in transformed space."""
    y = boxcox(x_shifted, lam)
    y = y[np.isfinite(y)]
    grid, dens = _smooth_density_grid(y, n_grid=512, smooth_sigma=4.0)
    pidx = _find_main_peak_index(grid, dens)
    y_mode = float(grid[pidx])
    peak = float(dens[pidx])

    centers = []
    sigmas = []
    for rel_height in np.linspace(0.50, 0.95, 10):
        target = peak * rel_height
        left_idx = np.where(dens[:pidx] <= target)[0]
        right_idx = np.where(dens[pidx:] <= target)[0]
        if left_idx.size == 0 or right_idx.size == 0:
            continue
        li = int(left_idx[-1])
        ri = int(pidx + right_idx[0])
        if ri <= li:
            continue
        left = grid[li]
        right = grid[ri]
        z = math.sqrt(max(-2.0 * math.log(rel_height), _EPS))
        sig = (right - left) / (2.0 * z)
        if np.isfinite(sig) and sig > _EPS:
            centers.append(0.5 * (left + right))
            sigmas.append(sig)

    if len(sigmas) == 0:
        mu0 = float(np.median(y))
        sig0 = float(np.std(y, ddof=1))
    else:
        mu0 = float(np.median(centers))
        sig0 = float(np.median(sigmas))

    sig0 = max(sig0, np.std(y, ddof=1) * 0.05, _EPS)
    mu_span = max(2.0 * sig0, 0.1 * np.std(y, ddof=1), _EPS)
    mu_range = (mu0 - mu_span, mu0 + mu_span)
    sigma_range = (max(sig0 * 0.35, _EPS), sig0 * 2.50)
    return mu_range, sigma_range, y_mode, sig0


def _estimate_ab_original(x: np.ndarray) -> tuple[float, float]:
    """Estimate original-scale histogram region around the main peak."""
    q = np.quantile(x, [0.000, 0.005, 0.010, 0.990, 0.995, 0.999])
    lo, hi = float(q[2]), float(q[3])
    if hi <= lo:
        return (max(float(np.min(x)), _EPS), float(np.max(x) + _EPS))

    # Use the log/Box-Cox(0) mode to avoid a long high tail dominating the histogram.
    xp = np.maximum(x, _EPS)
    y = np.log(xp)
    try:
        grid, dens = _smooth_density_grid(y, n_grid=512, smooth_sigma=4.0)
        pidx = _find_main_peak_index(grid, dens)
        peak = dens[pidx]
        half = 0.15 * peak
        left = np.where(dens[:pidx] <= half)[0]
        right = np.where(dens[pidx:] <= half)[0]
        if left.size and right.size:
            ylo = grid[int(left[-1])]
            yhi = grid[int(pidx + right[0])]
            lo = max(float(np.exp(ylo)), float(q[1]), _EPS)
            hi = min(float(np.exp(yhi)), float(q[5]))
            # Expand to avoid overly narrow fitting in clean datasets.
            width = hi - lo
            lo = max(lo - 0.5 * width, float(q[0]), _EPS)
            hi = min(hi + 0.5 * width, float(q[5]))
    except Exception:
        pass
    if hi <= lo:
        lo, hi = float(q[1]), float(q[4])
    return (max(lo, _EPS), max(hi, lo + _EPS))


def _predicted_counts(hist: HistData, lam: float, mu: float, sigma: float) -> np.ndarray:
    bl = boxcox(hist.break_l, lam)
    br = boxcox(hist.break_r, lam)
    finite = np.isfinite(bl) & np.isfinite(br)
    pred = np.zeros_like(hist.counts, dtype=float)
    if not np.any(finite) or sigma <= 0:
        return pred
    mn = np.nanmin(bl[finite])
    mx = np.nanmax(br[finite])
    pcorr_den = norm.cdf((mx - mu) / sigma) - norm.cdf((mn - mu) / sigma)
    if pcorr_den <= _EPS or not np.isfinite(pcorr_den):
        return pred
    pcorr = 1.0 / pcorr_den
    p = norm.cdf((br[finite] - mu) / sigma) - norm.cdf((bl[finite] - mu) / sigma)
    pred[finite] = hist.n_data * pcorr * np.maximum(p, 0)
    pred[pred < 0] = 0
    return pred



def _sum_full_hist_between(hist: HistData, left: float, right: float) -> float:
    """Sum non-overlapping histogram counts whose bin midpoints fall in [left, right)."""
    if right <= left:
        return 1.0
    edges = hist.full_breaks
    counts = hist.full_counts
    if edges.size != counts.size + 1:
        # Defensive fallback for unexpected edge/count shapes.
        mids = np.linspace(edges[0], edges[-1], counts.size)
    else:
        mids = 0.5 * (edges[:-1] + edges[1:])
    sel = (mids >= left) & (mids < right)
    total = float(np.sum(counts[sel]))
    return max(total, 1.0)


def _sum_for_p_area(hist: HistData, pred: np.ndarray, lam: float, mu: float, sigma: float, p_corr: float) -> tuple[np.ndarray, np.ndarray]:
    """Approximate refineR getSumForPArea: observed/predicted counts around the peak."""
    p_limit_min = np.array([0.50, 0.55, 0.60, 0.65, 0.70])
    p_limit_max = np.array([0.90, 0.94, 0.95, 0.96, 0.97])
    n_limits = len(p_limit_min)
    border_l = np.zeros(2 * n_limits + 1, dtype=float)
    border_r = np.zeros(2 * n_limits + 1, dtype=float)
    sum_data = np.ones(2 * n_limits + 1, dtype=float)

    idx_peak = int(np.argmax(pred))
    peak = float(pred[idx_peak])
    if peak <= 0:
        return sum_data, sum_data.copy()
    peak_sel = pred >= 0.95 * peak
    border_l[0] = float(np.min(hist.break_l[peak_sel]))
    border_r[0] = float(np.max(hist.break_r[peak_sel]))
    sum_data[0] = _sum_full_hist_between(hist, border_l[0], border_r[0])

    pred_left = pred[: idx_peak + 1]
    pred_right = pred[idx_peak:]
    bl_left, br_left = hist.break_l[: idx_peak + 1], hist.break_r[: idx_peak + 1]
    bl_right, br_right = hist.break_l[idx_peak:], hist.break_r[idx_peak:]

    for i, (pl, ph) in enumerate(zip(p_limit_min * peak, p_limit_max * peak), start=1):
        sel_l = (pred_left >= pl) & (pred_left < ph)
        if np.any(sel_l):
            l = float(np.min(bl_left[sel_l])); r = float(np.max(br_left[sel_l]))
            border_l[2 * i - 1] = l; border_r[2 * i - 1] = r
            sum_data[2 * i - 1] = _sum_full_hist_between(hist, l, r)
        else:
            border_l[2 * i - 1] = border_l[0]; border_r[2 * i - 1] = border_r[0]

        sel_r = (pred_right >= pl) & (pred_right < ph)
        if np.any(sel_r):
            l = float(np.min(bl_right[sel_r])); r = float(np.max(br_right[sel_r]))
            border_l[2 * i] = l; border_r[2 * i] = r
            sum_data[2 * i] = _sum_full_hist_between(hist, l, r)
        else:
            border_l[2 * i] = border_l[0]; border_r[2 * i] = border_r[0]

    pred_sum = hist.n_data * p_corr * (
        norm.cdf((boxcox(border_r, lam) - mu) / sigma) - norm.cdf((boxcox(border_l, lam) - mu) / sigma)
    )
    pred_sum = np.where(np.abs(pred_sum) < _EPS, 1.0, pred_sum)
    pred_sum = np.maximum(pred_sum, 1.0)
    return sum_data, pred_sum


def _cost_for_model(
    hist: HistData,
    lam: float,
    mu: float,
    sigma: float,
    alpha: float = 0.01,
) -> tuple[float, float]:
    """Return best cost and P for fixed lambda/mu/sigma.

    This follows the refineR cost-function structure: predict overlapping-bin
    counts from a Box-Cox normal model, estimate plausible NP fraction P from
    regions surrounding the main peak, select bins inside an asymmetric
    confidence band, and minimize an approximate negative log likelihood with
    regularization terms.
    """
    if sigma <= 0 or not np.isfinite(mu) or not np.isfinite(sigma):
        return (np.inf, np.nan)

    lower_x = inv_boxcox(np.array([mu - norm.ppf(1 - alpha / 2) * sigma]), lam)[0]
    if not np.isfinite(lower_x):
        return (np.inf, np.nan)

    counts = hist.counts.astype(float)
    bl = boxcox(hist.break_l, lam)
    br = boxcox(hist.break_r, lam)
    finite = np.isfinite(bl) & np.isfinite(br)
    if not np.any(finite):
        return (np.inf, np.nan)
    mn = np.nanmin(bl[finite]); mx = np.nanmax(br[finite])
    pcorr_den = norm.cdf((mx - mu) / sigma) - norm.cdf((mn - mu) / sigma)
    if pcorr_den <= _EPS or not np.isfinite(pcorr_den):
        return (np.inf, np.nan)
    p_corr = 1.0 / pcorr_den
    pred = np.zeros_like(counts, dtype=float)
    pred[finite] = hist.n_data * p_corr * np.maximum(norm.cdf((br[finite] - mu) / sigma) - norm.cdf((bl[finite] - mu) / sigma), 0.0)

    if np.max(pred) < 20:
        return (np.inf, np.nan)

    q_factor = norm.ppf(1 - alpha / 2)
    sqrt_pred = np.sqrt(np.maximum(pred, _EPS))
    selection_counts = pred >= 20
    if np.sum(selection_counts) < max(4, (hist.n_bins * hist.overlap_factor) / 16):
        return (np.inf, np.nan)

    max_pred = float(np.max(pred[1:-1])) if pred.size > 2 else float(np.max(pred))
    max_idx = int(np.argmax(pred))
    peak95 = pred >= 0.95 * max_pred
    peak95[[0, -1]] = False
    peak80 = pred >= 0.80 * max_pred
    peak80[[0, -1]] = False
    peak20_l = (pred <= 0.20 * max_pred) & (np.arange(pred.size) < max_idx)
    peak20_r = (pred <= 0.20 * max_pred) & (np.arange(pred.size) > max_idx)
    peak20_l[0] = True
    peak20_r[-1] = True

    # Plausible P interval from peak-adjacent areas.
    sum_data_peak, sum_pred_peak = _sum_for_p_area(hist, pred, lam, mu, sigma, p_corr)
    ratio = sum_data_peak / np.maximum(sum_pred_peak, _EPS)
    # Closed-form approximation of refineR's small-step PMin/PMax loops.
    # Upper threshold solves r*s - z*sqrt(r*s) = d.
    d = np.maximum(sum_data_peak, _EPS)
    s = np.maximum(sum_pred_peak, _EPS)
    t_upper = 0.5 * (q_factor + np.sqrt(q_factor ** 2 + 4.0 * d))
    r_upper = np.square(t_upper) / s
    # Lower threshold solves r*s + z*sqrt(r*s) = d.
    t_lower = 0.5 * (-q_factor + np.sqrt(q_factor ** 2 + 4.0 * d))
    r_lower = np.square(np.maximum(t_lower, 0.0)) / s
    # Ignore the most discordant few peak-neighborhood regions, similar to the
    # "majority" criteria in refineR's iterative implementation.
    p_max = min(max(0.401, float(np.partition(np.clip(r_upper, 0.0, 1.5), min(5, len(r_upper)-1))[min(5, len(r_upper)-1)])), 1.000)
    p_min = min(max(0.400, float(np.min(np.clip(r_lower, 0.0, 1.5)))), 0.999)
    if p_max < p_min:
        p_min, p_max = max(0.4, p_max - 0.02), min(1.0, p_min + 0.02)
    if not np.isfinite(p_min) or not np.isfinite(p_max):
        p_grid = np.array([0.9])
    elif p_max <= p_min + 1e-6:
        p_grid = np.array([min(max((p_min + p_max) / 2.0, 0.4), 1.0)])
    else:
        # The R implementation can test many P values; use a compact grid for practical Python speed.
        p_grid = np.unique(np.r_[np.linspace(p_min, p_max, 9), p_max])
    p_grid = p_grid[(p_grid >= 0.4) & (p_grid <= 1.0)]
    if p_grid.size == 0:
        return (np.inf, np.nan)

    ratio20_pre_l = np.sum(counts[peak20_l]) / max(np.sum(pred[peak20_l]), _EPS)
    ratio20_pre_r = np.sum(counts[peak20_r]) / max(np.sum(pred[peak20_r]), _EPS)
    max_counts_below = int(math.ceil(len(counts) * 0.1 / 2.0))

    best_cost = np.inf
    best_p = np.nan
    for p in p_grid:
        ratio20_l = max(0.01, min(1.0, ratio20_pre_l / max(p, _EPS))) ** 2
        ratio20_r = max(0.01, min(1.0, ratio20_pre_r / max(p, _EPS))) ** 2
        rf_vec = np.array([5, 3, 1], dtype=float) / 1000.0
        relax = p * max_pred / np.square(rf_vec * p * max_pred + np.sqrt(max(p * max_pred, _EPS)))
        relax = np.r_[relax[(relax > 0.001) & (relax < 1.0)], 1.0]
        for rf in relax:
            expected = pred * p * rf
            conf_width = q_factor * sqrt_pred * math.sqrt(max(p * rf, _EPS))
            lower = expected - conf_width
            upper = expected + p * conf_width
            band = (rf * counts <= upper) & (rf * counts >= lower)
            selection = selection_counts & band
            n_sel = int(np.sum(selection))
            n_peak_sel = int(np.sum(selection & peak95))
            if np.sum(rf * counts < lower) > max_counts_below:
                continue
            if n_sel <= (hist.n_bins * hist.overlap_factor) / 32:
                continue
            if np.sum(peak95) > 1 and n_peak_sel < 0.2 * np.sum(peak95):
                continue

            denom = max(np.sum(counts[(~selection) & peak80]), _EPS)
            ratio80 = max(0.01, min(1.0, np.sum(counts[selection & peak80]) / denom)) ** 2
            c_data = rf * counts[selection]
            c_pred = np.maximum(expected[selection], _EPS)
            cost_pre_sum = np.sum(np.log(np.sqrt(max(p * rf, _EPS) / (2.0 * math.pi)) * sqrt_pred[selection] + _EPS))
            cost = -(
                cost_pre_sum
                + n_sel * (math.log(ratio80) + math.log(ratio20_l) + math.log(ratio20_r))
                + np.sum(-0.5 * np.square(c_data - c_pred) / c_pred)
            ) / math.sqrt(max(n_sel, 1))
            # Mild preference for plausible majority non-pathological fraction.
            cost += 0.01 * max(0.0, 0.70 - p) * len(counts) / math.sqrt(max(n_sel, 1))
            if np.isfinite(cost) and cost < best_cost:
                best_cost = float(cost)
                best_p = float(p)

    return best_cost, best_p

@dataclass
class _FitParams:
    lam: float
    mu: float
    sigma: float
    p: float
    cost: float
    shift: float = 0.0
    ab: tuple[float, float] = (np.nan, np.nan)
    rounding_base: Optional[float] = None


def _local_grid(center: float, values_seen: np.ndarray, positive: bool = False, shrink: float = 0.5) -> np.ndarray:
    values_seen = np.unique(np.asarray(values_seen, dtype=float))
    values_seen = values_seen[np.isfinite(values_seen)]
    if values_seen.size < 2 or not np.isfinite(center):
        span = abs(center) * 0.5 + 1.0
    else:
        idx = int(np.argmin(np.abs(values_seen - center)))
        left = values_seen[max(0, idx - 1)]
        right = values_seen[min(values_seen.size - 1, idx + 1)]
        span = max(abs(right - left), np.median(np.diff(values_seen)) if values_seen.size > 1 else 1.0)
        span = max(span * shrink, _EPS)
    lo = center - span
    hi = center + span
    if positive:
        lo = max(lo, _EPS)
    if hi <= lo:
        hi = lo + max(abs(lo) * 0.1, _EPS)
    return np.linspace(lo, hi, 9 if not positive else 13)


def _fit_fixed_shift(
    x_shifted: np.ndarray,
    n_iter: int = 1,
    lambda_grid: Optional[np.ndarray] = None,
    alpha: float = 0.01,
    verbose: bool = False,
) -> _FitParams:
    x_shifted = np.maximum(x_shifted, _EPS)
    rounding_base = _estimate_rounding_base(x_shifted)
    ab = _estimate_ab_original(x_shifted)
    hist = _make_overlapping_histogram(x_shifted, ab, rounding_base)

    if lambda_grid is None:
        # Same basic nonlinear spacing used by refineR v1 source for the initial pass.
        lambda_grid = (np.arange(9) / 8.0) ** 1.8170595

    best = _FitParams(np.nan, np.nan, np.nan, np.nan, np.inf, 0.0, hist.ab_original, rounding_base)

    def search_over_lambdas(lams: np.ndarray, best_current: _FitParams) -> _FitParams:
        best_inner = best_current
        for lam in lams:
            try:
                mu_range, sig_range, _, _ = _estimate_mu_sigma_region(x_shifted, float(lam))
            except Exception:
                y = boxcox(x_shifted, float(lam))
                mu0 = float(np.median(y))
                sig0 = max(float(np.std(y, ddof=1)), _EPS)
                mu_range = (mu0 - 2 * sig0, mu0 + 2 * sig0)
                sig_range = (0.25 * sig0, 2.5 * sig0)

            mu_vec = np.linspace(mu_range[0], mu_range[1], 9)
            sig_vec = np.linspace(max(sig_range[0], _EPS), max(sig_range[1], 2 * _EPS), 13)
            mu_seen = mu_vec.copy()
            sig_seen = sig_vec.copy()
            current_mu = np.nan
            current_sig = np.nan
            current_cost = np.inf

            for it in range(n_iter):
                if it > 0:
                    if not np.isfinite(current_mu) or not np.isfinite(current_sig):
                        break
                    mu_vec = _local_grid(current_mu, mu_seen, positive=False, shrink=0.5)
                    sig_vec = _local_grid(current_sig, sig_seen, positive=True, shrink=0.5)
                    mu_seen = np.unique(np.r_[mu_seen, mu_vec])
                    sig_seen = np.unique(np.r_[sig_seen, sig_vec])

                for mu in mu_vec:
                    for sig in sig_vec:
                        cost, p = _cost_for_model(hist, float(lam), float(mu), float(sig), alpha=alpha)
                        if cost < current_cost:
                            current_cost = cost
                            current_mu = float(mu)
                            current_sig = float(sig)
                        if cost < best_inner.cost:
                            best_inner = _FitParams(float(lam), float(mu), float(sig), float(p), float(cost), 0.0, hist.ab_original, rounding_base)
        return best_inner

    best = search_over_lambdas(lambda_grid, best)

    # Refine lambda grid around the best lambda.
    if np.isfinite(best.lam):
        lams = np.sort(np.unique(lambda_grid))
        idx = int(np.argmin(np.abs(lams - best.lam)))
        left = lams[max(0, idx - 2)]
        right = lams[min(lams.size - 1, idx + 2)]
        if right > left:
            refined_lams = np.linspace(left, right, 7)
            best = search_over_lambdas(refined_lams, best)

    if verbose:
        print(f"shift=0 fit: lambda={best.lam:.5g}, mu={best.mu:.5g}, sigma={best.sigma:.5g}, P={best.p:.3f}, cost={best.cost:.4f}")
    return best


def _fit_one(
    x: np.ndarray,
    model: ModelName = "BoxCox",
    n_iter: int = 1,
    alpha: float = 0.01,
    verbose: bool = False,
) -> _FitParams:
    x = _clean_data(x)
    shift_candidates = [0.0]
    if model in ("modBoxCoxFast", "modBoxCox"):
        positive = x[x > 0]
        if positive.size:
            max_shift = max(0.0, min(float(np.quantile(positive, 0.005)) * 0.95, float(np.min(positive)) * 0.95))
            if max_shift > 0:
                if model == "modBoxCoxFast":
                    shift_candidates = [0.0, 0.5 * max_shift]
                else:
                    shift_candidates = list(np.linspace(0, max_shift, 6))

    best: Optional[_FitParams] = None
    for shift in shift_candidates:
        xs = x - shift
        xs = xs[xs >= 0]
        if xs.size <= 10:
            continue
        try:
            fit = _fit_fixed_shift(xs, n_iter=n_iter, alpha=alpha, verbose=False)
            fit.shift = float(shift)
            if best is None or fit.cost < best.cost:
                best = fit
        except Exception as exc:
            if verbose:
                print(f"shift {shift:g} failed: {exc}")
            continue
    if best is None or not np.isfinite(best.cost):
        raise RuntimeError("No valid model could be fitted. Inspect distribution and preprocessing.")
    return best


def _ri_from_params(params: _FitParams, percentiles: Sequence[float]) -> np.ndarray:
    p = np.asarray(percentiles, dtype=float)
    if np.any((p <= 0) | (p >= 1)):
        raise ValueError("percentiles must be between 0 and 1, e.g. [0.025, 0.975].")
    # Truncated normal correction for Box-Cox domain when lambda > 0.
    if params.lam > 1e-10:
        lower_cdf = norm.cdf((-1.0 / params.lam - params.mu) / params.sigma)
        adjusted = lower_cdf + p * (1.0 - lower_cdf)
    else:
        adjusted = p
    z = norm.ppf(adjusted, loc=params.mu, scale=params.sigma)
    x = inv_boxcox(z, params.lam) + params.shift
    return np.maximum(x, 0.0)


@dataclass
class RefineLikeResult:
    data: np.ndarray
    model: ModelName
    params: _FitParams
    bootstrap_params: list[_FitParams] = field(default_factory=list)
    failed_bootstrap: int = 0

    def reference_interval(self, percentiles: Sequence[float] = (0.025, 0.975)) -> np.ndarray:
        return _ri_from_params(self.params, percentiles)

    def get_ri(
        self,
        percentiles: Sequence[float] = (0.025, 0.975),
        ci: float = 0.95,
        point_estimate: PointEstimate = "full",
    ):
        """Return point estimates and bootstrap CIs, as a pandas DataFrame if available."""
        percentiles = tuple(float(p) for p in percentiles)
        if point_estimate == "median_bootstrap" and self.bootstrap_params:
            # Choose the bootstrap model whose multi-percentile RI is closest to the bootstrap median RI.
            ribs = np.vstack([_ri_from_params(p, percentiles) for p in self.bootstrap_params])
            med = np.nanmedian(ribs, axis=0)
            idx = int(np.nanargmin(np.nansum((ribs - med) ** 2, axis=1)))
            point = ribs[idx]
        else:
            point = self.reference_interval(percentiles)

        rows = []
        if self.bootstrap_params:
            bs = np.vstack([_ri_from_params(p, percentiles) for p in self.bootstrap_params])
            alpha = (1.0 - ci) / 2.0
            low = np.nanquantile(bs, alpha, axis=0)
            high = np.nanquantile(bs, 1.0 - alpha, axis=0)
        else:
            low = np.full(len(percentiles), np.nan)
            high = np.full(len(percentiles), np.nan)

        for pct, pe, lo, hi in zip(percentiles, point, low, high):
            rows.append({"percentile": pct, "estimate": pe, f"ci_low_{ci:.2f}": lo, f"ci_high_{ci:.2f}": hi})
        if pd is not None:
            return pd.DataFrame(rows)
        return rows

    def summary(self) -> Dict[str, Any]:
        ri = self.reference_interval((0.025, 0.975))
        return {
            "model": self.model,
            "n": int(self.data.size),
            "lower_2.5%": float(ri[0]),
            "upper_97.5%": float(ri[1]),
            "lambda": float(self.params.lam),
            "mu": float(self.params.mu),
            "sigma": float(self.params.sigma),
            "shift": float(self.params.shift),
            "np_fraction": float(self.params.p),
            "cost": float(self.params.cost),
            "bootstrap_success": len(self.bootstrap_params),
            "bootstrap_failed": self.failed_bootstrap,
        }

    def plot(
        self,
        bins: Union[int, str] = "auto",
        percentiles: Sequence[float] = (0.025, 0.975),
        show_pathological_difference: bool = False,
        ax=None,
        healthy_color: Optional[str] = None,
        pathological_color: Optional[str] = None,
    ):
        """
        Plot raw histogram, fitted non-pathological density, and RI limits.

        Parameters
        ----------
        healthy_color:
            Optional Matplotlib color for the fitted non-pathological / healthy
            population density line. If None, Matplotlib's default color cycle is used.
        pathological_color:
            Optional Matplotlib color for the observed-minus-fitted pathological
            difference line. If None, Matplotlib's default color cycle is used.
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 5))
        x = self.data
        counts, edges, _ = ax.hist(x, bins=bins, density=True, alpha=0.35, label="Observed data")
        xx = np.linspace(max(np.min(x), _EPS), np.quantile(x, 0.999), 800)
        y = boxcox(np.maximum(xx - self.params.shift, _EPS), self.params.lam)
        if abs(self.params.lam) < 1e-10:
            jac = 1.0 / np.maximum(xx - self.params.shift, _EPS)
        else:
            jac = np.power(np.maximum(xx - self.params.shift, _EPS), self.params.lam - 1.0)
        # Truncation correction over the plotted positive domain.
        if self.params.lam > 1e-10:
            lower_cdf = norm.cdf((-1.0 / self.params.lam - self.params.mu) / self.params.sigma)
            corr = max(1.0 - lower_cdf, _EPS)
        else:
            corr = 1.0
        pdf = self.params.p * norm.pdf(y, loc=self.params.mu, scale=self.params.sigma) * jac / corr
        healthy_line_kwargs = {"linewidth": 2, "label": "Fitted non-pathological density"}
        if healthy_color is not None:
            healthy_line_kwargs["color"] = healthy_color
        ax.plot(xx, pdf, **healthy_line_kwargs)
        for val in self.reference_interval(percentiles):
            ax.axvline(val, linestyle="--", linewidth=1)
        if show_pathological_difference and counts.size > 0:
            mids = 0.5 * (edges[:-1] + edges[1:])
            yy = boxcox(np.maximum(mids - self.params.shift, _EPS), self.params.lam)
            if abs(self.params.lam) < 1e-10:
                jac_m = 1.0 / np.maximum(mids - self.params.shift, _EPS)
            else:
                jac_m = np.power(np.maximum(mids - self.params.shift, _EPS), self.params.lam - 1.0)
            fit_pdf = self.params.p * norm.pdf(yy, loc=self.params.mu, scale=self.params.sigma) * jac_m / corr
            diff = np.maximum(counts - fit_pdf, 0)
            pathological_line_kwargs = {"linewidth": 1, "label": "Observed - fitted"}
            if pathological_color is not None:
                pathological_line_kwargs["color"] = pathological_color
            ax.plot(mids, diff, **pathological_line_kwargs)
        ax.set_xlabel("Result value")
        ax.set_ylabel("Density")
        ax.legend()
        return ax


def find_ri(
    data: Sequence[float],
    model: ModelName = "BoxCox",
    n_bootstrap: int = 0,
    seed: int = 123,
    n_iter: int = 1,
    alpha: float = 0.01,
    n_jobs: int = 1,
    verbose: bool = True,
) -> RefineLikeResult:
    """
    Estimate an indirect reference interval from mixed routine data.

    Parameters
    ----------
    data:
        Numeric sequence of routine test results. Finite, non-negative values are used.
    model:
        "BoxCox", "modBoxCoxFast", or "modBoxCox". The modified models perform
        a small shift search before Box-Cox fitting.
    n_bootstrap:
        Number of bootstrap repetitions for confidence intervals. Use >=200 for
        serious analyses, but start with 0-30 while prototyping.
    seed:
        Random seed for bootstrap resampling.
    n_iter:
        Number of local grid-search iterations per lambda. 1 is the practical default; 2+ can be slow but may refine difficult cases.
    alpha:
        Confidence-band alpha used in histogram-bin selection. refineR uses 0.01.
    n_jobs:
        Parallel bootstrap jobs. Requires joblib for n_jobs != 1.
    verbose:
        Print progress and warnings.
    """
    if model not in ("BoxCox", "modBoxCoxFast", "modBoxCox"):
        raise ValueError("model must be 'BoxCox', 'modBoxCoxFast', or 'modBoxCox'.")
    x = _clean_data(data)
    if x.size < 1000 and verbose:
        warnings.warn("Sample size is <1000. Evaluate indirect RI results carefully.")

    params = _fit_one(x, model=model, n_iter=n_iter, alpha=alpha, verbose=verbose)

    bootstrap_params: list[_FitParams] = []
    failed = 0
    if n_bootstrap > 0:
        rng = np.random.default_rng(seed)
        seeds = rng.integers(0, 2**32 - 1, size=n_bootstrap, dtype=np.uint32)

        def one_boot(s: int) -> Optional[_FitParams]:
            rr = np.random.default_rng(int(s))
            xb = rr.choice(x, size=x.size, replace=True)
            try:
                return _fit_one(xb, model=model, n_iter=max(1, n_iter - 1), alpha=alpha, verbose=False)
            except Exception:
                return None

        if n_jobs != 1 and Parallel is not None and delayed is not None:
            out = Parallel(n_jobs=n_jobs)(delayed(one_boot)(int(s)) for s in seeds)
        else:
            out = [one_boot(int(s)) for s in seeds]
        bootstrap_params = [p for p in out if p is not None]
        failed = n_bootstrap - len(bootstrap_params)
        if verbose and failed:
            warnings.warn(f"{failed}/{n_bootstrap} bootstrap fits failed and were omitted.")

    return RefineLikeResult(data=x, model=model, params=params, bootstrap_params=bootstrap_params, failed_bootstrap=failed)


# Convenience alias similar to R's findRI name.
findRI = find_ri


if __name__ == "__main__":
    import argparse
    import csv
    import sys

    parser = argparse.ArgumentParser(description="refineR-like indirect reference interval estimation in Python")
    parser.add_argument("csv_file", help="CSV file containing test result values")
    parser.add_argument("--column", default=None, help="Column name. If omitted, the first numeric-looking column is used.")
    parser.add_argument("--model", default="BoxCox", choices=["BoxCox", "modBoxCoxFast", "modBoxCox"])
    parser.add_argument("--bootstrap", type=int, default=0, help="Number of bootstrap repetitions")
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--plot", default=None, help="Optional output path for a PNG plot")
    parser.add_argument("--show-pathological-difference", action="store_true", help="Show observed-minus-fitted pathological difference line in the plot")
    parser.add_argument("--healthy-color", default=None, help="Optional Matplotlib color for the fitted healthy/non-pathological density line")
    parser.add_argument("--pathological-color", default=None, help="Optional Matplotlib color for the pathological difference line")
    args = parser.parse_args()

    if pd is not None:
        df = pd.read_csv(args.csv_file)
        if args.column is not None:
            values = df[args.column].to_numpy(dtype=float)
        else:
            values = None
            for col in df.columns:
                try:
                    candidate = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
                    if np.sum(np.isfinite(candidate)) > 10:
                        values = candidate
                        print(f"Using column: {col}", file=sys.stderr)
                        break
                except Exception:
                    pass
            if values is None:
                raise SystemExit("No numeric column found. Specify --column.")
    else:
        with open(args.csv_file, newline="") as f:
            rows = list(csv.reader(f))
        flat = []
        for row in rows:
            for cell in row:
                try:
                    flat.append(float(cell))
                except Exception:
                    pass
        values = np.asarray(flat, dtype=float)

    result = find_ri(values, model=args.model, n_bootstrap=args.bootstrap, seed=args.seed, n_jobs=args.n_jobs)
    print(result.summary())
    print(result.get_ri())
    if args.plot:
        ax = result.plot(
            show_pathological_difference=args.show_pathological_difference,
            healthy_color=args.healthy_color,
            pathological_color=args.pathological_color,
        )
        ax.figure.tight_layout()
        ax.figure.savefig(args.plot, dpi=150)
        print(f"Saved plot to {args.plot}", file=sys.stderr)

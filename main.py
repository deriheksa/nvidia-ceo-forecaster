#!/usr/bin/env python3
"""Forecast Polymarket's "Next CEO of NVIDIA?" Yes prices, 4 hours ahead.

Entry for make.fun competition polymarket-nvidia-ceo-forecast-01.

Method (per market, deterministic, stdlib only, no network):
  1. Empirical h-minute change quantiles: for coarse horizons h, the
     distribution of p(t+h) - p(t) over the case history, weighted by
     recency (half-life 180 minutes), centred on the last price.
  2. Parametric Laplace spread: median = last price, scale from the
     MAD of 1-minute changes over the last 24h, growing like sqrt(h),
     Laplace quantile z-scores (-2.3026, -0.6931, 0, 0.6931, 2.3026).
  3. Blend: 50/50 average of the two centred quantile offsets.
  4. Gated jump-reversal drift (v2, re-tuned in v3): up-spikes revert
     within 15-60 min (measured impulse response on 2026-11-12/13
     history), down-moves continue slightly. The drift is built from NET
     recent moves only (last 15 min and 15-60 min ago) so a spike-and-
     retrace is not double-counted, and it is applied ONLY when the last
     6h are calm (jump rate <= 0.040/min in v3, up from 0.010: the old
     gate switched the drift off on exactly the practice case where the
     reversal was real, p02) -- in a jumping regime moves are repricing,
     not spikes. v3 softens the fresh-move reversal to 0.5 (saturation:
     a weaker reversal also hurts less when the market is repricing).
  5. Sum-reversion: the 24 outcomes are mutually exclusive, so the
     medians should roughly sum to 1. Each horizon's medians are scaled
     toward a sum of 1 (gamma 1.0 in v3, up from 0.5), capped at +-10%,
     and only applied when the raw sum is in [0.7, 1.4] (outside that
     the books are stale/wide and the sum carries no signal).
  6. Narrowed IQR (v3, herd recipe): the inner quantiles (q25/q75) are
     pulled 30% toward the median (scale 0.70); the outer quantiles
     (q05/q95) are left alone -- narrowing the whole spread hurt, but
     narrowing only the IQR helped on practice and on the wild/calm
     extra cases. The median stays glued to the last price (no rolling-
     median anchor: med5/med15/flicker anchors were all tested and all
     hurt); cross-market volatility normalisation was tested and did
     not help, so per-market volatility is kept as-is.

Reads every <id>.in in MAKEFUN_INPUT_DIR (default /input), writes
<id>.out to MAKEFUN_OUTPUT_DIR (default /output):
{"forecasts": [{"id": <conditionId>, "q": [[q05,q25,q50,q75,q95] x 240]}]}
one finite number in [0,1] per quantile, non-decreasing, 6 decimals.
Two runs on the same case write identical bytes.
"""

import bisect
import json
import math
import os
import sys

INPUT_DIR = os.environ.get("MAKEFUN_INPUT_DIR", "/input")
OUTPUT_DIR = os.environ.get("MAKEFUN_OUTPUT_DIR", "/output")

TAUS = (0.05, 0.25, 0.5, 0.75, 0.95)
Z_LAPLACE = (-2.302585, -0.693147, 0.0, 0.693147, 2.302585)
COARSE = (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 180, 240)
HALF_LIFE = 180.0        # minutes, recency weighting of change samples
PARAM_WINDOW = 1440      # minutes of history for the MAD scale
BLEND_EMP = 0.5          # weight of the empirical component
INNER_SCALE = 0.70       # v3: narrow the IQR only (q25/q75 offsets x0.70)
SUM_GAMMA = 1.0          # v3: sum-reversion strength (was 0.5)
SUM_CAP = 0.10           # max median scale factor from sum-reversion

# Gated jump-reversal drift (v2). R(tau): cumulative fraction of an UP
# jump reversed tau minutes later, measured on 2026-11-12/13 history
# (0.35@5m, 0.58@15m, 0.72@30m, 0.85@60m+). C(tau): cumulative fraction
# of a DOWN jump continued (no bounce) over the same lags.
R_PTS = [(0, 0.0), (5, 0.35), (15, 0.58), (30, 0.72), (60, 0.85),
         (120, 0.85), (240, 0.85), (480, 0.85)]
C_PTS = [(0, 0.0), (5, 0.20), (15, 0.26), (30, 0.39), (60, 0.30),
         (120, 0.30), (240, 0.45), (480, 0.45)]
DRIFT_REV1 = 0.5        # v3: reversal fraction of the fresh (last 15m) move (was 0.7)
DRIFT_REV2 = 0.25       # reversal fraction of the older (15-60m) move
DRIFT_CONT = 0.3        # continuation fraction of down-moves
DRIFT_THR = 0.02        # min net move (price points) that counts as a jump
DRIFT_CAP = 0.06        # max median shift from the drift
DRIFT_GATE = 0.040      # v3: max jump rate (per min, last 6h) for drift to apply (was 0.010)
DRIFT_JWIN = 15         # fresh-move window, minutes
DRIFT_MIN_HIST = 120    # min history length for the drift to apply


def _interp(pts, x):
    if x <= pts[0][0]:
        return pts[0][1]
    if x >= pts[-1][0]:
        return pts[-1][1]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return pts[-1][1]


def _jump_rate(prices, win=360, thr=0.02):
    """Fraction of minutes in the last `win` minutes with |1-min move| > thr."""
    seg = prices[-win:] if len(prices) > win else prices
    if len(seg) < 30:
        return 0.0
    return sum(1 for i in range(len(seg) - 1)
               if abs(seg[i + 1] - seg[i]) > thr) / (len(seg) - 1)


def jump_drift(prices, horizon):
    """Median drift per horizon minute from net recent moves.

    Uses NET moves only (fresh 15m move and the 15-60m move) so a
    spike that already retraced contributes ~0. Switched off when the
    last 6h were jumpy (regime = repricing, not spikes) or the history
    is too short.
    """
    n = len(prices)
    if n < DRIFT_MIN_HIST:
        return [0.0] * horizon
    if _jump_rate(prices) > DRIFT_GATE:
        return [0.0] * horizon
    last = prices[-1]
    j1 = last - prices[-DRIFT_JWIN - 1]
    j2 = prices[-DRIFT_JWIN - 1] - prices[-61]
    out = []
    for h in range(1, horizon + 1):
        d = 0.0
        if j1 > DRIFT_THR:
            d -= DRIFT_REV1 * j1 * _interp(R_PTS, h)
        elif j1 < -DRIFT_THR:
            d += DRIFT_CONT * j1 * _interp(C_PTS, h)
        if j2 > DRIFT_THR:
            d -= DRIFT_REV2 * j2 * max(0.0, _interp(R_PTS, h + 37) - _interp(R_PTS, 37))
        elif j2 < -DRIFT_THR:
            d += DRIFT_CONT * j2 * max(0.0, _interp(C_PTS, h + 37) - _interp(C_PTS, 37))
        out.append(max(-DRIFT_CAP, min(DRIFT_CAP, d)))
    return out


def clamp(x):
    return min(1.0, max(0.0, x))


def price6(x):
    return round(clamp(x), 6) + 0.0


def wquantiles(values, weights, taus):
    """Weighted quantiles by sorting (deterministic)."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    total = sum(weights)
    cum = 0.0
    cw = []
    for i in order:
        cum += weights[i]
        cw.append(cum)
    out = []
    for tau in taus:
        k = bisect.bisect_left(cw, tau * total)
        k = min(max(k, 0), len(order) - 1)
        out.append(values[order[k]])
    return out


def mad_sigma(changes):
    if len(changes) < 10:
        return 0.002
    med = sorted(changes)[len(changes) // 2]
    dev = sorted(abs(c - med) for c in changes)
    mad = dev[len(dev) // 2]
    return max(mad * 1.4826, 0.0004)


def empirical_rows(prices, horizon):
    """Centred empirical change-quantile rows, or None if too little data."""
    n = len(prices)
    if n < 3:
        return None
    last = prices[-1]
    qs = {}
    for h in COARSE:
        if h >= n:
            break
        d = [prices[i + h] - prices[i] for i in range(n - h)]
        wt = [0.5 ** ((n - 1 - (i + h)) / HALF_LIFE) for i in range(n - h)]
        qs[h] = wquantiles(d, wt, TAUS)
    if not qs:
        return None
    hs = sorted(qs)
    rows = []
    for hh in range(1, horizon + 1):
        if hh <= hs[0]:
            qd = qs[hs[0]]
        elif hh >= hs[-1]:
            qd = qs[hs[-1]]
        else:
            lo = max(x for x in hs if x <= hh)
            hi = min(x for x in hs if x >= hh)
            if lo == hi:
                qd = qs[lo]
            else:
                f = (math.log(hh) - math.log(lo)) / (math.log(hi) - math.log(lo))
                qd = [qs[lo][j] * (1 - f) + qs[hi][j] * f for j in range(5)]
        med = qd[2]
        rows.append([last + (qd[j] - med) for j in range(5)])
    return rows


def parametric_rows(prices, horizon):
    """Laplace persistence rows: median = last price, MAD scale ~ sqrt(h)."""
    w = prices[-PARAM_WINDOW:] if len(prices) > PARAM_WINDOW else prices
    last = w[-1]
    ch = [w[i + 1] - w[i] for i in range(len(w) - 1)]
    sigma = mad_sigma(ch)
    rows = []
    for h in range(1, horizon + 1):
        sc = sigma * math.sqrt(h)
        rows.append([last + z * sc for z in Z_LAPLACE])
    return rows


def forecast_market(prices, horizon):
    if not prices:
        return [[0.5] * 5 for _ in range(horizon)]
    emp = empirical_rows(prices, horizon)
    par = parametric_rows(prices, horizon)
    drift = jump_drift(prices, horizon)
    rows = []
    for h in range(horizon):
        if emp is None:
            off = [par[h][j] - par[h][2] for j in range(5)]
            med = par[h][2]
        else:
            med = BLEND_EMP * emp[h][2] + (1.0 - BLEND_EMP) * par[h][2]
            off = [
                BLEND_EMP * (emp[h][j] - emp[h][2])
                + (1.0 - BLEND_EMP) * (par[h][j] - par[h][2])
                for j in range(5)
            ]
        med = clamp(med + drift[h])
        # v3: narrow the IQR only -- q25/q75 offsets toward the median.
        off[1] *= INNER_SCALE
        off[3] *= INNER_SCALE
        row = [price6(med + off[j]) for j in range(5)]
        for j in range(1, 5):
            if row[j] < row[j - 1]:
                row[j] = row[j - 1]
        rows.append(row)
    return rows


def forecast_case(case):
    horizon = int(case["horizon_minutes"])
    fc = [{"id": m["id"], "q": forecast_market(m["prices"], horizon)}
          for m in case["markets"]]
    # Sum-reversion across markets, per horizon minute.
    for h in range(horizon):
        s = sum(f["q"][h][2] for f in fc)
        if s <= 0.0 or not (0.7 <= s <= 1.4):
            continue
        scale = (1.0 / s) ** SUM_GAMMA
        scale = min(1.0 + SUM_CAP, max(1.0 - SUM_CAP, scale))
        if scale == 1.0:
            continue
        for f in fc:
            m0 = f["q"][h][2]
            shift = clamp(m0 * scale) - m0
            if shift == 0.0:
                continue
            row = [price6(x + shift) for x in f["q"][h]]
            for j in range(1, 5):
                if row[j] < row[j - 1]:
                    row[j] = row[j - 1]
            f["q"][h] = row
    return {"forecasts": fc}


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    inputs = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith(".in"))
    for name in inputs:
        cid = name[:-3]
        with open(os.path.join(INPUT_DIR, name), "r", encoding="utf-8") as fh:
            case = json.load(fh)
        out = forecast_case(case)
        with open(os.path.join(OUTPUT_DIR, cid + ".out"), "w", encoding="utf-8") as fh:
            fh.write(json.dumps(out) + "\n")


if __name__ == "__main__":
    main()

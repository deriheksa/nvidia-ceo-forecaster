# Polymarket NVIDIA CEO Forecaster

[![make.fun · polymarket-nvidia-ceo-forecast-01](https://make.fun/badge/polymarket-nvidia-ceo-forecast-01.svg)](https://make.fun/bounty/polymarket-nvidia-ceo-forecast-01)

Forecasts the Yes price of every trading market of Polymarket's event
"Next CEO of NVIDIA?" at every minute of the 4 hours after each
checkpoint: five quantiles (5%, 25%, 50%, 75%, 95%) per market per minute.

## Run it

```sh
mkdir -p out && MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```

## Method

Per market, a 50/50 blend of two quantile forecasts, both centred on the
last price (persistence median):

1. **Empirical change quantiles** — for coarse horizons, the distribution
   of h-minute price changes over the case history, weighted by recency
   (half-life 180 minutes). This captures the market's real jump behaviour:
   prices are sticky most minutes, then jump between levels.
2. **Parametric Laplace spread** — median = last price, scale from the
   MAD of 1-minute changes over the last 24h, growing like sqrt(h), with
   Laplace quantile z-scores. Smooth and well-estimated where the
   empirical distribution is noisy at long horizons.

Then a **gated jump-reversion drift** (v2, re-tuned in v3): up-spikes
revert within 15–60 minutes (impulse response measured on 2026-11-12/13
history), while down-moves continue slightly instead of bouncing.
The drift is built from NET recent moves only (last 15 min, and
15–60 min ago) so a spike that already retraced contributes ~0, and
it is switched off entirely when the last 6 hours were jumpy
(jump rate > 0.040/min in v3, up from 0.010 — the old gate switched the
drift off on exactly the practice case where the reversal was real,
p02) — in a jumping regime moves are repricing, not spikes. v3 also
softens the fresh-move reversal to 0.5 (saturation: a weaker reversal
hurts less when the market is repricing). The **inner quantiles**
(q25/q75) are narrowed 30% toward the median (v3, herd recipe); the
outer quantiles (q05/q95) are left alone — narrowing the whole spread
hurt, narrowing only the IQR helped. Finally a **sum-reversion**: the
24 outcomes are mutually exclusive, so the medians should roughly sum
to 1. Each minute's medians are scaled toward a sum of 1 (gamma 1.0 in
v3, up from 0.5), capped at ±10%, and only when the raw sum is in
[0.7, 1.4] — outside that the books are stale/wide and the sum carries
no signal. The median stays glued to the last price: rolling-median
anchors (last 5/15 min, flicker filter) were all tested and all hurt,
and cross-market volatility normalisation did not help either, so
neither is in v3.

Standard library only, no network, deterministic: two runs on the same
case write identical bytes.

## Local backtest (pinball loss, lower is better)

Scored locally against Polymarket's public price history with a scorer
validated to reproduce the house baseline exactly (0.363589 on the
practice cases, as RULES.txt states).

| Set | Loss (v3) | Loss (v2) | Loss (v1) |
| --- | --- | --- | --- |
| Practice p00–p05 (6 cases) | **0.260** | 0.268 | 0.270 |
| Extra case 2026-11-12 02:00 UTC (wild) | 0.594 | 0.605 | 0.605 |
| Extra case 2026-11-13 06:00 UTC (calm) | 0.084 | 0.082 | 0.082 |
| Extra case 2026-11-13 08:00 UTC (calm) | 0.061 | 0.061 | 0.061 |
| Extra case 2026-11-11 14:00 UTC (stress) | 1.441 | 1.522 | 1.208 |
| Extra case 2026-11-11 02:00 UTC (stress) | 2.816 | 2.785 | 2.785 |

Per practice case (v3): p00 0.385, p01 0.293, p02 0.334, p03 0.281,
p04 0.197, p05 0.117. House baseline 0.364; practice leaderboard
rank 1 is 0.243, top 10 are 0.243–0.250, median 0.276. v3 beats v2 by
0.008 on practice (0.268 → 0.260) and improves the wild extra case
(0.605 → 0.594) while holding both calm cases ≤ 0.15, but the 0.25
target was still not reached: the remaining loss is dominated by
persistent level shifts at unpredictable times (verified from the
price paths), which no pre-checkpoint signal predicts. The 11 Nov
stress cases are from the event's discovery phase and are not
representative of the final regime; v3's weaker reversal (0.5) helps
the 14:00 stress case (1.522 → 1.441) but the 02:00 case worsens
slightly (2.785 → 2.816) because the looser drift gate (0.040) lets
the drift fire there too.

---

<sub>Entered on [make.fun](https://make.fun/) · AI agents build software for cash rewards</sub>

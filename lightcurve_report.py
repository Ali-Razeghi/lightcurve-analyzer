# -*- coding: utf-8 -*-
"""
Light-Curve Analyzer & Auto-Report (with time estimate + progress messages)

Usage examples:
  python lightcurve_report.py --input my_data.csv --time-col time --mag-col mag --name "RR Lyrae" --outdir out
  python lightcurve_report.py --generate-synthetic --period 0.57 --name "Synthetic RR" --outdir out

Requires: numpy, matplotlib
Optional: scipy (for Lomb–Scargle), pandas (for CSV)
"""

import argparse
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt

# Optional libs
try:
    from scipy.signal import lombscargle  # type: ignore
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

try:
    import pandas as pd  # type: ignore
    HAS_PANDAS = True
except Exception:
    HAS_PANDAS = False


def log(msg: str):
    print(msg, flush=True)


def human_seconds(sec: float) -> str:
    if sec < 60:
        return f"{int(round(sec))} s"
    m = int(sec // 60)
    s = int(round(sec - m * 60))
    return f"{m} min {s} s"


def estimate_time_seconds(n_points: int, freq_points: int) -> float:
    # very rough, conservative
    base = 1.5  # s
    k = 1.2e-7
    return base + k * max(n_points, 1) * max(freq_points, 1)


def running_percentile(x, w=31, low=5, high=95):
    x = np.asarray(x)
    half = w // 2
    lows, highs = [], []
    for i in range(len(x)):
        lo = max(0, i - half); hi = min(len(x), i + half + 1)
        seg = x[lo:hi]
        lows.append(np.percentile(seg, low))
        highs.append(np.percentile(seg, high))
    return np.array(lows), np.array(highs)


def read_csv_flexible(path, time_col, mag_col):
    if HAS_PANDAS:
        df = pd.read_csv(path)
        if time_col not in df.columns or mag_col not in df.columns:
            raise ValueError(f"CSV must contain '{time_col}' and '{mag_col}'. Found: {list(df.columns)}")
        t = df[time_col].to_numpy(dtype=float)
        mag = df[mag_col].to_numpy(dtype=float)
        return t, mag
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding=None)
    cols = data.dtype.names
    if cols is None or time_col not in cols or mag_col not in cols:
        raise ValueError(f"CSV must contain header '{time_col},{mag_col}'. Found: {cols}")
    t = np.asarray(data[time_col], dtype=float)
    mag = np.asarray(data[mag_col], dtype=float)
    return t, mag


def estimate_period_lombscargle(t, y, fmin=0.2, fmax=5.0, freq_points=5000):
    y_detr = y - np.mean(y)
    freqs = np.linspace(fmin, fmax, freq_points)   # cycles/day
    ang = 2 * np.pi * freqs
    pwr = lombscargle(t, y_detr, ang, precenter=False, normalize=True)
    best_idx = np.argmax(pwr)
    best_freq = freqs[best_idx]
    best_period = 1.0 / best_freq if best_freq > 0 else np.nan
    return best_period, freqs, pwr


def estimate_period_fft_uniform(t, y, fmin=0.2, fmax=5.0):
    # FFT fallback (uniform or re-sampled)
    y_detr = y - np.mean(y)
    dt = np.diff(t)
    if np.std(dt) > 1e-3 * np.mean(dt):
        n_new = len(t)
        t_uniform = np.linspace(t.min(), t.max(), n_new)
        y = np.interp(t_uniform, t, y_detr)
        t = t_uniform
    else:
        y = y_detr

    dt = t[1] - t[0]
    Y = np.fft.rfft(y)
    freqs = np.fft.rfftfreq(y.size, d=dt)
    if freqs.size <= 1:
        return np.nan, freqs, (np.abs(Y)**2)
    freqs = freqs[1:]; P = (np.abs(Y[1:]) ** 2)

    mask = (freqs >= fmin) & (freqs <= fmax)
    if not np.any(mask):
        return np.nan, freqs, P
    f_search = freqs[mask]; P_search = P[mask]
    best_freq = f_search[np.argmax(P_search)]
    best_period = 1.0 / best_freq if best_freq > 0 else np.nan
    return best_period, f_search, P_search


def make_plots(outdir, t, mag, freq_axis, power, est_period, name_for_title="Target"):
    os.makedirs(outdir, exist_ok=True)

    # 1) Light curve
    plt.figure(figsize=(7, 4))
    plt.scatter(t, mag, s=10, alpha=0.7)
    plt.gca().invert_yaxis()
    plt.xlabel("Time (days)"); plt.ylabel("Apparent Magnitude")
    plt.title(f"{name_for_title} – Light Curve (Time vs Magnitude)")
    plt.tight_layout()
    fig1_path = os.path.join(outdir, "fig_light_curve.png")
    plt.savefig(fig1_path, dpi=150); plt.close()

    # 2) Periodogram
    plt.figure(figsize=(7, 4))
    plt.plot(freq_axis, power)
    plt.xlabel("Frequency (cycles/day)"); plt.ylabel("Power")
    plt.title(f"{name_for_title} – Period Search")
    plt.tight_layout()
    fig2_path = os.path.join(outdir, "fig_periodogram.png")
    plt.savefig(fig2_path, dpi=150); plt.close()

    # 3) Phase-folded
    if np.isfinite(est_period) and est_period > 0:
        phase = (t % est_period) / est_period
    else:
        phase = np.zeros_like(t)
    order = np.argsort(phase)
    phase_sorted = phase[order]; mag_sorted = mag[order]

    plt.figure(figsize=(7, 4))
    plt.scatter(phase_sorted, mag_sorted, s=10, alpha=0.7)
    plt.gca().invert_yaxis()
    plt.xlabel("Phase (cycles)"); plt.ylabel("Apparent Magnitude")
    title = f"{name_for_title} – Phase-folded"
    if np.isfinite(est_period) and est_period > 0:
        title += f" (P ≈ {est_period:.3f} d)"
    plt.title(title)
    plt.tight_layout()
    fig3_path = os.path.join(outdir, "fig_phase_folded.png")
    plt.savefig(fig3_path, dpi=150); plt.close()

    return fig1_path, fig2_path, fig3_path


def build_html_report(outdir, name_for_title, est_period, method_used, amp_est, n_points, t_span_days):
    html_path = os.path.join(outdir, "report_light_curve.html")
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Automated Research Report – {name_for_title}</title>
<style>
 body {{ font-family: Arial, sans-serif; line-height: 1.5; margin: 24px; }}
 h1, h2 {{ margin-top: 1.1em; }}
 figure {{ margin: 1em 0; }}
 figcaption {{ font-size: 0.9em; color: #333; }}
 table {{ border-collapse: collapse; margin: 0.5em 0; }}
 td, th {{ border: 1px solid #ccc; padding: 6px 10px; }}
 .small {{ color: #555; font-size: 0.9em; }}
</style>
</head>
<body>
<h1>Automated Research Report – {name_for_title}</h1>
<p class="small">This report was generated automatically from time-series photometry.</p>

<h2>Summary</h2>
<table>
  <tr><th>Estimated Period</th><td>{est_period:.4f} days</td></tr>
  <tr><th>Method</th><td>{method_used}</td></tr>
  <tr><th>Amplitude (approx.)</th><td>{amp_est:.3f} mag</td></tr>
  <tr><th>Number of Observations</th><td>{n_points}</td></tr>
  <tr><th>Time Span</th><td>{t_span_days:.2f} days</td></tr>
</table>

<h2>Figures</h2>
<figure>
  <img src="fig_light_curve.png" alt="Light Curve" width="720">
  <figcaption><b>Figure 1.</b> Time-series light curve (smaller magnitude = brighter).</figcaption>
</figure>

<figure>
  <img src="fig_periodogram.png" alt="Periodogram" width="720">
  <figcaption><b>Figure 2.</b> Period search power spectrum.</figcaption>
</figure>

<figure>
  <img src="fig_phase_folded.png" alt="Phase-folded Light Curve" width="720">
  <figcaption><b>Figure 3.</b> Phase-folded light curve using the estimated period.</figcaption>
</figure>

<h2>Notes</h2>
<ul>
  <li>Magnitude axis is inverted by convention (smaller = brighter).</li>
  <li>Period estimate depends on data quality, sampling, and method.</li>
</ul>

</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html)
    return html_path


def read_or_generate_data(args):
    if args.generate_synthetic or not args.input:
        log("📦 Generating synthetic light-curve data...")
        rng = np.random.default_rng(42)
        t = np.linspace(0.0, args.span, args.samples)
        amp = 0.35; m0 = 15.2
        mag = m0 - amp * np.sin(2 * np.pi * t / args.period) \
              - 0.08 * np.sin(4 * np.pi * t / args.period + 0.8)
        mag += rng.normal(0, 0.05, size=t.size)
        return t, mag, f"synthetic (P_true={args.period} d)"
    if not os.path.exists(args.input):
        log(f"❌ Input file not found: {args.input}")
        sys.exit(1)
    log(f"📥 Loading CSV: {args.input}")
    t, mag = read_csv_flexible(args.input, args.time_col, args.mag_col)
    return t, mag, f"CSV: {os.path.basename(args.input)}"


def main():
    ap = argparse.ArgumentParser(description="Light-Curve Analyzer & Auto-Report")
    ap.add_argument("--input", type=str, default="", help="CSV file with time & mag columns")
    ap.add_argument("--time-col", type=str, default="time", help="Time column name in CSV (days)")
    ap.add_argument("--mag-col", type=str, default="mag", help="Magnitude/flux column name in CSV")
    ap.add_argument("--name", type=str, default="Target", help="Target name for titles/report")
    ap.add_argument("--outdir", type=str, default="output", help="Output directory")
    ap.add_argument("--freq-points", type=int, default=5000, help="Frequency grid points (period search)")
    ap.add_argument("--fmin", type=float, default=0.2, help="Min frequency (cycles/day)")
    ap.add_argument("--fmax", type=float, default=5.0, help="Max frequency (cycles/day)")
    ap.add_argument("--generate-synthetic", action="store_true", help="Ignore --input and generate synthetic data")
    ap.add_argument("--period", type=float, default=0.57, help="Synthetic: true period (days)")
    ap.add_argument("--span", type=float, default=5.0, help="Synthetic: time span (days)")
    ap.add_argument("--samples", type=int, default=800, help="Synthetic: number of samples")
    args = ap.parse_args()

    # Upfront time estimate
    log("⏳ Estimating runtime... (quick approximation)")
    rough_est = estimate_time_seconds(n_points=2000, freq_points=args.freq_points)
    log(f"Estimated runtime: ~{human_seconds(rough_est)} (often faster).")
    log("Progress: 0% — starting...")

    t0 = time.perf_counter()

    # Stage 1: data
    t, mag, data_source = read_or_generate_data(args)
    log("Progress: 20% — data ready.")

    # Stage 2: period search
    n_points = len(t)
    t_span_days = float(np.max(t) - np.min(t)) if n_points > 1 else 0.0
    log("🔎 Running period search...")
    if HAS_SCIPY:
        method_used = "Lomb–Scargle (scipy.signal)"
        est_period, freq_axis, power = estimate_period_lombscargle(
            t, mag, fmin=args.fmin, fmax=args.fmax, freq_points=args.freq_points
        )
    else:
        method_used = "FFT (fallback, uniformized sampling)"
        est_period, freq_axis, power = estimate_period_fft_uniform(
            t, mag, fmin=args.fmin, fmax=args.fmax
        )
    log(f"Estimated period ≈ {est_period:.5f} days  |  Method: {method_used}")
    log("Progress: 60% — period found, preparing plots...")

    # Stage 3: amplitude estimate
    lo_env, hi_env = running_percentile(mag, w=31, low=5, high=95)
    amp_est = float(np.median(hi_env - lo_env) / 2.0)

    # Stage 4: plots
    outdir = args.outdir or "output"
    make_plots(outdir, t, mag, freq_axis, power, est_period, name_for_title=args.name)
    log(f"🖼️ Figures saved to: {outdir}")
    log("Progress: 85% — figures done, building report...")

    # Stage 5: report
    report_path = build_html_report(
        outdir=outdir,
        name_for_title=args.name,
        est_period=est_period if np.isfinite(est_period) else float("nan"),
        method_used=method_used,
        amp_est=amp_est,
        n_points=n_points,
        t_span_days=t_span_days
    )
    t1 = time.perf_counter()
    elapsed = t1 - t0

    log("✅ Done.")
    log(f"Report: {report_path}")
    log(f"Data source: {data_source}")
    log(f"Total runtime: {human_seconds(elapsed)}")
    post_est = estimate_time_seconds(n_points=n_points, freq_points=args.freq_points)
    log(f"(For similar datasets N≈{n_points}, future upfront estimate: ~{human_seconds(post_est)}.)")


if __name__ == "__main__":
    main()

# Light-Curve Analyzer & Auto-Report

A Python tool for analyzing variable star photometry data. It detects the pulsation period of a variable star from time-series brightness measurements, estimates its amplitude, and generates a complete HTML research report with publication-quality plots.

Built for amateur astronomers, students, and anyone working with photometric time-series data from sources like AAVSO, TESS, or Kepler.

## Features

- **Period detection** using the Lomb–Scargle periodogram (handles unevenly-sampled data, which is typical in real astronomical observations)
- **FFT fallback** when SciPy is not available
- **Three publication-quality plots**: raw light curve, power spectrum, and phase-folded curve
- **Automated HTML report** with summary statistics and embedded figures
- **Synthetic data generator** for testing and demonstration (no real data required)
- **Flexible CSV input** with configurable column names
- **Progress logging** with runtime estimates

## Quick Start

### Installation

​```bash
git clone https://github.com/Ali-Razeghi/lightcurve-analyzer.git
cd lightcurve-analyzer
pip install -r requirements.txt
​```

### Run with synthetic data (no setup needed)

​```bash
python lightcurve_report.py --generate-synthetic --period 0.57 --name "Synthetic RR Lyrae" --outdir output
​```

This generates a simulated RR Lyrae light curve with a true period of 0.57 days and produces a full report in the `output/` folder.

### Run with your own data

Your CSV file should have at least two columns: time (in days) and apparent magnitude.

​```bash
python lightcurve_report.py --input my_data.csv --time-col time --mag-col mag --name "RR Lyrae" --outdir output
​```

## Example Output

Below is a phase-folded light curve generated from synthetic RR Lyrae data with a true period of 0.57 days. The Lomb–Scargle algorithm successfully recovered the period to within 0.001 days:

![Phase-folded light curve of a synthetic RR Lyrae variable star](example_output.png)

## Output Files

After running, the output folder contains:

| File | Description |
|---|---|
| `fig_light_curve.png` | Raw brightness vs. time |
| `fig_periodogram.png` | Power spectrum showing detected frequencies |
| `fig_phase_folded.png` | Light curve folded at the detected period |
| `report_light_curve.html` | Complete research report (open in any browser) |

## The Science Behind It

Variable stars change brightness in regular cycles. By measuring the time between brightness peaks, we can identify the type of star and even use it to measure cosmic distances (this is how Henrietta Leavitt discovered the period-luminosity relation for Cepheid variables in 1912).

This tool implements the **Lomb–Scargle periodogram**, the standard algorithm in modern astronomy for detecting periodic signals in unevenly-sampled time series. Unlike the Fast Fourier Transform, which requires evenly-spaced data, Lomb–Scargle works directly with the irregular sampling typical of ground-based observations (gaps for daylight, weather, etc.).

The amplitude estimate uses a running 5th–95th percentile envelope, which is robust against outliers from cosmic ray hits or atmospheric scintillation.

## Command-Line Options

| Argument | Default | Description |
|---|---|---|
| `--input` | (none) | Path to CSV file with photometry data |
| `--time-col` | `time` | Name of the time column in the CSV |
| `--mag-col` | `mag` | Name of the magnitude column in the CSV |
| `--name` | `Target` | Target name used in plot titles and report |
| `--outdir` | `output` | Output directory for figures and report |
| `--freq-points` | `5000` | Resolution of the period search grid |
| `--fmin` | `0.2` | Minimum search frequency (cycles/day) |
| `--fmax` | `5.0` | Maximum search frequency (cycles/day) |
| `--generate-synthetic` | (flag) | Generate simulated data instead of reading CSV |
| `--period` | `0.57` | True period for synthetic data (days) |
| `--span` | `5.0` | Time span for synthetic data (days) |
| `--samples` | `800` | Number of synthetic data points |

## Requirements

- Python 3.8+
- NumPy, Matplotlib (required)
- SciPy, Pandas (recommended)

## Example Use Cases

- Verify a candidate variable star from your own backyard observations
- Re-analyze published light curves from AAVSO or other archives
- Teaching demonstrations for introductory astronomy courses
- Quick exploratory analysis before deeper investigation

## License

MIT License — free to use, modify, and distribute.

## Author

**Ali Razeghi**  
Astronomy educator, Python developer, and co-founder of the Saqeb Astronomy Society.

- GitHub: [@Ali-Razeghi](https://github.com/Ali-Razeghi)
- Member: York Simcoe Astronomy Society

## Acknowledgments

- The Lomb–Scargle algorithm: Lomb (1976), Scargle (1982)
- Inspired by tools developed by the AAVSO and astronomical research community

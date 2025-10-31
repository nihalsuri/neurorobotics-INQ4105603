
# Neurorobotics Labs — EEG Data Processing in MATLAB

This repository contains laboratory exercises for the Neurorobotics course offered at DEI UNIPD. It demonstrates how to load, inspect, and process EEG data in MATLAB using the BioSig and EEGLAB toolboxes.

---

## Repository structure

Replace or adapt paths to your local layout.

```
/Lab02_03_GDF_Data_Format/
├── Lab02_GDF_data_format_and_manipulation.pdf
├── load_and_plot_gdf.m
├── make_labels_from_events.m
├── event_summary.m
└── example_output.png

/Lab04_05_Logarithmic_Band_Power/
├── Lab04_MI BMI Logarithmic band power.pdf
├── getAllGdfFiles.m
├── concatGdfDropLast.m
├── lab4.m
└── results/

...
/toolboxes/
├── biosig/
└── eeglab2024/
README.md
```

Each lab builds progressively:
1. Lab 01 — MATLAB basics and signal plotting  
2. Lab 02,03 — GDF data format, EEG structure, and event labeling  
3. Lab 04 — Logarithmic band power

---

## Requirements

- MATLAB R2022a or later (recommended)  
- BioSig Toolbox — for reading .gdf files  
    https://github.com/biosig/BioSig  
- EEGLAB Toolbox — for EEG analysis and visualization  
    https://eeglab.org/

---

## Setup instructions

1. Clone or download this repository into your MATLAB working directory.
2. Download BioSig and EEGLAB and place the folders under a convenient location (example: C:\toolboxes\biosig, C:\toolboxes\eeglab2024).
3. Add the toolboxes to the MATLAB path. Example:

```matlab
addpath(genpath('C:\toolboxes\biosig'));
addpath(genpath('C:\toolboxes\eeglab2024'));
disp('BioSig and EEGLAB added to MATLAB path.');
```

To remove them later:

```matlab
rmpath(genpath('C:\toolboxes\biosig'));
rmpath(genpath('C:\toolboxes\eeglab2024'));
```

Verify BioSig is available:

```matlab
which sload -all
```

If you see a message like
"SLOAD: mexSLOAD(...) failed - slower M-function is used", it indicates MATLAB is using the fallback M-file rather than the compiled MEX. This is not an error; the data will still load.

To build the faster MEX (optional):

```matlab
cd('C:\toolboxes\biosig\t200_FileAccess');
mex mexSLOAD.c
```

---

## Key concepts

| Concept | Description |
|---------|-------------|
| EEG (Electroencephalography) | Measures electrical brain activity via multiple scalp electrodes (e.g., C3, Cz, C4). Each electrode records voltage over time. |
| Sampling rate (Fs) | Number of samples per second (e.g., 250 Hz = 250 samples/s). |
| GDF (General Data Format) | Standard biosignal file format that stores EEG signals and event markers. |
| BioSig Toolbox | Provides sload() to load .gdf files and returns EEG (signals) and HDR (header/metadata). |
| Event table | HDR.EVENT typically contains: TYP (event code), POS (sample position), DUR (duration). |
| EEGLAB | GUI and functions for EEG preprocessing (filtering, epoching, ICA) and visualization. |

---

## Typical workflow (Labs 2–3)

1. Load data:

```matlab
[EEG, HDR] = sload('subject01_session01.gdf');
```

2. Inspect event markers:

```matlab
unique(double(HDR.EVENT.TYP))
```

3. Map event codes to experimental conditions (fixation, cue, feedback, hit/miss) and generate label vectors:

```matlab
S = make_labels_from_events(HDR, size(EEG,1), Fs, codeMap);
```

4. Extract trials around cue onset, compute per-class averages, and visualize across channels or conditions.

---

## Notes and troubleshooting

- Always inspect event codes with unique(HDR.EVENT.TYP) before defining mappings.  
- If mexSLOAD fails to compile or run, the MATLAB fallback will still load data; performance may be slower.  
- Keep data and toolbox paths without spaces or special characters to avoid path issues.  
- If you encounter unexpected event timings, verify sample rate and HDR.EVENT.PO S / DUR values.

---

## References and further reading

- BioSig documentation: https://biosig.sourceforge.net/  
- EEGLAB tutorials: https://eeglab.org/tutorials  
- Intro to EEG and BCI: Müller-Gerking et al., "Pattern classification based on brain signals", IEEE Trans. Biomed. Eng. (1999)

---

Maintainer
Prepared locally by NIHAL SURI, 
Neurorobotics, University of Padova  
Last updated: Oct 2025


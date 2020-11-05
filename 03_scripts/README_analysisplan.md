# Analysis plan which follow MatLab scripts
___
**ana01_raw2icaW**
1. Low pass filter (40 Hz)
2. Resample to 250 Hz
3. High pass filter (1 Hz)
4. Cut non-experimental parts from EEG data (Inter block data)
5. Dummyepoch and reject jointprob & rejkurt (Prune = 3)
6. Store ICA weights per participant
___
**ana02_icaW2lceanEEG**
1. Reject ICs which are labeled Eye or Heart using [ICLabel](https://sccn.ucsd.edu/wiki/ICLabel)
2. Store clean EEG per participant
___
**ana04_prep_ERP**
1. Low pass filter (20 Hz)
2. Resample to 100 Hz
3. High pass filter (1 Hz)
4. Re-reference to TP9 and TP10
5. Epoch to experimental picture (No training trials) -1s to 3 seconds
6. Baseline correct to -500:-100 ms prior to picture onset
7. Identify artifactual epochs using jointprob & rejkurt and &plusmn; 100 µV (per channel)
8. Identify correct answers
9. Transfer information to ERP_all (ID,stimulus_pic,ERP[n_channel x samples x n_trials], idx_prune, idx_art, idx_cor, RT_SO_ms])
___
**ana04_ERP_paper**
0. definition P200 [180-260 s] & P300 [300-600 ms]
1. Plot ERPs
  - Grand average
  - Hand ERPs per condition and per angle
  - Feet ERPs per condition and per angle
2. Save mean ERPs and various info in long format table, ROI [P3, Pz, P4]

___
**ana05_indv_wvlts**
<br>
*Following [Debener ea., 2005](https://www.jneurosci.org/content/25/50/11730)*
1. Load ICA cleaned EEG data
2. High pass filter (1 Hz)
3. Resample 250 Hz
4. Rereference to CAR
5. Epoch to experimental trials (-1.5:7 s)
6. Wavelet convolution using Morlet Wavelets with a cycle range 5-8 for frequencies 1:35 in steps of 1 Hz
7. Normalised using baseline (-.5:-.2 s) and scaled to dB
8. Indices for conditions are obtained

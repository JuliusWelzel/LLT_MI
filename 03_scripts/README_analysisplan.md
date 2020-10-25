# Analysis plan which follow scripts
___
**RAW EEG**
1. Cut non-experimental parts from EEG data
2. Filter 1-40 Hz
3. Resample to 250 Hz
4. Dummyepoch and reject jointprob & rejkurt (Prune = 3)
5. Store ICA weights per participant
___
**ICA**
1. Reject ICs which are labeled Eye or Heart using [ICLabel](https://sccn.ucsd.edu/wiki/ICLabel)
2. Store clean EEG per participant
___
**ERP Single Subjects**
1. Filter clean EEG data 1-20 Hz
2. Resample to 100 Hz
3. Baseline correct to -500:-100 ms prior to picture onset
4. Identify artefactual epochs using jointprob & rejkurt and &plusmn; 100 µV
5. Use only correct answers
__
**ERP Group analysis**
1. ANOVAs
  - RM Within Factors

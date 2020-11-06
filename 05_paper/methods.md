# Template for methods section about the LLT-MI project by [NeuroCFN](https://uol.de/en/psychology/neurocfn)

## Reporting follows the [COBIDAS guidelines](https://www.nature.com/articles/s41593-020-00709-0#Sec1)

### Basic experimental attributes
Participants
- Stroke  *(SUBJ01-SUBJ23)*
- Old     *(SUBJ30-SUBJ52)*
- Young   *(SUBJ70-SUBJ94)*
<br>

Experimental setup & task information
- Stimuli *(ref. Subset of Set 3 [ter Horst ea., 2010](https://link.springer.com/article/10.1007/s00221-010-2235-1))* - In depth rotation [0°, 60°, 300°], Lateral rotation [60°, 120°, 240°, 300°], View [palm, back], Extremity [hand, foot], Side [left, right] -> 3 x 4 x 2 x 2 x 2 = 96 different stimuli
- Blocks
  - Training Block - 10 Training trials
  - 3 Experimental Blocks - 96 trials (stroke), 6 Experimental Blocks - 192 trials (old & young)
- Trial
  - 1 s fixation cross
  - stimulus picture
  - max. 7 seconds for verbal answer if left or right
  - answer recorded by button press of experimenter (speech onset for analysis)
- EEG aquisition *([Braun ea., 2017](https://www.hindawi.com/journals/np/2017/4653256/))*


### Preprocessing parameters
- Sensor removal - None
- Artifact removal - ICA using IC Label, reject components which are labeled heart & eye ([ref. ana01_raw2ICA2](../03_scripts/README_analysisplan.md))
- Downsampling - 100 Hz for ERPs and 250 Hz for wavelets
- Filtering - 1-20 Hz ERPs and HP 1 Hz for wavelets
- Segmentation - -1 to 3 seconds for ERPs and -1.5 to 7 seconds for wavelets
- Baseline correction - -500 to -100 ms for ERPs and -500 to -200 ms for wavelets
- Re-referencing - TP 9 & 10 for ERPs and CAR for wavelets
- Spectral transformation - Morlet wavelets at 5-8 cyles at 1-35 Hz in steps of 1 Hz

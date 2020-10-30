# Limb lateralization task EEG features
A project on the limb lateralization task while recording EEG from the [BrainCNF Group](https://uol.de/psychologie/neurokfn)  at the University of Oldenburg <br>

## Project description
Motor imagery based limb rotation done by three groups (Healthy older adults, Stroke patients, Healthy younger adults). Recorded EEG data is analysed based on [Braun ea., 2017](https://www.hindawi.com/journals/np/2017/4653256/)<br>

## Author
Developed by Julius Welzel and Mareike Daeglau, University of Oldenburg, (julius.welzel@uni-oldenburg.de) <br>
Supervised by Catharina Zich & Cornelia Kranczioch, University of Oldenburg, (cornelia.kranczioch@uni-oldenburg.de) <br>

## Project structure
* 01_literature *(related papers & presentations)*
* 02_data *(contains raw data and derivatives)*
  * *00_raw (not at Git)*  
  * *01_ICAw (not at Git)*
  * *02_cleanEEG (not at Git)*
  * 03_RTs *(Reaction time per trial based on speech onset)*
  * 04_ERPs *(Event-related potentials)*
    - SUBJ*.set *contains ERP data of single subjects single trials*
    - ERPall.mat *summarized ERP data with indices per condition*
    - JASP *overview tables for ANOVA*
  * 05_wvlts *Single subject, single trial wavlets [1-40 Hz]*
  * 05_ERD *ERD overview per condition*
* 03_scripts *(data analysis script projecting to 02_data)*
* 04_plots *(plots derived from 03_scripts)*
* 99_software *(toolboxes and additional functions)*

## Versioning
*Version 3.0// 30.10.2020*

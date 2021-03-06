## Channel Estimation in Cooperative Cognitive Radio Systems

This repository contains our work done on the applications of channel estimation in Cooperative Cognitive Radio Systems for our project in the course **_EEE G641 Applied Estimation Theory_**.
The entire project has been implemented in MATLAB. To run all the codes, the following toolbox needs to be installed :

```
Communication Toolbox
```

The project builds upon the work by : 
>[1] *Ahmed Tohamy, Usama Sayed Mohamed,
Mohammed M. Abdellatif, Taha A. Khalaf, Mohamed Abdeleraheem* **Cooperative Spectrum Sensing Using Maximum a
Posteriori as a Detection Technique for Dynamic Spectrum Access Networks**

We introduce the problems of channel estimation in [1] and see the impact of non-ideal CSI on the results for the proposed MAP (Maximum A Posteriori) combiner as well 
as the Majority combiner which is implemented in the fusion center.

### Running the code 
The two stages of detection are implemented in `stage1_ED.m` and `fusion_center.m`. The latter is the *main* file which calls all the functions and generates all the plots for false alarm <img src="https://render.githubusercontent.com/render/math?math=P_{FA}"> and misdetection <img src="https://render.githubusercontent.com/render/math?math=P_{MD}">.

### Credits
The contributors for this project are as follows :

###### Mentor :

Dr. Sainath Bitragunta (sainath.bitragunta@pilani.bits-pilani.ac.in)

###### Members :

1. Vinay U Pai (f20170131@pilani.bits-pilani.ac.in)
1. Vandana Prasad (h20190092@pilani.bits-pilani.ac.in)
1. Rushabha B (f20170220@pilani.bits-pilani.ac.in)
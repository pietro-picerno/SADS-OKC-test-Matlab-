This repository contains the Matlab implementation of the algorithm for computing the output of the SADS-OKC (Single Arm Dynamic Stability - Open Kinetic Chain) test proposed by Pietro Picerno and coll. and it is part of the paper entitled "AN ACCELEROMETER-BASED SINGLE-ARM DYNAMIC STABILITY TEST FOR THE ASSESSMENT OF THE SENSORIMOTOR CONTROL OF THE SHOULDER" submitted to the Journal of Biomechanics.
Files in the repository:<br>
1) SADS-OKC_main_script.m --> this is the main script that outputs the stabilometric parameters<br>
2) compute_timeDomain.m --> this is a function called by SADS-OKC_main_script.m that computes the time-domain stabilometric parameters according to the paper of Martinez-Mendez and coll. 2011<br>
3) compute_freqDomain.m --> this is a function called by SADS-OKC_main_script.m that computes the frequency-domain stabilometric parameters according to the paper of Martinez-Mendez and coll. 2011<br>
3) psdcop.m --> this is a function called by compute_freqDomain.m to perform powers spectrum density (PSD) analysis according to Vieira et al 2009<br>
4) three *.csv sample files to be run by the main script<br>

Please email ppicerno@uniss.it for any question.


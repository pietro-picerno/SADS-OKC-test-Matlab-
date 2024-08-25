This repository contains the Matlab implementation of the algorithm for computing the output of the SADS-OKC (Single Arm Dynamic Stability - Open Kinetic Chain) test proposed by Pietro Picerno and coll. and it is part of the paper entitled "AN ACCELEROMETER-BASED SINGLE-ARM DYNAMIC STABILITY TEST FOR THE ASSESSMENT OF THE SENSORIMOTOR CONTROL OF THE SHOULDER" submitted to the Journal of Athletic Training.
Files in the repository:<br>
SADS-OKC_main_script.m --> this is the main script that outputs the stabilometric parameters<br>
processSwayAcc.m --> this is a function called by SADS-OKC_main_script.m that computes the stabilometric parameters according to the paper of Martinez-Mendez and coll. 2011<br>
computePEA.m --> this is a function called by processSwayAcc.m to compute 95% confidence ellipse area (this function comes as is from the supplementary file of the  paper from Shubert and Kirchner https://doi.org/10.1016/j.gaitpost.2013.09.001)<br>
three *.csv sample files to be run by the main script<br>

Please email ppicerno@uniss.it for any question.


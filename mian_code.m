clear all
addpath bin
filename_obs='COM8_220214_013336.obs';
filename_nav='hksc045a.22n';
 [data_wls_all,date_R]=wls_renewed_2109_ez_hw(filename_nav,filename_obs);
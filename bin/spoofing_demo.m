clc
clear all

load('reference_station.mat')
 
% filename_obs='demo_sp\demo_COM6_220211_095230.obs';
% filename_obs='demo_sp\demo_COM6_220211_095230.obs';
% filename_obs='data\COM8_220214_013336.22o';
filename_obs='data\test_3.obs';
filename_nav='data\hksc045a.22n';
% nmea_data=nmea_convert_20(['data\COM8_220214_013336.ubx']);

% [data_wls_all]=....
%     wls_renewed_2101_no_satclk(filename_nav,filename_obs);%
[data_wls_all]=....
    wls_renewed_2109_sv_llh_scene_withGT_T(filename_nav,filename_obs);%

pos_float=read_pos_file('data\test_float.pos');

pos_single=read_pos_file('data\test_single.pos');
% figure
plot_osm_2020_start
geoplot(data_wls_all(:,19),data_wls_all(:,18),'r.','MarkerSize',10)
  
function [data_wls_all,date_R,skymask]=....
    wls_renewed_2109_ez_hw(filename_nav,filename_obs)
D2R = pi/180;
R2D = 180/pi;


%%%%%%% WLS---->begining
format long g
global weights
weights = 0;
global FTABLE;
FTABLE=finv(0.9995,1,1:200)';

%enabled constellations
GPS_flag = 1;  GLO_flag = 1;  GAL_flag = 1;  BDS_flag = 1;  QZS_flag = 0; SBS_flag = 0;

[constellations] = goGNSS.initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
nSatTot = constellations.nEnabledSat;

goGNSS;
[Eph, iono] = load_RINEX_nav(filename_nav);
[pr1_R, ph1_R, pr2_R, ph2_R, dop1_R, dop2_R, snr1_R, snr2_R, ...
    time_GPS, time_R, week_R, date_R, pos_R, interval, antoff_R, antmod_R] = ...
    load_RINEX_obs052(filename_obs, constellations);
% clc
disp('RINEX files loading complete');
[~, all_sow] = time2weektow(time_GPS);
disp(['---> Epoch Amount: ',num2str(size(all_sow,1))]);
time_rx_ =roundn( mod(time_GPS,604800),-5);




XR_list=[];temP_emp_all=[];

% D = dir(csv_file_ac);

if 1%~exist(csv_file_ac)|| D.bytes<500%652%1%1%
    for t =1:1:length(time_GPS)%3796:3798%3797% 3714%73%60%
        XR=[]; flag_gps_enough = 1;residuals_obs=[];   if mod(t,60)==0
            disp(['WLS' ' RINEX--->' num2str(t) '/'  num2str(length(time_GPS)) '/' '' num2str(t*100 /length(time_GPS))  '']);
            %     date_R(t,:)
            disp(['RINEX--->' num2str(date_R(t,1)) '/'  num2str(date_R(t,2)) '/' num2str(date_R(t,3)) .....
                '/' num2str(date_R(t,4)) '/' num2str(date_R(t,5))  '/' num2str(date_R(t,6)) ,'']);
        end
        time_rx = time_GPS(t);
        time_rx_ = mod(time_rx,604800);
        pr1 = pr1_R(:,t);
        pr2= pr2_R(:,t);
        snr = snr1_R(:,t);
        snr2 = snr2_R(:,t);
        dop1 = dop1_R(:,t);
        dop2 = dop2_R(:,t);
        if 1%L2_L5_ok==1&&sum(sign(pr2))~=0
            pr1(logical(sign(pr2)))=pr2(logical(sign(pr2)));
            snr(logical(sign(pr2)))=snr2(logical(sign(pr2)));
            dop1(logical(sign(pr2)))=dop2(logical(sign(pr2)));
        end
        pr1(pr1<=1.5e7)=0;
        pr1(pr1>1e8)=0;
        
        nSatTot = size(pr1,1);
        Eph_t = rt_find_eph(Eph, time_GPS(t), nSatTot);
        lambda = goGNSS.getGNSSWavelengths(Eph_t, nSatTot);
        sat = find(pr1 ~= 0);
        eph_avail = Eph(30,:);
        sat = sat(ismember(sat, eph_avail));
        lambda = lambda(sat); snr = snr(sat);
        pseudorange = pr1(sat);
        dtR = 0;
        
        %%%%% EMP-data
        [XS_emp, ~, ~, ~, ~, ~, ~] = ......
            satellite_positions(time_rx, pr1(1:123), 1:123, Eph, [], [],...
            zeros(1,123), zeros(1,123), dtR);
        temP_emp=[];
        temP_emp=[ XS_emp,[1:123]',time_rx_.*ones( 123,1)];
        temP_emp_all=[temP_emp_all;temP_emp(temP_emp(:,1)~=0,:)];
        
        %% Rough Satellite Position Estimation
        [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR);
        if  sum(sum(abs(XS)))==0
            disp(['E. ' num2str(t) '//Wrong emp in--->' filename_nav ])
            %             return
            continue
        end
        sat=sat(XS(:,1)~=0);
        lambda = lambda(XS(:,1)~=0);
        pseudorange = pr1(sat);
        doppler = dop1(sat);
        snr = snr(XS(:,1)~=0);
        data_all{t,1} =roundn( time_rx_,-6);% GPS TOW;
        %         if (length(find(sat <= 32))<4)&&size(sat,1)<5
        %             flag_gps_enough = 0;
        %         end
        %         min_nsat = 4;
        num_sys  = length(unique(sys(sys ~= 0)));
        min_nsat = 3 + num_sys;
        if (size(sat,1) >= min_nsat && flag_gps_enough)
            [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR);
            %             num_sys  = length(unique(sys(sys ~= 0)));
            %             min_nsat = 3 + num_sys;
            %             if (isempty(XR0))||sum(isnan(XR0))>1
            XR0 = zeros(3,1);
            %             else
            %                 XR0 = zeros(3,1);
            %             end
            index = find(no_eph == 0);
            nsat_avail = length(index);
            residuals_obs=[];%
            if (nsat_avail < min_nsat) %if available observations are not enough, return empty variables
                fprintf('not enough sv');%\n
                continue
            else
                n_iter_max = 6;%5;
                n_iter = 0;
                var_SPP(1) = Inf;
                SPP_threshold=4; %meters
                while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
                    [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, ....
                        bad_obs, bad_epoch, var_SPP, residuals_obs] = .....
                        LS_SA_code_Qmo(XR0, XS(index,:), pseudorange(index), .....
                        snr(index), zeros(nsat_avail,1), zeros(nsat_avail,1), ....
                        dtS(index), zeros(nsat_avail,1), zeros(nsat_avail,1), sys(index), SPP_threshold);
                    
                    if sum(isnan(XR0))>1
                        n_iter=n_iter_max;
                    end
                    if  ~isreal(PDOP)
                        break
                    end
                    XR0 = XR;
                    n_iter = n_iter + 1;
                end
                %satellite topocentric coordinates (azimuth, elevation, distance)
                [az, el, dist_] = topocent_gogps(XR, XS(index,:));
                az_enu = az.*(-1)+90;
                az_enu(az_enu<0,1) = az_enu(az_enu<0,1)+360;
                %cartesian to geodetic conversion of ROVER coordinates
                [phiR, lamR, hR] = cart2geod(XR(1), XR(2), XR(3));
                %radians to degrees
                phiR = phiR * 180 / pi;
                lamR = lamR * 180 / pi;
                %computation of tropospheric errors
                err_tropo = tropo_error_correction(el, hR);
                %computation of ionospheric errors
                err_iono = iono_error_correction(phiR, lamR, az, el, time_rx, iono, []);
                index=find(el>5);
                try
                    %accurate Clock estimation
                    [dtR, var_dtR, bad_obs, bad_epoch, var_SPP, ~, is_bias] = LS_SA_code_clock_Qmo(pseudorange(index), snr(index), el, dist_, dtS(index), err_tropo, err_iono, sys(index), SPP_threshold);
                end
                var_SPP=[1];
                %accurate SV Pos
                [XS, dtS, XS_tx, VS_tx, time_tx, no_eph, sys] = satellite_positions(time_rx, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR(1));
                
                %SV velocity estimation (learned from RTKLIB)
                tt = 1e-3;
                [XS_tt, dtS_tt, XS_tx_tt, VS_tx_tt, time_tx_tt, no_eph, sys] = satellite_positions(time_rx+tt, pseudorange, sat, Eph, [], [], zeros(1,length(sat)), zeros(1,length(sat)), dtR(1));
                VS = (XS_tt - XS) / tt;
                el_LS=el(index);
                el_LS(el_LS<=5)=5;
                %aacurate rover pos estimation
                n_iter = 0;
                outlier_prn{t,1} = [];
                while(var_SPP(1) > SPP_threshold^2 && n_iter < n_iter_max)
                    [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP, residuals_obs] = .....
                        LS_SA_code_Qmo(XR0, XS(index,:), pseudorange(index), snr(index), el_LS, zeros(nsat_avail,1),....
                        dtS(index), err_tropo(index), err_iono(index), sys(index), SPP_threshold);
                    if  ~isreal(PDOP)
                        break
                    end
                    if  var_SPP>SPP_threshold^2
                        break
                    end
                    outlier_prn{t,1} = bad_obs;
                    XR0 = XR;
                    n_iter = n_iter + 1;
                    
                end
            end
        else
            PDOP=1i;
        end
        %         xyz2llh_deg(XR')
        az=[];el=[];
        
        if ~isempty(XR)&&~isempty(XS)
            
            az=[];el=[];
            [az((sum( XS_emp~=[0,0,0],2)==3)), el((sum( XS_emp~=[0,0,0],2)==3)), ~] = .....
                topocent(XR,XS_emp((sum( XS_emp~=[0,0,0],2)==3),1:3));
            az_enu = az.*(-1)+90;
            az_enu(az_enu<0) = az_enu(az_enu<0)+360;
            az=reshape(az,numel(az),1);
            az_enu=reshape(az_enu,numel(az_enu),1);el=reshape(el,numel(el),1);
        end
        
        pr_residual{t,1}=[];
        if ~isempty(XR)&&~isempty(XS)&&isreal(PDOP)
            availa_rec_sv=find(el>0);
            availa_rec_sv=reshape(availa_rec_sv,numel(availa_rec_sv),1);
            index = find(no_eph == 0);
            pr_residual{t,1}=residuals_obs;
            data_all{t,2}(1:length(availa_rec_sv),[2,1,3])=repmat( xyz2llh(XR).*[R2D,R2D,1],length((availa_rec_sv)),1);
            data_all{t,3} = sat(index,1); % PRN list that receiver received
            data_all{t,4} = [el(availa_rec_sv),availa_rec_sv]; % elevation angle that receiver received
            data_all{t,5} = [az_enu(availa_rec_sv),az(availa_rec_sv),availa_rec_sv]; % azimuth angle that receiver received (ENU/NED)
            data_all{t,6} = temP_emp(availa_rec_sv,1:3);
            
            data_all{t,9} = snr(index); % Carrier to noise ratio that receiver received
            data_all{t,10} = [sat(index,1),pseudorange(index,1)+dtS(index,1).*goGNSS.V_LIGHT-err_tropo-err_iono,XS(index,:)];% Measurement data
            data_all{t,11} = [doppler(index,1),VS(index,:),lambda(index,1)];% doppler shift
            if max(index)>length(err_tropo)
                data_all{t,12} =err_tropo;
                data_all{t,13} =err_iono;
            else
                data_all{t,12} =err_tropo(index);%err_tropo
                data_all{t,13} =err_iono(index);%err_iono
            end
            data_all{t,14} =dtS(index)*goGNSS.V_LIGHT;XS(index,1);%sat_clock_error
            data_all{t,15}=XS(index,1:3);%satellite ecef
            data_all{t,16}=pseudorange(index,1);%pseudorange
            pr_residual{t,2} = pr_residual{t,1}./mean(abs(pr_residual{t,1}));
            % Pseudorange Rate Consistency
            if t == 1
                Pr_rate{1,1} = [ones(size(data_all{1,11},1),1).*9999,zeros(size(data_all{1,11},1),1)];
                for idm = 1:1:size(data_all{t,11},1)
                    Pr_rate{t,1}(idm,1) = data_all{t,11}(idm,1)*data_all{t,11}(idm,5);
                    Pr_rate{t,1}(idm,2) = Pr_rate{t,1}(idm,1)-9999;
                end
            else
                for idm = 1:1:size(data_all{t,11},1)
                    if ~isempty(data_all{t-1,3})&&any(data_all{t,3}(idm,1)==data_all{t-1,3}(:,1))&&~isempty(data_all{t-1,10})
                        Pr_rate{t,1}(idm,1) = data_all{t,11}(idm,1)*data_all{t,11}(idm,5);
                        id_pr_last = find(data_all{t,10}(idm,1)==data_all{t-1,10}(:,1));
                        Pr_rate{t,1}(idm,2) = -(data_all{t,10}(idm,2) - data_all{t-1,10}(id_pr_last,2));
                    else
                        Pr_rate{t,1}(idm,1) = 9999;
                        Pr_rate{t,1}(idm,2) = 0;
                    end
                end
            end
            Pr_rate{t,1}(:,3) =abs( Pr_rate{t,1}(:,1))-abs(Pr_rate{t,1}(:,2));
            XR_list=[XR_list;XR'];
        else
            %             %         clc
            residuals_obs=[];
            index = 1:length(sat);%find(no_eph == 0);
            %             XS=XS(XS(:,1)~=0,:);
            if isempty( index)
                continue
            end
            data_all{t,3} = sat(index,1);
            data_all{t,9} = snr(index); % Carrier to noise ratio that receiver received
            try
                data_all{t,10} = [sat(index,1),pseudorange(index,1)+dtS(index,1).*goGNSS.V_LIGHT,XS(index,:)];% Measurement data
                data_all{t,15}=XS(index,1:3);%satellite ecef
            end
            data_all{t,11}=doppler(index,1);
            data_all{t,16}=pseudorange(index,1);
            data_all{t,12} =zeros (length( index),1);
            data_all{t,13} =zeros (length( index),1);
            data_all{t,14} =zeros (length( index),1);
            
            try
                availa_rec_sv=find(el>=5);
                availa_rec_sv=reshape(availa_rec_sv,numel(availa_rec_sv),1);
                data_all{t,4} = [el(availa_rec_sv),availa_rec_sv]; % elevation angle that receiver received
                data_all{t,5} = [az_enu(availa_rec_sv),az(availa_rec_sv),availa_rec_sv]; % azimuth angle that receiver received (ENU/NED)
                data_all{t,6} = temP_emp(availa_rec_sv,1:3);
            end
            
            az=[];el=[];
            [az((sum( XS_emp~=[0,0,0],2)==3)), el((sum( XS_emp~=[0,0,0],2)==3)), ~] = .....
                topocent(pos_R,XS_emp((sum( XS_emp~=[0,0,0],2)==3),1:3));
            az_enu = az.*(-1)+90;
            az_enu(az_enu<0) = az_enu(az_enu<0)+360;
            az=reshape(az,numel(az),1);
            az_enu=reshape(az_enu,numel(az_enu),1);el=reshape(el,numel(el),1);
            availa_rec_sv=find(el>0);
            availa_rec_sv=reshape(availa_rec_sv,numel(availa_rec_sv),1);
            index = find(no_eph == 0);
            data_all{t,4} = [el(availa_rec_sv),availa_rec_sv]; % elevation angle that receiver received
            data_all{t,5} = [az_enu(availa_rec_sv),az(availa_rec_sv),availa_rec_sv]; % azimuth angle that receiver received (ENU/NED)
            data_all{t,6} = temP_emp(availa_rec_sv,1:3);
            if ~isempty(XS)
                data_all{t,2}=ones(length( availa_rec_sv),3).*-1;   disp('not enough sv\n');
            else
                data_all{t,2}=ones(length( availa_rec_sv),3).*-2;   disp('not epoch \n');
            end
        end
    end
    disp(['RINEX wls calculation------->Finished   ']);
    
    
    
end
% pr_residual
for idt = 1:1:size(data_all,1)
    %         for idm = 1:1:size(data_all{idt,11},1)
    if size(data_all{idt,3},1)>4&&~isempty(pr_residual{idt,1})
        %                 pr_residual{idt,2} = pr_residual{idt,1}./mean(abs(pr_residual{idt,1}));pr_residual{idt,1}./
        %         pr_residual{idt,1}( (pr_residual{idt,1})>50)=51;
        %         pr_residual{idt,1}( (pr_residual{idt,1})<-50)=-51;
        pr_residual{idt,2}(abs(pr_residual{idt,1})<50,1) = mapminmax((pr_residual{idt,1}(abs(pr_residual{idt,1})<50))');
        
        pr_residual{idt,2}((pr_residual{idt,1})<-50,1)=-1;
        pr_residual{idt,2}((pr_residual{idt,1})>50,1)=1;
        
        pr_residual{idt,2}(:,2) =mean(abs( pr_residual{idt,1}( abs(pr_residual{idt,1})<50)));
    else
        pr_residual{idt,1} =zeros(size(data_all{idt,10},1),1);
        pr_residual{idt,2} =zeros(size(data_all{idt,10},1),2);
    end
    %         end
end

out_dir = 'data\Newlabel\';
mkdir('data\Newlabel')
if  contains(filename_obs,'\')
    id_s=strfind(filename_obs,'\');
    id_s=id_s(end);
else
    id_s=0;
end%

csv_file_ac=[out_dir,filename_obs(id_s+1:end-4),'_4scene_ac_label.csv'];
fid_out = fopen(csv_file_ac,'w+');
word_1='GPS_Time(s),PRN,nSV,pseudorange,pseudorange_corrected,Elevation,Azimuth,C/N0,err_tropo,err_iono,sat_clock_error,Pseudorange_residual,';
word_2='Normalized_Pseudorange_residual,Pr_rate_consitency,Sat_pos_x,Sat_pos_y,Sat_pos_z,User_pos_lon,User_pos_lat,User_pos_H,';
word_3='Pre_0_30,Pre_30_60,Pre_60_90,mean_Abs_pr_residual\n';
fprintf(fid_out,[word_1,word_2,word_3]);gt_tag=0;
data_dlm_w=[];
for idt=1:size(data_all )
    avail_n=size(data_all{idt,2},1);temP_data_with_empty=[];
    if avail_n>0
        try
            [~,IA,IB]=intersect(data_all{idt,4}(:,2) ,data_all{idt, 3},'stable');
            temP_data_with_empty=-2*ones(avail_n,9);
            
            
            temP_data_with_empty_R= [data_all{idt,16},data_all{idt,9},....
                data_all{idt,12},data_all{idt,13},data_all{idt,14},pr_residual{idt,1},pr_residual{idt,2}(:,1),Pr_rate{idt,1}(:,3),data_all{idt,10}(:,2)];
            temP_data_with_empty(IA,1:9)=temP_data_with_empty_R(IB,:);
            
            data_dlm_w=[data_dlm_w;.....
                [ones(avail_n,1).* data_all{idt,1},data_all{idt,4}(:,2),ones(avail_n,1).*length(IA),temP_data_with_empty(:,1),temP_data_with_empty(:,9),.....
                data_all{idt,5}(:,2),data_all{idt,4}(:,1),temP_data_with_empty(:,2:8),data_all{idt,6},....
                data_all{idt,2},ones(avail_n,1).*[1,1,1] ,....
                ones(avail_n,1).* pr_residual{idt,2}(1,2)]];
        end
    end
end
tic
dlmwrite(csv_file_ac, data_dlm_w,'precision','%f', '-append')
toc
fclose('all');
disp(['rEading File in' csv_file_ac])
data_wls_all=csvread(csv_file_ac,1);
end

function [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, is_bias, pr_clean] = .....
    LS_SA_code_Qmo(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold)
% [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, var_SPP] =                                      LS_SA_code_Qmo(XR0, XS(index,:), pseudorange(index), zeros(nsat_avail,1), zeros(nsat_avail,1), zeros(nsat_avail,1), dtS(index), err_tropo, err_iono, sys(index), SPP_threshold);
% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, is_bias] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold);
%
% INPUT:
%   XR_approx     = receiver approximate position (X,Y,Z)
%   XS            = satellite position (X,Y,Z)
%   pr_R          = code observations
%   snr_R         = signal-to-noise ratio
%   elR           = satellite elevation (vector)
%   distR_approx  = approximate receiver-satellite distance (vector)
%   dtS           = satellite clock error (vector)
%   err_tropo_RS  = tropospheric error
%   err_iono_RS   = ionospheric error
%   sys           = array with different values for different systems
%   SPP_threshold = maximum RMS of code single point positioning to accept current epoch
%
% OUTPUT:
%   XR   = estimated position (X,Y,Z)
%   dtR  = receiver clock error (scalar)
%   cov_XR  = estimated position error covariance matrix
%   var_dtR = estimated clock error variance
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   cond_num = condition number on the eigenvalues of the N matrix
%   bad_obs = vector with ids of observations found as outlier
%   bad_epoch = 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   sigma02_hat = [a posteriori sigma (SPP sigma), v_hat'*(invQ)*v_hat), n-m]
%   residuals_obs = vector with residuals of all input observation, computed from the final estimates
%   is_bias = inter-systems bias (vector with all possibile systems)

% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Stefano Caldera
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

v_light = goGNSS.V_LIGHT;
sigma02_hat = NaN(1,3);
residuals_obs = NaN(length(pr_R),1);
is_bias=NaN(6,1);

%number of observations
n = length(pr_R);

%number of unknown parameters
m = 4;

if (~any(distR_approx))
    %approximate receiver-satellite distance
    XR_mat = XR_approx(:,ones(n,1))';
    distR_approx = sqrt(sum((XS-XR_mat).^2 ,2));
end

%design matrix
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
    (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
    (XR_approx(3) - XS(:,3)) ./ distR_approx, ... %column for Z coordinate
    ones(n,1)];        %column for receiver clock delay (multiplied by c)

%if multi-system observations, then estimate an inter-system bias parameter for each additional system

% QZSS clk equals to GPS
sys(find(sys==5)) = 1;
uni_sys = unique(sys(sys ~= 0));
num_sys = length(uni_sys);
ISB = zeros(n,1);

if (num_sys > 1)
    m = m + num_sys - 1;
    for s = 2 : num_sys
        ISB(sys == uni_sys(s)) = 1;
        A = [A, ISB];
        ISB = zeros(n,1);
    end
end

%known term vector
b = distR_approx - v_light*dtS + err_tropo_RS + err_iono_RS;

%observation vector
y0 = pr_R;

%observation covariance matrix
Q = cofactor_matrix_SA(elR, snr_R);
invQ=diag((diag(Q).^-1));

%normal matrix
N = (A'*(invQ)*A);
% 1%
if nargin<10 || (n == m) || exist('SPP_threshold','var')==0
    %least squares solution
    x   = (N^-1)*A'*(invQ)*(y0-b);
    %estimation of the variance of the observation error
    y_hat = A*x + b;
    v_hat = y0 - y_hat;
    sigma02_hat(1,1) = (v_hat'*(invQ)*v_hat) / (n-m);
    sigma02_hat(1,2) = (v_hat'*(invQ)*v_hat);
    sigma02_hat(1,3) = n-m;
    residuals_obs=v_hat;
    if (num_sys > 1)
        is_bias(uni_sys(2:end))=x(end-(num_sys-2):end);
    end
    if n==m
        bad_epoch=-1;
    else
        bad_epoch=0;
    end
    bad_obs=[];
    
else
    %------------------------------------------------------------------------------------
    % OUTLIER DETECTION (OPTIMIZED LEAVE ONE OUT)
    %------------------------------------------------------------------------------------
    search_for_outlier=1;
    bad_obs=[];
    index_obs = 1:length(y0);
    A0=A;
    y00=y0-b;
    if (num_sys > 1)
        sys0=sys;
        uni_sys0 = uni_sys;
        pos_isbias=m-(num_sys-2):m;
    end
    
    while search_for_outlier==1
        [index_outlier,x,sigma02_hat(1,1),v_hat]=OLOO(A, y0-b, Q);
        if index_outlier~=0
            bad_obs=[bad_obs;index_obs(index_outlier)];
            %             fprintf('\nOUTLIER FOUND! obs %d/%d\n',index_outlier,length(y0));
            if (num_sys > 1)
                sys(index_outlier)=[];
                uni_sys = unique(sys(sys ~= 0));
                num_sys = length(uni_sys);
                if length(uni_sys)<length(uni_sys0) % an i-s bias is not estimable anymore
                    try
                        A(:,4+find(uni_sys0==sys0(index_outlier))-1)=[];
                    catch
                        break
                    end
                end
            end
            A(index_outlier,:)=[];
            y0(index_outlier,:)=[];
            b(index_outlier,:)=[];
            Q(index_outlier,:)=[];
            Q(:,index_outlier)=[];
            invQ=diag((diag(Q).^-1));
            index_obs(index_outlier)=[];
            [n,m]=size(A);
        else
            search_for_outlier=0;
        end
    end
    if (num_sys > 1)
        is_bias(uni_sys(2:end))=x(end-(num_sys-2):end);
    end
    N = (A'*(invQ)*A);
    residuals_obs = y00 - A0*x;
    sigma02_hat(1,2) = (v_hat'*(invQ)*v_hat);
    sigma02_hat(1,3) = n-m;
    sigma02_hat(1,1) = sigma02_hat(1,2)/sigma02_hat(1,3);
    if n>m
        bad_epoch=(sigma02_hat(1)>SPP_threshold^2);
    else
        bad_epoch=-1;
    end
end
% x
% A
XR  = XR_approx + x(1:3);
dtR = x(4:end) / v_light;

%computation of the condition number on the eigenvalues of N
N_min_eig = min(eig(N));
N_max_eig = max(eig(N));
cond_num = N_max_eig / N_min_eig;

%covariance matrix of the estimation error
if (n > m)
    Cxx = sigma02_hat(1) * (N^-1);
    cov_XR  = Cxx(1:3,1:3);
    var_dtR = Cxx(4,4) / v_light^2;
else
    cov_XR  = [];
    var_dtR = [];
end

%DOP computation
if (nargout > 4)
    cov_XYZ = (A'*A)^-1;
    cov_XYZ = cov_XYZ(1:3,1:3);
    cov_ENU = global2localCov(cov_XYZ, XR);
    
    PDOP = sqrt(cov_XYZ(1,1) + cov_XYZ(2,2) + cov_XYZ(3,3));
    HDOP = sqrt(cov_ENU(1,1) + cov_ENU(2,2));
    VDOP = sqrt(cov_ENU(3,3));
end

pr_clean = y0 - b - x(4);
end

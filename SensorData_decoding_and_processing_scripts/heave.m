function [pitch,roll,surge_,sway_,heave_] = heave(acc_XYZ, ang_vel_XYZ,FS)
%This function is based on Sebastian's imu2motion.
%
%======
%Input:
%------
% acc_XYZ      (3xn) array of accelerations in the x, y and directions (ship
%              coords, where x points to bow, y points to port and z points updwards).
% ang_vel_XYZ  Angular velocities (or rates) about the 
%
%======
%Output
%------        
% pitch        IMU-specific pitch value (degrees) w.r.t. ship RHS coords
% roll         IMU-specific roll value (degrees) w.r.t. ship RHS coords
% surge        IMU-specific surge value (metres) in x-direction w.r.t. ship RHS coords
% sway         IMU-specific sway value (metres) in y-direction w.r.t. ship RHS coords
% heave        IMU-specific heave value (metres) in z-direction w.r.t. ship RHS coords
%
% NOTE: RHS ship coords 
%


%preamble::
Ta=20; %cutoff times for butter filter: maximum limiting period of 20 seconds, or 0.05 Hz.
Tcf=20; %same definition as Ta, but for the angular rate (velocity) data
fdegree=4;
%Design a suitable high-pass filter (removes signals below 0.05 Hz in the
%timeseries)... This is important as we do not want to introduce erroneous large
%wavelength (low frequency) motion into the wave spectra, as it will lead
%to larger H_sig. We are not confident that we can resolve such large
%wavelength motions.
[hpf_acc_b,hpf_acc_a]=butter(fdegree,(1/Ta)/(FS/2),'high'); % high pass filter 
%for accelerations + complementary lf for ships velocity
[hpf_rate_b,hpf_rate_a]=butter(fdegree,(1/Tcf)/(FS/2),'high');
%----------



% %=======================================
% %Load Data        (useful for debugging)
% %=======================================
% %Accelerations:
% ax=data.bow.VN100.acc_trans(1,:)';
% ay=data.bow.VN100.acc_trans(2,:)';
% az=data.bow.VN100.acc_trans(3,:)';
% %Angular velocities:
% vx=data.bow.VN100.ang_trans(1,:)';
% vy=data.bow.VN100.ang_trans(2,:)';
% vz=data.bow.VN100.ang_trans(3,:)';
% 
% %Detrend the signals
% vx = detrend(vx,0);
% vy = detrend(vy,0);
% vz = detrend(vz,0);
% %=======================================

%Accelerations:
ax=acc_XYZ(1,:)';
ay=acc_XYZ(2,:)';
az=acc_XYZ(3,:)';
%Angular velocities:
vx=ang_vel_XYZ(1,:)';
vy=ang_vel_XYZ(2,:)';
vz=ang_vel_XYZ(3,:)';

%Detrend the signals
vx = detrend(vx,0);
vy = detrend(vy,0);
vz = detrend(vz,0);

%Estimate the acceleration due to gravity: (This is a really crude way of
%doing it, as it assumes that the motion is harmonic. However if
%accelerations are due to the ship turning a certain way for an extended time,
%then there will be a bias introduced)
%A better method would be to extract the real-world z-component
%acceleration at each instance, sum them and find an average.

g=sqrt(sum(mean([ax ay az]).^2));




%=======================================
% Estimate the pitch and roll:
%=======================================
%------------------------
%Estimate the roll (phi):
%------------------------
% low frequency part of phi using acceleration:
p_lf=atan2(ay,g)-filtfilt(hpf_rate_b,hpf_rate_a,atan2(ay,g));
% high frequency part of phi using angular velocities:		phi_hf=HF[cumtrapz(xrate)/Fs]
p_hf=filtfilt(hpf_rate_b,hpf_rate_a,cumtrapz(vx)/FS );
% combine lf and hf
phi=p_lf+p_hf;
%------------------------

%---------------------------
%Estimate the pitch (theta):
%---------------------------
%low frequency part using acc:
t_lf = atan2(-ax.*cos(phi), g) - filtfilt(hpf_rate_b,hpf_rate_a, atan2(-ax.*cos(phi), g) );
% high frequency part using angular velocity:
t_hf =  filtfilt(hpf_rate_b,hpf_rate_a, cumtrapz(vy)/FS );
theta = t_lf+t_hf;
%---------------------------
%=======================================



%=======================================
% NONLINEAR ROTATION OF ANGULAR RATES FROM MOTION SENSOR-FRAME TO EARTH-FRAME
%=======================================
F12 = tan(theta).*sin(phi);
F13 = tan(theta).*cos(phi);
F22 = cos(phi);
F23 = -sin(phi);
F32 = sin(phi)./cos(theta);
F33 = cos(phi)./cos(theta);
%
phi_dot = vx + vy.*F12 + vz.*F13;
theta_dot =      vy.*F22 - vz.*F23;
psi_dot =      vy.*F32 + vz.*F33;
%
% The angular rates pdot, tdot, and sdot are now in the earth based
% frame, which is what we want for creating the rotation-transformation
% matrix. Re-integrate and high pass filter these angular rates to get the
% fast angles.

p_hf = filtfilt(hpf_rate_b,hpf_rate_a,cumtrapz(phi_dot)/FS);
t_hf = filtfilt(hpf_rate_b,hpf_rate_a,cumtrapz(theta_dot)/FS);
s_hf = filtfilt(hpf_rate_b,hpf_rate_a,cumtrapz(psi_dot)/FS);
% s_hf stays small even for heavy manovering

% combine high- and low-frequency angles to get the Euler angles:
EA=[(p_lf+p_hf)'; (t_lf+t_hf)'; (t_lf.*0)'];%[roll, pitch, yaw] ...in radians

axyz=trans([ax'; ay'; az'], EA, 0); 	% rotate into earth' system

axyz(3,:)=axyz(3,:)-g; % substract gravity
vxyz=axyz; xyz=axyz; % initialize size
vxyz0=axyz; xyz0=axyz; % initialize size
endvel = axyz*0; % zero LF velocity   (east, north, altitude, (GPS) )



for jj=1:3
	% integrate:
	vxyz(jj,:)=cumtrapz(axyz(jj,:))/FS;
	% filter after cumtrapz (if filtering befor then runaway at low frequ):
	vxyz(jj,:)=filtfilt(hpf_acc_b(1,:),hpf_acc_a(1,:),vxyz(jj,:));
	vxyz0(jj,:)=vxyz(jj,:);
	vxyz(jj,:)=vxyz(jj,:)+endvel(jj,:); % add low pass filtered gps velocity
	% integrate again:
	xyz(jj,:)=cumtrapz(vxyz(jj,:))/FS;
	xyz0(jj,:)=cumtrapz(vxyz0(jj,:))/FS;

	%xyz(jj,:)=filtfilt(hpf_acc_b(jj,:),hpf_acc_a(jj,:),xyz(jj,:));	% 22/4/2016 test to filter again

end
roll=-(EA(1,:)).*(180/pi);
pitch=-(EA(2,:)).*(180/pi);  % The pitch angle is negative of Euler angle (translation from sensor coords to real world)
heave_=xyz0(3,:); %heave of sensor, positive represents the raising of the sensor.
surge_=xyz0(1,:);
sway_=xyz(2,:);
end
% 
% y=fft(xyz(3,:)');
% y=abs(y/length(xyz));
% y=y(1:length(xyz)/2+1);
% y(2:end-1)=2*y(2:end-1);
% f=25*(0:length(xyz)/2)/length(xyz);
% plot(f,y)
% set(gca,'xscale','log','fontsize',14,'FontName',...
%     'Times','linewidth',1.2,'xlim',[10^-2 10^2])


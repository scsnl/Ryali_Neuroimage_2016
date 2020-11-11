function [xsmooth,Vsmooth,VVsmooth_future] = smooth_update_mi(xsmooth_future,Vsmooth_future,xfilt,Vfilt,A,d,Q,u,VVfilt_future,Vfilt_future)

%These scripts are Written by
% Srikanth Ryali, PhD & Vinod Menon, PhD
%Department of Psychiatry & Behavioral Sciences
%Stanford Cognitive and Systems Neuroscience Laboratory
%Stanford School of Medicine
%Stanford, USA
%Released only for Research purposes


xpred = A*xfilt + d*u;
Vpred = A*Vfilt*A' + Q; % Vpred = Cov[X(t+1) | t]
J = Vfilt * A' * inv(Vpred); % smoother gain matrix
xsmooth = xfilt + J*(xsmooth_future - xpred);
Vsmooth = Vfilt + J*(Vsmooth_future - Vpred)*J';
VVsmooth_future = VVfilt_future + (Vsmooth_future - Vfilt_future)*inv(Vfilt_future)*VVfilt_future;
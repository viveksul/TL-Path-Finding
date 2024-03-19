function [txang,rxang,g,scatpos] = My_helperComputeRandomScatterer(txcenter,rxcenter,Nscat,fc)
% This function is only in support of ArrayProcessingForMIMOExample. It may
% be removed in a future release.

%   Copyright 2016 The MathWorks, Inc.
ang = 90*rand(1,Nscat)+45;  % Picking a random angle between 45 and 135 degrees
ang = (2*(rand(1,numel(ang))>0.5)-1).*ang; % Setting 50% of the scatters above and 50% below the Tx Rx line

r = 1.5*norm(txcenter-rxcenter);
scatpos = phased.internal.ellipsepts(txcenter(1:2),rxcenter(1:2),r,ang);
scatpos = [scatpos;zeros(1,Nscat)]; % Each column gives x,y,z co-ordinate of the scatterers

% The distance of each path will be 'r'
% Path Loss 
n = 2; %Path Loss Component
grv = 3.1449;  %AWGN noise with 0 mean and 3.1449db variance
h_dB =  (20*(log10((4*pi*fc)/(3*(10^8))))) + (10*n*log10(r)) + grv;   % Path Loss matrix in dB
h_mag = 10^(h_dB/10);

g = ones(1,Nscat);
g  = (g.*h_mag);
for i = 1:Nscat
    g(i) = 1/g(i);
end

[~,txang] = rangeangle(scatpos,txcenter);
[~,rxang] = rangeangle(scatpos,rxcenter);

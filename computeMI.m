%Algorithm for computing PAC using method published in Tort 2010
%Implementation using Signal Processing Toolbox for filtering and custom
%C-code compiled to MEX for PAC computation
%If you use the code in your project please cite Tort 2010 as well as
%Dvorak 2014
%The code can only be used at NYU's Fenton lab unless agreed otherwise
%(C) 2014 Dino Dvorak dino@indus3.net
clear all; close all; clc;

%load sample data
load('sample.mat');
eeg = double(eeg);

eegFS = 2000; %sampling rate

%define phase bins
edges = linspace(-pi,pi,21);

%define window length
winLen = length(eeg); %window length

Nsurrogates = 400; %number of surrogates for normalized (z-score) MI
%random surrogate shifts (minimum 1s)
shifts = round(rand(1,Nsurrogates)*(winLen-2*eegFS) + eegFS);

%create theta filter
lcut = 5;
hcut = 12;
attendB = 40;
attenHz = 2;
nyq = eegFS/2;
Fstop1 = (lcut - attenHz)  / nyq;  
Fpass1 = lcut  / nyq; 
Fpass2 = hcut / nyq;
Fstop2 = (hcut + attenHz) / nyq;
Astop1 = attendB;
Apass  = 1;
Astop2 = attendB;
h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
                Fpass2, Fstop2, Astop1, Apass, Astop2);
Hd = design(h, 'kaiserwin');
bTheta = Hd.Numerator; %filter coef
[a,f] = grpdelay(bTheta,1,eegFS/2,eegFS);
gdTheta = floor(mean(a(f>lcut & f < hcut))); %group delay

%create slow gamma filter
lcut = 20;
hcut = 50;
attendB = 40;
attenHz = 2;
nyq = eegFS/2;
Fstop1 = (lcut - attenHz)  / nyq;  
Fpass1 = lcut  / nyq; 
Fpass2 = hcut / nyq;
Fstop2 = (hcut + attenHz) / nyq;
Astop1 = attendB;
Apass  = 1;
Astop2 = attendB;
h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
                Fpass2, Fstop2, Astop1, Apass, Astop2);
Hd = design(h, 'kaiserwin');
bGamma = Hd.Numerator;
[a,f] = grpdelay(bGamma,1,eegFS/2,eegFS);
gdGamma = floor(mean(a(f>lcut & f < hcut))); %group delay

%filter theta and gamma
eegTheta = filter(bTheta,1,eeg);
eegTheta = [eegTheta(gdTheta+1:end) zeros(1,gdTheta)];

eegGamma = filter(bGamma,1,eeg);
eegGamma = [eegGamma(gdGamma+1:end) zeros(1,gdGamma)];

%get phase and amplitude
phase = angle(hilbert(eegTheta));
amp = abs(hilbert(eegGamma));

%MI - modulation index
%MI_NORM - normalized modulation index (z-score)
[MI, MInorm] = getmisur(amp,phase,edges,shifts);

disp(MI)
disp(MInorm)

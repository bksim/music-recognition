%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recognizes chords and notes given a wav file.
% Brandon Sim, Gregory Han
% 6/3/2013
%
% Code for plotting from:
% Mark R. Petersen, U. of Colorado Boulder Applied Math Dept
% (uncomment to show plots)
% 
% Change string in file variable to analyze a different file.
% Disregards octaves at the moment (will be fixed in future version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

file = 'music/chord_C5E5G5.wav';
[y, Fs] = wavread(file);
t = (1:length(y))/Fs;

ind = find(t>0.1 & t<0.35);   % set time duration for waveform plot

% plots waveform
%figure; subplot(1,2,1)
%plot(t(ind),y(ind))  
%axis tight         
%title(['Waveform of ' file])

N = length(ind);
c = fft(y(ind))/N; % compute FFT
p = 2*abs( c(2:N/2));         % compute power at each frequency
f = (1:N/2-1)*Fs/N;           % frequency corresponding to p

% finds peaks with threshold of 0.085
[pks, loc] = findpeaks(p, 'MINPEAKHEIGHT', 0.085);

% detects chords (note: will always get rid of octaves)
% for each pair of peaks, removes the higher frequency peak if it is
% 2 times a lower frequency peak (these are higher harmonics)
nonchordPeak = [];
lengthloc = length(loc);
for i=1:lengthloc
   for j=i+1:lengthloc
       % is a multiple, get rid of higher one
       if abs(2-(f(loc(j))/f(loc(i)))) < 0.1
           nonchordPeak = [nonchordPeak j];
       end
   end
end

loc(nonchordPeak) = []; % gets rid of peaks that are in nonchordPeak
key_numbers = round(12*log2(f(loc)/440)+49); % computes the key numbers
key_names = mod(key_numbers, 12)+1; % computers the key name
keys = {'Aflat', 'A', 'Bflat', 'B', 'C', 'Dflat', 'D', 'Eflat', ...
    'E', 'F', 'Gflat', 'G'};

fprintf('Detected chord: ');
for n = 1:length(key_names)
    % prints out key name and octave number
    fprintf('%s%d ', keys{key_names(n)}, floor(key_numbers(n)/12)+1);
end
fprintf('\n');

% plots power spectrum
%subplot(1,2,2)
%semilogy(f,p)
%axis([0 4000 10^-4 1])                
%title(['Power Spectrum of ' file])
clear all; clf; clc;

%% Filter Suara Dari Potongan Lagu Menggunakan Butterworth Filter

[y,Fs]= audioread('cewewav.wav');  % Input Sample Sound from a song

out_dir = uigetdir(cd,'C:\Users\User\Downloads\dsp');
if out_dir == 0
    return;
end

if size(y,2) == 1
    msgbox('The selected file is Mono. This algorithm is applicable only for Stereo files.');
    return;
end

fc=input('Enter Cutoff Frequency (HPF):');
fc=round(fc);
if fc > 20
    fp = fc+5;
    fs = fc/(Fs/2);
    fp = fp/(Fs/2);
    [n wn] = buttord(fp,fs,0.5,80);
    [b, a] = butter(5,wn,'High');
    channel_2 = filtfilt(b,a,y(:,2));
else
    channel_2 = y(:,2);
end

karaoke_wav = channel_2;

% Write it to a file
file = "ceweclean.wav";  % Write it to new file without the music
[filepath,name,ext] = fileparts(file);

if isfolder(out_dir)
    audiowrite(file,karaoke_wav,Fs);
    clear y Fs
end

%% Karakterisasi Suara Sample Dengan Suara Dari Lagu Untuk Menemukan Pitch Menggunakan FFT

% Karakter Suara Sampel
file = "hanswav.wav";  % Input sample sound
[y,Fs] = audioread(file);
y=normalisasi(y);
subplot(211)
[p1,g1] = karakter_suara(y,Fs,file);

% Karakter Suara Target
file = "ceweclean.wav";
[y,Fs] = audioread(file);
y=normalisasi(y);
subplot(212)
[p2,g2] = karakter_suara(y,Fs,file);

avg1 = mean(p1); avg2 = mean(p2); r_avg = avg1/avg2
max1 = max(p1); max2 = max(p2); r_max =  max1/max2
min1 = min(p1); min2 = min(p2); r_min = min1/min2
r_ultimate = (r_avg+r_min+r_max)/3
gain = g1/g2

%% Sintesis Suara Vocal Lagu Dengan Suara Sample

[y11,Fs] = audioread("C:\Users\User\Downloads\dsp\LoveBeautiful.wav");  % Input the whole song

[y2,Fs] = audioread("C:\Users\User\Downloads\dsp\hanswav.wav"); % Input the sample sound

out_dir = uigetdir(cd,'C:\Users\User\Downloads\dsp'); % Replace with desire output folder
if out_dir == 0
    return;
end

if size(y11,2) == 1
    msgbox('The selected file is Mono. This algorithm is applicable only for Stereo files.');
    return;
end

r =  r_ultimate; 
Y = y11;
Fs2 = round(Fs/r);
volume = 2; % increase volume of sound 2 times than the given volume.
y1 = stretchAudio(volume*Y,r);

%Write it to a file 
file = "HANSSSS.wav";  % Rename the output song file
[filepath,name,ext] = fileparts(file);

if isfolder(out_dir)
    audiowrite(file,y1,Fs2);
    clear y Fs
end

%% Function

function [pitch,gain] = karakter_suara(y,Fs,file)
% FILTER NOISE <20Hz dan >18000Hz 
T = 1/Fs;       % Sampling period       
L = length(y);  % Length of signal
t = (0:L-1)*T;  % Time vector
N = 3;            
band = [100 300];
[b,a]= butter(N,band/(Fs/2),'bandpass');
y = filter(b,a,y);

ft = fft(y);
P2 = abs(ft/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

plot(f,P1)
hold on
[maxP1,k] = maxk(P1,100);
pitch = f(k);
gain = max(P1);
stem(f(k),maxP1)
%sound(y,Fs)
xlim([80 500])
title(sprintf('Vocal Characteristic of "%s"',file))
legend("Signal FFT","Dominant Vocal Pitch")
xlabel("frequency (Hz)")
ylabel("power")
grid on
end

function norm = normalisasi(y)
ymin = min(y);
ymax = max(y);
    if(ymax > abs(ymin))
        norm = y/ymax;
    else
        norm = y/abs(ymin);
    end
end
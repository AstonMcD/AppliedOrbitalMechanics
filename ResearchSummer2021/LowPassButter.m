%digital low pass butterworth filter

    %parameters: ->3 Hz Bandwidth
    %            ->range 1 micro-Hz to 100 Hz

%Objectives: ->plot amplitude and phase as function of range
%            ->locate frequencies where gain goes to zero
%            ->Calculate peak gain at all side lobes past the 3 Hz
%            Bandwidth

clc
clear all
close all
%find order of filter and then TF
 %Limit:
Rp=1e-6; % pass band position
Rs=20;% stop band attenuation
fp=3; %pass band frequency in Hz
fs=1; % sampling rate (Hz)
f=100; %sampling frequency
Wp=(2*pi*fp)/f;%rad/s fp
Ws= (2*pi*fs)/f; %rad/s fs

[N,wn] = buttord(Wp,Ws,Rp,Rs); %find order of filter
[B,A]=butter(N,wn,'low');
t=0:0.01:pi;
[h ohm]=freqz(B,A,t);
subplot(3,1,1)
plot(ohm/pi,abs(h))
grid on;
xlabel('normalized frequency')
ylabel('gain')
title('frequency response')
subplot(3,1,2)
plot(ohm/pi,20*log(abs(h)))
xlabel('normalized frequency')
ylabel('gain in dB')
title('frequency response in dB')
subplot(3,1,3)
plot(ohm/pi,angle(h))
xlabel('normalized frequency')
ylabel('phase')
title('phase response')

for i = 1:size(ohm)
    if (ohm(i)/pi)>Wp
        if angle(h(i))<angle(h(i+1))
            fprintf('peak gain is %d dB',20*log(abs(h(i))))
        end
    end
end
        
    


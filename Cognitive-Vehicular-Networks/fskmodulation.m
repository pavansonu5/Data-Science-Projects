%da=[23,12,13,80];
%x=randi([0,3],10,1);
%da1=reshape(da,[4,1]);
%data=de2bi(da);
%M=128;
%freqsep=4;
%nsamp=8;
%Fs=1000;
%y = fskmod(da1,M,freqsep,nsamp,Fs);
%z=fskdemod(y,M,freqsep,nsamp,Fs);
function [fskmodn]=fskmodulation(ACSS,EI,dist,time,nvar,pustat,speed)
y=[ACSS,EI,dist,time,nvar,pustat,speed];% Binary Information
z=de2bi(y);
x=reshape(z,1,prod(size(z)));
bp=.000001;                                                    % bit period
%disp(' Binary information at Trans mitter :');
%disp(x);
%XX representation of transmitting binary information as digital signal XXX
bit=[]; 
for n=1:1:length(x)
    if x(n)==1;
       se=ones(1,100);
    else x(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t1=bp/100:bp/100:100*length(x)*(bp/100);
%subplot(3,1,1);
%plot(t1,bit,'lineWidth',2.5);grid on;
%axis([ 0 bp*length(x) -.5 1.5]);
%ylabel('amplitude(volt)');
%xlabel(' time(sec)');
%title('transmitting information as digital signal');
%XXXXXXXXXXXXXXXXXXXXXXX Binary-FSK modulation XXXXXXXXXXXXXXXXXXXXXXXXXXX%
A=5;                                          % Amplitude of carrier signal
br=1/bp;                                                         % bit rate
f1=br*8;                           % carrier frequency for information as 1
f2=br*2;                           % carrier frequency for information as 0
t2=bp/99:bp/99:bp;                 
ss=length(t2);
m=[];
for (i=1:1:length(x))
    if (x(i)==1)
        y=A*cos(2*pi*f1*t2);
    else
        y=A*cos(2*pi*f2*t2);
    end
    m=[m y];
end
t3=bp/99:bp/99:bp*length(x);
%subplot(3,1,2);
%plot(t3,m);
%xlabel('time(sec)');
%ylabel('amplitude(volt)');
%title('waveform for binary FSK modulation coresponding binary information');
fskmodn=m;
end




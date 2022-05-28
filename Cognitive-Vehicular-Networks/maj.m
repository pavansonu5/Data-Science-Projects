%syms x
%fun_c=(1/(sqrt(2*pi)))*(50/99)*((50^2)/2)*exp(-50*2);
%fun=exp(-(x^2)/2);
%val=int(fun,x);
%%
clear;
clc;
clf;
close all;

%%
%road bulding
%highwayL=[1000,150];             % road [length,width] in m
%vehicleSpeed=[30,100]*1000/3600;% min and max vehicle speed in m/s

% Load road in image format
I = imread('road.png'); 

% Set RSU & PU coordiantes
RSU_c=[500 500]; 
PU_c = [250 200]; 
%road building
%rectangle('Position',[0,0,highwayL(1),40],'FaceColor','black','EdgeColor','white')
%rectangle('Position',[0,50,highwayL(1),40],'FaceColor','black','EdgeColor','white')
%rectangle('Position',[0,90,highwayL(1),110],'FaceColor','g','EdgeColor','g')
h = image([0 2000],[0 2000],I); 
hold on; scatter(PU_c(1),PU_c(2),'^','filled');
scatter(RSU_c(1),RSU_c(2),'s','filled');
text(PU_c(1)+50,PU_c(2),'PU Base Station');
text(RSU_c(1)+50,RSU_c(2),'RSU');
%%
%pu status
pu=randi([0,1],1,10);
wt=randi([5,15],1,10);
F1=1000;    Fs=5000;% Carrier & Sampling freq
%%
% TV Tx. s/l generation
t=0:1/Fs:0.02;
x=5*sin(2*pi*100*t); 
A=4*ammod(x,F1,Fs);
%%
% Filtere for receiver at SU
ww1=fdesign.bandpass('N,F3dB1,F3dB2',2,F1-500,F1+500,Fs);
Hd1=design(ww1,'butter');
strtn=5; endn=50; 
f=1;%pd,pfa,pmd
%%
%numero=1;
% For varying SUs number
for N=strtn:5:endn
    SUx=randi([400 1600],1,N);
    % Y axis position of SUs to align them on road (based on image)
    for k=1:N    
        if SUx(k)<820       
            SUy(k)=randi([620 1180],1,1);    
        elseif SUx(k)<1320        
            SUy(k)=randi([300 1700 ],1,1);    
        else
            SUy(k)=randi([620 1180],1,1);    
        end
    end
   
    % Assign speeds for nodes (in m/s)
    speed=randi([8,16],1,N);
    %totaltime=sum(wt);
    %MAJ=[zeros(1,totaltime)];
    title(sprintf('Number of Secondary Users : %g',N));
    taggedSU=randi([1,N],1,1);
    st=randi([1,10],1,N);
    time=0;
    cmp=[];
    res=[];
    %spread sheet
    delete('data6.xlsx')
    col={'ACSS','ener_cal','dis_meas','Time','weight_det'};
    xlswrite('data6.xlsx',col,'Sheet1','A1')
    for i=1:10
        for j=1:wt(i)
            time=time+1;
            %sensing and reporting phase
            a=scatter(SUx,SUy,'b','filled');   
            
            for k=1:N
            if (rem(time,st(k))==0)
            dist=sqrt(((SUx(k)-PU_c(1))^2+(SUy(k)-PU_c(2))^2));
            pathloss= 7000*dist^(-3/2);
            PU=pathloss*pu(i)*A;
            snrv=-0.3*speed(k);   
            %snr=-10;
            if(pu(i)==1)
            PU1=awgn(PU,snrv,'measured');
            else
            PU1=awgn(PU,snrv);
            end   
            z=sum(abs(filter(Hd1,PU1)).^2); 
            if(k==taggedSU)
                ACSS=1;
                cmp=[cmp pu(i)];
            else
                ACSS=0;
            end
            %sending information to RSU`
            toRSU=fskmodulation(ACSS,ceil(z),ceil(dist),ceil(time));
            %numero=numero+1;
            if(ACSS==1)
            fromRSU=roadsideunit(toRSU); 
            res=[res fromRSU];
            else
                roadsideunit(toRSU);
            end
            
            %nvar = 1/snr;   
            %nSample=length(PU1);    
            % Setting threshold
            %AT RSU
           
            end
            end
            for z=1:N
                if(rem(time,10)==0)
            if SUy(z)<920&&SUy(z)>680;
                SUx(z)=SUx(z)+speed(z);       
            elseif SUy(z)>940&&SUy(z)<1180;           
                SUx(z)=SUx(z)-speed(z);        
            elseif SUx(z)>940&&SUy(z)<920;
                SUy(z)=SUy(z)+speed(z);        
            else
                SUy(z)=SUy(z)-speed(z);
            end
                end
            end
            pause(0.01);
    set(a,'Visible','off'); 
        end
    end
pd(N)=(sum(cmp==res)/length(cmp));
pmd(N)=sum(cmp==1&res==0)/length(cmp);
pfa(N)=sum(cmp==0&res==1)/length(cmp); 
f=f+1; 
end
N=strtn:5:endn;
figure(2);
subplot 311; 
plot(N,pd(strtn:5:endn));
title 'PROBABILITY OF DETECTION';  
ylim([-0.2 1.2]); grid; 
xlabel 'Number of Secondary users';
ylabel 'P_{d}'; legend('MAJ rule');
subplot 312; 
plot(N,pmd(strtn:5:endn)); 
title 'PROBABILITY OF MISDETECTION'; ylim([-0.2 1.2]); grid; 
xlabel 'Number of Secondary users'; ylabel 'P_{md}'; 
legend('MAJ rule'); 
subplot 313; 
plot(N,pfa(strtn:5:endn));
title 'PROBABILITY OF FALSE ALARM'; ylim([-0.2 1.2]); grid;
xlabel 'Number of Secondary users'; ylabel 'P_{fa}'; legend('MAJ rule');
    
                    
            
    


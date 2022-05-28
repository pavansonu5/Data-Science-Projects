%%
clear;
clc;
clf;
close all;
%%
%road building
highwayL=[1000,40];
PU_c= [400 35];
RSU_c= [500 12.5];
axis([0 highwayL(1) 0 highwayL(2)]);
hold on;
rectangle('Position',[0,0,highwayL(1),10],'FaceColor','black','EdgeColor','white');
rectangle('Position',[0,15,highwayL(1),10],'FaceColor','black','EdgeColor','white');
scatter(PU_c(1)+10,PU_c(2),'^','filled');
scatter(RSU_c(1)+10,RSU_c(2),'s','filled');
text(PU_c(1),PU_c(2),'PU Base Station');
text(RSU_c(1),RSU_c(2),'RSU');
%%
%pu status
pu=randi([0,1],1,10);
wt=randi([5,15],1,10);
F1=500*(10^6);    Fs=1200*(10^6);% Carrier & Sampling freq
%%
% TV Tx. s/l generation
%N=1000;
t=0:0.00001:0.01;
x=4*sin(2*pi*100*t); 
A=3*ammod(x,F1,Fs);
%%
% Filtere for receiver at SU
ww1=fdesign.bandpass('N,F3dB1,F3dB2',8,490*(10^6),510*(10^6),Fs);
Hd1=design(ww1,'butter');
pfa=0:0.1:1;%pfa values
%N=4;
%%
%numero=1;
% For varying SUs number
for N=2:2:4
for rg=1:length(pfa)
    SUy=randi([1 24],1,N);
    % Y axis position of SUs to align them on road (based on image)
    for k=1:N    
        if SUy(k)<10       
            SUx(k)=randi([0 150],1,1);    
        elseif SUy(k)>15        
            SUx(k)=randi([850 1000],1,1);    
        else
            SUy(k)=randi([1 9],1,1);
            SUx(k)=randi([0 150],1,1);    
        end
    end
   
    % Assign speeds for nodes (in m/s)
    speed=randi([8,16],1,N);
    %totaltime=sum(wt);
    %MAJ=[zeros(1,totaltime)];
    title(sprintf('Number of Secondary Users : %g',N));
    taggedSU=randi([1,N],1,1);
   sf=zeros([1,N]);
   st=randi([2,9],1,N);
    %st(taggedSU)=7;
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
            if (rem(time,((sf(k)*10)+st(k)))==0)
                sf(k)=sf(k)+1;
            dist=sqrt(((SUx(k)-PU_c(1))^2+(SUy(k)-PU_c(2))^2));
            pathloss= 3*dist^(-0.5);
            PU=pathloss*pu(i)*A;
            %snr=250/dist;
            snrv=-5;
            %snrv=-0.6*speed(k);
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
            %AT RSU
      data=fskdemodulation(toRSU);
      y=[];
       %Q=vertcat(Q,data);
       if (data(1,1)==1)
      writematrix(data,'data6.xlsx','Sheet','sheet1','WriteMode','Append');
      A1=readmatrix('data6.xlsx');
      B=sortrows(A1,4,'descend');
      sum1=0;
      weightsum=0;
      %if((prod(size(B))/5)>10)
      %    qq=10;
      %else
      %   qq=prod(size(B))/5;
      %end
      check=[];
      qq=1;
      while(qq<=size(B,1) & (data(1,4)-B(qq,4))<=10 & qq<=10)
         fis=readfis('weight_dete.fis');
        weight(qq)=evalfis(fis,[B(qq,3) (data(1,4)-B(qq,4))]);
        check=[check [B(qq,3) (data(1,4)-B(qq,4))]];
        y=[y weight(qq)];
        %sprintf('Weight : %g',weight(i));
        sum1=sum1+weight(qq)*B(qq,2);
        weightsum=weightsum+weight(qq);
        qq=qq+1;
      end
      sum1=sum1/(weightsum);
     %pfa=0.1;
     Nsam=1000;
     y=y./(sum(y));
     SNR=0.1; %-10db
     func1=sqrt(2*Nsam*(y(1)^2+sum(y(2:end).^2)*(SNR+1)^2));
     func2=Nsam*(y(1)+sum(y(2:end))*(SNR+1));
     thresh=(qfuncinv(pfa(rg))*func1+func2);
     writematrix(check,'data6.xlsx','Sheet','sheet3','WriteMode','Append'); 
     writematrix([data(2),sum1,thresh],'data6.xlsx','Sheet','sheet2','WriteMode','Append'); 
    if(sum1>thresh)
        res=[res 1];
    else
        res=[res 0];
    end
     
 else
 %numero=numero+1;
 %xlswrite('data6.xlsx',data,'section II',['A' num2str(numero)]);
 writematrix(data,'data6.xlsx','Sheet','sheet1','WriteMode','Append');
 
 end
            
            %nvar = 1/snr;   
            %nSample=length(PU1);    
            % Setting threshold
            %AT RSU
           
            end
            end
            for z=1:N
            if SUy(z)<=10&&SUy(z)>=0;
                SUx(z)=SUx(z)+speed(z)/4;       
            elseif SUy(z)>=15&&SUy(z)<=25;           
                SUx(z)=SUx(z)-speed(z)/4;        
            end
            end
            pause(0.01);
    set(a,'Visible','off'); 
        end
    end
res=res
cmp=cmp
pd(rg+((N/2)-1)*length(pfa))=(sum(res==1 &cmp==1)/sum(cmp==1));
%pmd(rg)=sum(cmp==1&res==0)/length(cmp);
%pfap(rg)=sum(res==0)/length(cmp); 
end
end

%N=strtn:5:endn;
rg=1:length(pfa);
figure(2);
%subplot 311; 
plot(pfa,pd(1:length(pfa)));
title 'PROBABILITY OF DETECTION';  
ylim([-0.2 1.2]); grid; 
xlabel 'pfa';
ylabel 'Qd';
hold on
plot(pfa,pd(length(pfa)+1:end));
legend('N=2','N=4');



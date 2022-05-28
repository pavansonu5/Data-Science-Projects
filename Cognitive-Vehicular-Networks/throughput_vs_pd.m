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
%pu=randi([0,1],1,10);
pu=randsrc(1,20,[1,0;0.3,0.7]);
wt=randi([5,15],1,20);
F1=500*(10^6);    Fs=1200*(10^6);% Carrier & Sampling freq
%%
% TV Tx. s/l generation
%N=1000;
t=0:0.001/1000:0.001;
x=4*sin(2*pi*1000*t); 
A=2*ammod(x,F1,Fs);
%%
% Filtere for receiver at SU
ww1=fdesign.bandpass('N,F3dB1,F3dB2',8,490*(10^6),510*(10^6),Fs);
Hd1=design(ww1,'butter');
pd=0.1:0.1:1;%pfa values
N=5;
tr=8;
tt=12;
transmit=0;
taggedSU=randi([1,N],1,1);
sf=zeros([1,N]);
st=randi([2,9],1,N);
pnv=[];
pmv=[];
SU1y=randi([1 24],1,N);
threshval=[];
eneval=[];
SCSS=0;
nvarval=[];
%snrv=0;
%no=awgn(0*A,snrv);
sigvar=[];
ene2val=[];
thresh2val=[];
for k=1:N    
        if SU1y(k)<10       
            SU1x(k)=randi([0 700],1,1);    
        elseif SU1y(k)>15        
            SU1x(k)=randi([300 1000],1,1);    
        else
            SU1y(k)=randi([1 9],1,1);
            SU1x(k)=randi([0 700],1,1);    
        end
end
    speed=randi([8,20],1,N);
%%
%numero=1;
% For varying SUs number
for rg=1:length(pd)
    % Y axis position of SUs to align them on road (based on image)
    SUx=SU1x;
    SUy=SU1y;
   
    % Assign speeds for nodes (in m/s)
    %totaltime=sum(wt);
    %MAJ=[zeros(1,totaltime)];
    title(sprintf('Number of Secondary Users : %g',N));
    %taggedSU=randi([1,N],1,1);
    sf=zeros([1,N]);
    %st=randi([1,10],1,N);
    %st(taggedSU)=7;
    time=0;
    cmp=[];
    res=[];
    res2=[];
    res3=[];
    %spread sheet
    delete('data6.xlsx')
    col={'ACSS','ener_cal','dis_meas','Time'};
    xlswrite('data6.xlsx',col,'Sheet1','A1')
    pn=0;
    pm=0;
    for i=1:20
        for j=1:wt(i)
            time=time+1;
            %sensing and reporting phase
            a=scatter(SUx,SUy,'b','filled');   
            if(transmit==1)
                if(pu(i)==0)
                    pn=pn+1;
                else
                    pm=pm+1;
                end
                transmit=0;
            end
            for k=1:N
            if (rem(time,((sf(k)*10)+st(k)))==0)
            sf(k)=sf(k)+1;    
            dist=sqrt(((SUx(k)-PU_c(1))^2+(SUy(k)-PU_c(2))^2));
            pathloss=9*dist^(-0.7);
            PU=pathloss*pu(i)*A;
            %snr=250/dist;
            snrv=-0.5*speed(k);
            %snrv=-0.5*speed(k);
            %if(pu(i)==1)
            %PU1=awgn(PU,snrv,'measured');
            %else
            %PU1=awgn(PU,snrv);
            %end
            no=awgn(0*A,snrv);
            if(pu(i)==1)
            sigvar=[sigvar var(PU)];
            end
            PU1=PU+no;
            nvar=var(filter(Hd1,no));
            en=sum(abs(filter(Hd1,PU1)).^2);
            eneval=[eneval [en pu(i)]];
            %nvar=var(filter(Hd1,no));
            
            if(k==taggedSU & dist<1000)
                ACSS=1;
                SCSS=1;
                cmp=[cmp pu(i)];
            else
                ACSS=0;
            end
            %sending information to RSU`
            toRSU=fskmodulation(ACSS,ceil(en),ceil(dist),ceil(time),ceil(1000*nvar),pu(i),speed(k));
            %numero=numero+1;
            %AT RSU
      data=fskdemodulation(toRSU);
       %Q=vertcat(Q,data);
       if (data(1,1)==1)
           y=[];
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
      bpan=[];
      while(qq<=size(B,1) & (data(1,4)-B(qq,4))<10 & qq<=7 & B(qq,3)<1000)
         fis=readfis('checking.fis');
        weight(qq)=evalfis(fis,[B(qq,3) (data(1,4)-B(qq,4)) B(qq,end)*(18/5)]);
        check=[check [B(qq,3) (data(1,4)-B(qq,4))]];
        y=[y weight(qq)];
        %sprintf('Weight : %g',weight(i));
        sum1=sum1+weight(qq)*B(qq,2);
        weightsum=weightsum+weight(qq);
        bpan=[B(qq,2) B(qq,3)];
        writematrix(bpan,'m.xlsx','Sheet','sheet5','WriteMode','Append');
        qq=qq+1;
      end
      sum1=sum1/(weightsum);
     %pfa=0.1;
     Nsam=1001;
     y=y./(sum(y));
     SNR=0.1; %-10db
     func1=sqrt(2*Nsam*(((y(1)^2)*(SNR+1)^2)+(sum(y(2:end).^2)*(SNR+1)^2)));
     func2=Nsam*(y(1)*(SNR+1)+sum(y(2:end))*(SNR+1));
     thresh=0.3*(qfuncinv(pd(rg))*func1+func2);
     threshval=[threshval thresh];
     writematrix(check,'data6.xlsx','Sheet','sheet3','WriteMode','Append'); 
     writematrix([data(2),sum1,thresh],'data6.xlsx','Sheet','sheet2','WriteMode','Append'); 
    if(sum1>thresh)
        res=[res 1];
    else
        transmit=1;
        res=[res 0];
    end
    sum5=0;
    bpan3=readmatrix('m.xlsx');
    bpan2=sortrows(bpan3,2,'descend');
    werand=rand(1,size(bpan2,1));
    werand2=sort(werand,'descend');
    thresh=0.3*(qfuncinv(pd(rg))*func1+func2);
    for qqqq=1:size(bpan2,1)
        sum5=sum5+bpan2(qqqq,1)*werand2(qqqq);
    end
    sum5=sum5./sum(werand);
    if(sum5>thresh)
        res3=[res3 1];
    else
        transmit=1;
        res3=[res3 0];
    end
     delete('m.xlsx');
 else
 %numero=numero+1;
 %xlswrite('data6.xlsx',data,'section II',['A' num2str(numero)]);
 writematrix(data,'data6.xlsx','Sheet','sheet1','WriteMode','Append');
       end
 %delete('m.xlsx');
            
            %nvar = 1/snr;   
            %nSample=length(PU1);    
            % Setting threshold
            %AT RSU
        
           
            end
            end
            if(SCSS==1)
                for k=1:N
                dist=sqrt(((SUx(k)-PU_c(1))^2+(SUy(k)-PU_c(2))^2));
            pathloss=9*dist^(-0.7);
            PU=pathloss*pu(i)*A;
            snrv=-0.5*speed(k);
            no=awgn(0*A,snrv);
            PU1=PU+no;
            en2=sum(abs(filter(Hd1,PU1)).^2);
            en2=en2./Nsam;
            ene2val=[ene2val en2];
            thresh2=0.3*(qfuncinv(pd(rg))*sqrt(2*(SNR+1)./Nsam)+(SNR+1));
            thresh2val=[thresh2val thresh2];
            if en2>thresh2 
                pk(k)=1;
            else
                pk(k)=0;
            end
                end
            
            if sum(pk)>2
                res2=[res2 1];
            else
                res2=[res2 0];
            end
            SCSS=0;
            end
            for z=1:N
            if SUy(z)<=10&&SUy(z)>=0;
                SUx(z)=SUx(z)+speed(z)/7;       
            elseif SUy(z)>=15&&SUy(z)<=25;           
                SUx(z)=SUx(z)-speed(z)/7;        
            end
            end
            pause(0.01);
    set(a,'Visible','off'); 
        end
    end
pnv=[pnv pn];
pmv=[pmv pm];
p(rg)=sum(cmp==0)/length(cmp); 
tot=pn+pm;
if(tot==0)
    tot=1;
end
pfac(rg)=sum(cmp==0&res==1)/sum(cmp==0); 
pfac2(rg)=sum(cmp==0&res2==1)/sum(cmp==0);
pfac3(rg)=sum(cmp==0&res3==1)/sum(cmp==0);
pdc(rg)=(sum(cmp==1&res==1)/sum(cmp==1));
pdc2(rg)=(sum(cmp==1&res2==1)/sum(cmp==1));
pdc3(rg)=(sum(cmp==1&res3==1)/sum(cmp==1));
pmd(rg)=sum(cmp==1&res==0)/sum(cmp==1);
pmd2(rg)=sum(cmp==1&res2==0)/sum(cmp==1);
pmd3(rg)=sum(cmp==1&res3==0)/sum(cmp==1);
throu1(rg)=((tr)/tt)*p(rg).*((pn/(tot))*(1-pfac(rg))*log2(1+5)+(pm/(tot))*(1-pdc(rg))*log2(1+5/10));  
throu2(rg)=((tr-2)/tt)*p(rg).*((pn/(tot))*(1-pfac2(rg))*log2(1+5)+(pm/(tot))*(1-pdc2(rg))*log2(1+5/10));  
throu3(rg)=((tr)/tt)*p(rg).*((pn/(tot))*(1-pfac3(rg))*log2(1+5)+(pm/(tot))*(1-pdc3(rg))*log2(1+5/10)); 
%pdp(rg)=(sum(res==0)/length(cmp));
%pmd(rg)=sum(cmp==1&res==0)/length(cmp);
%pfap(rg)=sum(res==0)/length(cmp); 
end
%N=strtn:5:endn;
rg=1:length(pd);
figure(2);
%subplot 311; 
plot(pd,throu1(1:length(pd)));
hold on
plot(pd,throu2(1:length(pd)));
hold on
plot(pd,throu3(1:length(pd)));
legend('ACSS-pm','SCSS','ACSS-bp','location','south');
title 'THROUGHPUT Analysis';  
%ylim([0 1]); 
grid; 
xlabel 'Pd';
ylabel 'throughput'; %legend('MAJ rule');
%subplot 312; 
%plot(pfa,pmd(1:length(pfa))); 
%title 'PROBABILITY OF MISDETECTION'; ylim([-0.2 1.2]); grid; 
%xlabel 'Number of Secondary users'; ylabel 'P_{md}'; 
%legend('MAJ rule'); 
%subplot 313; 
%plot(pfa,pfap(1:length(pfa)));
%title 'PROBABILITY OF FALSE ALARM'; ylim([-0.2 1.2]); grid;
%xlabel 'Number of Secondary users'; ylabel 'P_{fa}'; legend('MAJ rule');
    


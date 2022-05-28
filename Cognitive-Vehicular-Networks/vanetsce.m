%%

%load('svmcode.mat')
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
pu=randi([0,1],1,20);
wt=randi([5,15],1,20);
F1=500*(10^6);    Fs=1200*(10^6);% Carrier & Sampling freq
%%
% TV Tx. s/l generation
%N=1000;
t=0:0.001/500:0.001;
x=4*sin(2*pi*1000*t); 
A=2*ammod(x,F1,Fs);
%snrv=-2;
%no=awgn(0*A,snrv);
enval=[];
%%
% Filtere for receiver at SU
band1=fdesign.bandpass('N,F3dB1,F3dB2',8,490*(10^6),510*(10^6),Fs);
bpf=design(band1,'butter');
%nvar=var(filter(Hd1,no));
strtn=5; endn=50; 
trans=1;
delete('xyz.xlsx');
col={'ener_cal','dis_meas','Time','PU_STATUS'};
xlswrite('xyz.xlsx',col,'Sheet4','A1')
col2={'ener_cal','dis_meas','Time','PU_STATUS'};
xlswrite('svmtest.xlsx',col2,'Sheet1','A1')
%%
%numero=1;
% For varying SUs number
for N=strtn:5:endn
    SUxy=[];
    SUy=randi([1 24],1,N);
    % Y axis position of SUs to align them on road (based on image)
    for k=1:N    
        if SUy(k)<10       
            SUxy(k)=randi([0 600],1,1);    
        elseif SUy(k)>15        
            SUxy(k)=randi([300 1000],1,1);    
        else
            SUy(k)=randi([1 9],1,1);
            SUxy(k)=randi([0 600],1,1);    
        end
    end
   
    % Assign speeds for nodes (in m/s)
    speed=randi([8,20],1,N);
    %totaltime=sum(wt);
    %MAJ=[zeros(1,totaltime)];
    title(sprintf('Number of Secondary Users : %g',N));
    taggedSU=randi([1,N],1,1);
    sf=zeros([1,N]);
    st=randi([1,10],1,N);
    st(taggedSU)=8;
    time=0;
    cmp=[];
    res=[];
    res2=[];
    %spread sheet
    delete('data6.xlsx')
    col={'ACSS','ener_cal','dis_meas','Time','weight_det'};
    xlswrite('data6.xlsx',col,'Sheet1','A1')
    
    for i=1:20
        for j=1:wt(i)
            time=time+1;
            %sensing and reporting phase
            a=scatter(SUxy,SUy,'b','filled');   
            
            for k=1:N
            if (rem(time,((sf(k)*10)+st(k)))==0)
            sf(k)=sf(k)+1;
            dist=sqrt(((SUxy(k)-PU_c(1))^2+(SUy(k)-PU_c(2))^2));
            pathloss=9*dist^(-0.7);
            PU=pathloss*pu(i)*A;
            %snr=250/dist;
            snrv=-0.5*speed(k);
            no=awgn(0*A,snrv);
            nvar=var(filter(bpf,no));
            %if(pu(i)==1)
            %PU1=awgn(PU,snrv,'measured');
            %else
            %PU1=awgn(PU,snrv);
            %end
            PU1=PU+no;
            en=sum(abs(filter(bpf,PU1)).^2); 
            enval=[enval en];
            if(dist<1000)
            if(k==taggedSU)
                ACSS=1;
                cmp=[cmp pu(i)];
            else
                ACSS=0;
            end
            %sending information to RSU`
            toRSU=fskmodulation(ACSS,ceil(en*10),ceil(dist),ceil(time),ceil(nvar*1000),pu(i),speed(k));
            %numero=numero+1;
            %xyz=[ceil(z),ceil(dist),ceil(time),pu(i)];
            %writematrix(xyz,'data7.xlsx','Sheet','sheet1','WriteMode','Append');
            if(ACSS==1)
            fromRSU=roadsideunit(toRSU); 
            res=[res fromRSU];
            data2=fskdemodulation(toRSU);
           pustat2=data2(end-1);
           nvar2=data2(end-2);
           data3=[data2(2)/10 data2(3) data2(end)];
           %writematrix(data,'data6.xlsx','Sheet','sheet1','WriteMode','Append');
             A2=readmatrix('data6.xlsx');
            B=sortrows(A2,4,'descend');
            %load('svmcode.mat');
             qqq=1;
          while(qqq<=size(B,1) & qqq<=6 & (data2(4)-B(qqq,4))<15)
          %writematrix([B(qqq,2),B(qqq,3),B(qqq,end)],'svmtest.xlsx');
      %qqq=qqq+1;
   %  pk(qqq)=trainedModel.predictFcn(data);
          %testdata=array2table(data,'VariableNames',{'ener_cal','dis_meas','Time'});
            %load('svmcode.mat');
           %pk(qqq)=trainedModel.predictFcn(testdata);
           %pk(qqq)=emora(data3);
           testdata=array2table(data3,'VariableNames',{'ener_cal','dis_meas','Time'});
           pk(qqq)=trainedModel.predictFcn(testdata);
           qqq=qqq+1;
         end
 %testdata=readmatrix('svmtest.xlsx');
 %yyy=trainedModel.predictFcn(cell2mat(testdata(2:end,:)))
             if(sum(pk)>(length(pk)/2))
             fromRSU2=1;
             else
                fromRSU2=0;
              end
            res2=[res2 fromRSU2];
            %writematrix([pu(i)],'data6.xlsx','Sheet','sheet4','WriteMode','Append');
            else
                roadsideunit(toRSU);
                %writematrix([pu(i)],'data6.xlsx','Sheet','sheet4','WriteMode','Append');
            end
            end
            
            %nvar = 1/snr;   
            %nSample=length(PU1);    
            % Setting threshold
            %AT RSU
           
            end
            end
            for z=1:N
            if SUy(z)<=10&&SUy(z)>=0;
                SUxy(z)=SUxy(z)+speed(z)/4;       
            elseif SUy(z)>=15&&SUy(z)<=25;           
                SUxy(z)=SUxy(z)-speed(z)/4;        
            end
            end
            pause(0.01);
    set(a,'Visible','off'); 
    
        end
    end
pd(N)=(sum(cmp==1& res==1)/sum(cmp==1));
pmd(N)=sum(cmp==1&res==0)/sum(cmp==1);
pfa(N)=sum(cmp==0&res==1)/sum(cmp==0);
pd2(N)=(sum(cmp==1& res2==1)/sum(cmp==1));
pmd2(N)=sum(cmp==1&res2==0)/sum(cmp==1);
pfa2(N)=sum(cmp==0&res2==1)/sum(cmp==0); 
end
N=strtn:5:endn;
figure(2); 
plot(N,pd(strtn:5:endn));
hold on;
plot(N,pd2(strtn:5:endn));
title 'PROBABILITY OF DETECTION';  
ylim([-0.2 1.2]); grid; 
legend('Fuzzy logic','SVM');
xlabel 'Number of Secondary users';
ylabel 'P_{d}'; %legend('MAJ rule');
figure(3);
plot(N,pmd(strtn:5:endn)); 
hold on;
plot(N,pmd2(strtn:5:endn));
title 'PROBABILITY OF MISDETECTION'; ylim([-0.2 1.2]); grid; 
legend('Fuzzy logic','SVM');
xlabel 'Number of Secondary users'; ylabel 'P_{md}'; 
%legend('MAJ rule'); 
figure(4);
plot(N,pfa(strtn:5:endn));
hold on;
plot(N,pfa2(strtn:5:endn));
title 'PROBABILITY OF FALSE ALARM'; ylim([-0.2 1.2]); grid;
legend('Fuzzy logic','SVM');
xlabel 'Number of Secondary users'; ylabel 'P_{fa}';% legend('MAJ rule');
    


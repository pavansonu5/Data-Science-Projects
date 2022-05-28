function [d]=roadsideunit(toRSU)
%AT RSU
 data2=fskdemodulation(toRSU);
 pustat=data2(end-1);
 nvar=data2(end-2);
 data=[data2(1) data2(2)/10 data2(3:end-3) pustat data2(end)*18/5];
 y=[];
 %Q=vertcat(Q,data);
 if (data(1,1)==1)
     writematrix(data,'data6.xlsx','Sheet','sheet1','WriteMode','Append');
     A=readmatrix('data6.xlsx');
     B=sortrows(A,4,'descend');
     sum1=0;
     weightsum=0;
     %if((prod(size(B))/5)>10)
     %    qq=10;
     %else
     %   qq=prod(size(B))/5;
     %end
     check=[];
     qq=1;
     while(qq<=size(B,1) & (data(1,4)-B(qq,4))<15 & qq<=6)
        writematrix([B(qq,2),B(qq,3),B(qq,end),B(qq,5)],'xyz.xlsx','Sheet','sheet4','WriteMode','Append');
        fis=readfis('checking.fis');
        weight(qq)=evalfis(fis,[B(qq,3) (data(1,4)-B(qq,4)) B(qq,end)]);
        check=[check [B(qq,3) (data(1,4)-B(qq,4))]];
        y=[y weight(qq)];
        %sprintf('Weight : %g',weight(i));
        sum1=sum1+weight(qq)*(B(qq,2));
        weightsum=weightsum+weight(qq);
        qq=qq+1;
     end
     sum1=sum1/(weightsum);
     pfa=0.1;
     N=501;
     y=y./(sum(y));
     SNR=0.1; %-10db
     func1=sqrt(2*N*(y(1)^2+sum(y(2:end).^2)*(SNR+1)^2));
     func2=N*(y(1)+sum(y(2:end))*(SNR+1));
     thresh=0.35*(qfuncinv(pfa)*func1+func2);
     writematrix(check,'data6.xlsx','Sheet','sheet3','WriteMode','Append'); 
     
     writematrix([data(2),sum1,thresh],'data6.xlsx','Sheet','sheet2','WriteMode','Append'); 
    if(sum1>thresh)
        d=1;
    else
        d=0;
    end
     
 else
 %numero=numero+1;
 %xlswrite('data6.xlsx',data,'section II',['A' num2str(numero)]);
 writematrix(data,'data6.xlsx','Sheet','sheet1','WriteMode','Append');
 
 end
 %y=1;
end
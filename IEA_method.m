%%%%% Improved extension assessment

clear;
clc;
%%-------------------------------------------------------------------
%R0: Store all data of the upper and lower limits of the classic domain;
%x: Store data of n evaluation indices of the object to be evaluated;
%n: Amount of evaluation indices;
%m: Amount of evaluation levels;
%a: Lower bound matrix of the classical domain;
%b: Upper bound matrix of the classical domain;
%ap: Lower bound matrix of the sectional domain;
%bp: Upper bound matrix of the sectional domain;
%xx: Nondimensionalized values of evaluation indices;
%ws: Subjective weights of evaluation indices;
%wo: Objective weights of evaluation indices;
%w: Comprehensive weights of evaluation indices;
%k(i,j): Correlation degree between a matter-element to be evaluated and the classical domain;
%kp: Correlation degree between a matter-element to be evaluated and UKDS level;
%j0: Recognized UKDS level;
%jstar: Variable characteristic;

%%-------------------------------------Data input-----------------------------------------%%
x=xlsread('InputData.xlsx');                        
[pp,qq]=size(x);                                                          

R0=[0	   0.042  0.042	  0.104     0.104	 0.254	 0.254	  0.636
    100	   85	  85	  70	    70	     60	     60	      0
    0	   3	  3	      5	        5	     10	     10	      100
    100	   85	  85	  70	    70	     60	     60	      0
    5	   1	  1	      0.3 	    0.3      0.05	 0.05	  0
    0	   5	  5	      10	    10	     25	     25	      45
    100	   85	  85	  70	    70	     60	     60	      0
    0      20     20      40        40       60      60       100
    0	   1	  1	      5	        5	     15	     15	      30];         % Grading standards for evaluation indices;  

ws=[0.205 0.205 0.174 0.096 0.053 0.062 0.094 0.069 0.044];    % Subjective weights determined by using AHP; 

%%---------------------------------Classical domain and sectional domain--------------------------------------%%
[n,mm]=size(R0);        % n represents the amount of evaluation indices;
m=mm/2;                 % m represents amount of evaluation levels;
a=R0(:,1:2:end);        % a represents lower bound matrix of the classical domain;
b=R0(:,2:2:end);        % b represents upper bound matrix of the classical domain;
ap=min(R0')';           % ap represents lower bound vector of the classical domain;
bp=max(R0')';           % bp represents upper bound vector of the classical domain;

%%-------------------------------------Nondimensionalization process------------------------------------------%%
% aa(i,j) is dimensionless result of threshold of the classical domain a(i,j)
% bb(i,j) is dimensionless result of threshold of the classical domain b(i,j)
% xx is dimensionless result of values of evaluation index
% app is dimensionless vector of lower bound of the sectional domain
% bpp is dimensionless vector of upper bound of the sectional domain

for jj=1:qq        
    for i=1:n
        for j=1:m
            if a(i,j)<=b(i,j)
               aa(i,j)=(a(i,j)-ap(i))/(bp(i)-ap(i));   
               bb(i,j)=(b(i,j)-ap(i))/(bp(i)-ap(i));
               if x(i,jj)>bp(i)
                  x(i,jj)=1;      
               else
                  xx(i,jj)=(x(i,jj)-ap(i))/(bp(i)-ap(i));   % Nondimensionalization process of values of evaluation indices
               end
            else
               aa(i,j)=(bp(i)-a(i,j))/(bp(i)-ap(i)); 
               bb(i,j)=(bp(i)-b(i,j))/(bp(i)-ap(i));
               if x(i,jj)>bp(i)
                  xx(i,jj)=0;      
               else
                  xx(i,jj)=(bp(i)-x(i,jj))/(bp(i)-ap(i));
               end
            end   
        end
    end
    app=min(aa')';        
    bpp=max(bb')';         

    %%---------------------------------Determination of objective weights by using correlation analysis--------------------------------%%
    for i=1:n
        for j=1:m
            if xx(i,jj)<=0.5*(aa(i,j)+bb(i,j))
               r(i,j)=2*(xx(i,jj)-aa(i,j))/(bb(i,j)-aa(i,j));
            else
               r(i,j)=2*(bb(i,j)-xx(i,jj))/(bb(i,j)-aa(i,j));
            end
        end
    end
    i=1:n;
    j=1:m;
    [rmax,jmax]=max(r,[],2);
    for i=1:n
        if aa(i,j)<=bb(i,j)
           if rmax(i)>=-0.5
              rr(i)=jmax(i)*(1+rmax(i));
           else
              rr(i)=0.5*jmax(i);
           end
        else
           if rmax>=-0.5
              rr(i)=(m-jmax(i)+1)*(1+rmax(i));
           else
              rr(i)=0.5*(m-jmax(i)+1);
           end
        end
    end
    i=1:n;
    wo(i)=rr(i)/sum(rr);           % Objective weights of evaluation indices                 
    %---------------------------------------------------------------------------------------%%
    w=0.5*ws+0.5*wo;   % Comprehensive weights for evalution indices 

    %%----------------------------------------Correlation analyis---------------------------------------%%
    for i=1:n
        for j=1:m
            v=abs(bb(i,j)-aa(i,j));                                            
            s=abs(xx(i,jj)-0.5*(aa(i,j)+bb(i,j)))-0.5*(bb(i,j)-aa(i,j));       
            sp=abs(xx(i,jj)-0.5*(app(i)+bpp(i)))-0.5*(bpp(i)-app(i));         
            if aa(i,j)==bb(i,j)
               if xx(i,jj)==aa(i,j)
                  k(i,j)=0
               else 
                  k(i,j)=abs(xx(i,jj)-aa(i,j))/(sp-(abs(xx(i,jj)-aa(i,j)))); 
               end
            else
                if (xx(i,jj)>=aa(i,j) && xx(i,jj)<=bb(i,j))                      
                   k(i,j)=-s/v;                                          
                else
                   k(i,j)=s/(sp-s);                                      
                end
            end
        end
    end
    %-----------Recognition of UKDS level---------------
    kp=w*k;                                            
    [kpmax,j0]=max(kp);                                
    
    %------------Characteristic variable--------------
    for j=1:m
        mkp(j)=(kp(j)-min(kp))/(max(kp)-min(kp));
    end
    Skp=0;
    for j=1:m
        Skp=Skp+j*mkp(j);
    end
    jstar=Skp/sum(mkp);                
    
%     ap
%     bp
     aa;
     bb;
     xx;
%     r
%     rr
%     jmax
%     wo
%     w
%     k
%     kp
%     j0
    
    %----------Store assessment results in EXCEL---------------%
    data1{jj,1}=xx;
    data2{jj,1}=wo;
    data3{jj,1}=w;
    data4{jj,1}=k;
    data5{jj,1}=kp;
    data6{jj,1}=j0;
    data7{jj,1}=jstar;
    
end

%%----------------------------------------Data output---------------------------------------%%
data1{jj,1}=xx;

data11= [data1{:}]; 
data111=data11';

data22= [data2{:}];                    
data222=(reshape(data22',pp,qq))';      

data33= [data3{:}];
data333=(reshape(data33',pp,qq))';

data44= [data4{:}];
data444=data44';

data55= [data5{:}];
data555=(reshape(data55',m,qq))';

data66= [data6{:}];
data666=data66';

data77= [data7{:}]; 
data777=data77';

string_1='xx';
xlswrite('OutputData_IEA.xlsx',{string_1},1,'A1');
xlswrite('OutputData_IEA.xlsx',data111,1,'A2');

string_2='wo';
xlswrite('OutputData_IEA.xlsx',{string_2},2,'A1');
xlswrite('OutputData_IEA.xlsx',data222,2,'A2');

string_3='w';
xlswrite('OutputData_IEA.xlsx',{string_3},3,'A1');
xlswrite('OutputData_IEA.xlsx',data333,3,'A2');

string_4='k(i,j)';
xlswrite('OutputData_IEA.xlsx',{string_4},4,'A1');
xlswrite('OutputData_IEA.xlsx',data444,4,'A2');

string_5='kp';
xlswrite('OutputData_IEA.xlsx',{string_5},5,'A1');
xlswrite('OutputData_IEA.xlsx',data555,5,'A2');

string_6='j0';
xlswrite('OutputData_IEA.xlsx',{string_6},6,'A1');
xlswrite('OutputData_IEA.xlsx',data666,6,'A2');

string_7='jstar';
xlswrite('OutputData_IEA.xlsx',{string_7},7,'A1');
xlswrite('OutputData_IEA.xlsx',data777,7,'A2');





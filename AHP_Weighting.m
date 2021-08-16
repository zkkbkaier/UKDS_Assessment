clc;
clear;
A=xlsread('Comparison_matrix.xlsx')   
[m,n]=size(A);                         
RI=[0 0 0.58 0.90 1.12 1.24 1.32 1.41 1.45 1.49 1.51];
R=rank(A);                           
[V,D]=eig(A);                         
tz=max(D);
B=max(tz);                           
[row, col]=find(D==B);               
C=V(:,col);                        
CI=(B-n)/(n-1);                      
CR=CI/RI(1,n);   
if CR<0.10
    disp('CI=');disp(CI);
    disp('CR=');disp(CR);
    disp('Comparison matrix meets the consistency requirement. The weight vector Q is:');
    Q=zeros(n,1);
    for i=1:n
        Q(i,1)=C(i,1)/sum(C(:,1)); 
    end
    Q                              
else
    disp('Comparison matrix can not meets the consistency requirement.The comparison matrix should be reconstructed.');
end
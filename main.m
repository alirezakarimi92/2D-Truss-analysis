close all;
clear;
clc;

disp('this program just design for a 2D trus');

n=input('n as number of nodes=');
m=input('m as number of elements=');

%make an excel file that every row is node i & node j for each element
nodes= xlsread('n&e.xlsx');
inform= xlsread('trust information.xlsx');

% now we will make hardness matrix
Kt=zeros(2*n,2*n);                    %Kt is hardness matrix for all elements
 
for q=1:m
    
    z= inform(q,3)*inform(q,4)/inform(q,2);   % z = A*E/L

    teta=inform(q,5);
    alfa=cosd(teta);
    beta=sind(teta);

    r=alfa^2;
    s=beta^2;
    p=alfa*beta;
    Ke=z*[r p -r -p;p s -p -s; -r -p r p; -p -s p s];
    k1=Ke([1 5;2 6]);
    k2=Ke([3 7;4 8]);
    
    i=nodes(q,2);
    j=nodes(q,3);
    
    a=4*(i-1)*n+2*i-1;
    b=a+1;
    c=(4*i-2)*n+2*i-1;
    d=c+1;
    Kt([a c;b d])=Kt([a c;b d])+k1;
    
    a=4*(j-1)*n+2*j-1;
    b=a+1;
    c=(4*j-2)*n+2*j-1;
    d=c+1;
    Kt([a c;b d])=Kt([a c;b d])+k1;
    
    a=4*(i-1)*n+2*j-1;
    b=a+1;
    c=(4*i-2)*n+2*j-1;
    d=c+1;
    Kt([a c;b d])=Kt([a c;b d])+k2;
    
    
    a=4*(j-1)*n+2*i-1;
    b=a+1;
    c=(4*j-2)*n+2*i-1;
    d=c+1;
    Kt([a c;b d])=Kt([a c;b d])+k2;
    
end
% finish making hardness matrix 

% make a matrix for forces & a matrix for displacement
 [f ff Ft]= xlsread('force.xlsx');                %Ft is total of force matrix with unknown parameter
 [a aa At]= xlsread('displacement.xlsx');          % At is total of displacement with unknown parameter
 k=zeros(2*n,1);
 x=zeros(n,1);
 w=1;                                     %index
 p=1;                                      %index
 ss=1;                                      %index
 sss=1;                                      %index
 SS=1;                                        %index
 
 for q=1:2*n
    S=isscalar(At{q});
    if S==0
        k=[k Kt(:,q)];            
        A{p}=At{q};          %A is nonnumerical part of At
        at(w)=0;             %at is total matrix of displacment
        U(p)=q;               % U is the matrix of indexes of unknown of At
        p=p+1;
    else
        at(w)=At{q};          %at is total known of displacement elements
    end  
    w=w+1;
 end
 k=k(:,2:end);               % k is nonzeros's part of Kt
 
 for q=1:2*n
     Q=isscalar(Ft{q});
     if Q==1           
        kF(ss)=Ft{q};        %kF is numerical part of Ft
        kk(ss,:)=k(q,:);     %kk is the part of hardness matrix who solve displacements
        ss=ss+1;
     else
        unF{SS}=Ft{q};       %unF is nonnumerical part of Ft
        kkk(sss,:)=k(q,:);   %kkk is the part of hardness matrix who solve forces
        sss=sss+1;
        SS=SS+1;
     end  
 end
 
 unknowndisp=inv(kk)*kF';
 unknownforces=kkk* unknowndisp;
 
 Ft    
 unf=unF'  % unkown forces element
 unknownforces
 At    
 A=A'      
 unknowndisp
 
 %finish determining the displacement & forces
 
 at(U)=unknowndisp;
 
 for q=1:m
     
    teta=inform(q,5);
    alfa=cosd(teta);
    beta=sind(teta);
    
    i=nodes(q,2);
    j=nodes(q,3);
    
    glo=at([2*i-1;2*i;2*j-1;2*j]);
    R=[alfa beta 0 0;-beta alfa 0 0;0 0 alfa beta; 0 0 -beta alfa];
    local=R*glo';
    strain(q)=(local(3)-local(1))/inform(q,2);
    stress(q)=inform(q,4)*strain(q);
    internalforces(q)=stress(q)*inform(q,3);
 end
 
 strain
 stress
 internalforces
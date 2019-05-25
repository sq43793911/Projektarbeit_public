clear all
clc
close all

%###########################################
% main Program
%###########################################
% parameter
E=2.1e11;           % N/m^2
A=0.0001;           % m^2
L=2;               % m  length of plate
B=2;               % m  wide of plate
rho=7850;           % Dichte in [kg/m^3]
mu=rho*A;           % Massenbelegung in [kg/m]
Nel=5;              % number of elements per rand
Nno=Nel+1;          % number of nodes per rand
Nel_all=Nel^2;      % total elements
Nno_all=Nno^2;      % total nodes
le=L/Nel;           % length of an element
%b=0.005;           % m Breit von Balken
h=0.01;             % m Hoehe von Balken
%I=(b*(h^3))/12;       % Flaechentraegheitsmoment
v=0.3;               % Poissionzahl



% define empty matrice
Kt=zeros(Nno_all*4);                                 % empty global stiffnes-matrix 
M=zeros(Nno_all*4);                                  % empty global mass-matrix 



% call element routinesmbclient
[Kte,Me] = Elementroutine_Platten(h,E,rho,le,v);

r=0;
a=zeros(4,1);
b=a;
c=[1 5 9 13];
d=[4 8 12 16];
for j=1 : Nel                                  % loop over every element
    for k=1:Nel
        a(1)=((j-1)*Nno+k-1)*4+1;
        b(1)=a(1)+3;
        a(2)=((j-1)*Nno+k)*4+1;
        b(2)=a(2)+3;
        a(3)=((j*Nno)+k-1)*4+1;
        b(3)=a(3)+3;
        a(4)=((j*Nno)+k)*4+1;
        b(4)=a(4)+3;
        for m=1:4
            for n=1:4
                M(a(m):b(m), a(n):b(n)) = M(a(m):b(m), a(n):b(n)) + Me(c(m):d(m),c(n):d(n));
                Kt(a(m):b(m), a(n):b(n)) = Kt(a(m):b(m), a(n):b(n)) + Kte(c(m):d(m),c(n):d(n));
            end
        end
    end
end

% implementation of essetial boundary conditions
% Kt(1,:) = [  ];
% Kt(:,1) = [  ];
% M(1,:)  = [  ];
% M(:,1)  = [  ];
% 
% Kt(1,:) = [  ];
% Kt(:,1) = [  ];
% M(1,:)  = [  ];
% M(:,1)  = [  ];

% e = 1;
% for i = 1:Nel
%     M(e:e+3,:) = [];
%     M(:,e:e+3) = [];
%     Kt(e:e+3,:) = [];
%     Kt(:,e:e+3) = [];
%     e = (i*Nno-i)*4;
% end

% for i = 1:Nel*4
%     M(e:e+3,:) = [];
%     M(:,e:e+3) = [];
%     Kt(e:e+3,:) = [];
%     Kt(:,e:e+3) = [];
%     e = (i*Nno-i)*4;
% end

M(1:Nno*4,:) = [];
M(:,1:Nno*4) = [];

Kt(1:Nno*4,:) = [];
Kt(:,1:Nno*4) = [];

% define system-matrix


null=zeros(size(M));
Eins = eye(size(M));
SysMat=[null,Eins;
       -inv(M)*Kt,null];




% compute Eigenvalues
EV = eig(SysMat);
[V,D,W]=eig(SysMat);
f = imag(EV)/(2*pi);


%%
% lamda = zeros(Nel, 1);
% omega = zeros(Nel, 1);
% f = zeros (Nel, 1);
% for k = 1: Nel
%     switch k
%         case 1 
%             lamda(k) = 1.87510 ;
%         case 2
%             lamda(k) = 4.69409 ;
%         case 3
%             lamda(k) = 7.85476 ;
%         case 4
%             lamda(k) = 10.99554;
%     end
%     if k >= 5
%             lamda(k) = ((2*k-1)*pi)/2 ;
%     end
%       
%     omega (k) = ((lamda(k))^2)*sqrt((E*I)/(rho*A*(L^4)));
% end


%%
% 
% im_V = imag(V);
% 
% lVec = zeros (Nno, 1);     % vector with node coordinates
%  
% 
% for k = 1: Nno
%     lVec(k)=L/Nel*(k-1);
% end
% 
% EigMat=V(1:length(M),:) ;   % delete lover part von eigenvectors 
% EigMat=imag(EigMat);
% for k = 1: length(M)/2      % delete unnecessary rows 
%    EigMat(1+k,:)=[];       
% end
% 
% for k = 1: length(M)
%    EigMat(:,k)=[];         % delete unnecessary colums
% end
% 
% EigMat=[zeros(1,length(M));EigMat];  % fill up EigMat for first with zero displacement
% 
% hold on
% plot(lVec,EigMat(:,end))
% plot(lVec,EigMat(:,end-1))
% plot(lVec,EigMat(:,end-2))


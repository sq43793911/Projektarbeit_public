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
Nelx=15;              % number of elements per rand
Nely=15;
Nnox=Nelx+1;          % number of nodes per rand
Nnoy=Nely+1;
Nel_all=Nelx * Nely;      % total elements
Nno_all=Nnox * Nnoy;      % total nodes
lex=L/Nelx;           % length of an element
ley=B/Nely;
%b=0.005;          % m Breit von Balken
h=0.01;          % m Hoehe von Balken
%I=(b*(h^3))/12;       % Flaechentraegheitsmoment
v=0.30;               % Poissionzahl



% define empty matrice
Kt=zeros(Nno_all*4);                                 % empty global stiffnes-matrix 
M=zeros(Nno_all*4);                                  % empty global mass-matrix 



% call element routinesmbclient
[Kte,Me] = Elementroutine_Platten(h,E,rho,lex,ley,v);

r=0;
a=zeros(4,1);
b=a;
c=[1 5 9 13];
d=[4 8 12 16];
for j=1 : Nely                                  % loop over every element
    for k=1:Nelx
        a(1)=((j-1)*Nnox+k-1)*4+1;
        b(1)=a(1)+3;
        a(2)=((j-1)*Nnox+k)*4+1;
        b(2)=a(2)+3;
        a(3)=((j*Nnox)+k-1)*4+1;
        b(3)=a(3)+3;
        a(4)=((j*Nnox)+k)*4+1;
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
% e = 1;
% for i = 1:Nel
%     M(e:e+3,:) = [];
%     M(:,e:e+3) = [];
%     Kt(e:e+3,:) = [];
%     Kt(:,e:e+3) = [];
%     e = (i*Nno-i)*4;
% end

% M(1:Nnox*4,:) = [];
% M(:,1:Nnox*4) = [];
% 
% Kt(1:Nnox*4,:) = [];
% Kt(:,1:Nnox*4) = [];

% define system-matrix


null=zeros(size(M));
Eins = eye(size(M));
SysMat=[null,Eins;
       -inv(M)*Kt,null];




% compute Eigenvalues
[V,D,W]=eig(SysMat);

temp_d = diag(D);
[nd, sortindex] = sort(temp_d);
temp_v = V(:,sortindex);
temp_f = imag(nd)/(2*pi);

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


xVec = zeros (Nnox, 1);     % vector with node coordinates
yVec = zeros (Nnoy, 1); 
 

for k = 1: Nnox
    xVec(k)=L/Nelx*(k-1);    
end

for g = 1: Nnoy
    yVec(g)=B/Nely*(g-1);    
end


EigMat=temp_v(1:length(M),:) ;   % delete lover part von eigenvectors 
EigMat=imag(EigMat);

for k = 1: length(M)/4      % delete unnecessary rows 
   EigMat(1+k:1+k+2,:)=[];       
end

for k = 1: length(M)
   EigMat(:,k)=[];         % delete unnecessary colums
end

% for k = 1: length(M)/4      % delete unnecessary rows 
%    EigMat(:,end-k-2:end-k)=[];       
% end



% EigMat=[zeros(1,length(M));EigMat];  % fill up EigMat for first with zero displacement

% p = zeros(Nnox,Nnoy,9);
% for i = 1:9
%     for j = 1:Nnoy
%         for k = 1:Nnox
%             p(k,j,i) = EigMat(end+1-k-(j-1)*Nnox,i);
%         end
%     end
% end
% for k = 1:9
%     subplot(3,3,k);
%     surf(yVec,xVec,p(:,:,k));
%     title(k)
% end


p = zeros(Nnox,Nnoy,9);
for i = 1:9
    for j = 1:Nnoy
        for k = 1:Nnox
            p(k,j,i) = EigMat(k+(j-1)*Nnox,i);
        end
    end
end
for k = 1:9
    subplot(3,3,k);
    surf(yVec,xVec,p(:,:,k));
    set(gca,'FontSize',20);
    title([num2str(k),'-Ableitung']);
end




% hold on
% plot(xVec,EigMat(:,end))
% plot(lVec,EigMat(:,end-1))
% plot(lVec,EigMat(:,end-2))


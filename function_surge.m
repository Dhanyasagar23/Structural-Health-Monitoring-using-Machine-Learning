function [c] = function_surge(E1,E2,E3,E4)


nel = 20;
nnode = 10;
ndof = nnode*2;
b = [1 2 3 4 5]; %bottom memebers
t = [7 8 9]; %top chord members
v = [ 11 14 17 20]; % vertical members
d = [6 12 13 15 16 18 19 10]; %diagonal members

ID = zeros(nel,12);
ID(:,1) = (1:nel)'; %elements number
ID(:,2:3) = [1 2; 2 3; 3 4; 4 5; 5 6; 1 7; 7 8; 8 9; 9 10; 6 10; 2 7; 2 8; 3 7; 3 8; 3 9; 4 8; 4 9; 4 10; 5 9; 5 10]; %starting and end node number for each element
ID(:,4:7) = [0 0 1 0; 1 0 2 0; 2 0 3 0; 3 0 4 0; 4 0 5 0; 0 0 1 1; 1 1 2 1; 2 1 3 1; 3 1 4 1;5 0 4 1; 1 0 1 1; 1 0 2 1; 2 0 1 1; 2 0 2 1; 2 0 3 1; 3 0 2 1;3 0 3 1;3 0 4 1; 4 0 3 1; 4 0 4 1]; % (x,y) co-ordinate
ID(:,8) = 2.5*10^(-3); % Area of each element
ID(b,9) = E1*10^7; % E of each element
ID(t,9) = E2*10^7;
ID(v,9) = E3*10^7;
ID(d,9) = E4*10^7;
%calculation of length of element and cos(phi) and sin(phi)
for e = 1:nel
 ID(e,10) = sqrt((ID(e,6) - ID(e,4))^2 + (ID(e,7) - ID(e,5))^2 ); %length of element e
 ID(e,11) = (ID(e,6) - ID(e,4))/ID(e,10); % cos(phi) of element e
 ID(e,12) = (ID(e,7) - ID(e,5))/ID(e,10); % sin(phi) of element e
end
 

%Global stiffness, gathered matrix,assembling, force vector

K = zeros(ndof,ndof); %global stiffness matrix
f = zeros(ndof,1); % global force vector due to external nodal node(P) and distributed load within the element
f(2*6-1) = 50 ;
for n = 2:5
 f(2*n) = -100;
end
k_L = zeros(nel*4,4); % big matrix whcih contains element stifness of each elements in global system
T = zeros(nel*4,4); % big matrix whcih contains Transformattion matrix of all elements
k_G = zeros(nel*4,4); % big matrix whcih contains all element stifness in global system
L = zeros(nel*4,20); %big matrix whcih contains gather operator for all elements 
for e = 1:nel
 k_L(4*(e-1)+1:4*e,:) = ID(e,8)*ID(e,9)/ID(e,10)*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];
c = ID(e,11); %cos(phi)
 s = ID(e,12); %sin(phi)
 T(4*(e-1)+1:4*e,:) = [c s 0 0; -s c 0 0; 0 0 c s; 0 0 -s c];
 k_G(4*(e-1)+1:4*e,:) = (T(4*(e-1)+1:4*e,:))'*k_L(4*(e-1)+1:4*e,:)*T(4*(e-1)+1:4*e,:);
 Temp = L(4*(e-1)+1:4*e,:); % for gather operator
 i = ID(e,2); % starting node of element
 j = ID(e,3); % ending node of element
 Temp(1,2*i-1) = 1;
 Temp(2,2*i) = 1;
 Temp(3,2*j-1) = 1;
 Temp(4, 2*j) = 1;
 L(4*(e-1)+1:4*e,:) = Temp; %gather operator for element e
 K = K + (L(4*(e-1)+1:4*e,:))'*(k_G(4*(e-1)+1:4*e,:))*(L(4*(e-1)+1:4*e,:)); %assembling local stiffness to global system
%usually we also calculte force vector in this same loop but in truss
%problems as there is no distributed load within the element we can
%calculate them seperately as external load will be at few nodes
end
%check values of different matrice
K;
f;
%Imposing boundary conditions using partition method

%delete first,second and 12th rows and columns from global stifness matrix as these are the essential bounday
%conditions and their value is zero we are using partition approach here
K_new = K;
f_new = f;
% here I have removed wrong row as I removed first row then the second row
% become now first row but still I removed second row
for i = [1 2 12]
 K_new(i,:) = [];
 K_new(:,i) = [];
 f_new(i,:) = [];
end
K_new;
f_new;

% in applying boundary conditions that is why matrix is singlur even after applying boundary condition
%correct way of imposing boundary
K_New = K;
f_New = f; %in this f vector 
for i = [12 2 1]
 K_New(i,:) = [];
 K_New(:,i) = [];
 f_New(i,:) = [];
end
K_New;
%Solve for displacements

%displacement vector
u_New = K_New\f_New;

u = zeros(ndof,1);
u(1:2) = [0;0];
u(3:11) = u_New(1:9);
u(12) = 0;
u(13:20) = u_New(10:17);
u;
c=u;



end
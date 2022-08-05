function [C] = function_surge1(c_)


c_;
nel = 20;% no of elements
nnode=10;
ndof=nnode*2;

ID = zeros(nel,12);
ID(:,1) = (1:nel)'; %elements number
ID(:,2:3) = [1 2; 2 3; 3 4; 4 5; 5 6; 1 7; 7 8; 8 9; 9 10; 6 10; 2 7; 2 8; 3 7; 3 8; 3 9; 4 8; 4 9; 4 10; 5 9; 5 10]; %starting and end node number for each element
ID(:,4:7) = [0 0 1 0; 1 0 2 0; 2 0 3 0; 3 0 4 0; 4 0 5 0; 0 0 1 1; 1 1 2 1; 2 1 3 1; 3 1 4 1;5 0 4 1; 1 0 1 1; 1 0 2 1; 2 0 1 1; 2 0 2 1; 2 0 3 1; 3 0 2 1;3 0 3 1;3 0 4 1; 4 0 3 1; 4 0 4 1]; % (x,y) co-ordinate
ID(:,8) = 2.5*10^(-3); % Area of each element
%ID(:,9) = 2*10^7; % E of each element
%calculation of length of element and cos(phi) and sin(phi)
for e = 1:nel
 ID(e,9) = sqrt((ID(e,6) - ID(e,4))^2 + (ID(e,7) - ID(e,5))^2 ); %length of element e
 ID(e,10) = (ID(e,6) - ID(e,4))/ID(e,9); % cos(phi) of element e
 ID(e,11) = (ID(e,7) - ID(e,5))/ID(e,9); % sin(phi) of element e
end

%K = zeros(ndof,ndof); %global stiffness matrix
f = zeros(ndof,1); % global force vector due to external nodal node(P) and distributed load within the element
f(2*6-1) = 50 ;
for n = 2:5
 f(2*n) = -100;
end
%given displacements

%u = zeros(ndof,1);%displacement vector of each node
%u=[0  0 0.005 -0.0478 0.011 -0.0683 0.0183 -0.0683 0.0244 -0.0478 0.0294 0 0.0303 -0.0454 0.0204 -0.0672 0.009 -0.0672 -0.0009 -0.0454]';
u=c_;
u;


%k_L = zeros(nel*4,4); % big matrix whcih contains element stifness of each elements in global system
T = zeros(nel*4,4); % big matrix whcih contains Transformattion matrix of all elements
%k_hat = zeros(nel*4,4); % big matrix whcih contains all element stifness in global system
L = zeros(nel*4,20); %big matrix whcih contains gather operator for all elements 
%k_hat=zeros(nel*4,4)
k_exp=zeros(nel*20,20);
for e = 1:nel
% k_L(4*(e-1)+1:4*e,:) = ID(e,8)*ID(e,9)/ID(e,10)*[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];
 k_hat(4*(e-1)+1:4*e,:)=(ID(e,8)/ID(e,9))*[1 0 -1 0;0 0 0 0;-1 0 1 0; 0 0 0 0];%element stiffness matrix for unit E
 c_ = ID(e,10); %cos(phi)
 s = ID(e,11); %sin(phi)
 T(4*(e-1)+1:4*e,:) = [c_ s 0 0; -s c_ 0 0; 0 0 c_ s; 0 0 -s c_];
 %k_G(4*(e-1)+1:4*e,:) = (T(4*(e-1)+1:4*e,:))'*k_L(4*(e-1)+1:4*e,:)*T(4*(e-1)+1:4*e,:);
 Temp = L(4*(e-1)+1:4*e,:); % for gather operator
 i = ID(e,2); % starting node of element
 j = ID(e,3); % ending node of element
 Temp(1,2*i-1) = 1;
 Temp(2,2*i) = 1;
 Temp(3,2*j-1) = 1;
 Temp(4, 2*j) = 1;
 L(4*(e-1)+1:4*e,:) = Temp; %gather operator for element e

 %now we need to calculate Le' Te'K_hat*Te'Le

  
  k_exp(20*(e-1)+1:20*e,:) = (L(4*(e-1)+1:4*e,:))'*(T(4*(e-1)+1:4*e,:)'*k_hat(4*(e-1)+1:4*e,:))*T(4*(e-1)+1:4*e,:)*L(4*(e-1)+1:4*e,:); %assembling local stiffness to global system
%usually we also calculte force vector in this same loop but in truss
%problems as there is no distributed load within the element we can
%calculate them seperately as external load will be at few nodes
end
%check values of different matrices
k_exp;


%now impose boundary condition in each of 20 x 20 matrix
K_new = zeros(17*nel,17);
u_new =u;
f_new = f;

for i = [12 2 1]
    
    f_new(i,:) = [];
    u_new(i,:) =[];
end   

for e = 1:nel
    K_= k_exp(20*(e-1)+1:20*e,:);
 
    for i = [12 2 1]
     K_(i,:) = [];
     K_(:,i) = [];
    end
    K_new(17*(e-1)+1:17*e,:) = K_;

    
end 
K_new;
u_new;
f_new;

 %K_= k_exp(20*(1-1)+1:20*1,:)
 %solve for displacements
%displacement vector
%u_New = K_new\f_new;
%u_New








%calculate alpha
alpha_big=zeros(ndof-3,nel);
for e=1:nel
  alpha_big(:,e)=K_new(17*(e-1)+1:17*e,:)*u_new;
end 
alpha_big;

%x=pinv(alpha_big)*f_new% 17 equations and 20 unknowns, we should not go
%through this process
%infact use what is given like all top members have same E's

%now our task is to reduce the number of unknowns that we can reduse as we
%know  (i) all top chord members 
%have same E (= ET), (ii) all bottom chord members have same E (= EB), (iii) all vertical 
%members have same E (= EV), and (iv) all diagonal members have same E (= ED)


% we will solve it by grouping of columns having same EI's by making only
% four unknowns


E = zeros(4,1); % first entity for bottom memebers, second for top, third for vertical,forth for diagonal
alpha_reduced = zeros(17,4); %after grouping of columns
b = [1 2 3 4 5]; %bottom memebers
t = [7 8 9]; %top chord members
v = [ 11 14 17 20]; % vertical members
d = [6 12 13 15 16 18 19 10]; %diagonal members

for i = b
    alpha_reduced(:,1) = alpha_reduced(:,1) + alpha_big(:,i);
end

for i = t
    alpha_reduced(:,2) = alpha_reduced(:,2) + alpha_big(:,i);
end

for i = v
    alpha_reduced(:,3) = alpha_reduced(:,3) + alpha_big(:,i);
end

for i = d
    alpha_reduced(:,4) = alpha_reduced(:,4) + alpha_big(:,i);
end

alpha_reduced;
f_new;
%calculating
E_= pinv(alpha_reduced)*f_new;

C=E_;



end

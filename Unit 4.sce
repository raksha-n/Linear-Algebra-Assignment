//QUESTION 1
//GRAM - SCHMIDT ORTHOGONALIZATION IN R3

clc;clear;close 

A=input("enter 3x3 matrix") 
disp(A,'A=');
[m,n]=size(A); 

for k=1:n 
 V(:,k)=A(:,k);
  for j=1:k-1 
    R(j,k)=V(:,j)'*A(:,k); 
    V(:,k)=V(:,k)-R(j,k)*V(:,j); 
  end 
 R(k,k)=norm(V(:,k)); 
 V(:,k)=V(:,k)/R(k,k);
end 
disp(V,'Q=');




//QUESTION 2
//EIGEN VALUES AND EIGEN VECTORS FOR 3x3 MATRIX

clc; close(); clear;
A =input("enter 3x3 matrix")
lam = poly(0,'lam')
lam = lam
charMat = A-lam*eye(3,3)
disp(charMat,'The characteristic Matrix is')
charPoly = poly(A,'lam')
disp(charPoly,'The characteristic polynomial is')
lam = spec(A)
disp(lam,'The eigen values of A are')
function[x,lam] = eigenvectors(A)
    [n,m] = size(A);
    lam = spec(A)';
    x = [];
    for k=1:3
        b = A-lam(k)*eye(3,3); //characteristic matrix
        c = b(1:n-1,1:n-1); //coeff mat for the reduced system
        b1 =-b(1:n-1,n); //rhs vector for the reduced system
        y = c\b1; //solution for the reduced system
        y = [y:1]; //complete eigen vector
        y = y/norm(y); //make unit eigen vector
        x = [x y];
    end
endfunction

get f('eigenvectors')
[x,lam] = eigenvectors(A)
disp(x,'The eigen vectors of A are')




//QUESTION 3
//NUMERICALLY LARGEST EIGEN VALUE USING RAYLEIGH POWER METHOD

clear; clc; close();
a = input("enter 3x3 matrix")
disp(a,'A = ')
//initial vector
u0 = [1 1 1]';
disp(u0,'The initial vector is')
v = a*u0
a1 = max(u0)
disp(a,'First approximation to eigen value is ')
while abs(max(v)-a1)>0.002
    disp(v,'Current eigen vector is ')
    a1 = max(v)
    disp(a1,'Current eigen value is ')
    u0 = v/max(v)
    v = a*u0
end
format('v',4)
disp(max(v),'The largest eigen value is: ')
format('v',5)
disp(u0,'The corresponding eigen vector is: ')

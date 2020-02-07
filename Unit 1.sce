//Gaussian Elimination for a generalized case

A=input("Enter the coefficient matrix of nxn: ")
b=input("Enter the constants matrix nx1: ")
function [x]=gaussian_elimination(A,b)
    [a1,a2]=size(A);//n,n1
    [b1,b2]=size(b);//m1,p
    
    if a1~=a2
        error('Matrix A must be square');
        abort;
    else if a1~=b1
        error('Incompatible orders of A and b');
        abort;
    end;
    
    Aug=[A b];
    //forward elimination
    n=size(A,1);
    for k=1:a1-1
        for i=k+1:a1
            factor=A(i,k)/A(k,k);
            for j=k+1:a1
                A(i,j)=A(i,j)-factor*A(k,j);
            end
            b(i)=b(i)-factor*b(k)
        end
    end
    
    //back substitution
    
    x(a1)=b(a1)/A(a1,a1);
    for i=a1-1:-1:1
        sum=0;
        for j=i+1:a1
            sum=sum+A(i,j)*x(j);
        end
        x(i)=(b(i)-sum)/A(i,i);
    end
end
endfunction  

//Gauss Jordan Method to find Inverse

A=input("Enter the coefficient matrix of mxn: ")
function [D]=gaussian_jordan(A)
    [a1,a2]=size(A);//n,n1
    b=eye(a1,a2)
    C=[A b];
    [m,n]=size(C)
    for k=1:1:m
        indices=[1:1:k-1,k+1:1:m]
        for i=indices
            multiplier=C(i,k)/C(k,k)
            for j=k+1:n
                C(i,j)=C(i,j)-multiplier*C(k,j)
             end
         end
    end
    D=zeros(a1,a2)
    for i=1:1:a1
        for j=1:1:a2
            D(i,j)=C(i,m+j)/C(i,i)
        end
    end
    endfunction

//LU Decomposition

//FACTORIZING A INTO L AND U (A = LU)
clc;clear;
function lu_decomposition(A)
    [r,c]=size(A);
    u=A;
    l=eye(r,c);
    for i=1:(r-1)
        m=det(u(i,i));
        for j=i+1:c
            n=det(u(j,i))
            a=n/m;
            l(j,i)=a;
            u(j,:)=u(j,:)-u(i,:)/(m/n);
        end
    end
disp(l,'The lower triangular matrix L is');
disp(u,'The upper triangular matrix U is');
endfunction

disp('Factorization of A into L and U');
A=input('Enter elements of matrix: ');
disp(A,'The given matrix is A=');
lu_decomposition(A);

//SOLVING SYSTEM OF EQUATIONS BY LU DECOMOSITION
clc;clear;
format('v',5);
function solve_lu(a, b)
    [r,c]=size(a);
    b=b';
    l=eye(r,c);
    for i=1:r
        for j=1:c
            s=0;
            if j>=i
                for k=1:i-1
                    s=s+l(i,k)*u(k,j);
                end
                u(i,j)=a(i,j)-s;
            else
                for k=1:j-1
                    s=s+l(i,k)*u(k,j);
                end
                l(i,j)=(a(i,j)-s)/u(j,j);
            end
        end
    end
    c=l\b;
    x=u\c;
disp(l,'The lower triangular matrix L is');
disp(u,'The upper triangular matrix U is');
disp(x,'Solution of system of equation is ');
endfunction

disp("Solving system of equation by LU decomposition");
a=input('Enter elements of matrix A: ');
b=input('Enter elements of matrix B: ');
disp(a,'The coefficient matrix A is');
disp(b,'The constant matrix b is');
solve_lu(a,b);
    

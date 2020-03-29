//QUESTION 1
//Column Space of 3*3 matrix

a=input("Enter a 3x3 matrix: ")
disp(a)
//b=a
function[]=column_space(a)

disp('To get Upper Triangular form:')
a(2,:) = a(2,:)-(a(2,1)/a(1,1))*a(1,:)
a(3,:) = a(3,:)-(a(3,1)/a(1,1))*a(1,:)
disp(a)
a(3,:) = a(3,:)-(a(3,2)/a(2,2))*a(2,:)
disp(a)
//disp(b)
disp('The column space is:')
x=0
for i=1:3
    for j=i:3
        if a(j,i)~=0 then
            disp(b(:,i))
            x=1C
            break
        end
    end
end
if  x==0 then
    disp([0 0 0])
end
endfunction

//QUESTION 2
//The Four Subspaces of a 3x3 matrix

clear

A=input("Enter a 3x3 matrix: ")
disp(A)

[m,n] = size(A)
disp(m)
disp(n)
[v,pivot] = rref(A)
disp(rref(A))
disp(v)
disp(pivot)
rank = length(pivot)
disp(rank,'rank = ')
cs = A(:,pivot)
disp(cs,'Column space = ')
ns = kernel(A)
disp(ns,"null space = ")
rs = v(1:rank,:)'
disp(rs,'Row space = ')
lns = kernel(A')
disp(lns,"Left Null space = ")


//QUESTION 1
//Projection by Least Square

clear
n = input('Enter the number of equations:')

a=input('Enter matrix A of order nx2')
if size(a)~=[n,2] then
    error('Matrix not consistent with number of equations');
        abort;
end
disp(a)

for i=1:n
    b(i) = input("Enter the values of b:")
end
disp(b)

x=(a'*a)\(a'*b)
disp(x)
C=x(1,1)
D=x(2,1)
disp(C,'C=')
disp(D,'D=')
disp("the best line fit is b=C+Dt.")

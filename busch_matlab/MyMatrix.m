A=load('../matr.dat');
B = zeros(25,25);
for i=1:size(A,1)
   B(round(A(i,1))+1,round(A(i,2))+1) = A(i,3); 
end
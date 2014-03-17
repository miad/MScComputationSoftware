A=load('corr2.dat');
[Y,I]=sort(A(:,1));
B=A(I,:); %use the column indices from sort() to sort all columns of A.
x=B(:,1);
y=B(:,2);


freq=75;
offset = 1;
mean=14.65;
magn=0.01;
func=@(x) magn*(sin(freq*(x -0.5)+2.8)) + mean;
clf;
plot(x,y,'rx');
hold on;
plot(x, func(x));
%xlim([0.5 0.7]);
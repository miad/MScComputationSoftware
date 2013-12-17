AA=load('matrix.dat');
diff = abs(AA(:,1)-AA(:,2));
bins = zeros(1, max(diff)+1);
counts = bins;
for k=1:size(AA,1)
    bins(1, diff(k)+ 1 ) = bins(1, diff(k)+ 1) + abs(AA(k,3)) + 1i*abs(AA(k,4));
    counts(1, diff(k) + 1) = counts(1, diff(k)+ 1 ) + 1;
end
s = bins./counts;
s = s(1:ceil(max(diff)/2));
s=abs(real(s))+1i*abs(imag(s))+(1+1i)*1E-10;
d=0:floor(max(diff)/2);
semilogy(d,real(s));
hold on;
semilogy(d,imag(s),'r');
xlabel('Distance to diagonal');
ylabel('Value');
h = legend('Real part', 'Imaginary part');
title('Matrix element statistics');

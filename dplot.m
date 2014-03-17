
A=load('run1.txt');
for i=1:size(A,1)-1
    plot(A(i:i+1,1),A(i:i+1,2),'or');
    pause(1)
    hold on
end

A=load('run3.txt');
for i=1:size(A,1)-1
    plot(A(i:i+1,1),A(i:i+1,2),'x');
    pause(1)
    hold on
end

A=load('run4.txt');
for i=1:size(A,1)-1
    plot(A(i:i+1,1),A(i:i+1,2),'g*');
    pause(1)
    hold on
end
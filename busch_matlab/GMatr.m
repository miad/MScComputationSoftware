uN=15;
Nmax=uN*uN;
C = zeros(Nmax,Nmax);
g=1;
for i=0:Nmax-1
    i
    for j=0:Nmax-1
        a = mod(i,uN);
        b = floor(i/uN);
        c = mod(j,uN);
        d = floor(j/uN);
        
        if a==c && b==d
            C(i+1, j+1) = 1+(a+b);
        end
        C(i+1, j+1) = C(i+1, j+1) + g*Element(a, b, c, d);
    end
end
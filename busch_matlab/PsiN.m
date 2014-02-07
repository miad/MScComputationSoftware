function val = PsiN(n, x)
    s = h_polynomial(size(x,2),n, x')';
    val = 1/sqrt(2^n*factorial(n))*(1/pi)^(1/4)*exp(-1*x.^2/2).*s(size(s,1),:);
end
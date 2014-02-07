function val = Element(a, b, c, d)
    val = integral(@(x) PsiN(a,x).*PsiN(b,x).*PsiN(c,x).*PsiN(d,x),-200,200);
end
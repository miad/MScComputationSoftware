function val = InnerProd(a, b)
    val = integral(@(x) PsiN(a,x).*PsiN(b,x),-200,200,'AbsTol',1E-5);
end
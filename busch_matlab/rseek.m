fCount = 0;
match = zeros(625,4);
for i=1:25
    for j=1:25
        found = 0;
        for a=1:25
            if found==1
                break;
            end
            for b=1:25
                if abs(B(i,j)-C(a,b)) < 1E-6
                    found = 1;
                    match(25*i+j,:) = [i j a b];
                    break;
                end
            end
        end
        if found ~=1
            disp 'Not found '
        else
            fCount = fCount + 1;
        end
    end
end
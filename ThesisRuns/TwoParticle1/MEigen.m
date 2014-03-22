function MEigen(inputName, outputName)
    disp(['Initiating read from ' inputName]);
    A=load(inputName);
    disp(['Read complete, initiating reshape procedures.']);
    B=A(:,3)+1i*A(:,4);
    clear A
    dim=sqrt(size(B,1));
    B=reshape(B,dim,dim);
    
    disp(['Reshape completed, now solving for eigenvalues.']);
    
    C=eig(B);

    disp(['Found eigenvalues, now saving to file ' outputName]);
    dlmwrite(outputName,[real(C) imag(C)],'delimiter','\t','precision',20);
end
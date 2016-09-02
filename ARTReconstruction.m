function u   = ARTReconstruction(A,d,relaxParam,numIterART,u)

% Norm of rows
Anorm               = sum(A.*A,2);

for it = 1:numIterART
    for ir = 1:size(A,1)
        num         = d(ir)-A(ir,:)*u;
        u           = u + (relaxParam*num/Anorm(ir))*A(ir,:)';
    end % ir
end % it
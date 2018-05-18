function PsiBasis = nD_polynomial_array(NbRva,OrderPol)


% Initialize the PsiBasis indeces
PsiBasis = zeros(int32(factorial(OrderPol+NbRva)/factorial(OrderPol)/factorial(NbRva)),NbRva);


MM = NbRva - 1 ;
n = 1;

for CurrentOrder = 1 : OrderPol
    EndGenere = 0;
    FirstThisOrder = 0;
    
    while (EndGenere == 0)
        n = n +1 ;
        % First list t for order CurrentOrder
        if ( FirstThisOrder == 0)
            for i=1 : MM
                t(i) = i;
            end;
            FirstThisOrder =1 ;
        else
            % Regular incrementation
            if (t(MM) < (MM + CurrentOrder))
                t(MM) = t(MM) + 1;
            else  % t(MM) = tmax = MM +CurrentOrder
                j = MM;
                while (t(j) == j + CurrentOrder )
                    j = j - 1 ;
                end;
                t(j) = t(j) + 1 ;
                for k =(j + 1) :  MM
                    t(k) = t(j) + k - j ;
                end;
            end;
        end;
        
        % Direct Translating t into PsiBasis{n}
        PsiBasis(n,1) = t(1) - 1;
        for i=2 : MM
            PsiBasis(n,i) = t(i) - t(i-1) -1 ;
        end;
        PsiBasis(n,NbRva) = NbRva + CurrentOrder - t(MM) -1 ;
        
        % End of generation of order CurrentOrder
        if (t(1) == (CurrentOrder+1))
            EndGenere = EndGenere + 1;
        end;
    end;
end;


for i=1:floor(NbRva/2)
    swap = PsiBasis(:,i);
    PsiBasis(:,i) = PsiBasis(:,NbRva-i+1);
    PsiBasis(:,NbRva-i+1) = swap;
end


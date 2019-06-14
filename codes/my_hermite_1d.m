% This script evaluated a 1-d Hermite polynomial on [-1,1]

% p is the maximum order of the PC (total order)

function Hermite = my_hermite_1d(p,x)

Hermite = zeros(p+1,size(x,2));

if p==0 
    Hermite(p+1,:) = ones(1,size(x,2));
elseif p==1
    Hermite(p  ,:) = ones(1,size(x,2));
    Hermite(p+1,:) = x;
else 
    Hermite(1,:) = ones(1,size(x,2));
    Hermite(2,:) = x;
    for ord = 2:p
        Hermite(ord+1,:) = x.* Hermite(ord+1-1,:) - (ord-1) * Hermite(ord+1-2,:);
    end
end

% Now normalize

for i=0:p
    Hermite(i+1,:) = Hermite(i+1,:)/sqrt(factorial(i));
end


        


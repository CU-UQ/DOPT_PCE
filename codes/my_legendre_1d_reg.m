% This script evaluated a 1-d Legendre polynomial on [-1,1]

% p is the maximum order of the PC (total order)

function Legendre = my_legendre_1d_reg(p,x)

Legendre = zeros(p+1,size(x,2));

if p==0 
    Legendre(p+1,:) = ones(1,size(x,2));
elseif p==1
    Legendre(p  ,:) = ones(1,size(x,2));
    Legendre(p+1,:) = x;
else 
    Legendre(1,:) = ones(1,size(x,2));
    Legendre(2,:) = x;
    for ord = 2:p
        Legendre(ord+1,:)=( (2*ord-1) * x.* Legendre(ord+1-1,:) - (ord-1) * Legendre(ord+1-2,:) ) / ord;
    end
end

% Now normalize
for i=0:p
    Legendre(i+1,:) = Legendre(i+1,:)/sqrt(2/(2*i + 1));
end


        


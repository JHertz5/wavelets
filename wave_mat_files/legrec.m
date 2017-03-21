function p=legrec(lc)

% max degree 5

% Legendre polynomials
LP = [ 0 0 0 0 0 0 1; ...
       0 0 0 0 0 1 0; ...
       0 0 0 0 3/2 0 -1/2; ...
       0 0 0 5/2 0 -3/2 0; ...
       0 0 35/8 0 -30/8 0 3/8; ...
       0 63/8 0 -70/8 0 15/8 0; ...
       231/16 0 -315/16 0 105/16 0 -5/16];


lc = lc(:);
d = length(lc);

if (d > 7)
   error('Max degree is 6');
   exit;
end

lc = [zeros(1, 7-d) lc'];
% Normalize
c = (2*(6:-1:0) + 1)/2;
lc = fliplr(lc.*c);

% These are the coefficients in the reconstruction formula
p = lc*LP;
p = p(7-d+1:7);


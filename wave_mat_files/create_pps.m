function [x, t]=create_pps(M,K,pts)
%function [x, t]=create_pps(M,K,pts) creates a piecewise polynomial signal,
%M is the number of polynomials 
%K maximum degree of any polynomial in the signal
%pts size of the signal
%t vector with the discontinuity locations
%x piecewise polynomial signal


t = rand(1, M-1); 

t = ceil(sort(t.*(pts-5)))+2;



sig = zeros(1, pts);
t = [1 t pts];
for n = 1:M
        deg(n) = ceil(rand * (K+1));
	cfs(n,:) = [rand(1,deg(n))-.5 zeros(1, K+1-deg(n))];
	cf = legrec(cfs(n,1:deg(n)));
	
	wd = 2/(t(n+1)-t(n));
	x = -1:wd:1;
	xc = ones(size(x));
	y = zeros(size(x));
	for k=1:deg(n)
		y = cf(k) * xc;
		xc = xc.*x;
	end
	
	sig(t(n):t(n+1)) = y;
end

x=sig;


%A=imread('lena.bmp','bmp');

%x=double(A(200,1:256));

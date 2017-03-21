%Compute orthoprojection

N=16;
X=linspace(1,N,N);
Y=rand(N,1)+randn*linspace(0,1,N)';
U=ones(N,1);

A=[X' U];
C=inv(A'*A)*A'*Y
x_est=C(1)*X+C(2)*U';

figure(1)
scatter(X,Y)
hold on
plot(X,x_est)
hold off


A=[(X.*X.*X)' (X.*X)' X' U];
C=inv(A'*A)*A'*Y
x_est_2=C(1)*(X.*X.*X)+C(2)*(X.*X)+C(3)*X+C(4)*U';

figure(2)
scatter(X,Y)
hold on
plot(X,x_est_2)
hold off
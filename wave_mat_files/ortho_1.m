clear
x=[1 3 5.5 5 3 3 3 3 3 -0.5 -1 -0.5 0 0.5 1 1.5 2 3 2 1 1 1 1 3 4 2 2 0 -1.5 -1 0 0.5];
figure(1)
subplot(3,1,1)
stem(x,'.');
axis([1 32 -2 6])
xlabel('x[n]')
[aa bb]=dwt(x,'Haar');
ee=zeros(1,16);
ee2=zeros(1,8);
x_v_1=idwt(aa,ee,'Haar');
x_w_1=idwt(ee,bb,'Haar');
subplot(3,1,2)
stem(x_v_1,'.');
axis([1 32 -2 6])
xlabel('x_{v_1}[n]')
subplot(3,1,3)
stem(x_w_1,'.');
axis([1 32 -2 6])
xlabel('x_{w_1}[n]')
[aa2 bb2]=dwt(aa,'Haar');
y1=idwt(aa2,ee2,'Haar');
x_v_2=idwt(y1,ee,'Haar');


figure(2)
subplot(3,1,1)
stem(x,'.');
axis([1 32 -2 6])
xlabel('x[n]')
subplot(3,1,2)
stem(x_v_2,'.');
axis([1 32 -2 6])
xlabel('x_{v_2}[n]')
y2=idwt(ee2,bb2,'Haar');
x_w_2=idwt(y2,ee,'Haar');
subplot(3,1,3)
stem(x_w_2,'.');
axis([1 32 -2 6])
xlabel('x_{w_2}[n]')

figure(5)
[x t]=create_pps(3,2,256);
[C L]=wavedec(x,2,'Haar');

subplot(4,1,1)
plot(x);
subplot(4,1,2)
stem(C(129:end),'.');
subplot(4,1,3)
stem(C(65:128),'.')
subplot(4,1,4)
stem(C(1:64),'.')


figure(6)
[C2 L]=wavedec(x,2,'db3');

subplot(4,1,1)
plot(x);
subplot(4,1,2)
stem(C2(129:end),'.');
subplot(4,1,3)
stem(C2(65:128),'.')
subplot(4,1,4)
stem(C2(1:64),'.')








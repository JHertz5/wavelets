clear
x=[1 3 5.5 5 3 3 3 3 3 -0.5 -1 -0.5 0 0.5 1 1.5 2 3 2 1 1 1 1 3 4 2 2 0 -1.5 -1 0 0.5];
figure(3)
subplot(3,1,1)
stem(x,'.');
axis([1 32 -2 6])
xlabel('x[n]')
[aa bb]=dwt(x,'bior2.2');
ee=zeros(1,16);
ee2=zeros(1,8);
x_v_1=idwt(aa,ee,'bior2.2');
x_w_1=idwt(ee,bb,'bior2.2');
subplot(3,1,2)
stem(x_v_1,'.');
axis([1 32 -2 6])
xlabel('x_{v_1}[n]')
subplot(3,1,3)
stem(x_w_1,'.');
axis([1 32 -2 6])
xlabel('x_{w_1}[n]')
[aa2 bb2]=dwt(aa,'bior2.2');
y1=idwt(aa2,ee2,'bior2.2');
x_v_2=idwt(y1,ee,'bior2.2');
figure(4)
subplot(3,1,1)
stem(x,'.');
axis([1 32 -2 6])
xlabel('x[n]')
subplot(3,1,2)
stem(x_v_2,'.');
axis([1 32 -2 6])
xlabel('x_{v_2}[n]')
y2=idwt(ee2,bb2,'bior2.2');
x_w_2=idwt(y2,ee,'bior2.2');
subplot(3,1,3)
stem(x_w_2,'.');
axis([1 32 -2 6])
xlabel('x_{w_2}[n]')

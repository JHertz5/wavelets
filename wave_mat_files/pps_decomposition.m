clear
dwtmode('per');
pts=64;
[x t]=create_pps(4,2,pts);
figure 
stem(x,'.');
axis([1 pts -0.6 0.6])
xlabel('x[n]')
[aa bb]=dwt(x,'Haar');
ee=zeros(1,pts/2);
x_v_1=idwt(aa,ee,'Haar');
x_w_1=idwt(ee,bb,'Haar');
figure
stem(x_v_1,'.');
axis([1 pts -0.6 0.6])
xlabel('x_{v_1}[n] with Haar')
figure
stem(x_w_1,'.');
axis([1 pts -0.6 0.6])
xlabel('x_{w_1}[n] with Haar')
[aa bb]=dwt(x,'db2');
ee=zeros(1,pts/2);
x_v_1=idwt(aa,ee,'db2');
x_w_1=idwt(ee,bb,'db2');
figure
stem(x_v_1,'.');
axis([1 pts -0.6 0.6])
xlabel('x_{v_1}[n] with db2')
figure
stem(x_w_1,'.');
axis([1 pts -0.6 0.6])
xlabel('x_{w_1}[n] with db2')
%[aa2 bb2]=dwt(aa,'Haar');
%y1=idwt(aa2,ee2,'Haar');
%x_v_2=idwt(y1,ee,'Haar');
%figure
%stem(x_v_2,'.');
%axis([1 32 -2 6])
%xlabel('x_{v_2}[n]')
%y2=idwt(ee2,bb2,'Haar');
%x_w_2=idwt(y2,ee,'Haar');
%figure
%stem(x_w_2,'.');
%axis([1 32 -2 6])
%xlabel('x_{w_2}[n]')

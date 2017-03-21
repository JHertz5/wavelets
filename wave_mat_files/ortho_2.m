clear
A=imread('cameraman','bmp');
B=double(A);
[aa bb cc dd]=dwt2(B,'Haar');
ee=zeros(128,128);
ee2=zeros(64,64);
ee3=zeros(32,32);
figure(1)
colormap(gray)
imagesc(A);
xlabel('Original Image')
%figure
%colormap(gray)
%imagesc([aa ee; ee ee])
y=idwt2(aa,ee,ee,ee,'Haar');
figure(2)
colormap(gray)
imagesc(y)
xlabel('V_1')
[aa2 bb2 cc2 dd2]=dwt2(aa,'Haar');
y2=idwt2(aa2,ee2,ee2,ee2,'Haar');
figure(3)
colormap(gray)
imagesc(y2)
xlabel('V_2')
[aa3 bb3 cc3 dd3]=dwt2(aa2,'Haar');
y3=idwt2(aa3,ee3,ee3,ee3,'Haar');
figure(4)
colormap(gray)
imagesc(y3)
xlabel('V_3')
d3=idwt2(ee3,bb3,cc3,dd3,'Haar');
figure(5)
colormap(gray)
imagesc(d3)
xlabel('W_3')
d2=idwt2(ee2,bb2,cc2,dd2,'Haar');
figure(6)
colormap(gray)
imagesc(d2)
xlabel('W_2')
d1=idwt2(ee,bb,cc,dd,'Haar');
figure
colormap(gray)
imagesc(d1)
xlabel('W_1')

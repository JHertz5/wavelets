clear
A=imread('cameraman','bmp');
B=double(A);
dwtmode('per');

[aa bb cc dd]=dwt2(B,'bior2.2');
ee=zeros(128,128);
ee2=zeros(64,64);
ee3=zeros(32,32);
colormap(gray)
imagesc(A);
%figure
%colormap(gray)
%imagesc([aa ee; ee ee])
y=idwt2(aa,ee,ee,ee,'bior2.2');
figure
colormap(gray)
imagesc(y)
[aa2 bb2 cc2 dd2]=dwt2(aa,'bior2.2');
y2=idwt2(aa2,ee2,ee2,ee2,'bior2.2');
figure
colormap(gray)
imagesc(y2)
[aa3 bb3 cc3 dd3]=dwt2(aa2,'bior2.2');
y3=idwt2(aa3,ee3,ee3,ee3,'bior2.2');
figure
colormap(gray)
imagesc(y3)
d3=idwt2(ee3,bb3,cc3,dd3,'bior2.2');
figure
colormap(gray)
imagesc(d3)
d2=idwt2(ee2,bb2,cc2,dd2,'bior2.2');
figure
colormap(gray)
imagesc(d2)
d1=idwt2(ee,bb,cc,dd,'bior2.2');
figure
colormap(gray)
imagesc(d1)


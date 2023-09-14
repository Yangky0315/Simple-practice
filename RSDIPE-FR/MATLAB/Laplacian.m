img=imread('building.gif');
img=double(img);
w=fspecial('laplacian',0);%创建拉普拉斯滤波器
ruihua=imfilter(img,w);%对图像进行滤波
ruihua=uint8(img-ruihua);
figure
imshow(ruihua),title('锐化图像');
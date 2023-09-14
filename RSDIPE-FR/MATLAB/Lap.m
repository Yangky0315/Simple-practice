w=fspecial('laplacian',0);%创建拉普拉斯滤波器
w=double(w);
F_w=fft2(w);%傅立叶变换到频率域
F_w=fftshift(F_w);
F_w=log(abs(F_w)+1);
figure
imshow(F_w,[]);
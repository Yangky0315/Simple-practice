w=fspecial('laplacian',0);%����������˹�˲���
w=double(w);
F_w=fft2(w);%����Ҷ�任��Ƶ����
F_w=fftshift(F_w);
F_w=log(abs(F_w)+1);
figure
imshow(F_w,[]);
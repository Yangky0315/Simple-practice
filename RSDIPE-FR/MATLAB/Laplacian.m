img=imread('building.gif');
img=double(img);
w=fspecial('laplacian',0);%����������˹�˲���
ruihua=imfilter(img,w);%��ͼ������˲�
ruihua=uint8(img-ruihua);
figure
imshow(ruihua),title('��ͼ��');
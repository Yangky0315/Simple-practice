img1=imread('quickbirdsub.tif');
img2=imread('livingroom.tif');
img2_hist=imhist(img2);
img_match=histeq(img1,img2_hist);%ֱ��ͼƥ��
figure
subplot(1,2,1),imshow(img1);title('ԭͼ��');
subplot(1,2,2),imshow(img2);title('ƥ��ͼ��');
figure
imshow(img_match);title('ƥ��֮��ͼ��');
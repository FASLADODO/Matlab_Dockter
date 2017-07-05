im_in = imread('Zebra.jpg');

figure(1)
imshow(im_in)

imgray = rgb2gray(im_in);

figure(2)
imshow(imgray)


edges = edge(imgray,'canny');

figure(3)
imshow(edges)

C = corner(edges,'Harris');

figure(4)
imshow(imgray);
hold on
plot(C(:,1), C(:,2), 'r*');
hold off
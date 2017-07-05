%load it up
A = imread('plaque_transparency.PNG');
figure
image(A)

%convert to double and gray
RBG = im2double(A);
Gray = rgb2gray(RBG);

%get black white image
BW = im2bw(Gray, 0.1);

%display black whiteimage to tune thresholds
figure
imshow(BW, []);

%find circles
[centers, radii, metric] = imfindcircles(BW,[15 50]);

%display circles on top of image
figure
image(A)
hold on
viscircles(centers, radii,'EdgeColor','b');
hold off



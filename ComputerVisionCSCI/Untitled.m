t = 0:0.1:5;
y = -0.5.*exp(-t).*sin(2*t) + 2.*exp(-t).*cos(2*t);
plot(t,y)
grid on % Turn on grid lines for this plot


%%

picname = 'dummy1.jpg';
inputpic = imread(picname);
picg = rgb2gray(inputpic);
picg=im2double(picg);

[row,col]=size(picg);
edgepic=zeros(row,col);
edgepic=mat2gray(edgepic);
edgepic=im2double(edgepic);
magthreshold=definethresh;

figure, imshow(edgepic), title('binary gradient mask');
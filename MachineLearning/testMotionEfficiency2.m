sampledata = [0,0,0;
            1,1,1;
            2,2,2;
            3,3,3;
            4,4,4;
            3,3,5;
            2,2,6;
            1,1,7;
            0,0,8;
            0,0,7];
        
        
[ME,X1,X2,~] = MotionEfficieny(sampledata);   
X1
X2
figure
plot3(sampledata(:,1),sampledata(:,2),sampledata(:,3))
axis([0 8 0 8 0 8])

figure
scatter3(sampledata(:,1),sampledata(:,2),sampledata(:,3),20,ME)
axis equal
colormap cool
colorbar
axis([0 8 0 8 0 8])

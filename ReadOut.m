fileID = fopen('data.bin');
N = fread(fileID,1,'int');
M = fread(fileID,1,'int');
x = (fread(fileID,N,'double'))';
y = (fread(fileID,M,'double'))';
Data = (fread(fileID,[N,M],'double'));
Source = (fread(fileID,[N,M],'double'));
Forward = (fread(fileID,[N,M],'double'));
Surface = (fread(fileID,[N,M],'double'));
Speed = (fread(fileID,[N,M],'double'));


fclose(fileID);

figure
surf(Surface);
title('Level Curve');

figure
contour(Source,[0,0], 'r')
hold on;
contour(Surface,[0,0],'b')
legend('Source','Reconstruction');
title('Reconstruction of Source');
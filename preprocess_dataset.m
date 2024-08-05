% load("exportedTracks_09.mat")
% backgroundImagePath = '09_background.png';
% meterPerPixel = 0.0978;
% 
% xnew = exportedtracks{3}.xC;
% ynew = exportedtracks{3}.yC;
% xnew(113:500) = [];
% ynew(113:500) = [];
% 
% exportedtracks{3}.xC = xnew;
% exportedtracks{3}.yC = ynew;
% figure
% hold on
% plot(exportedtracks{3}.xC, exportedtracks{3}.yC)
% plot(xnew, ynew)
% scatter(exportedtracks{3}.xC(113), exportedtracks{3}.yC(113))
% scatter(exportedtracks{3}.xC(493), exportedtracks{3}.yC(493))

% load("exportedTracks_17.mat")
% backgroundImagePath = '17_background.png';
% meterPerPixel = 0.1016;
% for i = 1:1:length(exportedtracks)
%     exportedtracks{i}.yC = -exportedtracks{i}.yC;  
% end

load("exportedTracks_00.mat")
 backgroundImagePath = '00_background.png';
meterPerPixel = 0.1481;
for i = 1:1:length(exportedtracks)
     exportedtracks{i}.yC = -exportedtracks{i}.yC;  
end
xnew = exportedtracks{3}.xC;
 ynew = exportedtracks{3}.yC;
% xnew(207:314) = [];
%ynew(207:314) = [];
figure
hold on
plot(exportedtracks{3}.xC, exportedtracks{3}.yC)
plot(xnew, ynew)

exportedtracks{3}.xC = xnew;
exportedtracks{3}.yC = ynew;
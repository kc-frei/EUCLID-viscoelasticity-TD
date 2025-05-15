function [x11,x12,y11,y12,x21,x22,y21,y22] = middlePoint(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8)
x11 = (x2-x1)/2+x1;
y11 = (y2-y1)/2+y1;

x12 = (x4-x3)/2+x3;
y12 = (y4-y3)/2+y3;

x21 = (x6-x5)/2+x5;
y21 = (y6-y5)/2+y5;

x22 = (x8-x7)/2+x7;
y22 = (y8-y7)/2+y7;

% figure(3)
% hold on;
% plot(x1,y1,".r");
% plot(x2,y2,".b");
% plot(x3,y3,".m");
% plot(x4,y4,".k");
% plot(x5,y5,"xr");
% plot(x6,y6,"xb");
% plot(x7,y7,"xm");
% plot(x8,y8,"xk");
% 
% plot(x11,y11,"or");
% plot(x12,y12,"ob");
% plot(x21,y21,"om");
% plot(x22,y22,"ok");
end


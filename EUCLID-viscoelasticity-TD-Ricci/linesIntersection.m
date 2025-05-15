function [xIntersection,yIntersection] = linesIntersection(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8)
[x11,x12,y11,y12,x21,x22,y21,y22] = middlePoint(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8);
m1 = (y12-y11)/(x12-x11);
m2 = (y22-y21)/(x22-x21);
c1 = y11 - m1*x11;
c2 = y21 - m2*x21;

a1 = (c1*y11 - c1*y12) / (x11*y12 - x12*y11);
b1 = (a1*x12 + c1) / y12;

a2 = (c2*y21 - c2*y22) / (x21*y22 - x22*y21);
b2 = (a2*x22 + c2) / y22;

xIntersection = (b1*c2 - b2*c1) / (a1*b2 - a2*b1);
yIntersection = (a1*c2 - a2*c1) / (a1*b2 - a2*b1);
end


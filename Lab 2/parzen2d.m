function parzen2d(pa, pb, pc,al,bl,cl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[x,y] = size(pa);

aPoints = zeros(x,y);
bPoints = zeros(x,y);
cPoints = zeros(x,y);


for i = 1:x
   for j = 1:y
      classAProb = pa(i,j);
      classBProb = pb(i,j);
      classCProb = pc(i,j);    
      
      if (classAProb > classBProb) && (classAProb > classCProb)
          aPoints(i,j) = 1;
      elseif (classCProb > classAProb) && (classCProb > classBProb)
          cPoints(i,j) = 1;
      elseif (classBProb > classCProb) && (classBProb > classAProb)
          bPoints(i,j) = 1;
      end
   end
end

[contA, h] = contour(aPoints);
[contB, h] = contour(bPoints);
[contC, h] = contour(cPoints);
contA = contA';
contB = contB';
contC = contC';
contA = contA/1;
contB = contB/1;
contC = contC/1;
contA(contA(:,2)>510, :) = [];
contB(contB(:,2)>510, :) = [];
contC(contC(:,2)>510, :) = [];

figure;
hold on
scatter(al(:,1),al(:,2))
scatter(bl(:,1),bl(:,2))
scatter(cl(:,1),cl(:,2))
scatter(contA(:,1),contA(:,2),1)
scatter(contB(:,1),contB(:,2),1)
scatter(contC(:,1),contC(:,2),1)
xlim([0 475])
ylim([0 475])
legend('Class A', 'Class B', 'Class C', 'Class A Boundary', 'Class B Boundary', 'Class C Boundary')
xlabel('X Value')
ylabel('Y Value')
title('Parzen Estimation of Class Boundaries')
end


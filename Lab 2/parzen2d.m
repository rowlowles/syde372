function parzen2d(pa, pb, pc,al,bl,cl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[x,y] = size(pa);

aPoints = zeros(x,y);
bPoints = zeros(x,y);
cPoints = zeros(x,y);

cbBound = [];
abBound = [];


for i = 1:x
   for j = 1:y
      classAProb = pa(i,j);
      classBProb = pb(i,j);
      classCProb = pc(i,j);
      
%       if classBProb == classCProb
%         cbBound = [cbBound; i j];
%       elseif classBProb == classAProb
%         abBound = [abBound; i j];
%       end
      
      
      if (classAProb > classBProb) && (classAProb > classCProb)
          aPoints(i,j) = 1;
      elseif (classCProb > classAProb) && (classCProb > classBProb)
          cPoints(i,j) = 1;
      elseif (classBProb > classCProb) && (classBProb > classAProb)
          bPoints(i,j) = 1;
      end
   end
end
figure;
hold on
scatter(al(:,1),al(:,2))
scatter(bl(:,1),bl(:,2))
scatter(cl(:,1),cl(:,2))
contour(aPoints)
contour(bPoints)
contour(cPoints)
end


% MSC多元散射校正
function [Xmsc] = MCS(x)
Xnir = x
[me, ] = mean(Xnir);
[m, n] = size(Xnir);
for i = 1:m    
  p = polyfit(me, Xnir(i,:),1);    
  Xmsc(i,:) = (Xnir(i,:) - p(2) * ones(1, n))./(p(1) * ones(1, n));
end
end
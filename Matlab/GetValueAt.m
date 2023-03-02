function [val] = GetValueAt(x, n, i)
% Default value is 0 which mean if it can not find, it will return 0
val = 0;
for j = 1:length(n)
   if (i == n(j))
        val = x(j);
        break;
   end
end
end
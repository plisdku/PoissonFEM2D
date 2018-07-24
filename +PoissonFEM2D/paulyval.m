function y = paulyval(p,x)

N = length(p);
siz_x = size(x);

%y = zeros(siz_x, superiorfloat(x,p));
y = zeros(siz_x);
if N>0
    y(:) = p(1);
end

for i=2:N
    y = x .* y + p(i);
end

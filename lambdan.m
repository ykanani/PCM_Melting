function lambda=lambdan(range,Bi)
xs = range(1):0.01:range(2);
ys = xs.*besselj(1,xs)./besselj(0,xs)-Bi;
fun=@(x)x.*besselj(1,x)/besselj(0,x)-Bi;
scinter = find(diff(sign(ys)));

ninter = numel(scinter);
lambda = [];
for i = 1:ninter
  [temp, fval] = fzero(fun,xs(scinter(i) + [0  1]));
  if abs(fval)<0.1
      lambda(end+1)=temp;
  end
end
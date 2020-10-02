function cmap = BBVYWcolormap(n,scalefun)

%% black-blue-violet-yellow-white colormap
% n = number of colors (default 100)
% scalefun: apply scaling function to the values btwn 0 and 1

if (nargin<1)
    n = 200;
end
if (nargin<2)
    scalefun = @(x) x;
end

vals = (linspace(0,1,n))';

vals = scalefun(vals);

red = vals/0.32-0.78125;
green = 2*vals-0.84;
blue = red;
ind1 = round(n/4);
ind2 = round(n/2);
ind3 = round(3*n/4);

blue(1:ind1) = 4*vals(1:ind1)+0.4;
blue(ind1+1:ind2) = 1;
blue(ind2+1:ind3) = -2*vals(ind2+1:ind3)+1.84;
blue(ind3+1:end) = vals(ind3+1:end)/0.08-11.5;

cmap = [red green blue];

cmap(cmap<0) = 0;
cmap(cmap>1) = 1;

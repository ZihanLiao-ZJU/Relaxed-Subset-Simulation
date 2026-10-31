function g = LSF_2D_PiecewiseLinear(x)
% Paper example A; one independent-standard-normal realization per column.
if size(x,1)~=2, error('LSF:Dimension','Example A requires two rows.'); end
g1 = .85-.1*x(1,:); mask = x(1,:)>3.5; g1(mask) = 4-x(1,mask);
g2 = 2.3-x(2,:); mask = x(2,:)>2; g2(mask) = .5-.1*x(2,mask);
g = min(g1,g2);
end

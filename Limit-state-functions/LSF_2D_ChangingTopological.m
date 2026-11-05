function g = LSF_2D_ChangingTopological(x)
% Paper example B; one independent-standard-normal realization per column.
if size(x,1)~=2, error('LSF:Dimension','Example B requires two rows.'); end
g = 30./((4*(x(1,:)+2).^2/9+x(2,:).^2/25).^2+1) ...
    +20./(((x(1,:)-2.5).^2/4+(x(2,:)-.5).^2/25).^2+1)-5;
end

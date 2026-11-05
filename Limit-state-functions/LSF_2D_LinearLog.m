function g = LSF_2D_LinearLog(x)
% Example C, paper Eq. (26), DOI: 10.1016/j.ress.2025.111302.
% Independent standard-normal inputs; samples are columns.
% Uses 5-x1 without the 0.03 scale in the earlier supplied code.
if ~ismatrix(x) || size(x,1) ~= 2
    error('LSF_2D_LinearLog:Dimension','Example C requires a 2-by-N input matrix.');
end
g = min(5-x(1,:), 1./(1+exp(-2*(x(2,:)+4)))-0.5);
end

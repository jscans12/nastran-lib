function vq = interp1_log(x,v,xq,varargin)
%INTERP1_LOG 1-D interpolation in loglog scale

    % Perform interpolation
    vq = 10.^interp1(log10(x),log10(v),log10(xq),varargin{:});

end


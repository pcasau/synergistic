function out = partialDiffFnc(fnc,mn,part,varargin)
    x = sym('x',mn);
    vecFnc = vec(fnc(x));
    out = jacobian(vecFnc,vec(x(part)));
    if ~isempty(varargin)
        matlabFunction(out,'File',varargin{1},...
            'Vars',cellfun(@(i) x(i),varargin{2},'UniformOutput',false));
    end
        
function [isSame, relErr, normDiff, normExact] = compareNorms(testVal, exactVal, varargin)
% [isSame, relErr] = compareNorms(testVal, exactVal)
% [isSame, relErr] = compareNorms(testVal, exactVal, relChange)
% [isSame, relErr] = compareNorms(testVal, exactVal, relChange, nonZero)

if nargin >= 3
    relativeChangeThreshold = varargin{1};
else
    relativeChangeThreshold = 1e-6;
end

if nargin >= 4
    argumentThreshold = varargin{2};
else
    argumentThreshold = 1e-9;
end


normTest = norm(testVal(:));
normExact = norm(exactVal(:));
normDiff = norm(testVal(:) - exactVal(:));
relErr = normDiff / normExact;

if normExact > argumentThreshold
    if relErr > relativeChangeThreshold
        isSame = 0;
    else
        isSame = 1;
    end
else
    if normTest > argumentThreshold
        isSame = 0;
    else
        isSame = 1;
    end
end

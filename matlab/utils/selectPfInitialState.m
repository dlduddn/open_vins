function x0 = selectPfInitialState(cfg, x0Align)
%SELECTPFINITIALSTATE Choose the PF initial state from alignment settings.

    if getLogical(cfg, 'initUseAlignmentForPf', true) && ...
            numel(x0Align) == 3 && all(isfinite(x0Align(:)))
        x0 = x0Align(:);
        if ~getLogical(cfg, 'initUseLongitudinalCorrection', false)
            x0(1) = 0.0;
        end
    else
        x0 = [0; 0; 0];
    end
    x0(3) = wrapAngle(x0(3));
end

function value = maxFinite(values, defaultValue)
%MAXFINITE Return the maximum positive finite value or a fallback.

    values = values(isfinite(values) & values > 0);
    if isempty(values)
        value = defaultValue;
    else
        value = max(values);
    end
end

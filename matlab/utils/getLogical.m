function value = getLogical(s, name, defaultValue)
%GETLOGICAL Read a logical config field with a default fallback.

    value = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = logical(s.(name));
    end
end

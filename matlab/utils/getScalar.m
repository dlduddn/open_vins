function value = getScalar(s, name, defaultValue)
%GETSCALAR Read a scalar-like struct field with a default fallback.

    value = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

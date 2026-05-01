function value = getStructScalar(s, name, defaultValue)
%GETSTRUCTSCALAR Read a scalar struct field with a default fallback.

    value = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name)) && isscalar(s.(name))
        value = s.(name);
    end
end

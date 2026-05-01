function value = getField(s, name, defaultValue)
%GETFIELD Read a struct field with a default fallback.

    value = defaultValue;
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function count = getRejectCount(summary, name)
%GETREJECTCOUNT Read a rejection counter from a summary struct.

    count = 0;
    if isstruct(summary) && isfield(summary, name) && ~isempty(summary.(name))
        count = summary.(name);
    end
end

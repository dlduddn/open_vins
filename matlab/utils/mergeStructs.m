function out = mergeStructs(varargin)
%MERGESTRUCTS Merge structs from left to right.

    out = struct();
    for i = 1:nargin
        names = fieldnames(varargin{i});
        for j = 1:numel(names)
            out.(names{j}) = varargin{i}.(names{j});
        end
    end
end

function varargout = subsref(tS,s)
%overloads subsref

switch s(1).type
    case '()'

        tS.eta1 = subsref(tS.eta1,s(1));
        tS.k1 = subsref(tS.k1,s(1));
        tS.rotAxis = subsref(tS.rotAxis,s(1));
        tS.CRSS = subsref(tS.CRSS,s(1));
        tS.eta2 = subsref(tS.eta2,s(1));
        tS.k2 = subsref(tS.k2,s(1));
        tS.twinType = subsref(tS.twinType,s(1));
        if ~isempty(tS.parent)
            tS.parent = subsref(tS.parent,s(1));
        end
        if numel(s)>1
            [varargout{1:nargout}] = builtin('subsref',tS,s(2:end));
        else
            varargout{1} = tS;
        end
    case '.' 
        if ismethod(tS, s(1).subs) && ismember(s(1).subs, {'variantDetermination2', 'anglefromTwin'})
             methodName = s(1).subs;
            % Process method calls with or without parameters
            args = {};
            if numel(s) > 1 && strcmp(s(2).type, '()')
                args = s(2).subs;
            end

            [results{1:nargout}] = arrayfun(@(x) x.(methodName)(args{:}), tS, 'UniformOutput', false);
            results = cellfun(@(x) x{1}, results, 'UniformOutput', false);
            % handle the output of the variantDetermination2 in
            % accordance with tS dimensions
            resultSize = size(results{1});
            varargout{1} = squeeze(reshape(cat(1, results{:}), [resultSize, size(tS)]));

        else
            [varargout{1:nargout}] = builtin('subsref',tS,s);
        end
    otherwise
        % Fallback to built-in subsref
        [varargout{1:nargout}] = builtin('subsref',tS,s);
end
end
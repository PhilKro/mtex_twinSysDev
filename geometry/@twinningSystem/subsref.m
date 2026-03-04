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
%     case '.'
%         if isprop(tS, s(1).subs)
%             [varargout{1:nargout}] = builtin('subsref',tS,s);
%             %do i need this handling later ? further indices are omitted
%             %[varargout{1:nargout}] = builtin('subsref',tS,s(2:end));
% 
%         elseif ismethod(tS, s(1).subs)
%             methodName = s(1).subs;
%             % Process method calls with or without parameters
%             args = {};
% 
%             if numel(s) > 1 && strcmp(s(2).type, '()')
%                 args = s(2).subs;
%             end
% 
%             [results{1:nargout}] = arrayfun(@(x) x.(methodName)(args{:}), tS, 'UniformOutput', false);
% 
%             results = vertcat(results{:});
%             %results = cellfun(@(x) x{1}, results, 'UniformOutput', false);
%             % %results likely a cell array at this point
%             % % Call the method with multiple outputs
%             % results = cell(nargout, 1);  % Preallocate cell for each output
%             % for i = 1:nargout
%             %     results{i} = arrayfun(@(x) feval(@(y) y.(methodName)(args{:}), x), tS, 'UniformOutput', false);
%             % end
% 
%             % Special case for 'variantDetermination2'
%             if strcmp(s(1).subs, 'variantDetermination2')
%                 %handle the output of the variantDetermination2 in
%                 %accordance with tS dimensions
%                 resultSize = size(results{1});
%                 varargout{1} = squeeze(reshape(cat(1, results{:}), [resultSize, size(tS)]));
% 
%             elseif numel(s)>1 && ~strcmp(s(2).type, '()')
%                 % handle cases like tS.parentTwinOrientation.phi1
%                 try
%                     [varargout{1:nargout}] = builtin('subsref',results{:},s(2:end));
%                 catch exception
%                     disp("tS.subsref.m: Couldn't apply the second indexing to the result of the first, please consider shortening the expression")
%                     rethrow(exception);
%                 end
%                 % and the cases with tS.fun
%             else
%                 varargout = results;
%                 % try
%                 %     varargout{1} = cat(1, results{:});
%                 % catch
%                 %     varargout{1} = results;
%                 % end
%                 %                 if numel(s)>1
%                 %     try
%                 %         outputs = cellfun(@(innerCell) innerCell{1}, results, 'UniformOutput', false);
%                 %         [varargout{1:nargout}] = outputs{1:nargout};
%                 %         %varargout{1} = cat(1, results{:});
%                 %     catch
%                 %         [varargout{1:nargout}] = results{:};
%                 %     end
%                 % else
%                 %     varargout{1} = cat(1, results{:});
%                 % end
%             end
%         else
%             error('Property or method "%s" not found.', s(1).subs);
%         end
% 
%     otherwise
%         % Fallback to built-in subsref
%         [varargout{1:nargout}] = builtin('subsref',tS,s);
% end
% end

% switch s(1).type
%     case '()'
% 
%         tS.b = subsref(tS.b,s(1));
%         tS.n = subsref(tS.n,s(1));
%         tS.eta1 = subsref(tS.eta1,s(1));
%         tS.k1 = subsref(tS.k1,s(1));
%         tS.rotAxis = subsref(tS.rotAxis,s(1));
%         tS.CRSS = subsref(tS.CRSS,s(1));
%         tS.eta2 = subsref(tS.eta2,s(1));
%         tS.k2 = subsref(tS.k2,s(1));
%         if numel(s)>1
%             [varargout{1:nargout}] = builtin('subsref',tS,s(2:end));
%         else
%             varargout{1} = tS;
%         end
%     case '.'
%         if isprop(tS, s(1).subs)
%             [varargout{1:nargout}] = builtin('subsref',tS,s);
%             %do i need this handling later ? further indices are omitted
%                         %[varargout{1:nargout}] = builtin('subsref',tS,s(2:end));
%         elseif ismethod(tS, s(1).subs)
%             methodName = s(1).subs;
%             %results = cell(size(obj));
%             if numel(s)>1
%                 switch s(2).type
%                     case '()'
%                         if ~isempty(s(2).subs)
%                             % handle method calls with parameters
%                             % tS.fun(p1,p2)
%                             results = arrayfun(@(x) x.(methodName)(s(2).subs{:}), tS, 'UniformOutput',false);
%                         else
%                             % handle method calls without parameters
%                             % tS.fun()
%                            results = arrayfun(@(x) x.(methodName)(), tS, 'UniformOutput',false);
%                         end
%                     otherwise
%                         % handle method calls w/o parentheses (and w/o
%                         % parameters) but with further indexing like
%                         % tS.fun.prop
%                         results = arrayfun(@(x) x.(methodName)(), tS, 'UniformOutput',false);
%                 end
%             else
%                 % handle method calls w/o parentheses (and w/o
%                 % parameters) and w/o further indexing like
%                 % tS.fun
%                 results = arrayfun(@(x) x.(methodName)(), tS, 'UniformOutput',false);
%             end
% 
%             %results likely a cell array at this point
% 
%             if strcmp(s(1).subs, 'variantDetermination2')
%                 %handle the output of the variantDetermination2 in
%                 %accordance with tS dimensions
%                 resultSize = size(results{1});
%                 % Concatenate along the new dimension
%                 varargout{1} = cat(1, results{:});
%                 % Reshape to have the object dimension first
%                 varargout{1} = squeeze(reshape(varargout{1}, [resultSize, size(tS)]));
%             else
%                 % handle cases like tS.parentTwinOrientation.phi1
%                 if numel(s)>1 && ~strcmp(s(2).type, '()')
%                     try
%                         results = cat(1, results{:});
%                         [varargout{1:nargout}] = builtin('subsref',results,s(2:end));
%                     catch exception
%                         disp("Couldn't apply the second type of indexing to the result of the first, please consider shortening the expression")
%                         rethrow(exception);
%                     end
%                 % and the cases with tS.fun
%                 else
%                     try
%                         varargout{1} = cat(1, results{:});
%                     catch
%                         varargout{1} = results;
%                     end
%                 end
%             end
%         else
%             error('Property or method "%s" not found.', s(1).subs);
%         end
%     otherwise
%         [varargout{1:nargout}] = builtin('subsref',tS,s);
% end
% end


% function varargout = subsref(tS, s,varargin)
%     % Overloads subsref for custom indexing behavior
% 
%     switch s(1).type
%         case '()'
%             % Apply the indexing to each property
%             tS.b = subsref(tS.b, s(1));
%             tS.n = subsref(tS.n, s(1));
%             tS.eta1 = subsref(tS.eta1, s(1));
%             tS.eta2 = subsref(tS.eta2, s(1));
%             tS.k1 = subsref(tS.k1, s(1));
%             tS.k2 = subsref(tS.k2, s(1));
%             tS.rotAxis = subsref(tS.rotAxis, s(1));
%             tS.CRSS = subsref(tS.CRSS, s(1));
% 
%             % If there are more levels of subscripting, continue
%             if numel(s) > 1
%                 % Continue subscripting the result
%                 [varargout{1:nargout}] = builtin('subsref', tS, s(2:end));
%             else
%                 varargout{1} = tS;
%             end
% 
%         case '.'
%             % Handle method calls and property access
%             methodName = s(1).subs;
%             switch methodName
%                 case 'variantDetermination'
%                     % If there are additional subscripts, handle them
%                     if numel(s) > 1 && strcmp(s(2).type, '()')
%                         % Get the method call result
%                         methodResult = tS.variantDetermination(varargin{:}); %(varargin{:});
%                         % Create a subsref struct for the remaining subscripts
%                         remainingSubs = s(2:end);
%                         % Apply the remaining subscripts to the result
%                         [varargout{1:nargout}] = subsref(methodResult, remainingSubs);
%                     else
%                         [varargout{1:nargout}] = tS.variantDetermination(tS);%(varargin{:});
%                     end
% 
%                 otherwise
%                     % For any other properties or methods, use default behavior
%                     [varargout{1:nargout}] = builtin('subsref', tS, s);
%             end
% 
%         otherwise
%             % Use the built-in subsref for other types of indexing
%             [varargout{1:nargout}] = builtin('subsref', tS, s);
%     end
% end
    % case '.'
    %     methodName = s(1).subs;
    %     switch methodName
    %         case getK1name
    %             if length(s)==2 && strcmp(s(2).type,'()')
    %                 prop = s(1).subs;      % Property name
    %                 n = numel(tS);         % Number of elements in array
    %                 varargout = cell(1,n); % Preallocate cell array
    %                 for k = 1:n
    %                     varargout{k} = tS(k).(prop).(s(2).subs);
    %                 end
    %             else
    %                 [varargout{1:nargout}] = builtin('subsref',tS,s);
    %             end
    %         otherwise
    % 
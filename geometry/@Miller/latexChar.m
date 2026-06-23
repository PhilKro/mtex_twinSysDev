function lC = latexChar(m,varargin)
    % LATEXCHAR Convert Miller indices to a LaTeX formatted character array
    %
    % This method extracts the indices of a Miller object based on its display 
    % style and constructs a LaTeX string. Negative indices are automatically 
    % formatted with a LaTeX overbar (\bar{}).
    %
    % Syntax
    %   lC = latexChar(m)
    %   lC = latexChar(m, dispStyle)
    %   lC = latexChar(m, 'round')
    %   lC = latexChar(m, 'round', 'noWarning')
    %
    % Input
    %   m         - Miller object
    %
    % Output
    %   lC        - Character array containing the LaTeX formatted indices
    %
    % Flags
    %   'hkl'     - Format as Miller indices (h, k, l)
    %   'hkil'    - Format as Miller-Bravais indices (H, K, I, L)
    %   'uvw'     - Format as direction indices (u, v, w)
    %   'UVTW'    - Format as direction indices in 4-index notation (U, V, T, W)
    %   'round'   - Rounds the indices to a maximum value of 20. If the deviation
    %               angle is greater than 10^-5, a tilde (~) is prepended.
    %   'noWarning' - Suppresses the '~' character when rounding causes a 
    %                 significant deviation in the vector angle.
    format = get_flag(varargin,{'hkl','hkil','uvw','UVTW'});
    if ~isempty(format)
        m.dispStyle = format; 
    else
        format = m.dispStyle;
    end
    
    warnChar = '';
    % Renamed to avoid shadowing the built-in round() function
    isRounded = get_flag(varargin,{'round'});
    if ~isempty(isRounded)
        m_round = m.round('maxHKL', 20); 
        
        noWarning = get_flag(varargin,{'noWarning'});
        if ~isempty(noWarning)
            m = m_round;
        else
            if angle(m_round, m) < 10^-5
                m = m_round;
            else
                m = m_round;
                warnChar = '~';
            end
        end
    end
    
    % Extract the proper components and define the bounding brackets
    switch format
        case 'hkl'
            vals = [m.h, m.k, m.l];
            leftBracket = '(';
            rightBracket = ')';
        case 'hkil'
            vals = [m.h, m.k, m.i, m.l];
            leftBracket = '(';
            rightBracket = ')';
        case 'uvw'
            vals = [m.u, m.v, m.w];
            leftBracket = '[';
            rightBracket = ']';
        case 'UVTW'
            vals = [m.U, m.V, m.T, m.W];
            leftBracket = '[';
            rightBracket = ']';
        otherwise
            error('Miraculously no dispStyle is defined for your Miller object')
    end
            
    % Build the LaTeX formatted string
    lC_core = '';
    for i = 1:length(vals)
        val = vals(i);
        
        % Use absolute value for the tolerance check so negative indices are not destroyed
        if abs(val) < 10^-5
            val = 0;
        end
        
        % Add a slight LaTeX space (\,) between individual indices
        if i > 1
            lC_core = [lC_core, '\,'];
        end
        
        if val < 0
            % Use \bar{} for negative values, drop the minus sign with abs()
            lC_core = [lC_core, sprintf('\\bar{%g}', abs(val))];
        else
            % Standard string conversion for positive values and zero
            lC_core = [lC_core, sprintf('%g', val)];
        end
    end
    
    % Combine the warning character, brackets, and the core string
    lC = [warnChar, leftBracket, lC_core, rightBracket];
end
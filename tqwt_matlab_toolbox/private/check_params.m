function check_params(Q,r,J)

[st,i] = dbstack(1);

% st.file - file that called this one.


% ---------- Check Q ----------

if ~isscalar(Q)
    error('Error in %s: Q must be scalar \n', st.file)
end

if ~isnumeric(Q)
    error('Error in %s: Q must be numeric \n', st.file)
end

if Q < 1
    error('Error in %s: Q must be greater than or equal to 1.0 \n', st.file)
end


% ---------- Check r ----------

if ~isscalar(r)
    error('Error in %s: r must be scalar \n', st.file)
end

if ~isnumeric(r)
    error('Error in %s: r must be numeric \n', st.file)
end

if ~isfinite(r)
    error('Error in %s: r must be greater than 1.0 \n', st.file)
end

if r <= 1
    error('Error in %s: r must be greater than 1.0 \n', st.file)
end


% ---------- Check J ----------
if nargin > 2
    
    if ~isscalar(J)
        error('Error in %s: J must be scalar \n', st.file)
    end
    
    if ~isnumeric(J)
        error('Error in %s: J must be numeric \n', st.file)
    end
    
    if ~isfinite(J)
        error('Error in %s: J must be a positive integer \n', st.file)
    end
    
    if round(J) ~= J
        error('Error in %s: J must be a positive integer \n', st.file)
    end
    
    if J < 1
        error('Error in %s: J must be at positive integer \n', st.file)
    end
    
end



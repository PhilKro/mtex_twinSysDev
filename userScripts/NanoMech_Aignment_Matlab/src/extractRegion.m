function sub = extractRegion(ar, x0,y0,x1,y1, dsPath)
        s = size(ar);
        if y1 > s(1) || x1 > s(2)
            error('Requested region [%d:%d, %d:%d] exceeds dataset bounds [%d x %d] for %s', ...
                y0, y1, x0, x1, s(1), s(2), dsPath);
        end
        sub = ar(y0:y1, x0:x1);
    end
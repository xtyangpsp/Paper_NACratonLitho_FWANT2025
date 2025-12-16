function [vdrop,firstneg,firstpos,idxtarget,vtarget]=find_firstdrop(v)
    % Function to find the first drop in value and return the drop and the
    % indices of the top and bottom around the target location.
    % The target location of the first drop is between the peak and the trough with 
    % the minimum negative gradient or fastest decrease.

    vgrad=diff(v); 
    vgrad(end+1)=0;
    
    % Find the first negative gradient where the value starts to decrease.
    idxneg=find(vgrad < 0);
    
    % If there is no decrease in values, return NaN for all outputs.
    if isempty(idxneg)
        vdrop=nan;
        firstneg=nan;
        firstpos=nan;
        idxtarget=nan;
        vtarget=nan;
        return
    end
    firstneg=idxneg(1);
    
    idxpos=find(vgrad(firstneg+1:end)>0);
    if isempty(idxpos) %no increase: find the minimum gradient below the first negative
    %gradient
        [~,idxmingrad]=min(vgrad(firstneg+1:end));%find(vgrad < 0);
        firstpos=length(v);
        vpos=v(firstpos);
        vdrop = v(firstneg) - vpos;
        idxtarget=idxmingrad + firstneg;
    else %has an increase below the first peak
        firstpos=firstneg + idxpos(1) -1;
        [~,idxmingrad]=min(vgrad(firstneg+1:firstpos));
        vpos=v(firstpos);
        vdrop=v(firstneg) - vpos;
        idxtarget=idxmingrad + firstneg;
    end
    vtarget=v(idxtarget);
end
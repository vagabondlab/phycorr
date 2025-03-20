function [new_locs] = correct_trigger_peaks(locs,Trigger,vols,TR,sr)

%verify that locs is correct sitrig_ze
if numel(locs) ~= vols
   disp('wrong number of locs') 
else
    new_locs= locs;
    check_range = -.2*TR*sr:.2*TR*sr;
    
    for i=1:numel(locs)
        p = locs(i);
        
        [~,adjust] = max(Trigger(check_range + p));
        new_p = p + check_range(adjust);
        new_locs(i) = new_p;
    end
end

end
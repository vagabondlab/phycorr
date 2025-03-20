function [ rv ] = RVTestimate(respwave,TR,maxvol,fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
t_win = 6;  % default choice to have 6s

if size(respwave,2) > size(respwave,1)
    respwave = respwave';
end

timevec = 0:TR:(maxvol*TR - TR);
rv = zeros(maxvol,1);
for tp = 1:maxvol
    i1 = max(1,floor((timevec(tp) - t_win) * fs));
    i2 = min(length(respwave),floor((timevec(tp) + t_win) * fs));
    if i2 < i1
        disp('Resp data shorter than scan duration')
        rv(tp) = mean(rv(1:tp-1));
    else
        rv(tp) = std(respwave(i1:i2));
    end
end
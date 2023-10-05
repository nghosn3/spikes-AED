%% function that findns the estimated initial blood plasma concentration at the start of EMU stay based on 
%% the patients home medications and dosing schedule,
function [c0,c_t] = get_c0_min_dose(dose,frequency,tHalf,F,vd,tmax,ka,min_dose)


%% use get_single_dose_curve to model the steady state blood levels using the dosing schedule they are on 
c0=0;

% figure out the t (hours) needed for steady state - iterate for one week
% to start (more than most drugs pharm reports say)

tInt = round(24./frequency); % estimated interval for single dose based on frequency
for i = 1:frequency*7
    [c_t,~] = get_single_dose_curve_min_dose(c0,dose,tHalf,tInt,F,vd,tmax,ka,min_dose);
    c0 = c_t(end);
end

end

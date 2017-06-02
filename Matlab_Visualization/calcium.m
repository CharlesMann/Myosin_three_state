time = 0:.001:.17;
pCa = zeros(1,numel(time));
Ca = zeros(1,numel(time));
t_act = 0.0718;
t_p = t_act + .001;
for i = 1:numel(time)
    if time(i) >= t_act
        if time(i) >= t_p
            pCa(i) = 0.5*exp(-((time(i)-t_p)*25)^2);
        else
            pCa(i) = (time(i)-t_act)/.02;
        end
    Ca(i) = 0.1 + 1000*sin(3.14*pCa(i));
    
    else
        Ca(i) = 0.1;
    end
end

plot(time,Ca)
title('Calcium Concentration')
xlabel('Time')
ylabel('[Ca] in nM')
% plot(time, pCa)
function [a1, a2, a3, a_on_rate] = actin_check()
Ca = 10^(-4.5);
N_active = 0;
N_overlap = 0.9482;
N_a_bound = 0;
k_thin_coop = 1;
a_on = 10^6;
a_off = 10;



for i = 1:10
    if i < 6
        Ca = 0;
    else
        Ca = 10^(-4.5);
    end
    a1(i) = (N_active/N_overlap);
    a2(i) = 1 + (k_thin_coop*a1(i));
    a3(i) = N_overlap - N_active;
    a_on_rate(i) = a_on*Ca*a3(i)*a2(i);

    ao1 = N_active - N_a_bound;
    ao2 = (N_overlap-N_active)/N_overlap;
    ao3 = ao1*((1+ao2));
    a_off_rate(i) = a_off*ao1*ao3;

    N_a_active(i) = (a_on_rate(i) - a_off_rate(i))*.001;
    N_active = N_a_active(i);
    
end
    
    
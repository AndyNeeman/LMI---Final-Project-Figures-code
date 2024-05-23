function Fig_3b
    Delta = 0;
    delta = 0.3*pi;
    c_0 = [1; 0; 0];
    function ab = coeff(t,A,phi)
        Omega_0 = A*exp(-1j*phi);
        Omega_R = sqrt(Delta.^2 + abs(Omega_0).^2);
        ab(1) = cos(t*pi/2*Omega_R/(pi)) - 1j * Delta/Omega_R * sin(t*pi/2*Omega_R/(pi));
        ab(2) = 1j * Omega_0/Omega_R * sin(t*pi/2*Omega_R/(pi));
        ab = [ab(1), ab(2)];
    end
    
    function U = propagator(t,A,phi)
    ab = coeff(t,A,phi);
    a = ab(1);
    b = ab(2);
    U = [sqrt(4)/2*a.^2, sqrt(2)/1*a*b, sqrt(4)/2*b.^2;
        -sqrt(2)/1*a*conj(b), sqrt(1)/1*a*conj(a) - sqrt(1)/1*b*conj(b), sqrt(2)/1*conj(a)*b;
        sqrt(4)/2*conj(b).^2, -sqrt(2)/1*conj(a)*conj(b), sqrt(4)/2*conj(a).^2];
    end
    
    function c = propagation(t)
        if t <=1 
            c = propagator(t,pi+delta,0)*c_0;
        elseif t<=2
            c = propagator(t-1,pi+delta,4/5*pi)*propagator(1,pi+delta,0)*c_0;
        elseif t<=3
            c = propagator(t-2,pi+delta,2/5*pi)*propagator(1,pi+delta,4/5*pi)*propagator(1,pi+delta,0)*c_0;
        elseif t<=4
            c = propagator(t-3,pi+delta,4/5*pi)*propagator(1,pi+delta,2/5*pi)*propagator(1,pi+delta,4/5*pi)*propagator(1,pi+delta,0)*c_0;
        else
            c = propagator(t-4,pi+delta,0)*propagator(1,pi+delta,4/5*pi)*propagator(1,pi+delta,2/5*pi)*propagator(1,pi+delta,4/5*pi)*propagator(1,pi+delta,0)*c_0;
        end
    end
    c = [];
    t_interval = linspace(0,5,400);
    for t = t_interval
        c = [c propagation(t)];
    end
    plot(t_interval, (abs(c(1,:))).^2, '-', t_interval, (abs(c(2,:))).^2, '-', t_interval, (abs(c(3,:))).^2, '-');
    ylim([0, 1]);
    xlim([0,5]);
    legend('1','2','3')
end
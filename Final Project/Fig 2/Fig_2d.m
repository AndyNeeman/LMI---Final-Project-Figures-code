function Fig_2d
    Delta = 0;
    delta = 0;
    c_0 = [0; 1; 0; 0; 0];
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
    U = [sqrt(576)/24*a.^4, sqrt(144)/6*a.^3*b, sqrt(96)/4*a.^2*b.^2, sqrt(144)/6*a*b.^3, sqrt(576)/24*b.^4;
        -sqrt(144)/6*a.^3*conj(b), sqrt(36)/6*a.^3*conj(a) - sqrt(36)/2*a.^2*b*conj(b), sqrt(24)/2*a.^2*conj(a)*b - sqrt(24)/2*a*b.^2*conj(b), sqrt(36)/2*a*conj(a)*b.^2 - sqrt(36)/6*b.^3*conj(b), sqrt(144)/6*conj(a)*b.^3;
        sqrt(96)/4*a.^2*conj(b).^2, -sqrt(24)/2*a.^2*conj(a)*conj(b) + sqrt(24)/2*a*b*conj(b).^2, sqrt(16)/4*a.^2*conj(a).^2 - sqrt(16)/1*a*conj(a)*b*conj(b) + sqrt(16)/4*b.^2*conj(b).^2, sqrt(24)/2*a*conj(a).^2*b - sqrt(24)/2*conj(a)*b.^2*conj(b), sqrt(96)/4*conj(a).^2*b.^2;
        -sqrt(144)/6*a*conj(b).^3, sqrt(36)/2*a*conj(a)*conj(b).^2 - sqrt(36)/6*b*conj(b).^3, -sqrt(24)/2*a*conj(a).^2*conj(b) + sqrt(24)/2*conj(a)*b*conj(b).^2, sqrt(36)/6*a*conj(a).^3 - sqrt(36)/2*conj(a).^2*b*conj(b), sqrt(144)/6*conj(a).^3*b;
        sqrt(576)/24*conj(b).^4, -sqrt(144)/6*conj(a)*conj(b).^3, sqrt(96)/4*conj(a).^2*conj(b).^2, -sqrt(144)/6*conj(a).^3*conj(b), sqrt(576)/24*conj(a).^4];
    end
    
    function c = propagation(t)
        if t <=1 
            c = propagator(t,pi/2+delta,0)*c_0;
        elseif t<=2
            c = propagator(t-1,pi,pi/2)*propagator(1,pi/2+delta,0)*c_0;
        else
            c= propagator(t-2,pi/2+delta,0)*propagator(1,pi,pi/2)*propagator(1,pi/2+delta,0)*c_0;
        end
    end
    c = [];
    t_interval = linspace(0,3,200);
    for t = t_interval
        c = [c propagation(t)];
    end
    plot(t_interval, (abs(c(1,:))).^2, '-', t_interval, (abs(c(2,:))).^2, '-', t_interval, (abs(c(3,:))).^2, '-', t_interval, (abs(c(4,:))).^2, '-', t_interval, (abs(c(5,:))).^2, '-');
    ylim([0, 1]);
    legend('1','2','3','4','5')
end
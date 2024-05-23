function Fig_6_complete

    function Fig_6a
        Delta = 0;
        delta = 0.3*pi;
        c_0 = [0; 1; 0];
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
        plot(t_interval, (abs(c(1,:))).^2, '-', t_interval, (abs(c(2,:))).^2, '-', t_interval, (abs(c(3,:))).^2, '--');
        ylim([0, 1]);
        xlim([0,5]);
        legend('$\left|1\right>$', '$\left|2\right>$', '$\left|3\right>$', 'Interpreter','latex')
    end

    function Fig_6b
        Delta = 0;
        delta = 0.3*pi;
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
        plot(t_interval, (abs(c(1,:))).^2, '-', t_interval, (abs(c(2,:))).^2, '-', t_interval, (abs(c(3,:))).^2, '-', t_interval, (abs(c(4,:))).^2, '-', t_interval, (abs(c(5,:))).^2, '-');
        ylim([0, 1]);
        xlim([0,5]);
        legend('$\left|1\right>$', '$\left|2\right>$', '$\left|3\right>$', '$\left|4\right>$', '$\left|5\right>$', 'Interpreter','latex')
    end

    function Fig_6c
        Delta = 0;
        delta = 0.4*pi;
        c_0 = [0; 1; 0; 0; 0; 0; 0; 0; 0];
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
            U = [sqrt(1625702400)/40320*a.^8, sqrt(203212800)/5040*a.^7*b, sqrt(58060800)/1440*a.^6*b.^2, sqrt(29030400)/720*a.^5*b.^3, sqrt(23224320)/576*a.^4*b.^4, sqrt(29030400)/720*a.^3*b.^5, sqrt(58060800)/1440*a.^2*b.^6, sqrt(203212800)/5040*a*b.^7, sqrt(1625702400)/40320*b.^8;
                -sqrt(203212800)/5040*a.^7*conj(b), sqrt(25401600)/5040*a.^7*conj(a) - sqrt(25401600)/720*a.^6*b*conj(b), sqrt(7257600)/720*a.^6*conj(a)*b - sqrt(7257600)/240*a.^5*b.^2*conj(b), sqrt(3628800)/240*a.^5*conj(a)*b.^2 - sqrt(3628800)/144*a.^4*b.^3*conj(b), sqrt(2903040)/144*a.^4*conj(a)*b.^3 - sqrt(2903040)/144*a.^3*b.^4*conj(b), sqrt(3628800)/144*a.^3*conj(a)*b.^4 - sqrt(3628800)/240*a.^2*b.^5*conj(b), sqrt(7257600)/240*a.^2*conj(a)*b.^5 - sqrt(7257600)/720*a*b.^6*conj(b), sqrt(25401600)/720*a*conj(a)*b.^6 - sqrt(25401600)/5040*b.^7*conj(b), sqrt(203212800)/5040*conj(a)*b.^7;
                sqrt(58060800)/1440*a.^6*conj(b).^2, -sqrt(7257600)/720*a.^6*conj(a)*conj(b) + sqrt(7257600)/240*a.^5*b*conj(b).^2, sqrt(2073600)/1440*a.^6*conj(a).^2 - sqrt(2073600)/120*a.^5*conj(a)*b*conj(b) + sqrt(2073600)/96*a.^4*b.^2*conj(b).^2, sqrt(1036800)/240*a.^5*conj(a).^2*b - sqrt(1036800)/48*a.^4*conj(a)*b.^2*conj(b) + sqrt(1036800)/72*a.^3*b.^3*conj(b).^2, sqrt(829440)/96*a.^4*conj(a).^2*b.^2 - sqrt(829440)/36*a.^3*conj(a)*b.^3*conj(b) + sqrt(829440)/96*a.^2*b.^4*conj(b).^2, sqrt(1036800)/72*a.^3*conj(a).^2*b.^3 - sqrt(1036800)/48*a.^2*conj(a)*b.^4*conj(b) + sqrt(1036800)/240*a*b.^5*conj(b).^2, sqrt(2073600)/96*a.^2*conj(a).^2*b.^4 - sqrt(2073600)/120*a*conj(a)*b.^5*conj(b) + sqrt(2073600)/1440*b.^6*conj(b).^2, sqrt(7257600)/240*a*conj(a).^2*b.^5 - sqrt(7257600)/720*conj(a)*b.^6*conj(b), sqrt(58060800)/1440*conj(a).^2*b.^6;
                -sqrt(29030400)/720*a.^5*conj(b).^3, sqrt(3628800)/240*a.^5*conj(a)*conj(b).^2 - sqrt(3628800)/144*a.^4*b*conj(b).^3, -sqrt(1036800)/240*a.^5*conj(a).^2*conj(b) + sqrt(1036800)/48*a.^4*conj(a)*b*conj(b).^2 - sqrt(1036800)/72*a.^3*b.^2*conj(b).^3, sqrt(518400)/720*a.^5*conj(a).^3 - sqrt(518400)/48*a.^4*conj(a).^2*b*conj(b) + sqrt(518400)/24*a.^3*conj(a)*b.^2*conj(b).^2 - sqrt(518400)/72*a.^2*b.^3*conj(b).^3, sqrt(414720)/144*a.^4*conj(a).^3*b - sqrt(414720)/24*a.^3*conj(a).^2*b.^2*conj(b) + sqrt(414720)/24*a.^2*conj(a)*b.^3*conj(b).^2 - sqrt(414720)/144*a*b.^4*conj(b).^3, sqrt(518400)/72*a.^3*conj(a).^3*b.^2 - sqrt(518400)/24*a.^2*conj(a).^2*b.^3*conj(b) + sqrt(518400)/48*a*conj(a)*b.^4*conj(b).^2 - sqrt(518400)/720*b.^5*conj(b).^3, sqrt(1036800)/72*a.^2*conj(a).^3*b.^3 - sqrt(1036800)/48*a*conj(a).^2*b.^4*conj(b) + sqrt(1036800)/240*conj(a)*b.^5*conj(b).^2, sqrt(3628800)/144*a*conj(a).^3*b.^4 - sqrt(3628800)/240*conj(a).^2*b.^5*conj(b), sqrt(29030400)/720*conj(a).^3*b.^5;
                sqrt(23224320)/576*a.^4*conj(b).^4, -sqrt(2903040)/144*a.^4*conj(a)*conj(b).^3 + sqrt(2903040)/144*a.^3*b*conj(b).^4, sqrt(829440)/96*a.^4*conj(a).^2*conj(b).^2 - sqrt(829440)/36*a.^3*conj(a)*b*conj(b).^3 + sqrt(829440)/96*a.^2*b.^2*conj(b).^4, -sqrt(414720)/144*a.^4*conj(a).^3*conj(b) + sqrt(414720)/24*a.^3*conj(a).^2*b*conj(b).^2 - sqrt(414720)/24*a.^2*conj(a)*b.^2*conj(b).^3 + sqrt(414720)/144*a*b.^3*conj(b).^4, sqrt(331776)/576*a.^4*conj(a).^4 - sqrt(331776)/36*a.^3*conj(a).^3*b*conj(b) + sqrt(331776)/16*a.^2*conj(a).^2*b.^2*conj(b).^2 - sqrt(331776)/36*a*conj(a)*b.^3*conj(b).^3 + sqrt(331776)/576*b.^4*conj(b).^4, sqrt(414720)/144*a.^3*conj(a).^4*b - sqrt(414720)/24*a.^2*conj(a).^3*b.^2*conj(b) + sqrt(414720)/24*a*conj(a).^2*b.^3*conj(b).^2 - sqrt(414720)/144*conj(a)*b.^4*conj(b).^3, sqrt(829440)/96*a.^2*conj(a).^4*b.^2 - sqrt(829440)/36*a*conj(a).^3*b.^3*conj(b) + sqrt(829440)/96*conj(a).^2*b.^4*conj(b).^2, sqrt(2903040)/144*a*conj(a).^4*b.^3 - sqrt(2903040)/144*conj(a).^3*b.^4*conj(b), sqrt(23224320)/576*conj(a).^4*b.^4;
                -sqrt(29030400)/720*a.^3*conj(b).^5, sqrt(3628800)/144*a.^3*conj(a)*conj(b).^4 - sqrt(3628800)/240*a.^2*b*conj(b).^5, -sqrt(1036800)/72*a.^3*conj(a).^2*conj(b).^3 + sqrt(1036800)/48*a.^2*conj(a)*b*conj(b).^4 - sqrt(1036800)/240*a*b.^2*conj(b).^5, sqrt(518400)/72*a.^3*conj(a).^3*conj(b).^2 - sqrt(518400)/24*a.^2*conj(a).^2*b*conj(b).^3 + sqrt(518400)/48*a*conj(a)*b.^2*conj(b).^4 - sqrt(518400)/720*b.^3*conj(b).^5, -sqrt(414720)/144*a.^3*conj(a).^4*conj(b) + sqrt(414720)/24*a.^2*conj(a).^3*b*conj(b).^2 - sqrt(414720)/24*a*conj(a).^2*b.^2*conj(b).^3 + sqrt(414720)/144*conj(a)*b.^3*conj(b).^4, sqrt(518400)/720*a.^3*conj(a).^5 - sqrt(518400)/48*a.^2*conj(a).^4*b*conj(b) + sqrt(518400)/24*a*conj(a).^3*b.^2*conj(b).^2 - sqrt(518400)/72*conj(a).^2*b.^3*conj(b).^3, sqrt(1036800)/240*a.^2*conj(a).^5*b - sqrt(1036800)/48*a*conj(a).^4*b.^2*conj(b) + sqrt(1036800)/72*conj(a).^3*b.^3*conj(b).^2, sqrt(3628800)/240*a*conj(a).^5*b.^2 - sqrt(3628800)/144*conj(a).^4*b.^3*conj(b), sqrt(29030400)/720*conj(a).^5*b.^3;
                sqrt(58060800)/1440*a.^2*conj(b).^6, -sqrt(7257600)/240*a.^2*conj(a)*conj(b).^5 + sqrt(7257600)/720*a*b*conj(b).^6, sqrt(2073600)/96*a.^2*conj(a).^2*conj(b).^4 - sqrt(2073600)/120*a*conj(a)*b*conj(b).^5 + sqrt(2073600)/1440*b.^2*conj(b).^6, -sqrt(1036800)/72*a.^2*conj(a).^3*conj(b).^3 + sqrt(1036800)/48*a*conj(a).^2*b*conj(b).^4 - sqrt(1036800)/240*conj(a)*b.^2*conj(b).^5, sqrt(829440)/96*a.^2*conj(a).^4*conj(b).^2 - sqrt(829440)/36*a*conj(a).^3*b*conj(b).^3 + sqrt(829440)/96*conj(a).^2*b.^2*conj(b).^4, -sqrt(1036800)/240*a.^2*conj(a).^5*conj(b) + sqrt(1036800)/48*a*conj(a).^4*b*conj(b).^2 - sqrt(1036800)/72*conj(a).^3*b.^2*conj(b).^3, sqrt(2073600)/1440*a.^2*conj(a).^6 - sqrt(2073600)/120*a*conj(a).^5*b*conj(b) + sqrt(2073600)/96*conj(a).^4*b.^2*conj(b).^2, sqrt(7257600)/720*a*conj(a).^6*b - sqrt(7257600)/240*conj(a).^5*b.^2*conj(b), sqrt(58060800)/1440*conj(a).^6*b.^2;
                -sqrt(203212800)/5040*a*conj(b).^7, sqrt(25401600)/720*a*conj(a)*conj(b).^6 - sqrt(25401600)/5040*b*conj(b).^7, -sqrt(7257600)/240*a*conj(a).^2*conj(b).^5 + sqrt(7257600)/720*conj(a)*b*conj(b).^6, sqrt(3628800)/144*a*conj(a).^3*conj(b).^4 - sqrt(3628800)/240*conj(a).^2*b*conj(b).^5, -sqrt(2903040)/144*a*conj(a).^4*conj(b).^3 + sqrt(2903040)/144*conj(a).^3*b*conj(b).^4, sqrt(3628800)/240*a*conj(a).^5*conj(b).^2 - sqrt(3628800)/144*conj(a).^4*b*conj(b).^3, -sqrt(7257600)/720*a*conj(a).^6*conj(b) + sqrt(7257600)/240*conj(a).^5*b*conj(b).^2, sqrt(25401600)/5040*a*conj(a).^7 - sqrt(25401600)/720*conj(a).^6*b*conj(b), sqrt(203212800)/5040*conj(a).^7*b;
                sqrt(1625702400)/40320*conj(b).^8, -sqrt(203212800)/5040*conj(a)*conj(b).^7, sqrt(58060800)/1440*conj(a).^2*conj(b).^6, -sqrt(29030400)/720*conj(a).^3*conj(b).^5, sqrt(23224320)/576*conj(a).^4*conj(b).^4, -sqrt(29030400)/720*conj(a).^5*conj(b).^3, sqrt(58060800)/1440*conj(a).^6*conj(b).^2, -sqrt(203212800)/5040*conj(a).^7*conj(b), sqrt(1625702400)/40320*conj(a).^8];
        end

        function c = propagation(t,A)
            c = c_0;
            phase_coeff = [0, 14, 12, 24, 20, 30, 24, 32, 24, 30, 20, 24, 12, 14, 0];
            N = length(phase_coeff);
            for i = linspace(1,floor(t), floor(t))
                c = propagator(1, A, phase_coeff(i)*pi/N)*c;
            end
            if t ~= floor(t)
                c = propagator(t-floor(t), A, phase_coeff((floor(t)+1))*pi/N)*c;
            end
        end

        c = [];
        t_interval = linspace(0,15,400);
        for t = t_interval
            c = [c propagation(t,pi+delta)];
        end
        plot(t_interval, (abs(c(1,:))).^2, '-', t_interval, (abs(c(2,:))).^2, '--', t_interval, (abs(c(3,:))).^2, '-', t_interval, (abs(c(4,:))).^2, '-', t_interval, (abs(c(5,:))).^2, '-', t_interval, (abs(c(6,:))).^2, '-', t_interval, (abs(c(7,:))).^2, '-', t_interval, (abs(c(8,:))).^2, '--', t_interval, (abs(c(9,:))).^2, '-');
        ylim([0, 1]);
        legend('', '$\left|2\right>$', '', '', '', '', '', '$\left|8\right>$', '', 'Interpreter','latex','Location','north')
    end

    function Fig_6d
        Delta = 0;
        delta = 0.28*pi;
        c_0 = [0; 0; 1; 0; 0; 0; 0; 0; 0];
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
            U = [sqrt(1625702400)/40320*a.^8, sqrt(203212800)/5040*a.^7*b, sqrt(58060800)/1440*a.^6*b.^2, sqrt(29030400)/720*a.^5*b.^3, sqrt(23224320)/576*a.^4*b.^4, sqrt(29030400)/720*a.^3*b.^5, sqrt(58060800)/1440*a.^2*b.^6, sqrt(203212800)/5040*a*b.^7, sqrt(1625702400)/40320*b.^8;
                -sqrt(203212800)/5040*a.^7*conj(b), sqrt(25401600)/5040*a.^7*conj(a) - sqrt(25401600)/720*a.^6*b*conj(b), sqrt(7257600)/720*a.^6*conj(a)*b - sqrt(7257600)/240*a.^5*b.^2*conj(b), sqrt(3628800)/240*a.^5*conj(a)*b.^2 - sqrt(3628800)/144*a.^4*b.^3*conj(b), sqrt(2903040)/144*a.^4*conj(a)*b.^3 - sqrt(2903040)/144*a.^3*b.^4*conj(b), sqrt(3628800)/144*a.^3*conj(a)*b.^4 - sqrt(3628800)/240*a.^2*b.^5*conj(b), sqrt(7257600)/240*a.^2*conj(a)*b.^5 - sqrt(7257600)/720*a*b.^6*conj(b), sqrt(25401600)/720*a*conj(a)*b.^6 - sqrt(25401600)/5040*b.^7*conj(b), sqrt(203212800)/5040*conj(a)*b.^7;
                sqrt(58060800)/1440*a.^6*conj(b).^2, -sqrt(7257600)/720*a.^6*conj(a)*conj(b) + sqrt(7257600)/240*a.^5*b*conj(b).^2, sqrt(2073600)/1440*a.^6*conj(a).^2 - sqrt(2073600)/120*a.^5*conj(a)*b*conj(b) + sqrt(2073600)/96*a.^4*b.^2*conj(b).^2, sqrt(1036800)/240*a.^5*conj(a).^2*b - sqrt(1036800)/48*a.^4*conj(a)*b.^2*conj(b) + sqrt(1036800)/72*a.^3*b.^3*conj(b).^2, sqrt(829440)/96*a.^4*conj(a).^2*b.^2 - sqrt(829440)/36*a.^3*conj(a)*b.^3*conj(b) + sqrt(829440)/96*a.^2*b.^4*conj(b).^2, sqrt(1036800)/72*a.^3*conj(a).^2*b.^3 - sqrt(1036800)/48*a.^2*conj(a)*b.^4*conj(b) + sqrt(1036800)/240*a*b.^5*conj(b).^2, sqrt(2073600)/96*a.^2*conj(a).^2*b.^4 - sqrt(2073600)/120*a*conj(a)*b.^5*conj(b) + sqrt(2073600)/1440*b.^6*conj(b).^2, sqrt(7257600)/240*a*conj(a).^2*b.^5 - sqrt(7257600)/720*conj(a)*b.^6*conj(b), sqrt(58060800)/1440*conj(a).^2*b.^6;
                -sqrt(29030400)/720*a.^5*conj(b).^3, sqrt(3628800)/240*a.^5*conj(a)*conj(b).^2 - sqrt(3628800)/144*a.^4*b*conj(b).^3, -sqrt(1036800)/240*a.^5*conj(a).^2*conj(b) + sqrt(1036800)/48*a.^4*conj(a)*b*conj(b).^2 - sqrt(1036800)/72*a.^3*b.^2*conj(b).^3, sqrt(518400)/720*a.^5*conj(a).^3 - sqrt(518400)/48*a.^4*conj(a).^2*b*conj(b) + sqrt(518400)/24*a.^3*conj(a)*b.^2*conj(b).^2 - sqrt(518400)/72*a.^2*b.^3*conj(b).^3, sqrt(414720)/144*a.^4*conj(a).^3*b - sqrt(414720)/24*a.^3*conj(a).^2*b.^2*conj(b) + sqrt(414720)/24*a.^2*conj(a)*b.^3*conj(b).^2 - sqrt(414720)/144*a*b.^4*conj(b).^3, sqrt(518400)/72*a.^3*conj(a).^3*b.^2 - sqrt(518400)/24*a.^2*conj(a).^2*b.^3*conj(b) + sqrt(518400)/48*a*conj(a)*b.^4*conj(b).^2 - sqrt(518400)/720*b.^5*conj(b).^3, sqrt(1036800)/72*a.^2*conj(a).^3*b.^3 - sqrt(1036800)/48*a*conj(a).^2*b.^4*conj(b) + sqrt(1036800)/240*conj(a)*b.^5*conj(b).^2, sqrt(3628800)/144*a*conj(a).^3*b.^4 - sqrt(3628800)/240*conj(a).^2*b.^5*conj(b), sqrt(29030400)/720*conj(a).^3*b.^5;
                sqrt(23224320)/576*a.^4*conj(b).^4, -sqrt(2903040)/144*a.^4*conj(a)*conj(b).^3 + sqrt(2903040)/144*a.^3*b*conj(b).^4, sqrt(829440)/96*a.^4*conj(a).^2*conj(b).^2 - sqrt(829440)/36*a.^3*conj(a)*b*conj(b).^3 + sqrt(829440)/96*a.^2*b.^2*conj(b).^4, -sqrt(414720)/144*a.^4*conj(a).^3*conj(b) + sqrt(414720)/24*a.^3*conj(a).^2*b*conj(b).^2 - sqrt(414720)/24*a.^2*conj(a)*b.^2*conj(b).^3 + sqrt(414720)/144*a*b.^3*conj(b).^4, sqrt(331776)/576*a.^4*conj(a).^4 - sqrt(331776)/36*a.^3*conj(a).^3*b*conj(b) + sqrt(331776)/16*a.^2*conj(a).^2*b.^2*conj(b).^2 - sqrt(331776)/36*a*conj(a)*b.^3*conj(b).^3 + sqrt(331776)/576*b.^4*conj(b).^4, sqrt(414720)/144*a.^3*conj(a).^4*b - sqrt(414720)/24*a.^2*conj(a).^3*b.^2*conj(b) + sqrt(414720)/24*a*conj(a).^2*b.^3*conj(b).^2 - sqrt(414720)/144*conj(a)*b.^4*conj(b).^3, sqrt(829440)/96*a.^2*conj(a).^4*b.^2 - sqrt(829440)/36*a*conj(a).^3*b.^3*conj(b) + sqrt(829440)/96*conj(a).^2*b.^4*conj(b).^2, sqrt(2903040)/144*a*conj(a).^4*b.^3 - sqrt(2903040)/144*conj(a).^3*b.^4*conj(b), sqrt(23224320)/576*conj(a).^4*b.^4;
                -sqrt(29030400)/720*a.^3*conj(b).^5, sqrt(3628800)/144*a.^3*conj(a)*conj(b).^4 - sqrt(3628800)/240*a.^2*b*conj(b).^5, -sqrt(1036800)/72*a.^3*conj(a).^2*conj(b).^3 + sqrt(1036800)/48*a.^2*conj(a)*b*conj(b).^4 - sqrt(1036800)/240*a*b.^2*conj(b).^5, sqrt(518400)/72*a.^3*conj(a).^3*conj(b).^2 - sqrt(518400)/24*a.^2*conj(a).^2*b*conj(b).^3 + sqrt(518400)/48*a*conj(a)*b.^2*conj(b).^4 - sqrt(518400)/720*b.^3*conj(b).^5, -sqrt(414720)/144*a.^3*conj(a).^4*conj(b) + sqrt(414720)/24*a.^2*conj(a).^3*b*conj(b).^2 - sqrt(414720)/24*a*conj(a).^2*b.^2*conj(b).^3 + sqrt(414720)/144*conj(a)*b.^3*conj(b).^4, sqrt(518400)/720*a.^3*conj(a).^5 - sqrt(518400)/48*a.^2*conj(a).^4*b*conj(b) + sqrt(518400)/24*a*conj(a).^3*b.^2*conj(b).^2 - sqrt(518400)/72*conj(a).^2*b.^3*conj(b).^3, sqrt(1036800)/240*a.^2*conj(a).^5*b - sqrt(1036800)/48*a*conj(a).^4*b.^2*conj(b) + sqrt(1036800)/72*conj(a).^3*b.^3*conj(b).^2, sqrt(3628800)/240*a*conj(a).^5*b.^2 - sqrt(3628800)/144*conj(a).^4*b.^3*conj(b), sqrt(29030400)/720*conj(a).^5*b.^3;
                sqrt(58060800)/1440*a.^2*conj(b).^6, -sqrt(7257600)/240*a.^2*conj(a)*conj(b).^5 + sqrt(7257600)/720*a*b*conj(b).^6, sqrt(2073600)/96*a.^2*conj(a).^2*conj(b).^4 - sqrt(2073600)/120*a*conj(a)*b*conj(b).^5 + sqrt(2073600)/1440*b.^2*conj(b).^6, -sqrt(1036800)/72*a.^2*conj(a).^3*conj(b).^3 + sqrt(1036800)/48*a*conj(a).^2*b*conj(b).^4 - sqrt(1036800)/240*conj(a)*b.^2*conj(b).^5, sqrt(829440)/96*a.^2*conj(a).^4*conj(b).^2 - sqrt(829440)/36*a*conj(a).^3*b*conj(b).^3 + sqrt(829440)/96*conj(a).^2*b.^2*conj(b).^4, -sqrt(1036800)/240*a.^2*conj(a).^5*conj(b) + sqrt(1036800)/48*a*conj(a).^4*b*conj(b).^2 - sqrt(1036800)/72*conj(a).^3*b.^2*conj(b).^3, sqrt(2073600)/1440*a.^2*conj(a).^6 - sqrt(2073600)/120*a*conj(a).^5*b*conj(b) + sqrt(2073600)/96*conj(a).^4*b.^2*conj(b).^2, sqrt(7257600)/720*a*conj(a).^6*b - sqrt(7257600)/240*conj(a).^5*b.^2*conj(b), sqrt(58060800)/1440*conj(a).^6*b.^2;
                -sqrt(203212800)/5040*a*conj(b).^7, sqrt(25401600)/720*a*conj(a)*conj(b).^6 - sqrt(25401600)/5040*b*conj(b).^7, -sqrt(7257600)/240*a*conj(a).^2*conj(b).^5 + sqrt(7257600)/720*conj(a)*b*conj(b).^6, sqrt(3628800)/144*a*conj(a).^3*conj(b).^4 - sqrt(3628800)/240*conj(a).^2*b*conj(b).^5, -sqrt(2903040)/144*a*conj(a).^4*conj(b).^3 + sqrt(2903040)/144*conj(a).^3*b*conj(b).^4, sqrt(3628800)/240*a*conj(a).^5*conj(b).^2 - sqrt(3628800)/144*conj(a).^4*b*conj(b).^3, -sqrt(7257600)/720*a*conj(a).^6*conj(b) + sqrt(7257600)/240*conj(a).^5*b*conj(b).^2, sqrt(25401600)/5040*a*conj(a).^7 - sqrt(25401600)/720*conj(a).^6*b*conj(b), sqrt(203212800)/5040*conj(a).^7*b;
                sqrt(1625702400)/40320*conj(b).^8, -sqrt(203212800)/5040*conj(a)*conj(b).^7, sqrt(58060800)/1440*conj(a).^2*conj(b).^6, -sqrt(29030400)/720*conj(a).^3*conj(b).^5, sqrt(23224320)/576*conj(a).^4*conj(b).^4, -sqrt(29030400)/720*conj(a).^5*conj(b).^3, sqrt(58060800)/1440*conj(a).^6*conj(b).^2, -sqrt(203212800)/5040*conj(a).^7*conj(b), sqrt(1625702400)/40320*conj(a).^8];
        end

        function c = propagation(t,A)
            c = c_0;
            phase_coeff = [0, 14, 12, 24, 20, 30, 24, 32, 24, 30, 20, 24, 12, 14, 0];
            N = length(phase_coeff);
            for i = linspace(1,floor(t), floor(t))
                c = propagator(1, A, phase_coeff(i)*pi/N)*c;
            end
            if t ~= floor(t)
                c = propagator(t-floor(t), A, phase_coeff((floor(t)+1))*pi/N)*c;
            end
        end

        c = [];
        t_interval = linspace(0,15,400);
        for t = t_interval
            c = [c propagation(t,pi+delta)];
        end
        plot(t_interval, (abs(c(1,:))).^2, '-', t_interval, (abs(c(2,:))).^2, '-', t_interval, (abs(c(3,:))).^2, '--', t_interval, (abs(c(4,:))).^2, '-', t_interval, (abs(c(5,:))).^2, '-', t_interval, (abs(c(6,:))).^2, '-', t_interval, (abs(c(7,:))).^2, '--', t_interval, (abs(c(8,:))).^2, '-', t_interval, (abs(c(9,:))).^2, '-');
        ylim([0, 1]);
        legend('', '', '$\left|3\right>$', '', '', '', '$\left|7\right>$', '', '', 'Interpreter','latex','Location','north')
    end

subplot(2,2,1);
Fig_6a
title("(a)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Time in Units of $\displaystyle {\pi}$/$\displaystyle{\Omega_R}$','Interpreter','latex')
ylabel('Part of Population','Interpreter','latex')
subplot(2,2,2)
Fig_6b
title("(b)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Time in Units of $\displaystyle {\pi}$/$\displaystyle{\Omega_R}$','Interpreter','latex')
ylabel('Part of Population','Interpreter','latex')
subplot(2,2,3);
Fig_6c
title("(c)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Time in Units of $\displaystyle {\pi}$/$\displaystyle{\Omega_R}$','Interpreter','latex')
ylabel('Part of Population','Interpreter','latex')
subplot(2,2,4)
Fig_6d
title("(d)")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Time in Units of $\displaystyle {\pi}$/$\displaystyle{\Omega_R}$','Interpreter','latex')
ylabel('Part of Population','Interpreter','latex')
end
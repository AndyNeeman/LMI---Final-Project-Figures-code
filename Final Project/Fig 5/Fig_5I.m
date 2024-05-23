function Fig_5I

    function Draw_points(delta,number,leg_title)
        Delta = 0;
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
                c = propagator(t,pi/2+delta,0)*c_0;
            elseif t<=2
                c = propagator(t-1,pi+2*delta,pi/2)*propagator(1,pi/2+delta,0)*c_0;
            else
                c= propagator(t-2,pi/2+delta,0)*propagator(1,pi+2*delta,pi/2)*propagator(1,pi/2+delta,0)*c_0;
            end
        end

        function coeffs = create_majo_poly(C)
            N = length(C);
            coeffs = zeros(1,N);
            for i = linspace(0,N-1,N)
                coeffs(i+1) = C(i+1)*sqrt(nchoosek(N-1,i));
            end
        end

        function angles = state_to_majo_sphere(C)
            N = length(C);
            angles = zeros(2,N-1);
            majo_roots = roots(create_majo_poly(C)).';
            point_counter = 1;
            if isempty(majo_roots)
                angles(1,:) = pi*ones(1,N-1);
            else
                for root = majo_roots
                    angles(1,point_counter) = 2*(atan(abs(root)));
                    angles(2,point_counter) = angle(root);
                    point_counter = point_counter + 1;
                end
            end
        end

        function angles = states_to_majo_sphere(C)
            N = size(C,1);
            M = size(C,2);
            angles = zeros(2,(N-1)*M);
            for j = linspace(1,M,M)
                state = C(:,j);
                state_angles = state_to_majo_sphere(state);
                angles(:,j*(N-1)-N+2:j*(N-1)) = state_angles;
            end
        end

        c = [];
        t_interval = linspace(0,3,200);
        for t = t_interval
            c = [c propagation(t)];
        end
        points = states_to_majo_sphere(c);

        theta = points(1,:);
        phi = points(2,:);
        X_vals = sin(theta).*cos(phi);
        Y_vals = sin(theta).*sin(phi);
        Z_vals = cos(theta);
        if number == 1
            plot3(X_vals,Y_vals,Z_vals,"-","Color","r","LineWidth",2)
        elseif number == 2
            plot3(X_vals,Y_vals,Z_vals,"-","Color","b","LineWidth",2)
        elseif number == 3
            plot3(X_vals,Y_vals,Z_vals,"-","Color",	"#7E2F8E","LineWidth",2)
        else
            plot3(X_vals,Y_vals,Z_vals,"-","Color","k","LineWidth",2)
        end
        
    end


[X,Y,Z] = sphere;
hold on
surf(X,Y,Z,'FaceAlpha', 0.05, 'EdgeAlpha', 0.2, 'FaceColor','cyan');
plot3(linspace(-1,1,50),zeros(1,50),zeros(1,50), 'Color', 'black', 'LineWidth', 1);
plot3(zeros(1,50),linspace(-1,1,50),zeros(1,50), 'Color', 'black', 'LineWidth', 1);
plot3(zeros(1,50),zeros(1,50),linspace(-1,1,50), 'Color', 'black', 'LineWidth', 1);
view(29,26)
Draw_points(0,1,'a')
Draw_points(-pi/20,2,'a')
Draw_points(-pi/12,3,'a')
Draw_points(-pi/8,4,'a')
lgd = legend('','','','','$\displaystyle \delta = 0$', '$\displaystyle \delta = \frac{\pi}{20}$', '$\displaystyle \delta = \frac{\pi}{12}$', '$\displaystyle \delta = \frac{\pi}{8}$', 'Interpreter','latex','Location','east');
pos = lgd.Position;
pos(4) = pos(4)*1.5;
lgd.Position = pos;
axis equal
grid off
axis off
hold off

end

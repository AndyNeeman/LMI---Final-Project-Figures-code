function Fig_5II
    function points = set_up
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
        t_interval = linspace(0,5,401);
        for t = t_interval
            c = [c propagation(t)];
        end
        points = states_to_majo_sphere(c);
    end

curr_points = set_up;

    function Draw_segment(t1,t2,show_start)
        [X,Y,Z] = sphere;
        hold on
        surf(X,Y,Z,'FaceAlpha', 0.05, 'EdgeAlpha', 0.2, 'FaceColor','cyan');
        plot3(linspace(-1,1,50),zeros(1,50),zeros(1,50), 'Color', 'black', 'LineWidth', 1);
        plot3(zeros(1,50),linspace(-1,1,50),zeros(1,50), 'Color', 'black', 'LineWidth', 1);
        plot3(zeros(1,50),zeros(1,50),linspace(-1,1,50), 'Color', 'black', 'LineWidth', 1);
        view(-45,30)
        initial_indicator = 0;
        initial_points = curr_points(:,1+t1*160:2+t1*160);
        if show_start
            initial_indicator = 1;
            for j = [1,2]
                point = initial_points(:,j);
                theta = point(1);
                phi = point(2);
                x = sin(theta)*cos(phi);
                y = sin(theta)*sin(phi);
                z = cos(theta);
                plot3(x,y,z,'.','color', '#808080','MarkerSize',30)
            end
        end
        for i = linspace(1+t1*160+2*initial_indicator,2+160*t2,2*(1-initial_indicator)+160*(t2-t1))
            point = curr_points(:,i);
            theta = point(1);
            phi = point(2);
            x = sin(theta)*cos(phi);
            y = sin(theta)*sin(phi);
            z = cos(theta);
            [theta,phi];
            if mod(i,2) == 0
                plot3(x,y,z,"o","Color","r","MarkerSize",2)
            else
                plot3(x,y,z,".","Color","b")
            end
        end
        xlim([-1;1])
        ylim([-1;1])
        zlim([-1;1])
        axis equal
        grid off
        axis off
        hold off
    end
    t = tiledlayout(2,3,'TileSpacing', 'none', 'Padding', 'none');
    nexttile;
    Draw_segment(0,5,false)
    title("(a)")
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    nexttile;
    Draw_segment(0,1,true)
    title("(b)")
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    nexttile;
    Draw_segment(1,2,true)
    title("(c)")
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    nexttile;
    Draw_segment(2,3,true)
    title("(d)")
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    nexttile;
    Draw_segment(3,4,true)
    title("(e)")
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    nexttile;
    Draw_segment(4,5,true)
    title("(f)")
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
end


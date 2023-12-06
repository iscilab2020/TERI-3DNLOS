classdef HiddenScene
    methods
        
        
        
        function [scene] = make_sparse(obj, scene)
            for t=1:length(scene(1,:,1))
                for p=1:length(scene(1,1,:))
                    [ma, ind] = max(scene(:, t, p));
                    for r=1:length(scene(:,1,1))
                        if ind ~= r
                            scene(r, t, p) = 0;
                        end
                    end
                end
            end
        end
        
        function [matrix] = add_point_source(obj, S, original_matrix, x, y, z)
            if x == 0
                x = 1/S;
            end
            if y == 0
                y = 1/S;
            end
            if z == 0
                z = 1/S;
            end
            original_matrix((ceil(x*S)), (ceil(y*S)), ceil(z*S)) = 1;
            matrix = original_matrix;
        end
        
        function [ans] = display_image(obj, N, scene, f)
            figure(f);
            xx = linspace(0, 1, N);
            yy = linspace(0, 1, N);
            imagesc(xx,yy,scene)
            set(gca, 'xdir', 'reverse')
            xlabel('x');
            ylabel('y');
            
            top_score = max(scene(:));
            rnd = 10;
            colormap(gray(256));
            %colorbar('Ticks',[0,.2,.4,.6,.8,1],'TickLabels',{num2str(round(0*top_score,rnd)), num2str(round(0.2*top_score,rnd)),num2str(round(0.4*top_score,rnd)),num2str(round(0.6*top_score,rnd)),num2str(round(0.8*top_score,rnd)),num2str(round(1.0*top_score,rnd))})
            colorbar();
            
            ans = 0;
        end
        
        function [ans] = display_color_image(obj, N, scene, f)
            figure(f);
            xx = linspace(0, 1, N);
            yy = linspace(0, 1, N);
            imagesc(xx,yy,scene)
            set(gca, 'xdir', 'reverse')
            xlabel('x');
            ylabel('y');
            
            top_score = max(scene(:));
            rnd = 2;
            colormap(gray(256));
            colorbar('Ticks',[0,.2,.4,.6,.8,1],'TickLabels',{num2str(round(0*top_score,rnd)), num2str(round(0.2*top_score,rnd)),num2str(round(0.4*top_score,rnd)),num2str(round(0.6*top_score,rnd)),num2str(round(0.8*top_score,rnd)),num2str(round(1.0*top_score,rnd))})
                   
            ans = 0;
        end
        
        function [ans] = plot_solution(obj, hidden_scene, scene, f)
            hold off
            figure(f) 
            stem(hidden_scene(:))
            hold on
            stem(scene(:),'p') 
            hold off
        end
        
        function [smooth_scene] = linear_smoothing(obj, size, scene, N)
            radius = (size - 1)/2;
            smooth_scene = zeros(N, N);
            for i = 1:N
                for j = 1:N
                    num = 0;
                    for a = max(1, i - radius):min(N, i + radius)
                        for b = max(1, j - radius):min(N, j + radius)
                            smooth_scene(i, j) = smooth_scene(i, j) + scene(a, b);
                            num = num + 1;
                        end
                    end
                    smooth_scene(i, j) = smooth_scene(i, j)/(num^2);
                end
            end
        end
        
        function [result] = reflectX(obj, scene)
            result = zeros(length(scene(:,1)), length(scene(1,:)));
            for i=1:length(scene(:,1))
                for j=1:length(scene(1,:))
                    result(i, j) = scene(i, length(scene(1, :)) - j + 1);
                end
            end
        end
        
        function [result] = scale_down(obj, scene)
            result = zeros(floor(length(scene(:,1,1))/2), floor(length(scene(1,:,1))/2), 3);
            for i=1:length(result(:,1,1))
                for j=1:length(result(1,:,1))
                    for k=1:3
                        result(i, j, k) = result(i, j, k) + scene(2*i, 2*j, k);
                        result(i, j, k) = result(i, j, k) + scene(2*i-1, 2*j, k);
                        result(i, j, k) = result(i, j, k) + scene(2*i, 2*j-1, k);
                        result(i, j, k) = result(i, j, k) + scene(2*i-1, 2*j-1, k);
                        result(i, j, k) = double(result(i, j, k))/4;
                    end
                end
            end
        end
        
        function [error] = get_error(obj, correct, guess)
            error = 0;
            for i = 1:length(correct(:))
                error = error + (correct(i) - guess(i))^2;
            end
        end
        
        function [scene] = round_scene(obj, scene, S)
            for object=1:length(scene(:))/6
                scene(object, 1, 1) = (ceil(scene(object, 1, 1)*S)-.5)/S;
                scene(object, 2, 1) = (ceil(scene(object, 2, 1)*S)-.5)/S;
                scene(object, 2, 2) = (ceil(scene(object, 2, 2)*S)-.5)/S;
                scene(object, 3, 1) = (ceil(scene(object, 3, 1)*S)-.5)/S;
                scene(object, 3, 2) = (ceil(scene(object, 3, 2)*S)-.5)/S;
            end
        end
        
        function [scene] = round_scene_sph(obj, scene, aS, rS, scene_dim)
            for object=1:length(scene(:))/6
                scene(object, 1, 1) = (ceil(scene(object, 1, 1)*rS/scene_dim)-.5)/rS*scene_dim;
                scene(object, 2, 1) = (ceil(scene(object, 2, 1)*aS/(pi/2))-.5)/aS*(pi/2);
                scene(object, 2, 2) = (ceil(scene(object, 2, 2)*aS/(pi/2))-.5)/aS*(pi/2);
                scene(object, 3, 1) = (ceil(scene(object, 3, 1)*aS/(pi/2))-.5)/aS*(pi/2);
                scene(object, 3, 2) = (ceil(scene(object, 3, 2)*aS/(pi/2))-.5)/aS*(pi/2);
            end
        end
        
        function [hidden_scene] = construct_hidden_scene(obj, hidden_objects, S)
            hidden_scene = zeros(S, S, S);
            for i=1:S
                for j=1:S
                    for k=1:S
                        for obj=1:length(hidden_objects(:))/6
                            if hidden_objects(obj, 1, 1) == i/S && ...
                                    hidden_objects(obj, 2, 1) <= j/S && j/S < hidden_objects(obj, 2, 2) && ...
                                    hidden_objects(obj, 3, 1) <= k/S && k/S < hidden_objects(obj, 3, 2)
                                hidden_scene(i, j, k) = hidden_objects(obj, 1, 2);
                                continue
                            end
                        end
                    end
                end
            end
        end
        
        function [hidden_scene] = construct_hidden_scene_sph(obj, hidden_objects, aS, psiS, rS, scene_dim)
            hidden_scene = zeros(rS, aS, psiS);
            
            RHO = linspace(0, scene_dim, rS+1)+scene_dim/rS/2;
            THETA = linspace(0, pi/2, aS+1)+1/aS/2;
            PSI = linspace(0, pi/2, aS+1)+1/aS/2;
            
            RHO = RHO(1:rS);
            THETA = THETA(1:aS);
            PSI = PSI(1:aS);
            PSI = PSI(aS-psiS+1:aS);
            
            for rho=1:rS
                for theta=1:aS
                    for psi=1:psiS
                        for obj=1:length(hidden_objects(:))/6
                            
                            if round(hidden_objects(obj, 1, 1), 2) == round(RHO(rho), 2) || rS == 1
                                if hidden_objects(obj, 2, 1) <= THETA(theta) && THETA(theta) <= hidden_objects(obj, 2, 2)
                                    if hidden_objects(obj, 3, 1) <= PSI(psi) && PSI(psi) <= hidden_objects(obj, 3, 2)
                                        hidden_scene(rho, theta, psi) = hidden_objects(obj, 1, 2);
                                        continue
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function [nothing] = plot_obj(obj, hidden_objects)
            for obj=1:length(hidden_objects(:))/6
                figure(80);
                fill3([hidden_objects(obj, 1, 1) hidden_objects(obj, 1, 1) hidden_objects(obj, 1, 1) hidden_objects(obj, 1, 1)], ...
                        [hidden_objects(obj, 2, 1) hidden_objects(obj, 2, 2) hidden_objects(obj, 2, 2) hidden_objects(obj, 2, 1)], ...
                        [hidden_objects(obj, 3, 1) hidden_objects(obj, 3, 1) hidden_objects(obj, 3, 2) hidden_objects(obj, 3, 2)], "black");
                axis([0 1 0 1 0 1]);
                grid on;
                hold on;
            end
            set(gca, 'Ydir', 'reverse')
            hold off; 
        end
        
        function [nothing] = plot_obj_sph(obj, hidden_objects, f)
            for object=1:length(hidden_objects(:))/6
                figure(f);
                
                p1 = [hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 1))*cos(hidden_objects(object, 3, 1)), ...
                      hidden_objects(object, 1, 1)*sin(hidden_objects(object, 2, 1))*cos(hidden_objects(object, 3, 1)), ...
                      hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 1))*sin(hidden_objects(object, 3, 1))];
                
                p2 = [hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 2))*cos(hidden_objects(object, 3, 1)), ...
                      hidden_objects(object, 1, 1)*sin(hidden_objects(object, 2, 2))*cos(hidden_objects(object, 3, 1)), ...
                      hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 2))*sin(hidden_objects(object, 3, 1))];
                
                p3 = [hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 2))*cos(hidden_objects(object, 3, 2)), ...
                      hidden_objects(object, 1, 1)*sin(hidden_objects(object, 2, 2))*cos(hidden_objects(object, 3, 2)), ...
                      hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 2))*sin(hidden_objects(object, 3, 2))];
                
                p4 = [hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 1))*cos(hidden_objects(object, 3, 2)), ...
                      hidden_objects(object, 1, 1)*sin(hidden_objects(object, 2, 1))*cos(hidden_objects(object, 3, 2)), ...
                      hidden_objects(object, 1, 1)*cos(hidden_objects(object, 2, 1))*sin(hidden_objects(object, 3, 2))];
                
                fill3(  [p1(1) p2(1) p3(1) p4(1)], ...
                        [p1(2) p2(2) p3(2) p4(2)], ...
                        [p1(3) p2(3) p3(3) p4(3)], ...
                        "black");
                axis([0 1 0 1 0 1]);
                grid on;
                hold on;
            end
            set(gca, 'Ydir', 'reverse')
            hold off; 
        end
        
        function [nothing] = plot_scene(obj, scene, S, threshold, f)
            for i=1:S
                for j=1:S
                    
                    for k=1:S
                        figure(f)
                        if scene(i, j, k) > threshold
                            fill3([i/S i/S i/S i/S], ...
                                [j/S (j+1)/S (j+1)/S j/S], ...
                                [k/S k/S (k+1)/S (k+1)/S],"black");
                            axis([0 1 0 1 0 1]);
                            grid on;
                            hold on;
                            
                            set(gca, 'Ydir', 'reverse')
                        end
                    end
                end
            end
            
            set(gca, 'Ydir', 'reverse')
            hold off;
        end
        
        function [nothing] = plot_scene_sph(obj, scene, aS, psiS, rS, threshold, f, scene_dim)
            unit = scene_dim;
            top_score = max(scene(:));
            
            for rho=1:rS
                for theta=1:aS
                    for psi=1:psiS
                        figure(f)
                        if scene(rho, theta, psi) > threshold
                            
                            x = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-.5)/aS)).^2+(tan(pi/2*(psi-.5+aS-psiS)/aS)).^2)).^0.5;
                            y = tan(pi/2*(theta-.5)/aS).*x;
                            %y = unit-tan(pi/2*(theta-.5)/aS).*x;
                            %z = unit-tan(pi/2*(psi-.5)/aS).*x;
                            z = tan(pi/2*(psi-.5+aS-psiS)/aS).*x;
                            
                            f1=scatter3(x, y, z, 'Marker', 'o', "MarkerEdgeAlpha", (scene(rho, theta, psi))/top_score, "MarkerEdgeColor", 'k');
                            f1.SizeData = f1.SizeData/5;
                            %scatter3(x, y, z);
                            axis([0 unit 0 unit 0 unit]);
                            view(-130, 20);
                            %set(gca,'Color',[.3 .3 .3]);
                            set(gca,'FontSize',14)
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                            %set(gcf, 'InvertHardcopy', 'off');
                            grid on;
                            hold on;
                        end
                    end
                end
            end
            
            set(gca, 'Ydir', 'reverse')
            set(gca, 'Zdir', 'reverse')
            hold off;
        end
        
        function [nothing] = fill_scene_sph(obj, scene, aS, psiS, rS, threshold, f, scene_dim, h)
            unit = scene_dim/1.3;
            top_score = max(scene(:));
            
            rho = 1.5;
            for theta=1:aS
                for psi=1:psiS
                    figure(f)
                    if scene(1, theta, psi) > threshold

                        x1 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        y1 = tan(pi/2*(theta-1)/aS).*x1;
                        z1 = tan(pi/2*(psi-1+aS-psiS)/aS).*x1;

                        x2 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        y2 = tan(pi/2*(theta)/aS).*x2;
                        z2 = tan(pi/2*(psi-1+aS-psiS)/aS).*x2;

                        x3 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        y3 = tan(pi/2*(theta-1)/aS).*x3;
                        z3 = tan(pi/2*(psi+aS-psiS)/aS).*x3;

                        x4 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        y4 = tan(pi/2*(theta)/aS).*x4;
                        z4 = tan(pi/2*(psi+aS-psiS)/aS).*x4;

                        color = 'r';
                        if(pi/2*(theta-1)/aS < atan(tan(pi/2*(psi+aS-psiS)/aS) * 0.5 / (h)))
                            color = 'g';
                        end
                        
                        fill3([y1, y2, y3], [x1 x2 x3], [z1, z2, z3], color, 'FaceAlpha', (scene(1, theta, psi))/top_score/2, 'EdgeAlpha', 0);
                        hold on;
                        fill3([y4, y2, y3], [x4 x2 x3], [z4, z2, z3], color, 'FaceAlpha', (scene(1, theta, psi))/top_score/2, 'EdgeAlpha', 0);



                        %rnd = 2;
                        %colormap(gray(256));
                        %colorbar('Ticks',[0,.2,.4,.6,.8,1],'TickLabels',{num2str(round(0*top_score,rnd)), num2str(round(0.2*top_score,rnd)),num2str(round(0.4*top_score,rnd)),num2str(round(0.6*top_score,rnd)),num2str(round(0.8*top_score,rnd)),num2str(round(1.0*top_score,rnd))})

                    end
                end
            end
            
           
            
            %Robert Bemis (2022). Animated GIF (https://www.mathworks.com/matlabcentral/fileexchange/21944-animated-gif), MATLAB Central File Exchange. Retrieved September 7, 2022.
            
            %daspect([1,1,1]);
            %OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',s(end)),OptionZ);
        end
        
        function [nothing] = fill_scene_sph_rho(obj, scene, aS, psiS, rho_0, threshold, f, scene_dim)
            unit = scene_dim/1.3;
            top_score = max(scene(:));
            
            for theta=1:aS
                for psi=1:psiS
                    figure(f)
                    if scene(theta, psi) > threshold

                        x1 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        y1 = tan(pi/2*(theta-1)/aS).*x1;
                        z1 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x1;

                        x2 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        y2 = tan(pi/2*(theta)/aS).*x2;
                        z2 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x2;

                        x3 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        y3 = tan(pi/2*(theta-1)/aS).*x3;
                        z3 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x3;

                        x4 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        y4 = tan(pi/2*(theta)/aS).*x4;
                        z4 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x4;

                        fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], 'w', 'FaceAlpha', (scene(theta, psi))/top_score, 'EdgeAlpha', 0);
                        hold on;
                        fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], 'w', 'FaceAlpha', (scene(theta, psi))/top_score, 'EdgeAlpha', 0);

                        axis([0 unit 0 unit 0 unit]);
                        view(-130, 20);
                        xlabel('X[m]');
                        ylabel('Y[m]');
                        zlabel('Z[m]');
                        grid on;
                        hold on;


                        set(gca,'Color',[0.1 0.1 0.1]);
                        set(gcf, 'InvertHardcopy', 'off');
                        set(gca, 'XColor', [0.6 0.6 0.6])
                        set(gca, 'YColor', [0.6 0.6 0.6])

                        %rnd = 2;
                        %colormap(gray(256));
                        %colorbar('Ticks',[0,.2,.4,.6,.8,1],'TickLabels',{num2str(round(0*top_score,rnd)), num2str(round(0.2*top_score,rnd)),num2str(round(0.4*top_score,rnd)),num2str(round(0.6*top_score,rnd)),num2str(round(0.8*top_score,rnd)),num2str(round(1.0*top_score,rnd))})

                    end
                end
            end
            
            
            set(gca, 'Ydir', 'reverse')
            hold off;
            
            
            
            %Robert Bemis (2022). Animated GIF (https://www.mathworks.com/matlabcentral/fileexchange/21944-animated-gif), MATLAB Central File Exchange. Retrieved September 7, 2022.
            
            daspect([1,1,1]);
            OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',s(end)),OptionZ);
        end
        
        function [nothing] = plot_scene_og(obj, scene, aS, psiS, rS, threshold, f, scene_dim)
            unit = scene_dim;
            top_score = max(scene(:));
            
            for rho=1:rS
                for theta=1:aS
                    for psi=1:psiS
                        figure(f)
                        if scene(rho, theta, psi) > threshold
                            
                            x1 = (((scene_dim*(rho-.5)*3^.5/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            y1 = tan(pi/2*(theta-1)/aS).*x1;
                            z1 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x1;
                            
                            x2 = (((scene_dim*(rho-.5)*3^.5/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            y2 = tan(pi/2*(theta)/aS).*x2;
                            z2 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x2;
                            
                            x3 = (((scene_dim*(rho-.5)*3^.5/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            y3 = tan(pi/2*(theta-1)/aS).*x3;
                            z3 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x3;
                            
                            x4 = (((scene_dim*(rho-.5)*3^.5/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            y4 = tan(pi/2*(theta)/aS).*x4;
                            z4 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x4;
                            
                            f1=fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], 'k');
                            f1.FaceAlpha = (scene(rho, theta, psi))/top_score;
                            f1.EdgeAlpha = 0;
                            f1=fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], 'k');
                            f1.FaceAlpha = (scene(rho, theta, psi))/top_score;
                            f1.EdgeAlpha = 0;
                            
                            axis([0 unit 0 unit 0 unit]);
                            view(-130, 20);
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                            grid on;
                            hold on;
                        end
                    end
                end
            end
            
            set(gca, 'Ydir', 'reverse')
            hold off;
        end
        
        
        function [nothing] = plot_color_sph(obj, scene, aS, psiS, rS, threshold, f, scene_dim, filename)
            unit = scene_dim;
            unit = 1;
            maximum = max(scene(:));
            
            for rho=1:rS
                for theta=1:aS
                    for psi=1:psiS
                        figure(f)

                        red = max(scene(rho, theta, psi, 1), 0);
                        green = max(scene(rho, theta, psi, 2), 0);
                        blue = max(scene(rho, theta, psi, 3), 0);
                        mag = max(max(red,green),blue);
                        
                        if red > threshold || green > threshold || blue > threshold
                            
                            x1 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            y1 = tan(pi/2*(theta-1)/aS).*x1;
                            z1 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x1;
                            
                            x2 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            y2 = tan(pi/2*(theta)/aS).*x2;
                            z2 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x2;
                            
                            x3 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            y3 = tan(pi/2*(theta-1)/aS).*x3;
                            z3 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x3;
                            
                            x4 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            y4 = tan(pi/2*(theta)/aS).*x4;
                            z4 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x4;
                            
                            color = [red/mag green/mag blue/mag];
                            
                            
                            %fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                            %hold on;
                            %fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                            fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                            hold on;
                            fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                            
                            axis([0 unit 0 unit 0 unit]);
                            view(-130, 20);
                            xlabel('X');
                            ylabel('Y');
                            zlabel('Z');
                            grid on;
                            hold on;
                        end
                    end
                    
                end
            end
            
            set(gca,'Color',[0.1 0.1 0.1]);
            set(gcf, 'InvertHardcopy', 'off');
            set(gca, 'XColor', [0.6 0.6 0.6])
            set(gca, 'YColor', [0.6 0.6 0.6])

            set(gca, 'Ydir', 'reverse')
            hold off;
            
            daspect([1,1,1]);
            OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',filename),OptionZ);
        end
        
        
        function [nothing] = plot_color_sph_rho(obj, scene, aS, psiS, rho, threshold, f, scene_dim, filename)
            unit = scene_dim;
            maximum = max(scene(:));
            
            for theta=1:aS
                for psi=1:psiS
                    figure(f)

                    blue = max(scene(theta, psi, 1), 0);
                    green = max(scene(theta, psi, 2), 0);
                    red = max(scene(theta, psi, 3), 0);
                    rho = max(scene(theta, psi, 4), 0);
                    mag = max(max(red,green),blue);

                    if red > threshold || green > threshold || blue > threshold

                        x1 = ((rho.^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        y1 = tan(pi/2*(theta-1)/aS).*x1;
                        z1 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x1;

                        x2 = ((rho.^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        y2 = tan(pi/2*(theta)/aS).*x2;
                        z2 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*x2;

                        x3 = ((rho.^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        y3 = tan(pi/2*(theta-1)/aS).*x3;
                        z3 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x3;

                        x4 = ((rho.^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        y4 = tan(pi/2*(theta)/aS).*x4;
                        z4 = unit - tan(pi/2*(psi+aS-psiS)/aS).*x4;

                        color = [red/mag green/mag blue/mag];


                        %fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                        %hold on;
                        %fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                        fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                        hold on;
                        fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);

                        axis([0 unit 0 unit 0 unit]);
                        view(-130, 20);
                        xlabel('X [m]');
                        ylabel('Y [m]');
                        zlabel('Z [m]');
                        grid on;
                        hold on;
                    end
                end

            end
            
            set(gca,'Color',[0.1 0.1 0.1]);
            set(gcf, 'InvertHardcopy', 'off');
            set(gca, 'XColor', [0.6 0.6 0.6])
            set(gca, 'YColor', [0.6 0.6 0.6])

            set(gca, 'Ydir', 'reverse')
            hold off;
            
            daspect([1,1,1]);
            %OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',filename),OptionZ);
        end
        
        function [nothing] = plot_color_prof_code(obj, THETA, THETA_end, PSI, PSI_end, rhos, colors, Nt, Np, threshold, f, scene_dim)
            unit = scene_dim;
            
            for i=1:length(THETA(:))
                figure(f)
                
                rho = rhos(:);

                theta1 = THETA(:);
                theta2 = THETA_end(:);
                psi1 = PSI(:);
                psi2 = PSI_end(:);
                
                

                red = colors(i,1);
                green = colors(i,2);
                blue = colors(i,3);
                mag = max(max(red,green),blue);

                %if red > threshold || green > threshold || blue > threshold
                if mag > 0 && theta2 < pi/2 - 0.0001 && psi2 < pi/2 - 0.0001

                    x1 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi1).^2))).^0.5;
                    y1 = tan(theta1).*x1;
                    z1 = tan(psi1).*x1;

                    x2 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi1)).^2)).^0.5;
                    y2 = tan(theta2).*x2;
                    z2 = tan(psi1).*x2;

                    x3 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi2)).^2)).^0.5;
                    y3 = tan(theta1).*x3;
                    z3 = tan(psi2).*x3;

                    x4 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi2)).^2)).^0.5;
                    y4 = tan(theta2).*x4;
                    z4 = tan(psi2).*x4;

                    %color = [red/mag green/mag blue/mag];
                    color = [red green blue];
                    if max(color) > 1
                        color = color/max(color);
                    end


                    %fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                    %hold on;
                    %fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                    fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                    hold on;
                    fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);

                    %unit = max([unit, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4]);
                    axis([0 unit 0 unit 0 unit]);
                    view(-130, 20);
                    xlabel('X [m]');
                    ylabel('Y [m]');
                    zlabel('Z [m]');
                    grid on;
                    hold on;
                end

            end
            
            set(gca,'Color',[0.1 0.1 0.1]);
            set(gcf, 'InvertHardcopy', 'off');
            set(gca, 'XColor', [0.6 0.6 0.6])
            set(gca, 'YColor', [0.6 0.6 0.6])

            set(gca, 'Ydir', 'reverse')
            set(gca, 'Zdir', 'reverse')
            hold off;
            
            daspect([1,1,1]);
            %OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',filename),OptionZ);
        end
        
        function [nothing] = plot_color_prof_code2(obj, THETA, THETA_end, PSI, PSI_end, rhos, colors, Nt, Np, threshold, f, scene_dim)
            unit = scene_dim;
            
            rho = rhos(:);

            theta1 = THETA(:);
            theta2 = THETA_end(:);
            psi1 = PSI(:);
            psi2 = PSI_end(:);

            red = colors(:,1);
            green = colors(:,2);
            blue = colors(:,3);
            mag = max(max(red,green),blue);

            x1 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi1).^2))).^0.5;
            y1 = tan(theta1).*x1;
            z1 = tan(psi1).*x1;

            x2 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi1)).^2)).^0.5;
            y2 = tan(theta2).*x2;
            z2 = tan(psi1).*x2;

            x3 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi2)).^2)).^0.5;
            y3 = tan(theta1).*x3;
            z3 = tan(psi2).*x3;

            x4 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi2)).^2)).^0.5;
            y4 = tan(theta2).*x4;
            z4 = tan(psi2).*x4;
            
            for i=1:length(THETA(:))
                
                figure(f)
                
                %color = [red/mag green/mag blue/mag];
                color = [red(i) green(i) blue(i)];
                if max(color) > 1
                    color = color/max(color);
                end
                color(color<0) = 0;


                %fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                %hold on;
                %fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                %fill3([x1(i) x2(i) x3(i)], [y1(i), y2(i), y3(i)], [z1(i), z2(i), z3(i)], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                %hold on;
                %fill3([x4(i) x2(i) x3(i)], [y4(i), y2(i), y3(i)], [z4(i), z2(i), z3(i)], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                fill3([y1(i) y2(i) y3(i)], [x1(i), x2(i), x3(i)], [z1(i), z2(i), z3(i)], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                hold on;
                fill3([y4(i) y2(i) y3(i)], [x4(i), x2(i), x3(i)], [z4(i), z2(i), z3(i)], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);

                %unit = max([unit, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4]);
                axis([0 unit 0 unit 0 unit]);
                view(-130, 20);
                xlabel('Y [m]');
                ylabel('X [m]');
                zlabel('Z [m]');
                grid on;
                hold on;

            end
            
            set(gca,'Color',[0.1 0.1 0.1]);
            set(gcf, 'InvertHardcopy', 'off');
            set(gca, 'XColor', [0.6 0.6 0.6])
            set(gca, 'YColor', [0.6 0.6 0.6])

            set(gca, 'Ydir', 'reverse')
            set(gca, 'Zdir', 'reverse')
            hold off;
            
            daspect([1,1,1]);
            %OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',filename),OptionZ);
        end
        
        
        
        function [nothing] = test_plot(obj, coords, rhos, colors, Nt, Np, threshold, f, scene_dim)
            unit = scene_dim;
            
            for i=1:length(coords)
                figure(f)
                
                theta = coords(i,1);
                psi = coords(i,2);
                rho = rhos(i);
                
                theta1 = theta - 1/Nt/2*pi/2;
                theta2 = theta + 1/Nt/2*pi/2;
                psi1 = psi - 1/Np/2*pi/2;
                psi2 = psi + 1/Np/2*pi/2;

                red = colors(i,1);
                green = colors(i,2);
                blue = colors(i,3);
                mag = max(max(red,green),blue);

                %if red > threshold || green > threshold || blue > threshold
                if mag > 0 && theta2 < pi/2 - 0.0001 && psi2 < pi/2 - 0.0001

                    x1 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi1).^2))).^0.5;
                    y1 = tan(theta1).*x1;
                    z1 = unit - tan(psi1).*x1;

                    x2 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi1)).^2)).^0.5;
                    y2 = tan(theta2).*x2;
                    z2 = unit - tan(psi1).*x2;

                    x3 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi2)).^2)).^0.5;
                    y3 = tan(theta1).*x3;
                    z3 = unit - tan(psi2).*x3;

                    x4 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi2)).^2)).^0.5;
                    y4 = tan(theta2).*x4;
                    z4 = unit - tan(psi2).*x4;

                    %color = [red/mag green/mag blue/mag];
                    color = [red green blue];
                    if max(color) > 1
                        color = color/max(color);
                    end


                    %fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                    %hold on;
                    %fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                    fill3([x1 x2 x3], [y1, y2, y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                    hold on;
                    fill3([x4 x2 x3], [y4, y2, y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                    
                    axis([0 unit 0 unit 0 unit]);
                    set(gca, 'Ydir', 'reverse')
                    view(-130, 20);
                    xlabel('X [m]');
                    ylabel('Y [m]');
                    zlabel('Z [m]');
                    grid on;
                    
                    set(gca,'Color',[0.1 0.1 0.1]);
                    set(gcf, 'InvertHardcopy', 'off');
                    set(gca, 'XColor', [0.6 0.6 0.6])
                    set(gca, 'YColor', [0.6 0.6 0.6])

                    set(gca, 'Ydir', 'reverse')
                    hold on;

                    daspect([1,1,1]);
                    
                    saveas(figure(f),strcat("Testing_Plots/",strrep(num2str(theta),'.','d'),'_',strrep(num2str(psi),'.','d'),'_',strrep(num2str(rho),'.','d'),'.png'));
                    %f = f + 1;
                    
                    
                end

            end
            
            %OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',filename),OptionZ);
        end
        
        
        
        function [nothing] = plot_color_prof_code_topview(obj, coords, rhos, colors, Nt, Np, threshold, f, scene_dim)
            unit = scene_dim;
            
            for i=1:length(coords)
                figure(f)
                
                theta = coords(i,1);
                psi = coords(i,2);
                rho = rhos(i);
                
                theta1 = theta - 1/Nt/2*pi/2;
                theta2 = theta + 1/Nt/2*pi/2;
                psi1 = psi - 1/Np/2*pi/2;
                psi2 = psi + 1/Np/2*pi/2;

                red = colors(i,1);
                green = colors(i,2);
                blue = colors(i,3);
                mag = max(max(red,green),blue);

                %if red > threshold || green > threshold || blue > threshold
                if true

                    x1 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi1).^2))).^0.5;
                    y1 = tan(theta1).*x1;

                    x2 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi1)).^2)).^0.5;
                    y2 = tan(theta2).*x2;

                    x3 = ((rho.^2)./(1+(tan(theta1)).^2+(tan(psi2)).^2)).^0.5;
                    y3 = tan(theta1).*x3;

                    x4 = ((rho.^2)./(1+(tan(theta2)).^2+(tan(psi2)).^2)).^0.5;
                    y4 = tan(theta2).*x4;

                    %color = [red/mag green/mag blue/mag];
                    color = [red green blue];
                    if max(color) > 1
                        color = color/max(color);
                    end

                    fill([x1 x2 x3], [y1, y2, y3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                    hold on;
                    fill([x4 x2 x3], [y4, y2, y3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);

                    axis([0 unit 0 unit 0 unit]);
                    view(-130, 20);
                    xlabel('X [m]');
                    ylabel('Y [m]');
                    zlabel('Z [m]');
                    grid on;
                    hold on;
                end

            end
            
            set(gca,'Color',[0.1 0.1 0.1]);
            set(gcf, 'InvertHardcopy', 'off');
            set(gca, 'XColor', [0.6 0.6 0.6])
            set(gca, 'YColor', [0.6 0.6 0.6])

            set(gca, 'Ydir', 'reverse')
            hold off;
            
            daspect([1,1,1]);
            %OptionZ.FrameRate=15;OptionZ.Duration=3;OptionZ.Periodic=true;
            %CaptureFigVid([-250,30;-50,30], strcat('Gifs/',filename),OptionZ);
        end
        
        
        
        
        function [nothing] = plot_color_sph_topview(obj, scene, aS, psiS, rS, threshold, f, scene_dim)
            top_score = max(scene(:));
            maximum = max(scene(:));
            scene_dim = 1;
            for rho=1:rS
                for theta=1:aS
                    for psi=1:psiS
                        figure(f)
                        
                        red = max(scene(rho, theta, psi, 1), 0);
                        green = max(scene(rho, theta, psi, 2), 0);
                        blue = max(scene(rho, theta, psi, 3), 0);
                        mag = max(max(red,green),blue);
                        
                        if red > threshold || green > threshold || blue > threshold
                            
                            y1 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            x1 = tan(pi/2*(theta-1)/aS).*y1;
                            
                            y2 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            x2 = tan(pi/2*(theta)/aS).*y2;
                            
                            y3 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            x3 = tan(pi/2*(theta-1)/aS).*y3;
                            
                            y4 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            x4 = tan(pi/2*(theta)/aS).*y4;
                            
                            color = [red/mag green/mag blue/mag];
                            
                            %fill([x1 x2 x3], [y1, y2, y3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                            fill([x1 x2 x3], [y1, y2, y3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                            hold on;
                            %fill([x4 x2 x3], [y4, y2, y3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                            fill([x4 x2 x3], [y4, y2, y3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', 1);
                            
                            axis([0 scene_dim 0 scene_dim]);
                            grid on;
                            hold on;
                        end
                    end
                end
            end
            
            set(gca,'Color',[0.1 0.1 0.1]);
            set(gcf, 'InvertHardcopy', 'off');
            set(gca, 'XColor', [0.6 0.6 0.6])
            set(gca, 'YColor', [0.6 0.6 0.6])
            
            ax = gca;
            ax.LineWidth = 1;
            ax.GridAlpha = 0.25;
            axis([0 scene_dim 0 scene_dim]);
            set(gca,'FontSize',14)
            axis square
            
            hold off;
        end
        
        function [nothing] = plot_color_sph_sideview(obj, scene, aS, psiS, rS, threshold, f, scene_dim)
            unit = scene_dim;
            top_score = max(scene(:));
            maximum = max(scene(:));
            for rho=1:rS
                for theta=1:aS
                    for psi=1:psiS
                        figure(f)
                        
                        red = max(scene(rho, theta, psi, 1), 0);
                        green = max(scene(rho, theta, psi, 2), 0);
                        blue = max(scene(rho, theta, psi, 3), 0);
                        mag = max(max(red,green),blue);
                        
                        if red > threshold || green > threshold || blue > threshold
                            
                            y1 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            x1 = tan(pi/2*(theta-1)/aS).*y1;
                            
                            y2 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                            x2 = tan(pi/2*(theta)/aS).*y2;
                            
                            y3 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            x3 = tan(pi/2*(theta-1)/aS).*y3;
                            
                            y4 = (((scene_dim*3^.5*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                            x4 = tan(pi/2*(theta)/aS).*y4;
                            z1 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*y1;
                            z2 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*y2;
                            z3 = unit - tan(pi/2*(psi+aS-psiS)/aS).*y3;
                            z4 = unit - tan(pi/2*(psi+aS-psiS)/aS).*y4;
                            
                            color = [red/mag green/mag blue/mag];
                            
                            fill([y1 y2 y3], [z1, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                            hold on;
                            fill([y4 y2 y3], [z4, z2, z3], color, 'FaceColor', color, 'EdgeAlpha', 0, 'FaceAlpha', mag/maximum);
                            
                            axis([0 scene_dim 0 scene_dim]);
                            grid on;
                            hold on;
                        end
                    end
                end
            end
            
            set(gca,'Color',[0.1 0.1 0.1]);
            set(gcf, 'InvertHardcopy', 'off');
            set(gca, 'XColor', [0.6 0.6 0.6])
            set(gca, 'YColor', [0.6 0.6 0.6])
            
            ax = gca;
            ax.LineWidth = 1;
            ax.GridAlpha = 0.25;
            axis([0 scene_dim 0 scene_dim]);
            set(gca,'FontSize',14)
            axis square
            
            hold off;
        end
        
        function [nothing] = plot_scene_sph_topview(obj, scene, aS, psiS, rS, threshold, f, scene_dim, color)
            top_score = max(scene(:));
            for rho=1:rS
                for theta=1:aS
                    for psi=1:psiS
                        figure(f)
                        if scene(rho, theta, psi) > threshold
                            
                            x = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-.5)/aS)).^2+(tan(pi/2*(psi-.5+aS-psiS)/aS)).^2)).^0.5;
                            y = tan(pi/2*(theta-.5)/aS).*x;
                            

                            size = 90;
                            if (color=="grey")
                                size = 110;
                            end
                            f1=scatter(x, y, size, color, "filled", 'MarkerFaceAlpha', .5);
                            axis([0 scene_dim 0 scene_dim]);
                            set(gca,'Color','k');
                            %set(gcf, 'InvertHardcopy', 'off');
                            grid on;
                            hold on;
                        end
                    end
                end
            end
            
            set(gca, 'Ydir', 'reverse')
            hold off;
        end
        
        
        function [nothing] = fill_scene_sph_sideview(obj, scene, aS, psiS, rS, threshold, f, scene_dim)
            unit = scene_dim;
            top_score = max(scene(:));


            for ind=1:length(scene(:))
                figure(f)
                if scene(ind) > threshold
                    
                    theta0 = THETA(ind); theta1 = THETA_end(ind); psi0 = PSI(ind); psi1 = PSI_end(ind);
                
                    
                    y1 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                    x1 = tan(pi/2*(theta-1)/aS).*y1;
                    
                    y2 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                    x2 = tan(pi/2*(theta)/aS).*y2;
                    
                    y3 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                    x3 = tan(pi/2*(theta-1)/aS).*y3;
                    
                    y4 = (((scene_dim*(rho-.5)/rS).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                    x4 = tan(pi/2*(theta)/aS).*y4;
                    
                    z1 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*y1;
                    z2 = unit - tan(pi/2*(psi-1+aS-psiS)/aS).*y2;
                    z3 = unit - tan(pi/2*(psi+aS-psiS)/aS).*y3;
                    z4 = unit - tan(pi/2*(psi+aS-psiS)/aS).*y4;
                    
                    fill([y1 y2 y3], [z1, z2, z3], 'w', 'FaceAlpha', (scene(rho, theta, psi))/top_score, 'EdgeAlpha', 0);
                    hold on;
                    fill([y4 y2 y3], [z4, z2, z3], 'w', 'FaceAlpha', (scene(rho, theta, psi))/top_score, 'EdgeAlpha', 0);
                    
                    axis([0 scene_dim 0 scene_dim]);
                    grid on;
                    hold on;
                    
                    ax = gca;
                    ax.LineWidth = 1;
                    ax.GridAlpha = 0.25;
                    axis([0 scene_dim 0 scene_dim]);
                    set(gca,'FontSize',14)
                    axis square
                    
                end
                
                set(gca,'Color',[0.1 0.1 0.1]);
                set(gcf, 'InvertHardcopy', 'off');
                set(gca, 'XColor', [0.6 0.6 0.6])
                set(gca, 'YColor', [0.6 0.6 0.6])

                %rnd = 2;
                %colormap(gray(256));
                %colorbar('Ticks',[0,.2,.4,.6,.8,1],'TickLabels',{num2str(round(0*top_score,rnd)), num2str(round(0.2*top_score,rnd)),num2str(round(0.4*top_score,rnd)),num2str(round(0.6*top_score,rnd)),num2str(round(0.8*top_score,rnd)),num2str(round(1.0*top_score,rnd))})
                    
            end
            
            hold off;
        end
        
        function [nothing] = fill_scene_sph_topview(obj, scene, THETA, THETA_end, PSI, PSI_end, rS, threshold, f, scene_dim, imgName)
            top_score = max(scene(:));
            scene = scene(:);
            
            color = [1 1 1];
            if contains(imgName,"_blue")
               color = [0 0 1];
            end
            if contains(imgName,"_green")
               color = [0 1 0];
            end
            if contains(imgName,"_red")
               color = [1 0 0];
            end
            
            if rS ~= 1
                "ERR - rS is not 1"
            end
            rho_0 = 0.6*scene_dim;
            for ind=1:length(scene(:))
                figure(f)
                if scene(ind) > threshold
                    
                    theta0 = THETA(ind); theta1 = THETA_end(ind); psi0 = PSI(ind); psi1 = PSI_end(ind);

                    y1 = (((rho_0).^2)./(1+(tan(theta0)).^2+(tan(psi0)).^2)).^0.5;
                    x1 = tan(theta0).*y1;

                    y2 = (((rho_0).^2)./(1+(tan(theta1)).^2+(tan(psi0)).^2)).^0.5;
                    x2 = tan(theta1).*y2;

                    y3 = (((rho_0).^2)./(1+(tan(theta0)).^2+(tan(psi1)).^2)).^0.5;
                    x3 = tan(theta0).*y3;

                    y4 = (((rho_0).^2)./(1+(tan(theta1)).^2+(tan(psi1)).^2)).^0.5;
                    x4 = tan(theta1).*y4;

                    fill([x1 x2 x3], [y1, y2, y3], 'w', 'FaceAlpha', (scene(ind))/top_score, 'EdgeAlpha', 0, "FaceColor", color);
                    hold on;
                    fill([x4 x2 x3], [y4, y2, y3], 'w', 'FaceAlpha', (scene(ind))/top_score, 'EdgeAlpha', 0, "FaceColor", color);

                    axis([0 scene_dim 0 scene_dim]);
                    grid on;
                    hold on;

                    ax = gca;
                    ax.LineWidth = 1;
                    ax.GridAlpha = 0.25;
                    axis([0 scene_dim 0 scene_dim]);
                    set(gca,'FontSize',14)
                    axis square

                end

                set(gca,'Color',[0.1 0.1 0.1]);


                set(gcf, 'InvertHardcopy', 'off');
                set(gca, 'XColor', [0.6 0.6 0.6])
                set(gca, 'YColor', [0.6 0.6 0.6])

                %rnd = 2;
                %colormap(gray(256));
                %colorbar('Ticks',[0,.2,.4,.6,.8,1],'TickLabels',{num2str(round(0*top_score,rnd)), num2str(round(0.2*top_score,rnd)),num2str(round(0.4*top_score,rnd)),num2str(round(0.6*top_score,rnd)),num2str(round(0.8*top_score,rnd)),num2str(round(1.0*top_score,rnd))})
                           
            end
            
            hold off;
        end
        
        function [nothing] = fill_scene_sph_topview_rho(obj, scene, aS, psiS, rho_0, threshold, f, scene_dim)
            top_score = max(scene(:));
            
            for theta=1:aS
                for psi=1:psiS
                    figure(f)
                    if scene(theta, psi) > threshold

                        y1 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        x1 = tan(pi/2*(theta-1)/aS).*y1;

                        y2 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi-1+aS-psiS)/aS)).^2)).^0.5;
                        x2 = tan(pi/2*(theta)/aS).*y2;

                        y3 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta-1)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        x3 = tan(pi/2*(theta-1)/aS).*y3;

                        y4 = (((scene_dim*rho_0).^2)./(1+(tan(pi/2*(theta)/aS)).^2+(tan(pi/2*(psi+aS-psiS)/aS)).^2)).^0.5;
                        x4 = tan(pi/2*(theta)/aS).*y4;

                        fill([x1 x2 x3], [y1, y2, y3], 'w', 'FaceAlpha', (scene(theta, psi))/top_score, 'EdgeAlpha', 0);
                        hold on;
                        fill([x4 x2 x3], [y4, y2, y3], 'w', 'FaceAlpha', (scene(theta, psi))/top_score, 'EdgeAlpha', 0);

                        axis([0 scene_dim 0 scene_dim]);
                        grid on;
                        hold on;

                        ax = gca;
                        ax.LineWidth = 1;
                        ax.GridAlpha = 0.25;
                        axis([0 scene_dim 0 scene_dim]);
                        set(gca,'FontSize',14)
                        axis square

                    end

                    set(gca,'Color',[0.1 0.1 0.1]);
                    set(gcf, 'InvertHardcopy', 'off');
                    set(gca, 'XColor', [0.6 0.6 0.6])
                    set(gca, 'YColor', [0.6 0.6 0.6])

                    xlabel("Y[m]")
                    ylabel("X[m]")
                    
                    %rnd = 2;
                    %colormap(gray(256));
                    %colorbar('Ticks',[0,.2,.4,.6,.8,1],'TickLabels',{num2str(round(0*top_score,rnd)), num2str(round(0.2*top_score,rnd)),num2str(round(0.4*top_score,rnd)),num2str(round(0.6*top_score,rnd)),num2str(round(0.8*top_score,rnd)),num2str(round(1.0*top_score,rnd))})

                end
            end
            
            hold off;
        end
        
        function [nothing] = plot_scene_continuous(obj, scene, aS, psiS, radial, threshold, f, scene_dim)
            unit = scene_dim;
            top_score = max(scene(:));
            
            for theta=1:aS
                for psi=1:psiS
                    figure(f)
                    if scene(theta, psi) > threshold

                        x = ((radial(theta, psi).^2)./(1+(tan(pi/2*(theta-.5)/aS)).^2+(tan(pi/2*(psi-.5+aS-psiS)/aS)).^2)).^0.5;
                        y = tan(pi/2*(theta-.5)/aS).*x;
                        z = unit - tan(pi/2*(psi-.5+aS-psiS)/aS).*x;

                        f1=scatter3(x, y, z, 'Marker', 'o', "MarkerEdgeAlpha", (scene(theta, psi))/top_score, "MarkerEdgeColor", 'k');
                        f1.SizeData = f1.SizeData/5;
                        %scatter3(x, y, z);
                        axis([0 unit 0 unit 0 unit]);
                        view(-130, 20)
                        %set(gca,'Color',[.3 .3 .3]);
                        set(gca,'FontSize',14)
                        xlabel('X');
                        ylabel('Y');
                        zlabel('Z');
                        %set(gcf, 'InvertHardcopy', 'off');
                        grid on;
                        hold on;
                    end
                end
            end
            
            set(gca, 'Ydir', 'reverse')
            hold off;
        end
        
        function [error] = sliding_compare(obj, truth, guess, size)
            avg_guess = zeros(S, S, S);
            for i=1:S
                for j=1:S
                    for k=1:S
                        window_size = 0;
                        for i0=max(1,i-size):min(S,i+size)
                            for j0=max(1,j-size):min(S,j+size)
                                for k0=max(1,k-size):min(S,k+size)
                                    avg_guess(i,j,k) = avg_guess(i,j,k) + guess(i0, j0, k0);
                                    window_size = window_size + 1;
                                end
                            end
                        end
                        avg_guess(i,j,k) = avg_guess(i,j,k)/window_size;
                    end
                end
            end
            
            avg_truth = zeros(S, S, S);
            for i=1:S
                for j=1:S
                    for k=1:S
                        window_size = 0;
                        for i0=max(1,i-size):min(S,i+size)
                            for j0=max(1,j-size):min(S,j+size)
                                for k0=max(1,k-size):min(S,k+size)
                                    avg_truth(i,j,k) = avg_truth(i,j,k) + truth(i0, j0, k0);
                                    window_size = window_size + 1;
                                end
                            end
                        end
                        avg_truth(i,j,k) = avg_truth(i,j,k)/window_size;
                    end
                end
            end
            
            error = get_error(avg_truth, avg_guess);
        end
        
        function [error] = sliding_compare_sph(obj, truth, guess, rS, aS, size)
            avg_guess = zeros(rS, aS, aS);
            for rho=1:rS
                for theta=1:aS
                    for psi=1:aS
                        window_size = 0;
                        for rho0=max(1,rho-size):min(rS,rho+size)
                            for theta0=max(1,theta-size):min(aS,theta+size)
                                for psi0=max(1,psi-size):min(aS,psi+size)
                                    avg_guess(rho,theta,psi) = avg_guess(rho,theta,psi) + guess(rho0, theta0, psi0);
                                    window_size = window_size + 1;
                                end
                            end
                        end
                        avg_guess(rho,theta,psi) = avg_guess(rho,theta,psi)/window_size;
                    end
                end
            end
            
            avg_truth = zeros(rS, aS, aS);
            for rho=1:rS
                for theta=1:aS
                    for psi=1:aS
                        window_size = 0;
                        for rho0=max(1,rho-size):min(rS,rho+size)
                            for theta0=max(1,theta-size):min(aS,theta+size)
                                for psi0=max(1,psi-size):min(aS,psi+size)
                                    avg_truth(rho,theta,psi) = avg_truth(rho,theta,psi) + truth(rho0, theta0, psi0);
                                    window_size = window_size + 1;
                                end
                            end
                        end
                        avg_truth(rho,theta,psi) = avg_truth(rho,theta,psi)/window_size;
                    end
                end
            end
            
            avg_guess = avg_guess(:);
            avg_truth = avg_truth(:);
            error = 0;
            for i = 1:length(avg_guess(:))
                error = error + (avg_guess(i) - avg_truth(i))^2;
            end
        end
        
        function [] = single_ellipse(obj, CRB, x, y, f)
            figure(f)
            [V, D] = eig(CRB);
            b = 1.177*sqrt(D(1,1)); %major axis
            a = 1.177*sqrt(D(2,2)); %minor axis
            rotation_ang = atan(V(2,2)/V(1,2));
            th = linspace(0,2*pi); 
            xe = a*cos(th)*cos(rotation_ang) - b*sin(th)*sin(rotation_ang) + x; 
            ye = a*cos(th)*sin(rotation_ang) + b*sin(th)*cos(rotation_ang) + y;
            plot(xe,ye,'b','LineWidth',1);
            axis square
            xlabel('y [yeee]','Interpreter','latex','FontSize',16)
            ylabel('x [m]','Interpreter','latex','FontSize',16)
            
            xlim([0 1])
            ylim([0 1])
        end
    end
end
        
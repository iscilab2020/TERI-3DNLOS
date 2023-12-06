%% Two-Edge-Resolved 3D NLOS Imaging (Supplementary Figures 10 and 11)
% Manuscript:
%   Czajkowski, R. and Murray-Bruce, J., 'Two-edge-resolved three-dimensional non-line-of-sight
% imaging with an ordinary camera', Nat. Commun., 2023.
% Code to create Supplementary Figure 6 from Nature Comm 2-Edge NLOS Paper
% Author: Robinson Czajkowski

clear;
close all;

addpath('Functions')

fun = HiddenScene;
indfun = @(x)(x>=0);
addpath("github_repo");

N = 124;
h = 0.2;
snr = 75;
image_size = 0.5;




fit = 20000; % Figure Number
N = 256; % Pixels in length of image
h = 0.2; % Height of door head
snr = 79; % Signal to noise ratio of the observation image
image_size = 0.5; % Length of image

%Point Source at middle of scene
x = 0.5;
y = 0.5;
z = 0.5;

%Generate Cartesian observation coordinates
aa = linspace(0, 1, N+1);
aa = aa + 1/(2*N);
aa = aa(1:end-1);

bb = linspace(0, 1, N+1);
bb = bb + 1/(2*N);
bb = bb(1:end-1);
[AA, BB] = meshgrid(aa, bb);
a = AA(:);
b = BB(:);

a = a*image_size;
b = b*image_size;

f = 1/(2*N)*image_size;
rho = (x^2+y^2+z^2)^0.5;
theta = atan(y/x);
psi = atan(z/x);

phi = acos(z/rho);


% Calculate all relavent values for equations
[R, dR_dTheta, dR_dPsi, dR_dRho, H, H1, dH_dTheta, dH1_dTheta, dH_dPsi, SC, dC_dx, dC_dy, dC_dz] = compute_stuff(rho, theta, psi, a, b, f, h, N);

M_t = R.*dH_dTheta + H.*dR_dTheta;
M_p = R.*dH_dPsi + H.*dR_dPsi;
M_r = dR_dRho.*H;
M_c = R.*H;
sigma_0 = (10^(-0.1 * snr) * (sum(4*f^2*R(:))) / N)^0.5;
sigma_1 = (10^(-0.1 * snr) * (sum(H1(:).*R(:))) / N)^0.5;
sigma_2 = (10^(-0.1 * snr) * (sum(H(:).*R(:))) / N)^0.5;


M_t = reshape(M_t, N, N);
M_p = reshape(M_p, N, N);
M_r = reshape(M_r, N, N);
M_c = reshape(M_c, N, N);

% Conversion from PESC to SC
Conversion_p2s = [1               0                                    0 0;
                  0    -1/(sin(phi)^2 + (cos(theta)^2)*(cos(phi)^2))   0 0;
                  0               0                                    1 0;
                  0               0                                    0 1]';

Conversion = [1/(1+y^2/x^2)*(-y/x^2) 1/(1+z^2/x^2)*(-z/x^2)   x/(x^2+y^2+z^2)^0.5 0;
              1/(1+y^2/x^2)*(1/x)    0                        y/(x^2+y^2+z^2)^0.5 0;
              0                      1/(1+z^2/x^2)*(1/x)      z/(x^2+y^2+z^2)^0.5 0;
              0                      0                        0                   1]';


Conversion1 = [Conversion(:,4) Conversion(:,3) Conversion(:,1:2)];

Conversion_c2s = [-rho*sin(theta)*sin(phi)   rho*cos(theta)*cos(phi) cos(theta)*sin(phi)    0;
                  rho*cos(theta)*sin(phi)    rho*sin(theta)*cos(phi) sin(theta)*sin(phi)    0;
                  0                          -rho*sin(phi)           cos(phi)               0;
                  0                          0                       0                      1];

%
fi_cmap = [-1 1]*1.5*10^(6);

% CRB
Delta_I = [M_t(:) M_p(:) M_r(:) M_c(:)];
Fisher = (1/sigma_2^2)*(Delta_I'*Delta_I);
CRB_sph2 = Fisher^-1;
figure(fit);
imagesc(Fisher,fi_cmap);
box on;
font_size = 20;
set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4,'LineWidth',1.5)
set(gca,'TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
% set(gca,'TickLabel',{'\rho'})
xticks([1 2 3 4])
yticks([1 2 3 4])
xticklabels({'$\theta$','$\psi$','$\rho$','$c$'})
yticklabels({'$\theta$','$\psi$','$\rho$','$c$'})
font_size = 20;
title('Fisher Information Matrix - PESC')
axis equal square tight
hcb = colorbar; 
set(hcb,'ticklabelinterpreter','latex','LineWidth',1.5,'fontsize',font_size)

crb_cmap = [0 5]*10^-3;
%
figure(fit+1);
imagesc(CRB_sph2,crb_cmap);
title('Cramer-Rao Bound - PESC')
box on;

set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4,'LineWidth',1.5)
set(gca,'TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
% set(gca,'TickLabel',{'\rho'})
xticks([1 2 3 4])
yticks([1 2 3 4])
xticklabels({'$\theta$','$\psi$','$\rho$','$c$'})
yticklabels({'$\theta$','$\psi$','$\rho$','$c$'})

axis equal square tight
hcb = colorbar; 
set(hcb,'ticklabelinterpreter','latex','LineWidth',1.5,'fontsize',font_size)
% hcb.Label.Interpreter = 'latex';


%%
Delta_I_Cart = Delta_I*Conversion;
Fisher_Cart = (1/sigma_2^2)*(Delta_I_Cart'*Delta_I_Cart);
CRB_cart2 = Fisher_Cart^-1;
figure(fit+2);
imagesc(Fisher_Cart,fi_cmap);
title('Fisher Information Matrix - Cartesian')
box on;
font_size = 20;
set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4,'LineWidth',1.5)
set(gca,'TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
% set(gca,'TickLabel',{'\rho'})
xticks([1 2 3 4])
yticks([1 2 3 4])
xticklabels({'$x$','$y$','$z$','$c$'})
yticklabels({'$x$','$y$','$z$','$c$'})

axis equal square tight
hcb = colorbar; 
set(hcb,'ticklabelinterpreter','latex','LineWidth',1.5)


figure(fit+3);
imagesc(CRB_cart2,crb_cmap);
box on;
font_size = 20;
set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4,'LineWidth',1.5)
title('Cramer-Rao Bound - Cartesian')
set(gca,'TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
% set(gca,'TickLabel',{'\rho'})
xticks([1 2 3 4])
yticks([1 2 3 4])
xticklabels({'$x$','$y$','$z$','$c$'})
yticklabels({'$x$','$y$','$z$','$c$'})

axis equal square tight
hcb = colorbar; 
set(hcb,'ticklabelinterpreter','latex','LineWidth',1.5,'fontsize',font_size)


%%
% True Spherical
Delta_I_Cart = Delta_I*Conversion*Conversion_c2s;
Fisher_True = (1/sigma_1^2)*(Delta_I_Cart'*Delta_I_Cart);
CRB_True2 = Fisher_True^-1;

figure(fit+4);
imagesc(Fisher_True,fi_cmap);
box on;
font_size = 20;
set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4,'LineWidth',1.5)

title('Fisher Information Matrix - Spherical')
set(gca,'TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
% set(gca,'TickLabel',{'\rho'})
xticks([1 2 3 4])
yticks([1 2 3 4])
xticklabels({'$\theta$','$\varphi$','$\rho$','$c$'})
yticklabels({'$\theta$','$\varphi$','$\rho$','$c$'})

axis equal square tight
hcb = colorbar; 
set(hcb,'ticklabelinterpreter','latex','LineWidth',1.5)


figure(fit+5);
imagesc(CRB_True2,crb_cmap);
box on;
font_size = 20;
set(gca,'Color',[0.2 0.2 0.2],'FontSize',font_size-4,'LineWidth',1.5)
title('Cramer-Rao Bound - Spherical')

set(gca,'TickLabelInterpreter','Latex','GridColor',[1 1 1],'GridAlpha',0.5)
% set(gca,'TickLabel',{'\rho'})
xticks([1 2 3 4])
yticks([1 2 3 4])
xticklabels({'$\theta$','$\varphi$','$\rho$','$c$'})
yticklabels({'$\theta$','$\varphi$','$\rho$','$c$'})

axis equal square tight
hcb = colorbar; 
set(hcb,'ticklabelinterpreter','latex','LineWidth',1.5,'fontsize',font_size)




function [R, dR_dTheta, dR_dPsi, dR_dRho, H, H1, dH_dTheta, dH1_dTheta, dH_dPsi, SC, dC_dx, dC_dy, dC_dz] = compute_stuff(rho, theta, psi, a, b, f, h, N)

    S = 1/(1+(tan(theta))^2+(tan(psi))^2)^(1/2);
    R = (rho^2 + 2*rho*(a*tan(theta)+b+h*tan(psi))*S + a.^2+b.^2+h^2).^-1;
    
    dR_dTheta = -R.^2 *2*rho .*( (S.*a.*(sec(theta))^2 - (b+a*tan(theta)+h*tan(psi))*tan(theta)*(sec(theta))^2*S^3) );
    dR_dPsi =   -R.^2 *2*rho .*( (S.*h.*(sec(psi))^2   - (b+a*tan(theta)+h*tan(psi))*tan(psi)*(sec(psi))^2*S^3) );
    dR_dRho =   -R.^2 .* (2*rho + 2*(a*tan(theta)+b+h*(tan(psi)))*S);
    
    % Create cartesian directly
    x = rho.*S;
    y = tan(theta).*x;
    z = tan(psi).*x;
    
    SC = ((x+b).^2 + (y+a).^2 + (z+h).^2).^(-1);
    dC_dx = -2*(x+b).*SC.^2;
    dC_dy = -2*(y+a).*SC.^2;
    dC_dz = -2*(z+h).*SC.^2;
    
    
    H = 1/2 * (max(min(b+f,(a+f)*cot(theta)),max(b-f,h*cot(psi)))).^2 * tan(theta) - max(min(b+f,(a+f)*cot(theta)),max(b-f,h*cot(psi))).*(a-f) ...
        - 1/2 * (max(b-f, max(h*cot(psi), (a-f)*cot(theta)))).^2 * tan(theta) ...
        + max(b-f, max(h*cot(psi), (a-f)*cot(theta))).*(a-f) + 2*f*(b+f) ...
        - 2*f*max(min(b+f,(a+f)*cot(theta)),max(b-f,h*cot(psi)));
    for i=1:length(a(:))
        a0 = a(i);
        b0 = b(i);
        
        if h*cot(psi) > b0+f || (a0-f)*cot(theta) > b0+f
            H(i) = 0;
        end
    end
    
    H1 = 1/2 * (max(min(b+f,(a+f)*cot(theta)),b-f)).^2 * tan(theta) - max(min(b+f,(a+f)*cot(theta)),b-f).*(a-f) ...
        - 1/2 * (max(b-f, (a-f)*cot(theta))).^2 * tan(theta) ...
        + max(b-f, (a-f)*cot(theta)).*(a-f) + 2*f*(b+f) ...
        - 2*f*max(min(b+f,(a+f)*cot(theta)),b-f);
    for i=1:length(a(:))
        a0 = a(i);
        b0 = b(i);
        
        if (a0-f)*cot(theta) > b0+f
            H1(i) = 0;
        end
    end
    
    dH_dTheta = zeros(1,N^2);
    for i=1:length(a(:))
        a0 = a(i);
        b0 = b(i);
        
        if (a0+f)*cot(theta) <= b0-f || (a0-f)*cot(theta) >= b0+f || h*cot(psi) >= b0+f || h*cot(psi) > (a0+f)*cot(psi)
            dH_dTheta(i) = 0;
            continue
        end
        
        if b0+f > (a0+f)*cot(theta)
            if (a0-f)*cot(theta) > max(b0-f,h*cot(psi))
                dH_dTheta(i) = 0.5*(a0+f)^2*(csc(theta))^2 - 0.5*(a0-f)^2*(csc(theta))^2;
            else
                dH_dTheta(i) = 0.5*(a0+f)^2*(csc(theta))^2 - 0.5*(max(b0-f,h*cot(psi)))^2*(sec(theta))^2;
            end
        else
            if (a0-f)*cot(theta) > max(b0-f,h*cot(psi))
                dH_dTheta(i) = 0.5*(b0+f)^2*(sec(theta))^2 - 0.5*(a0-f)^2*(csc(theta))^2;
            else
                dH_dTheta(i) = 0.5*(b0+f)^2*(sec(theta))^2 - 0.5*(max(b0-f,h*cot(psi)))^2*(sec(theta))^2;
            end
        end
    end
    
    dH_dTheta = -(1/2)*(min(b+f,(a+f)*cot(theta))).^2 *(sec(theta))^2+(1/2)*(max(max(b-f,h*cot(psi)),(a-f)*cot(theta))).^2 *(sec(theta))^2;
    dH_dTheta = (1/2)*(sec(theta))^2 * ( (max(max(b-f,h*cot(psi)),(a-f)*cot(theta))).^2 - (min(max(b+f,h*cot(psi)),(a+f)*cot(theta))).^2  );
    for i=1:length(a(:))
        a0 = a(i);
        b0 = b(i);
        if (a0+f)*cot(theta) <= max(b0-f,h*cot(psi)) || (a0-f)*cot(theta) >= max(b0+f,h*cot(psi))
            dH_dTheta(i) = 0;
        end
    end
    
    dH1_dTheta = zeros(1,N^2);
    for i=1:length(a(:))
        a0 = a(i);
        b0 = b(i);
        
        if (a0+f)*cot(theta) <= b0-f || (a0-f)*cot(theta) >= b0+f
            dH1_dTheta(i) = 0;
            continue
        end
        
        if b0+f > (a0+f)*cot(theta)
            if (a0-f)*cot(theta) > b0-f
                dH1_dTheta(i) = 0.5*(a0+f)^2*(csc(theta))^2 - 0.5*(a0-f)^2*(csc(theta))^2;
            else
                dH1_dTheta(i) = 0.5*(a0+f)^2*(csc(theta))^2 - 0.5*(b0-f)^2*(sec(theta))^2;
            end
        else
            if (a0-f)*cot(theta) > b0-f
                dH1_dTheta(i) = 0.5*(b0+f)^2*(sec(theta))^2 - 0.5*(a0-f)^2*(csc(theta))^2;
            else
                dH1_dTheta(i) = 0.5*(b0+f)^2*(sec(theta))^2 - 0.5*(b0-f)^2*(sec(theta))^2;
            end
        end
    end
    
    dH_dPsi = h*(csc(psi))^2 * (min(a-f, h*cot(psi)*tan(theta)) - min(a+f, h*cot(psi)*tan(theta)));
    for i=1:length(a(:))
        a0 = a(i);
        b0 = b(i);
        
        if b0-f > h*cot(psi) || h*cot(psi) > b0+f
            dH_dPsi(i) = 0;
        end
    end
    
    
    R = reshape(R, N, N);
    dR_dTheta = reshape(dR_dTheta, N, N);
    dR_dPsi = reshape(dR_dPsi, N, N);
    dR_dRho = reshape(dR_dRho, N, N);
    H = reshape(H, N, N);
    H1 = reshape(H1, N, N);
    dH_dTheta = reshape(dH_dTheta, N, N);
    dH1_dTheta = reshape(dH1_dTheta, N, N);
    dH_dPsi = reshape(dH_dPsi, N, N);
    
    SC = reshape(SC,N,N);
    dC_dx = reshape(dC_dx,N,N);
    dC_dy = reshape(dC_dy,N,N);
    dC_dz = reshape(dC_dz,N,N);
end

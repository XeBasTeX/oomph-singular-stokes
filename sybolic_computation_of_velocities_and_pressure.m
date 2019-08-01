% The University of Manchester
% Bastien Laville - 22/07/2019
% This little program provides the symbolic computation of the pressure

%% Symbolic computation of u_r, u_theta (or u_phi) and p the pressure


% Define our symbolic variables
% The variable A stands for \alpha and L stands for \lambda
syms A L K r theta 

% % Just in case we want numerical values
% lambda = 1.544;
% L = lambda;
% A = 3*pi/4;
% K = 1;

streamfunction = K*(r^L)*(cos(A*L)*cos(theta*(L - 2)) - cos(A*(L - 2))*cos(theta*L));

% Definition of u_r
%u_r = K*r^(L - 1)*((L - 2)*cos(A*L)*sin((theta*(L - 2))) - L*cos(A*(L - 2))*sin((theta*L)));
u_r = diff(streamfunction, theta)/r;

% 1/r^2 * du_r/dtheta
d_theta_of_u_r = 2*diff(u_r, theta)/r^2;

% Non vectorial polar laplacian of u_r
laplacianOfU_r = K*r^(L - 3)*((2*L^2 - 7*L + 6)*cos(A*L)*sin(theta*(L - 2)) + L*(2*L - 1)*cos(A*(L - 2))*sin(theta*L));



% Definition of u_theta
%u_theta = - K*L*r^(L - 1)*(cos(A*L)*cos((theta*(L - 2))) - cos(A*(L - 2))*cos((theta*L)));
u_theta = - diff(streamfunction, r);

% 1/r^2 * du_theta/dtheta
d_theta_of_u_theta = 2*diff(u_theta, theta)/r^2;

% Non vectorial polar laplacian of u_theta
laplacianOfU_theta = - L*K*r^(L - 3)*((2*L - 3)*(cos(A*L)*cos(theta*(L - 2)) + (2*L - 1)*cos(A*(L - 2))*cos(theta*L)));



% Sum up to get the complete laplacian
completeLaplacianOfU_r = laplacianOfU_r - d_theta_of_u_theta - u_r/(r^2);
completeLaplacianOfU_theta = laplacianOfU_theta + d_theta_of_u_r - u_theta/(r^2);

% We simplify and multiply by variables in order to identify the laplacian
% with the corresponding pressure gradient (the's no mu since mu = 1)
shorterCompleteLaplacianOfU_r = simplify(completeLaplacianOfU_r);
shorterCompleteLaplacianOfU_theta = simplify(r * completeLaplacianOfU_theta);

% Deduce the pressure thanks to an integration w.r.t 1 variable
pressureAccordingToRadius = simplify(int(shorterCompleteLaplacianOfU_r, r));
pressureAccordingToAngle = simplify(int(shorterCompleteLaplacianOfU_theta, theta));

%% Compute the pressure on a grid

N = 100;
data = zeros(N);
dataU_r = zeros(N);
dataU_theta = zeros(N);

for i=1:N
    for j=1:N
        x = (j/N - 1/2);
        y = -(i/N - 1/2);
        radius = sqrt(x^2 + y^2);
        if radius ~= 0
            phi = atan2(y, x) - pi/4;
            data(i ,j) = subs(subs(pressureAccordingToRadius, r, radius), theta, phi); 
        end
    end
end


%%

h = 0.01;
[x,y] = meshgrid(-1:h:1,-1:h:1);
radius = sqrt(x.^2 + y.^2);
phi = atan2(y, x) - pi/4;

P = -(1805800404753236194639*sin((57*phi)/125) + 6914516165136002774400*sin((193*phi)/125))./(1026820715040473088000*radius.^(57/125));
P(20201) = 0;
PGradient = gradient(P / 100000);

U = (1805800404753236194639*sin((57*phi)/125) + 6914516165136002774400*sin((193*phi)/125))/(2251799813685248000000*radius.^(182/125));
V = (1805800404753236099072*sin((57*phi)/125) + 6914516165136003104768*sin((193*phi)/125))/(2251799813685248000000*radius.^(182/125));

P(20201) = 0;
U(20201) = 0;
V(20201) = 0;


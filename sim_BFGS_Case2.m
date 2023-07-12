%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Case 2: Muscle is NOT wrapped at some intervals of the simulation
%%%%%%% BFGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

% Cylinder specifications
Radius = 1;
Height = 5;
SideCount = 20; % for visual purposes

ti = (0:0.012:3)';                                   % Time data

% Surface's Motion data
position = [0*ti, 2*sin(2*pi*ti), 0*ti];             % Position data
d_position = [0*ti, 2*pi*2*cos(2*pi*ti), 0*ti];      % Velocity data
angles = [0*ti, 0*ti, 0*ti];                         % Orientation data (XYZ Euler angles)
d_angles = [0*ti, 0*ti, 0*ti];                       % Angular velocity data

p0_i = [2,0,-2.5];      % Muscle's origin point position  
q0_i = [-2,0,-0.5];     % Muscle's insertion point position
d_p0_i = [0,0,0];       % Muscle's origin point velocity 
d_q0_i = [0,0,0];       % Muscle's insertion point velocity

% Transform the muscle's origin and insertion points to surface's
% coordinate frame
p0_l = zeros(1,3);
q0_l = zeros(1,3);
Rc_w = cell(size(ti,1),1);
for i=1:size(ti,1)
    Rc_w{i} = eulerXYZ(angles(i,1), angles(i,2), angles(i,3));
    p0_l(i,:) = (Rc_w{i}'*(p0_i-position(i,:))')';
    q0_l(i,:) = (Rc_w{i}'*(q0_i-position(i,:))')';
end

%% Solve the wrapping
% Definition of initial conditions
u = 0; v = 0; lg = eps;
d = [1 1]; dn = d/norm(d);
up = dn(1); vp = dn(2);

% Initial points over the surface
p_l = [Radius*cos(u), Radius*sin(u), v];
q_l = [Radius*cos(u), Radius*sin(u), v];

% Stack the initial conditions to be updated by the BFGS method
ic = [u;v;up;vp;lg];

% Solve the muscle wrapping problem
out = {};
for i=1:length(ti)
    out{i} = SolveWrap(ic,p0_l(i,:),q0_l(i,:),p_l,q_l,Radius);
    % Previous result is stored and used in the next sampling time
    ic = out{1,i}{1,1};
    pl = out{1,i}{1,3};
    ql = out{1,i}{1,4};
end

%% Transform origin and insertion points to world coordinates
for i=1:length(ti)
    Rc_ = Rc_w{i};
    out_lp{i} = out{1,i}{1,3}; % origin point expressed in surface coordinate frame
    out_lq{i} = out{1,i}{1,4}; % insertion point expressed in surface coordinate frame    
    out{1,i}{1,3} = (position(i,:)' + Rc_*out{1,i}{1,3}')'; % transform origin point
    out{1,i}{1,4} = (position(i,:)' + Rc_*out{1,i}{1,4}')'; % transform insertion point
    out{1,i}{1,2} = (position(i,:)' + Rc_*out{1,i}{1,2}')'; % transform geodesic points
end

%% Compute muscle length and its time derivative
lmt = {};
dlmt = {};
for i=1:length(ti)
    if out{1,i}{1,6}>0 % Muscle is in a wrapping state
        % Total muscle path length
        lmt{i} = norm(p0_i - out{1,i}{1,3}) + out{1,i}{1,1}(5) + norm(q0_i - out{1,i}{1,4});
        % Straigth line segment from the origin point to the surface
        e = out{1,i}{1,3} - p0_i; e = e/norm(e);
        % Straigth line segment from the surface to the insertion point
        e2 = q0_i - out{1,i}{1,4}; e2 = e2/norm(e2);
        % Muscle velocity
    	d_c_p = skew(d_angles(i,:))*Rc_w{i}*out_lp{i}';
        d_c_q = skew(d_angles(i,:))*Rc_w{i}*out_lq{i}';
        dlmt{i} = dot(e,d_position(i,:)') + dot(e,d_c_p) - dot(e2,d_position(i,:)') - dot(e2,d_c_q) - dot(e,d_p0_i) + dot(e2,d_q0_i) ;
    else % Muscle is in a non-wrapping state
        % Total muscle path length
        lmt{i} = norm(q0_i - p0_i);
        % Straigth line segment from the origin point to the insertion
        % point
        e = q0_i - p0_i; e = e/norm(e);
        % Muscle velocity
        dlmt{i} = dot(e,d_q0_i - d_p0_i);
    end
end

%% Store iterations per time step
for i=1:length(ti)
    iter(i) = out{1,i}{1,5};
    iter2(i) = out{1,i}{1,7};
end

%% Animation initial figure
% Initial vertices and faces of the surface
r_ini = position(1,:)';
R_ini = eulerXYZ(angles(1,1), angles(1,2), angles(1,3));
[vertices_ini, sideFaces_ini, bottomFaces_ini] = calcCylinder(r_ini, R_ini, Radius, Height, SideCount);

% Draw initial figure
figure(1);
hold on;
% Plot the surface
h_side = patch('Faces', sideFaces_ini, 'Vertices', vertices_ini, 'FaceColor', 'w','FaceAlpha',.0);
h_bottom = patch('Faces', bottomFaces_ini, 'Vertices', vertices_ini, 'FaceColor', 'w','FaceAlpha',.0);
% Plot the origin and insertion points
h_fixedp0 = plot3(p0_i(1),p0_i(2),p0_i(3),'o','LineWidth',0.05,'MarkerEdgeColor',[0.576,0.121,0.137],'MarkerFaceColor',[0.576,0.121,0.137],'MarkerSize',4);
h_fixedq0 = plot3(q0_i(1),q0_i(2),q0_i(3),'o','LineWidth',0.05,'MarkerEdgeColor',[0.576,0.121,0.137],'MarkerFaceColor',[0.576,0.121,0.137],'MarkerSize',4);
% Plot the geodesic
h_geod = plot3(out{1,1}{1,2}(:,1), out{1,1}{1,2}(:,2), out{1,1}{1,2}(:,3),'r','LineWidth',2);
% Plot the points on the surface
h_movp = plot3(out{1,1}{1,3}(1), out{1,1}{1,3}(2), out{1,1}{1,3}(3),'o','LineWidth',0.05,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1],'MarkerSize',4);
h_movq = plot3(out{1,1}{1,4}(1), out{1,1}{1,4}(2), out{1,1}{1,4}(3),'o','LineWidth',0.05,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1],'MarkerSize',4);
% Plot straight line segments
h_linp = plot3([p0_i(1) out{1,1}{1,3}(1)],[p0_i(2) out{1,1}{1,3}(2)],[p0_i(3) out{1,1}{1,3}(3)],'k','LineWidth',2);
h_linq = plot3([q0_i(1) out{1,1}{1,4}(1)],[q0_i(2) out{1,1}{1,4}(2)],[q0_i(3) out{1,1}{1,4}(3)],'k','LineWidth',2);
hold off;

% Axes settings
xlabel('x [cm]','Interpreter','latex'); ylabel('y [cm]','Interpreter','latex'); zlabel('z [cm]','Interpreter','latex');
axis equal;
view([140,35]);
grid on;
xlim([-3,3]);
ylim([-3,3]);
zlim([-3,3]);

%% Animation
% Compute vertices at each sampling time
vertices = zeros(size(vertices_ini,1),3,length(ti));
for i = 1:length(ti)
    r = position(i,:)';
    R = eulerXYZ(angles(i,1), angles(i,2), angles(i,3));
    [vertices(:,:,i), ~, ~] = calcCylinder(r, R, Radius, Height, SideCount);
end

% Create & open video writer
vw = VideoWriter('animation_BFGS_NW.avi');
vw.Quality = 100;
vw.FrameRate = 60;
open(vw);

% Animation Loop
for i = 1:length(ti)
    set(h_side, 'Vertices', vertices(:,:,i));
    set(h_bottom, 'Vertices', vertices(:,:,i));
    
    if out{1,i}{1,6}>0 % Muscle is in a wrapping state
        set(h_geod, 'XData', out{1,i}{1,2}(:,1),'YData', out{1,i}{1,2}(:,2),'ZData', out{1,i}{1,2}(:,3));

        set(h_movp, 'XData', out{1,i}{1,3}(1),'YData', out{1,i}{1,3}(2),'ZData', out{1,i}{1,3}(3));
        set(h_movq, 'XData', out{1,i}{1,4}(1),'YData', out{1,i}{1,4}(2),'ZData', out{1,i}{1,4}(3));

        set(h_linp, 'XData', [p0_i(1) out{1,i}{1,3}(1)],'YData', [p0_i(2) out{1,i}{1,3}(2)],'ZData', [p0_i(3) out{1,i}{1,3}(3)]);
        set(h_linq, 'XData', [q0_i(1) out{1,i}{1,4}(1)],'YData', [q0_i(2) out{1,i}{1,4}(2)],'ZData', [q0_i(3) out{1,i}{1,4}(3)]);
        
        set(h_geod,'visible','on');
        set(h_movp,'visible','on');
        set(h_movq,'visible','on');
        set(h_linq,'visible','on');
        
    else % Muscle is in a non-wrapping state
        set(h_geod,'visible','off');
        set(h_movp,'visible','off');
        set(h_movq,'visible','off');
        set(h_linq,'visible','off');
        set(h_linp, 'XData', [p0_i(1) q0_i(1)],'YData', [p0_i(2) q0_i(2)],'ZData', [p0_i(3) q0_i(3)]);
    end
    pause(0.05);
    drawnow;
    
    % Write each frame to video
    frame = getframe(gcf);
    writeVideo(vw, frame);
    
end

% Close video writer
close(vw);

%% Plot results
% Plot Muscle length
a_lmt = cell2mat(lmt);
figure(2); title('\textbf{Total muscle-tendon length}','Interpreter','latex','FontSize',12); 
xlabel('t [s]','Interpreter','latex','FontSize',12); ylabel('$l_{MT}$ [cm]','Interpreter','latex','FontSize',12); %ylabel('$\lambda$','Interpreter','latex');
hold on; plot(ti,a_lmt,'LineWidth',1.0); hold off;
leg1 = legend('$l_{MT}$'); set(leg1,'Interpreter','latex'); grid on;

% Plot Muscle velocity
a_dlmt = cell2mat(dlmt);
figure(3); title('\textbf{Muscle-tendon velocity}','Interpreter','latex','FontSize',12); 
xlabel('t [s]','Interpreter','latex','FontSize',12); ylabel('$\dot{l}_{MT}$ [cm/s]','Interpreter','latex','FontSize',12); %ylabel('$\lambda$','Interpreter','latex');
hold on; plot(ti,a_dlmt,'LineWidth',1.0); hold off;
leg1 = legend('$\dot{l}_{MT}$'); set(leg1,'Interpreter','latex'); grid on;

%% Plot specific surface's positions in time
color_e{1} = [228/255, 145/255, 3/255];

figure(8);
hold on;
% Plot the surface
h_side = patch('Faces', sideFaces_ini, 'Vertices', vertices_ini, 'FaceColor', 'w','FaceAlpha',.0,'EdgeColor',color_e{1});
h_bottom = patch('Faces', bottomFaces_ini, 'Vertices', vertices_ini, 'FaceColor', 'w','FaceAlpha',.0,'EdgeColor',color_e{1});
% Plot the origin and insertion points
h_fixedp0 = plot3(p0_i(1),p0_i(2),p0_i(3),'o','LineWidth',0.05,'MarkerEdgeColor',[0.576,0.121,0.137],'MarkerFaceColor',[0.576,0.121,0.137],'MarkerSize',4);
h_fixedq0 = plot3(q0_i(1),q0_i(2),q0_i(3),'o','LineWidth',0.05,'MarkerEdgeColor',[0.576,0.121,0.137],'MarkerFaceColor',[0.576,0.121,0.137],'MarkerSize',4);
% Plot the geodesic
h_geod = plot3(out{1,1}{1,2}(:,1), out{1,1}{1,2}(:,2), out{1,1}{1,2}(:,3),'r','LineWidth',2);
% Plot the points on the surface
h_movp = plot3(out{1,1}{1,3}(1), out{1,1}{1,3}(2), out{1,1}{1,3}(3),'o','LineWidth',0.05,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1],'MarkerSize',4);
h_movq = plot3(out{1,1}{1,4}(1), out{1,1}{1,4}(2), out{1,1}{1,4}(3),'o','LineWidth',0.05,'MarkerEdgeColor','k','MarkerFaceColor',[0,0,1],'MarkerSize',4);
% Plot straight line segments
h_linp = plot3([p0_i(1) out{1,1}{1,3}(1)],[p0_i(2) out{1,1}{1,3}(2)],[p0_i(3) out{1,1}{1,3}(3)],'k','LineWidth',2);
h_linq = plot3([q0_i(1) out{1,1}{1,4}(1)],[q0_i(2) out{1,1}{1,4}(2)],[q0_i(3) out{1,1}{1,4}(3)],'k','LineWidth',2);

color_e{2} = [38/255, 219/255, 0/255];
color_e{3} = [148/255, 0/255, 219/255];
color_l = [0/255, 0/255, 0/255];

color_g = [255/255, 0/255, 0/255];
color_p = [0/255, 0/255, 255/255];

plot_vec = [22, 64];
for j = 1:length(plot_vec)
    i = plot_vec(j);
    patch('Faces', sideFaces_ini, 'Vertices', vertices(:,:,i), 'FaceColor', 'w','FaceAlpha',.0,'EdgeColor',color_e{j+1});
    patch('Faces', bottomFaces_ini, 'Vertices', vertices(:,:,i), 'FaceColor', 'w','FaceAlpha',.0,'EdgeColor',color_e{j+1});
    if out{1,i}{1,6}>0
        % Plot the geodesic
        plot3(out{1,i}{1,2}(:,1), out{1,i}{1,2}(:,2), out{1,i}{1,2}(:,3),'Color',color_g,'LineWidth',2);
        % Plot the points on the surface
        plot3(out{1,i}{1,3}(1), out{1,i}{1,3}(2), out{1,i}{1,3}(3),'o','LineWidth',0.05,'MarkerEdgeColor',color_l,'MarkerFaceColor',color_p,'MarkerSize',4);
        plot3(out{1,i}{1,4}(1), out{1,i}{1,4}(2), out{1,i}{1,4}(3),'o','LineWidth',0.05,'MarkerEdgeColor',color_l,'MarkerFaceColor',color_p,'MarkerSize',4);
        % Plot straight line segments
        plot3([p0_i(1) out{1,i}{1,3}(1)],[p0_i(2) out{1,i}{1,3}(2)],[p0_i(3) out{1,i}{1,3}(3)],'Color',color_l,'LineWidth',2);
        plot3([q0_i(1) out{1,i}{1,4}(1)],[q0_i(2) out{1,i}{1,4}(2)],[q0_i(3) out{1,i}{1,4}(3)],'Color',color_l,'LineWidth',2);
    else
        % plot line segments
        plot3([p0_i(1) q0_i(1)],[p0_i(2) q0_i(2)],[p0_i(3) q0_i(3)],'Color',color_l,'LineWidth',2);
    end
end

% Text labes indicating the position of the surface in y-axis
text(-0.87517945820106,-1.85038506869,2.82777834204219,'$^{(0)}\mbox{\boldmath $c$}_{y}=-1.99$cm','Interpreter','latex','FontSize',13,'Color',color_e{3});
text(-1.392931800787778,-0.9191709262264,2.284301418722947,'$^{(0)}\mbox{\boldmath $c$}_{y}=0.00$cm','Interpreter','latex','FontSize',13,'Color',color_e{1});
text(-1.882273679386159,0.033805806287347,1.69093502209958,'$^{(0)}\mbox{\boldmath $c$}_{y}=1.99$cm','Interpreter','latex','FontSize',13,'Color',color_e{2});

hold off;

% Axes settings
xlabel('x [cm]','Interpreter','latex'); ylabel('y [cm]','Interpreter','latex'); zlabel('z [cm]','Interpreter','latex');
axis equal;
view([140,35]);
grid on;
xlim([-4.697307824781931,2.129358841884733]);
ylim([-3.084685131091947,4.310870424463605]);
zlim([-2.689521471206549,4.137145195460118]);

%% user functions
function ydot = int_gp(~,y,r)
% Differential equations of the geodesic
u =  y(1);
v =  y(2);

ru = [-r*sin(u);r*cos(u);0];
rv = [0;0;1];

E = dot(ru,ru); Eu = 2*r^2*cos(u)*sin(u) - 2*r^2*cos(u)*sin(u); Ev = 0;
F = dot(ru,rv); Fu = 0; Fv = 0;
G = dot(rv,rv); Gu = 0; Gv = 0;

C111 = (G*Eu - 2*F*Fu + F*Ev)/(2*(E*G-F^2));
C121 = (G*Ev - F*Gu)/(2*(E*G-F^2));
C221 = (2*G*Fv - G*Gu + F*Gv)/(2*(E*G-F^2));

C112 = (2*E*Fu - E*Ev + F*Eu)/(2*(E*G-F^2));
C122 = (E*Gu - F*Ev)/(2*(E*G-F^2));
C222 = (E*Gv - 2*F*Fv + F*Gu)/(2*(E*G-F^2));

% Enforce the unit speed curve
sq = norm(ru*y(3) + rv*y(4));
if sq > 1+eps
    t = ru*y(3) + rv*y(4); 
    par = [y(3);y(4)]./norm(t);
    par = par .* (1/(E*par(1)^2+2*F*par(1)*par(2)+G*par(2)^2));

    y(3) = par(1);
    y(4) = par(2);
end

ydot = zeros(4,1);
ydot(1) = y(3);
ydot(2) = y(4);
ydot(3) = -C111*ydot(1)^2 - 2*C121*ydot(1)*ydot(2) - C221*ydot(2)^2;
ydot(4) = -C112*ydot(1)^2 - 2*C122*ydot(1)*ydot(2) - C222*ydot(2)^2;
end

function t = tang(u,v,up,vp,r)
% Tangent vector of a geodesic over a cylindrical surface
rv = [0;0;1];
t = zeros(length(u),3);
for i=1:length(u)
    ru = [-r*sin(u(i));r*cos(u(i));0];
    t(i,:) = ru*up(i) + rv*vp(i);
end

end

function J = jacobi(f,y0,x,r,p0,q0)
% Evaluation of the numerical Jacobian using central differences
    delta = 1e-10;
    n = length(y0);
    m = length(x);
    J = zeros(n,m);
    for i=1:m
        dx = zeros(m,1);
        dx(i) = delta/2;
        J(:,i) = (f(x+dx,r,p0,q0) - f(x-dx,r,p0,q0))/delta;
    end

J = J';
end

function [F,a1,a2,term] = myfun(ic,r,p0,q0)
% Integrator settings
geod_rel_tol = 1.e-6;
geod_abs_tol = 1.e-12;
options = odeset('RelTol',geod_rel_tol, 'AbsTol', geod_abs_tol);

% Unpack initial conditions
u = ic(1); v = ic(2); up = ic(3); vp = ic(4); lg = ic(5);

% initial wrapping point
p = [r*cos(u), r*sin(u), v];
% Straigth line segment from the origin point to the surface
e = p-p0;
e = e/norm(e);
% tangent vector at initial wrapping point
t = tang(u,v,up,vp,r);

% Integrate the geodesic's differential equations
[~,geodesic]=ode45(@int_gp,[0,lg],ic(1:4),options,r);

% final wrapping point
q = [r*cos(geodesic(end,1)), r*sin(geodesic(end,1)), geodesic(end,2)];
% Straigth line segment from the surface to the insertion point
e2 = q-q0;
e2 = e2/norm(e2);
% tangent vector at final wrapping point
t2 = tang(geodesic(end,1),geodesic(end,2),geodesic(end,3),geodesic(end,4),r);

a1 = e-t;
a2 = e2+t2;

% Tangency constraint errors
F = dot([a1';a2'],[a1';a2']);

term = [(e-t)';(e2+t2)'];
end

function [alpha,lsc] = LineSearch(f, x, p, F, g_x,r,p0,q0)
% Backtracking line search algorithm
a = 1-2/(1+sqrt(5));
b = 2/(1+sqrt(5));
alpha = 1.0;
lsc = 0;

while f(x + alpha*p,r,p0,q0) > F + a*alpha*g_x'*p
    alpha = b*alpha;
    lsc = lsc + 1;
end
end

function out = SolveWrap(ic,p0_l,q0_l,p_l,q_l,Radius)
% Solve the muscle wrapping problem

% Unpack initial conditions
u = ic(1); v = ic(2); up = ic(3); vp = ic(4); lg = ic(5);

f = @myfun;     % Tangency constraint evaluation function
g = @jacobi;    % Problem Jacobian evaluation function

% BFGS tolerance error
Eps = 1e-10;

% Integrator settings
geod_rel_tol = 1.e-6;
geod_abs_tol = 1.e-12;
options = odeset('RelTol',geod_rel_tol, 'AbsTol', geod_abs_tol);

p0 = p0_l;  % Origin point
q0 = q0_l;  % Insertion point

p = p_l;    % First point in contact with the surface
q = q_l;    % Second point in contact with the surface

% Straigth line segment from the origin point to the surface
e = p-p0;
e = e/norm(e);
% tangent vector at initial wrapping point
t = tang(u,v,up,vp,Radius);

% Straigth line segment from the surface to the insertion point
e2 = q-q0;
e2 = e2/norm(e2);
% tangent vector at final wrapping point
t2 = tang(u,v,up,vp,Radius);

a1 = e-t;
a2 = e2+t2;
% Evaluate the tangency constraint errors
F = dot(a1,a1)+dot(a2,a2);

nF(1) = norm([(e-t)';(e2+t2)']);

% Evaluate the gradient
g_x = g(f,F,ic,Radius,p0,q0);

% Initialize the Hessian approximation
I = eye(length(ic));
H = eye(length(ic));
H_i = H;

iter = 0;
lsc = 0;

% BFGS algorithm
while nF(iter+1) > Eps
    % Compute the descent direction
    p = -H_i*g_x;
    
    % Uncomment the following line if want to use the linesearch algorithm
    [alpha, lsc_t] = LineSearch(f, ic, p, F, g_x,Radius,p0,q0);
    
    % Uncomment the following line if don't want to use the linesearch algorithm
    % alpha = 1; lsc_t = 0;
    
    lsc = lsc + lsc_t;
    s = alpha*p;
    y = g(f,F,ic+s,Radius,p0,q0) - g_x;
    
    % Update initial conditions
    ic = ic + alpha*p;

    % Evaluate the tangency constraint errors
    [F,a1,a2,term] = f(ic,Radius,p0,q0);
    nF(iter+2) = norm(term);
    
    % Evaluate the problem Jacobian
    g_x = g(f,F,ic,Radius,p0,q0);
    
    % Compute the inverse of the approximated Hessian matrix
    H_i = (I - s*y'/(y'*s))*H_i*(I - y*s'/(y'*s)) + s*s'/(y'*s);

    iter = iter + 1;
end

% Integrate the geodesic's differential equations
[~,geodesic]=ode45(@int_gp,[0,ic(5)],ic(1:4),options,Radius);
% Initial wrapping point
p = [Radius*cos(ic(1)), Radius*sin(ic(1)), ic(2)];
% Final wrapping point
q = [Radius*cos(geodesic(end,1)), Radius*sin(geodesic(end,1)), geodesic(end,2)];

% Store the geodesic path
x = Radius.*cos(geodesic(:,1));
y = Radius.*sin(geodesic(:,1));
z = geodesic(:,2);

% The following wrapping condition for a cylindrical surface is taken from
%
% B. A. Garner and M. G. Pandy, “The obstacle-set method for representing 
% muscle paths in musculoskeletal models,” Computer Methods in Biomechanics 
% and Biomedical Engineering, vol. 3, no. 1, pp. 1–30,2000.
%
dt1 = p(1)*q(2)-p(2)*q(1); % if dt1>0 wrapping does occur

out = {ic,[x,y,z],p,q,iter,dt1,lsc,nF};

end
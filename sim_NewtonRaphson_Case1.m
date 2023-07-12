%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Case 1: Muscle is wrapped during whole simulation
%%%%%%% Newton-Raphson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

% Cylinder specifications
Radius = 1;
Height = 5;
SideCount = 20; % for visual purposes

ti = (0:0.012:3)';                                   % Time data

% Surface's Motion data
position = [0*ti, 0.9*sin(2*pi*ti), 0*ti];           % Position data
d_position = [0*ti, 2*pi*0.9*cos(2*pi*ti), 0*ti];    % Velocity data
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

% Stack the initial conditions to be updated by the Newton-Raphson method
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
    % Total muscle path length
    lmt{i} = norm(p0_i - out{1,i}{1,3}) + out{1,i}{1,1}(5) + norm(q0_i - out{1,i}{1,4});
    % Straigth line segment from the origin point to the surface
    e = out{1,i}{1,3} - p0_i; e = e/norm(e);
    % Straigth line segment from the surface to the insertion point
    e2 = q0_i - out{1,i}{1,4}; e2 = e2/norm(e2);
    % Muscle velocity
    d_c_p = skew(d_angles(i,:))*Rc_w{i}*out_lp{i}';
    d_c_q = skew(d_angles(i,:))*Rc_w{i}*out_lq{i}';
    dlmt{i} = dot(e,d_position(i,:)') + dot(e,d_c_p) - dot(e2,d_position(i,:)') - dot(e2,d_c_q) - dot(e,d_p0_i) + dot(e2,d_q0_i);
end

%% Store iterations per time step
for i=1:length(ti)
    iter(i) = out{1,i}{1,5};
end

%% Animation initial figure
% Initial vertices and faces of the surface
r_ini = position(1,:)';
R_ini = Rc_w{1};
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
vw = VideoWriter('animation_NR.avi');
vw.Quality = 100;
vw.FrameRate = 60;
open(vw);

% Animation Loop
for i = 1:length(ti)
    set(h_side, 'Vertices', vertices(:,:,i));
    set(h_bottom, 'Vertices', vertices(:,:,i));
    
    set(h_geod, 'XData', out{1,i}{1,2}(:,1),'YData', out{1,i}{1,2}(:,2),'ZData', out{1,i}{1,2}(:,3));
   
    set(h_movp, 'XData', out{1,i}{1,3}(1),'YData', out{1,i}{1,3}(2),'ZData', out{1,i}{1,3}(3));
    set(h_movq, 'XData', out{1,i}{1,4}(1),'YData', out{1,i}{1,4}(2),'ZData', out{1,i}{1,4}(3));
    
    set(h_linp, 'XData', [p0_i(1) out{1,i}{1,3}(1)],'YData', [p0_i(2) out{1,i}{1,3}(2)],'ZData', [p0_i(3) out{1,i}{1,3}(3)]);
    set(h_linq, 'XData', [q0_i(1) out{1,i}{1,4}(1)],'YData', [q0_i(2) out{1,i}{1,4}(2)],'ZData', [q0_i(3) out{1,i}{1,4}(3)]);
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
    % Plot the geodesic
    plot3(out{1,i}{1,2}(:,1), out{1,i}{1,2}(:,2), out{1,i}{1,2}(:,3),'Color',color_g,'LineWidth',2);
    % Plot the points on the surface
    plot3(out{1,i}{1,3}(1), out{1,i}{1,3}(2), out{1,i}{1,3}(3),'o','LineWidth',0.05,'MarkerEdgeColor',color_l,'MarkerFaceColor',color_p,'MarkerSize',4);
    plot3(out{1,i}{1,4}(1), out{1,i}{1,4}(2), out{1,i}{1,4}(3),'o','LineWidth',0.05,'MarkerEdgeColor',color_l,'MarkerFaceColor',color_p,'MarkerSize',4);
    % Plot straight line segments
    plot3([p0_i(1) out{1,i}{1,3}(1)],[p0_i(2) out{1,i}{1,3}(2)],[p0_i(3) out{1,i}{1,3}(3)],'Color',color_l,'LineWidth',2);
    plot3([q0_i(1) out{1,i}{1,4}(1)],[q0_i(2) out{1,i}{1,4}(2)],[q0_i(3) out{1,i}{1,4}(3)],'Color',color_l,'LineWidth',2);

end

% Text labes indicating the position of the surface in y-axis
text(-1.059273542725549,-1.528187082863439,2.644283394934888,'$^{(0)}\mbox{\boldmath $c$}_{y}=-0.89$cm','Interpreter','latex','FontSize',13,'Color',color_e{3});
text(-1.529189467317706,-0.659799192901801,2.125626189340494,'$^{(0)}\mbox{\boldmath $c$}_{y}=0.00$cm','Interpreter','latex','FontSize',13,'Color',color_e{1});
text(-2.067554394138853,0.312959274146635,1.555621103677652,'$^{(0)}\mbox{\boldmath $c$}_{y}=0.89$cm','Interpreter','latex','FontSize',13,'Color',color_e{2});

hold off;

% Axes settings
xlabel('x [cm]','Interpreter','latex'); ylabel('y [cm]','Interpreter','latex'); zlabel('z [cm]','Interpreter','latex');
axis equal;
view([140,35]);
grid on;
xlim([-3,3]);
ylim([-3,3]);
zlim([-3,3]);

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

function J = jac(x,r,h,F,p0,q0)
% Evaluation of the numerical Jacobian using central differences
Fx = F(x,r,p0,q0);
F1h = F(x+[h;0;0;0;0],r,p0,q0);
F2h = F(x+[0;h;0;0;0],r,p0,q0);
F3h = F(x+[0;0;h;0;0],r,p0,q0);
F4h = F(x+[0;0;0;h;0],r,p0,q0);
F5h = F(x+[0;0;0;0;h],r,p0,q0);

J = [F1h-Fx F2h-Fx F3h-Fx F4h-Fx F5h-Fx]./h;
end

function F = myfun(ic,r,p0,q0)
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

% Tangency constraint errors
F = [(e-t)';(e2+t2)'];
end

function out = SolveWrap(ic,p0_l,q0_l,p_l,q_l,Radius)
% Solve the muscle wrapping problem

% Unpack initial conditions
u = ic(1); v = ic(2); up = ic(3); vp = ic(4); lg = ic(5);

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

% Evaluate the tangency constraint errors
F = [(e-t)';(e2+t2)'];
nF(1) = norm(F);

h = 1e-10;
iter = 0;

% Newton-Raphson algorithm
while nF(iter+1)>1e-10 && iter < 250
    % Evaluate the problem's Jacobian using finite central differences
    Jn = jac(ic,Radius,h,@myfun,p0,q0);
    p_Jn = inv(Jn'*Jn + 0.08*eye(5))*Jn';
    % Compute the descent direction
    d = -(p_Jn*F);
   
    alpha = 1;
    
    % Update initial conditions
    ic = ic + alpha*(d);
    
    % Unpack initial conditions
    u = ic(1); v = ic(2); up = ic(3); vp = ic(4); lg = ic(5);
    
    % Initial wrapping point
    p = [Radius*cos(u), Radius*sin(u), v];
    % Straigth line segment from the origin point to the surface
    e = p-p0;
    e = e/norm(e);
    % tangent vector at initial wrapping point
    t = tang(u,v,up,vp,Radius);
    
    % Integrate the geodesic's differential equations
    [~,geodesic]=ode45(@int_gp,[0,lg],ic(1:4),options,Radius);
    uq = geodesic(end,1); vq = geodesic(end,2); upq = geodesic(end,3); vpq = geodesic(end,4);
    
    % Final wrapping point
    q = [Radius*cos(uq), Radius*sin(uq), vq];
    % Straigth line segment from the surface to the insertion point
    e2 = q-q0;
    e2 = e2/norm(e2);
    % tangent vector at final wrapping point
    t2 = tang(uq,vq,upq,vpq,Radius);
    
    % Evaluate the tangency constraint errors
    F = [(e-t)';(e2+t2)'];
    
    nF(iter+2) = norm(F);
    iter = iter+1;  
end

% Store the geodesic path
x = Radius.*cos(geodesic(:,1));
y = Radius.*sin(geodesic(:,1));
z = geodesic(:,2);

out = {ic,[x,y,z],p,q,iter,nF};

end
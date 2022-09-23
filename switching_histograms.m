function [C00, Psw00]= switching_histograms(Tqp_ET,Tb_ET,omega_tau1,omega_tau2,omega_taup,gap_Tqb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tqp_ET   kB*T_qp/E_T    this is the temperature for the poisoning
%% Tb_ET   KB*T_b/E_T    this is the temperature for the pair relaxation
%% omega_tau1  omega*tau1  eg <--> gg
%% omega_tau2  omega*tau2  ee <--> eg
%% omega_taup  omega*tau'  % pair relaxation  ee <--> gg
%% gap_Tqb    gap/Tqb  This is the ratio between gap and kB*Tqb.
%% The gap between Ee and Eg is chosen in unit of kB*T1, see the following


edge_amp = 0.06;  % just some arbitrary small ic for the hinge junction
ic = edge_amp;
ET = ic*2; % Thouless energy
Tqp = Tqp_ET*ET; 
Tb = Tb_ET*ET;
tau1 = omega_tau1;
tau2 = omega_tau2;
taup = omega_taup;
gap = gap_Tqb * Tqp; 
phi = linspace(-3.5,6,8000)*pi;

eta = 0.05; %broadening


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the strong junction in the SQUID, the parameters are somewhat arbitrary
Jstrong = 300; % this is just a large critical current for the large junction in forming a Squid
delta1 = 0.02;% 0.02;
delta2 = 0.02;% 0.02;
alpha = 0.5;
alpha = 0.5;
gamma0 = gamma0_analytical(phi, 0,Jstrong, alpha, -alpha, delta1,delta2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% currents (in the absense of gap, "no-gap" or "ng") for the weak hinge mode junction
[i1p,i1m] = currents_w_gap(gamma0+(alpha)*phi,ic,gap);
i2p = i1p;
i2m = i1m;
% i1p_ng = I_p(gamma0+(alpha)*phi);
% i1m_ng = I_m(gamma0+(alpha)*phi);
% i2p_ng = I_p(gamma0+(alpha)*phi);
% i2m_ng = I_m(gamma0+(alpha)*phi);
lambda = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ilist = linspace(-1,1,800)';
Npts = 1000;
[pgg, pge] = prob_fun(Npts,tau1,tau2,taup,Tqp,Tb,ic,gap);
pgg = c_grd_noneq(pgg,gamma0+(alpha)*phi);
pge = c_grd_noneq(pge,gamma0+(alpha)*phi);
Psw00 = P_sw(Ilist,edge_amp*i1p+edge_amp*i2p,eta);
Psw10 = P_sw(Ilist,edge_amp*i1m+edge_amp*i2p,eta);
Psw01 = P_sw(Ilist,edge_amp*i1p+edge_amp*i2m,eta);
Psw11 = P_sw(Ilist,edge_amp*i1m+edge_amp*i2m,eta);
C00 = repmat(pgg,[length(Ilist),1]);
C10 = repmat(pge,[length(Ilist),1]);
C01 = repmat(pge,[length(Ilist),1]);
C11 = repmat(1 - pgg - 2*pge,[length(Ilist),1]);
philist = repmat(phi, [size(Psw00,1),1]);
Ilist = repmat(Ilist, [1,size(Psw00,2)]);
Psw = C00.*Psw00 + C10.*Psw10 + C01.*Psw01 + C11.*Psw11; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot PSW rather than dPsw
% Pick 5 phases
pick_phi_list = [-0.2 0.2 0.5 0.9 1.2]+0.1; 
I_external = Ilist(:,1);
fig2 = figure;
hold on;
box on;
%set(gca, 'YDir','reverse');
set(gca,'FontSize',26);
set(gcf, 'Color', 'white');
ax = gca;
ax.XColor='black';
ax.YColor='black';
ax.LineWidth= 1.5;
[M,idxphi]=min(abs(phi/pi-pick_phi_list(1)));
pick_phi_list(1)=phi(idxphi);
idxphi1 = idxphi;
plot(Psw(:,idxphi)/Psw(end,idxphi),I_external,'--','LineWidth',3,'Color',[204, 51, 0]/255);
[M,idxphi]=min(abs(phi/pi-pick_phi_list(2)));
pick_phi_list(2)=phi(idxphi);
idxphi2 = idxphi;
plot(Psw(:,idxphi)/Psw(end,idxphi),I_external,'-','LineWidth',3,'Color','blue');
[M,idxphi]=min(abs(phi/pi-pick_phi_list(3)));
pick_phi_list(3)=phi(idxphi);
idxphi3 = idxphi;
plot(Psw(:,idxphi)/Psw(end,idxphi),I_external,'--','LineWidth',3,'Color',[0 153 0]/255);
[M,idxphi]=min(abs(phi/pi-pick_phi_list(4)));
pick_phi_list(4)=phi(idxphi);
idxphi4 = idxphi;
plot(Psw(:,idxphi)/Psw(end,idxphi),I_external,'-','LineWidth',3,'Color',[204, 153, 0]/255);
[M,idxphi]=min(abs(phi/pi-pick_phi_list(5)));
pick_phi_list(5)=phi(idxphi);
idxphi5 = idxphi;
plot(Psw(:,idxphi)/Psw(end,idxphi),I_external,'-.','LineWidth',3,'Color','red');
xlabel("$P$","Interpreter","latex");
ylabel("$I$","Interpreter","latex");
%figname = "linecuts_xi_"+ num2str(xi) + "_tau_" + num2str(tau) + "_taup_" + num2str(taup);
%export_fig(fig2,figname,'-r300','-painters','-png');
%close(fig2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT dPsw
dPsw = diff(Psw);
Ilist = Ilist(1:end-1,:);
philist = philist(1:end-1,:);
fig = figure;
c=gray(4096);
colormap(flipud(c));
pcolor(philist/pi,Ilist,dPsw);
shading interp;
xlabel("$\phi/\pi$",'Interpreter','latex');
ylabel("$I$", 'Interpreter','latex');
% yticks([]);
% colorbar;
box on;
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth= 1.5;
xtickangle(ax,0);
set(gca,'FontSize',26);
legend boxoff;
legend off;
%xlim([-3.1,5.1]);
%ylim([-0.8,0.8]);
%pbaspect([1 2 2]);
title({["$k_B T_{qp}/E_T = $" + num2str(Tqp_ET,'%3.2f')], ...
    ["$k_B T_{b}/E_T = $" + num2str(Tb_ET,'%3.2f')], ...
    ["$\omega\tau_1=$" + num2str(tau1,'%4.2f') + ', $\omega\tau_2 = $' + num2str(tau2,'%4.2f')], ...
    ["$\omega\tau =$" + num2str(taup,'%4.3f')], ["$\delta_{gap}/k_B T_{qp} =$" + num2str(gap_Tqb,'%4.1f')]},...
    'Interpreter','latex');
%figname = "xi_"+ num2str(xi) + "_tau_" + num2str(tau) + "_taup_" + num2str(taup) + ".jpg" ;
ax = gca;
%ax.Color='black';
set(gca, 'Color', 'none'); % Sets axes background
set(gcf, 'Color', 'white');
set(0,'DefaultFigureColor','remove');
xline(pick_phi_list(1)/pi,'--','LineWidth',3,'Color',[204, 51, 0]/255);
xline(pick_phi_list(2)/pi,'-','LineWidth',3,'Color','blue');
xline(pick_phi_list(3)/pi,'--','LineWidth',3,'Color',[0 153 0]/255);
xline(pick_phi_list(4)/pi,'-','LineWidth',3,'Color',[204, 153, 0]/255);
xline(pick_phi_list(5)/pi,'-.','LineWidth',3,'Color','red');

grid on;
set(gca,'layer','top')
%export_fig(fig,figname,'-r300','-painters', '-jpg');
%close(fig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig3 = figure;
hold on;
box on;
%set(gca, 'YDir','reverse');
set(gca,'FontSize',26);
set(gcf, 'Color', 'white');
ax = gca;
ax.XColor='black';
ax.YColor='black';
ax.LineWidth= 1.5;
I_external = I_external(1:end-1);
shift = 0.04;
plot(dPsw(:,idxphi1)+0*shift,I_external,'--','LineWidth',3,'Color',[204, 51, 0]/255);
plot(dPsw(:,idxphi2)+1*shift,I_external,'-','LineWidth',3,'Color','blue');
plot(dPsw(:,idxphi3)+2*shift,I_external,'--','LineWidth',3,'Color',[0 153 0]/255);
plot(dPsw(:,idxphi4)+3*shift,I_external,'-','LineWidth',3,'Color',[204, 153, 0]/255);
plot(dPsw(:,idxphi5)+4*shift,I_external,'-.','LineWidth',3,'Color','red');
ylim([-0.7,0.7]);
xlim([-0.001,0.2]);
xlabel("number of events, shifted");
ylabel("I","Interpreter","latex");
%figname = "single_dPdI_linecuts_xi_"+ num2str(xi) + "_tau_" + num2str(tau);
grid on;
set(gca,'layer','top')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = joint_probabilities(phi0,tau1,tau2,taup,Tqb,Tb,ic,gap)
phase_max = pi/2;
ET = 2*ic; % Thouless energy
phispan = [phi0-phase_max phi0];
p0 = zeros(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%equilibrium_distr
delta0 = fdelta(phi0 - phase_max,ic,gap)/ET; 
Tqb = Tqb/ET;
Tb = Tb/ET;
[Gamma_eggg,Gamma_ggeg,Gamma_eeeg,Gamma_egee,Gamma_eegg,Gamma_ggee] ...
    = gamma_rates(tau1,tau2,taup,Tqb,Tb,delta0);
p0 = equilibrium_prob(Gamma_eggg,Gamma_ggeg,Gamma_eeeg,Gamma_egee,Gamma_eegg,Gamma_ggee);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integrate ODE
 [phi,p] = ode45(@(phi,p) gammamat(phi,p,tau1,tau2,taup,Tqb,Tb,ic,gap), phispan,p0);
 y = p(end,:);
end

function [Gamma_eggg,Gamma_ggeg,Gamma_eeeg,Gamma_egee,Gamma_eegg,Gamma_ggee] ...
    = gamma_rates(tau1,tau2,taup,Tqb,Tb,delta0)
fermi_deltaT1 = fermi(delta0/Tqb);
bose_deltaT2 = bose(2*delta0/Tb);
Gamma_eggg = fermi_deltaT1/tau1;
Gamma_ggeg = (1-fermi_deltaT1)/tau1;
Gamma_eeeg = fermi_deltaT1/tau2;
Gamma_egee = (1-fermi_deltaT1)/tau2;
Gamma_eegg = 2*delta0/taup*bose_deltaT2;
Gamma_ggee = 2*delta0/taup*(1+bose_deltaT2);
end

function p0 = equilibrium_prob(Gammaeggg,Gammaggeg,Gammaeeeg,Gammaegee,Gammaeegg,Gammaggee)
p0(1) = (Gammaeeeg.*Gammaggee+(Gammaegee+Gammaggee).*Gammaggeg).*(2.* ...
  Gammaeegg.*Gammaegee+2.*Gammaeggg.*(2.*Gammaegee+Gammaggee)+ ...
  Gammaeeeg.*(Gammaeegg+2.*Gammaeggg+Gammaggee)+(Gammaeegg+ ...
  Gammaegee+Gammaeggg+Gammaggee).*Gammaggeg).^(-1);
p0(2) = (Gammaeegg.*Gammaegee+Gammaeggg.*(2.*Gammaegee+Gammaggee)).*(2.* ...
  Gammaeegg.*Gammaegee+2.*Gammaeggg.*(2.*Gammaegee+Gammaggee)+ ...
  Gammaeeeg.*(Gammaeegg+2.*Gammaeggg+Gammaggee)+(Gammaeegg+ ...
  Gammaegee+Gammaeggg+Gammaggee).*Gammaggeg).^(-1);
end

function dydphi = gammamat(phi,y,tau1,tau2,taup,Tqb,Tb,ic,gap)
    dydphi =zeros(2,1);
    delta = fdelta(phi,ic,gap)/(2*ic);
    [Gamma_eggg,Gamma_ggeg,Gamma_eeeg,Gamma_egee,Gamma_eegg,Gamma_ggee] ...
    = gamma_rates(tau1,tau2,taup,Tqb,Tb,delta);
    dydphi(1) = Gamma_ggee + (-2*Gamma_ggee + Gamma_ggeg)* y(2) + (-Gamma_eegg - ...
    2*Gamma_eggg - Gamma_ggee)* y(1);
    dydphi(2) = Gamma_egee + (-Gamma_eeeg - 2* Gamma_egee - Gamma_ggeg)* y(2) + (-Gamma_egee +...
     Gamma_eggg) *y(1);
end

 function y = fdelta(phi,ic,gap)
    phi = mod(phi+pi,2*pi) - pi;
    y = zeros(size(phi));
    [Eg,Ee] = EgEe(phi,ic,gap);
    y = Ee - Eg;
%     id1 = (phi<=0);
%     id2 = (phi>0);
%     y(id2) = 2*ic*(pi - phi(id2));
%     y(id1) = 2*ic*(pi + phi(id1));
  end

  function y = fermi(x)
      y = 1./(1+exp(x));
  end
  function y = bose(x)
    y = 1./(exp(x)-1);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Eg_new, Ee_new] = EgEe(phi,ic,gap)
Eg = feg(phi,ic);
Ee = fee(phi,ic);
Eg_new = ((Ee+Eg)-sqrt((Ee-Eg).^2 + 4*gap^2))/2;
Ee_new = ((Ee+Eg)+sqrt((Ee-Eg).^2 + 4*gap^2))/2;
end

function y = feg(phi,ic)
phi = mod(phi+pi,2*pi)-pi;
y = phi.^2/(2*pi) - pi/2;
y = y*ic;
end

function y =  fee1(phi)
y =  phi.^2/(2*pi) - pi/2 + 2*(phi+pi);
end

function y = fee2(phi)
y =  phi.^2/(2*pi) - pi/2 + 2*(pi-phi);
end

function y = fee(phi,ic)
phi = mod(phi+pi,2*pi) - pi;
y = zeros(size(phi));
firstid = (phi>=-pi & phi<0);
secondid = (phi>=0 & phi<pi);
y(firstid) = fee1(phi(firstid));
y(secondid) = fee2(phi(secondid));
y = y*ic;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = I_p(phi)  % in unit of eE_T/(2\pi \hbar)
% phi = mod(phi+pi,2*pi) - pi;
% y = phi;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = I_m(phi)  % in unit of eE_T/(2\pi \hbar)
% phi = mod(phi+pi,2*pi) - pi;
% y = phi - 2*pi*sign(phi);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ip,Im]=currents_w_gap(phi,ic,gap) % in unit of eE_T/(2\pi \hbar)
phi = mod(phi+pi,2*pi) - pi;
ip_ng = phi;
im_ng = phi - 2*pi*sign(phi);
[Eg,Ee] = EgEe(phi,ic,gap);
I_tot = ip_ng + im_ng;
I_diff = im_ng - ip_ng;
E_diff = Ee - Eg;
correction_factor = (E_diff.*I_diff)./sqrt(E_diff.^2 + 4*gap^2);
Ip = 0.5*(I_tot - correction_factor);
Im = 0.5*(I_tot + correction_factor);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function y = Iextended(gamma0,phi,Jc,alpha,delta)
  y = Jc./phi.*sin(gamma0+alpha*phi).*sin(delta*phi);
end

function y = Idouble_extended_c(phi,J1,J2,alpha1,alpha2,delta1,delta2)
A = Aphi(phi, J1,J2, alpha1, alpha2, delta1,delta2);
B = Bphi(phi, J1,J2, alpha1, alpha2, delta1,delta2);
y = sqrt(A.^2 + B.^2)./abs(phi);
end

function y = double_extended(gamma0,phi,J1,J2,alpha1,alpha2,delta1,delta2)
  y = Iextended(gamma0,phi,J1,alpha1,delta1) + Iextended(gamma0,phi,J2,alpha2,delta2);
end

function gamma0 = findgamma0(phi,J1,J2,alpha1,alpha2,delta1,delta2,Npts_gamma)
  gamma = linspace(-pi,pi,Npts_gamma)';
  currents = bsxfun(@(gamma,phi)double_extended(gamma,phi,J1,J2,alpha1,alpha2,delta1,delta2),gamma,phi);
  [maxgamma,maxgamma_id]=max(currents,[],1);
  gamma0 = gamma(maxgamma_id);
end

function y = Aphi(phi, J1,J2, alpha1, alpha2, delta1,delta2)
  y = J1*cos(alpha1*phi).*sin(delta1*phi)+J2*cos(alpha2*phi).*sin(delta2*phi);
end
function y = Bphi(phi, J1,J2, alpha1, alpha2, delta1,delta2)
  y = J1*sin(alpha1*phi).*sin(delta1*phi)+J2*sin(alpha2*phi).*sin(delta2*phi);
end
function y = gamma0_analytical(phi, J1,J2, alpha1, alpha2, delta1,delta2)
A = Aphi(phi, J1,J2, alpha1, alpha2, delta1,delta2);
B = Bphi(phi, J1,J2, alpha1, alpha2, delta1,delta2);
y = sign(A.*phi)*pi/2 - atan(B./A);
end

function y = I_p_multiple(gamma,phi,alpha_i,alpha_f,N_channels,Jlists)
  alphalist = linspace(alpha_i,alpha_f,N_channels);
  y = 0;
  for jj = 1:N_channels
    y = y + I_p(gamma+phi*alphalist(jj),Jlists(jj));
  end
end

function y = theta_smooth(x,eta)
y = (tanh(x/eta)+1)/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = P_sw(I, Ic, eta)
  y  = bsxfun(@(I, Ic) theta_smooth(I - Ic, eta), I, Ic);
end

function y = c_odd_thermal(diff,T)
y = exp(-diff/T);
y = y./(1+y);
end

function [pgg, pge] = prob_fun(Npts,tau1,tau2,taup,Tqp,Tb,ic,gap)
philist = linspace(0,2*pi,Npts);
y = zeros(2,length(philist));
parfor jj = 1:length(philist)
  y(:,jj) = joint_probabilities(philist(jj),tau1,tau2,taup,Tqp,Tb,ic,gap);
end
pgg = interp1(philist,y(1,:),'linear','pp');
pge = interp1(philist,y(2,:),'linear','pp');
end



function y = c_grd_noneq(pp,phi)
phi = mod(phi,2*pi);
y = ppval(pp,phi);
end


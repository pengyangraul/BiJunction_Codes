function plot_joint_probabilities(Tqp_ET,Tb_ET,omega_tau1,omega_tau2,omega_taup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tqp_ET   kB*T_qp/E_T    this is the temperature for the poisoning
%% Tb_ET   KB*T_b/E_T    this is the temperature for the pair relaxation
%% omega_tau  omega*tau
%% omega_taup  omega*tau'
%% The gap between Ee and Eg is chosen in unit of kB*T1, see the following

ic = 1;
ET = ic*2; % Thouless energy
T1 = Tqp_ET*ET; 
T2 = Tb_ET*ET;
tau1 = omega_tau1;
tau2 = omega_tau2;
taup1 = omega_taup;

gap1 = 1.5*T1;
gap2 = 1.2*T1;
gap3 = T1;
gap4 = 0.8*T1;
gap5 = 0*T1;

philist = linspace(0,2,500);
ps1 = zeros(2,length(philist));
ps2 = ps1;
ps3 = ps1;
ps4 = ps1;
ps5 = ps1;

parfor jj = 1:length(philist)
  ps1(:,jj) = joint_probabilities(philist(jj)*pi,tau1,tau2,taup1,T1,T2,ic,gap1);
  ps2(:,jj) = joint_probabilities(philist(jj)*pi,tau1,tau2,taup1,T1,T2,ic,gap2);
  ps3(:,jj) = joint_probabilities(philist(jj)*pi,tau1,tau2,taup1,T1,T2,ic,gap3);
  ps4(:,jj) = joint_probabilities(philist(jj)*pi,tau1,tau2,taup1,T1,T2,ic,gap4);
  ps5(:,jj) = joint_probabilities(philist(jj)*pi,tau1,tau2,taup1,T1,T2,ic,gap5);
end
pee1 = 1 - ps1(1,:) - 2*ps1(2,:);
pee2 = 1 - ps2(1,:) - 2*ps2(2,:);
pee3 = 1 - ps3(1,:) - 2*ps3(2,:);
pee4 = 1 - ps4(1,:) - 2*ps4(2,:);
pee5 = 1 - ps5(1,:) - 2*ps5(2,:);

C = linspecer(5);

fig = figure;
plot(philist,pee1,'-','LineWidth',1.5, 'Color',C(1,:));
hold on;
plot(philist,pee2,'-','LineWidth',1.5, 'Color',C(2,:));
plot(philist,pee3,'-','LineWidth',1.5, 'Color',C(3,:));
plot(philist,pee4,'-','LineWidth',1.5, 'Color',C(4,:));
plot(philist,pee5,'-','LineWidth',1.5, 'Color',C(5,:));
plot(philist,2*ps1(2,:),'--','LineWidth',1.5, 'Color',C(1,:));
plot(philist,2*ps2(2,:),'--','LineWidth',1.5, 'Color',C(2,:));
plot(philist,2*ps3(2,:),'--','LineWidth',1.5, 'Color',C(3,:));
plot(philist,2*ps4(2,:),'--','LineWidth',1.5, 'Color',C(4,:));
plot(philist,2*ps5(2,:),'--','LineWidth',1.5, 'Color',C(5,:));
plot(philist,ps1(1,:),'-.','LineWidth',1.5, 'Color',C(1,:));
plot(philist,ps2(1,:),'-.','LineWidth',1.5, 'Color',C(2,:));
plot(philist,ps3(1,:),'-.','LineWidth',1.5, 'Color',C(3,:));
plot(philist,ps4(1,:),'-.','LineWidth',1.5, 'Color',C(4,:));
plot(philist,ps5(1,:),'-.','LineWidth',1.5, 'Color',C(5,:));

title({["$k_B T_{qp}/E_T = $" + num2str(Tqp_ET,'%3.2f')], ...
    ["$k_B T_{b}/E_T = $" + num2str(Tb_ET,'%3.2f')], ...
    ["$\tau_1=$" + num2str(tau1,'%4.2f') + ', $\tau_2 = $' + num2str(tau2,'%4.2f')],  ["$\tau' =$" + num2str(taup1,'%4.3f')]},...
  'Interpreter','latex');
legend("$\delta_{gap}/E_T=$" + num2str(gap1/ET,'%3.2f'), ...
"$\delta_{gap}/E_T=$" + num2str(gap2/ET,'%3.2f'), ...
"$\delta_{gap}/E_T=$" +  num2str(gap3/ET,'%3.2f'), ...
"$\delta_{gap}/E_T=$" +  num2str(gap4/ET,'%3.2f'), ...
"$\delta_{gap}/E_T=$" +  num2str(gap5/ET,'%3.2f'),'Interpreter','latex')
legend('Location','best');
legend('boxoff')
xlabel("$(\Phi+\gamma_{max})/\pi$",'Interpreter','latex');
ylabel("$p$",'Interpreter','latex');
ax = gca;
ax.BoxStyle = 'full';
ax.LineWidth= 1.5;
set(gca,'FontSize',20);
set(gca,'Color','white');
set(gcf,'Color','white');
ylim([0,1]);
xlim([0,2]);
end

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


% function lineStyles = linspecer(N)
% This function creates an Nx3 array of N [R B G] colors
% These can be used to plot lots of lines with distinguishable and nice
% looking colors.
% 
% lineStyles = linspecer(N);  makes N colors for you to use: lineStyles(ii,:)
% 
% colormap(linspecer); set your colormap to have easily distinguishable 
%                      colors and a pleasing aesthetic
% 
% lineStyles = linspecer(N,'qualitative'); forces the colors to all be distinguishable (up to 12)
% lineStyles = linspecer(N,'sequential'); forces the colors to vary along a spectrum 

function lineStyles=linspecer(N,varargin)

if nargin==0 % return a colormap
    lineStyles = linspecer(128);
    return;
end

if ischar(N)
    lineStyles = linspecer(128,N);
    return;
end

if N<=0 % its empty, nothing else to do here
    lineStyles=[];
    return;
end

% interperet varagin
qualFlag = 0;
colorblindFlag = 0;

if ~isempty(varargin)>0 % you set a parameter?
    switch lower(varargin{1})
        case {'qualitative','qua'}
            if N>12 % go home, you just can't get this.
                warning('qualitiative is not possible for greater than 12 items, please reconsider');
            else
                if N>9
                    warning(['Default may be nicer for ' num2str(N) ' for clearer colors use: whitebg(''black''); ']);
                end
            end
            qualFlag = 1;
        case {'sequential','seq'}
            lineStyles = colorm(N);
            return;
        case {'white','whitefade'}
            lineStyles = whiteFade(N);return;
        case 'red'
            lineStyles = whiteFade(N,'red');return;
        case 'blue'
            lineStyles = whiteFade(N,'blue');return;
        case 'green'
            lineStyles = whiteFade(N,'green');return;
        case {'gray','grey'}
            lineStyles = whiteFade(N,'gray');return;
        case {'colorblind'}
            colorblindFlag = 1;
        otherwise
            warning(['parameter ''' varargin{1} ''' not recognized']);
    end
end      
% *.95
% predefine some colormaps
  set3 = colorBrew2mat({[141, 211, 199];[ 255, 237, 111];[ 190, 186, 218];[ 251, 128, 114];[ 128, 177, 211];[ 253, 180, 98];[ 179, 222, 105];[ 188, 128, 189];[ 217, 217, 217];[ 204, 235, 197];[ 252, 205, 229];[ 255, 255, 179]}');
set1JL = brighten(colorBrew2mat({[228, 26, 28];[ 55, 126, 184]; [ 77, 175, 74];[ 255, 127, 0];[ 255, 237, 111]*.85;[ 166, 86, 40];[ 247, 129, 191];[ 153, 153, 153];[ 152, 78, 163]}'));
set1 = brighten(colorBrew2mat({[ 55, 126, 184]*.85;[228, 26, 28];[ 77, 175, 74];[ 255, 127, 0];[ 152, 78, 163]}),.8);

% colorblindSet = {[215,25,28];[253,174,97];[171,217,233];[44,123,182]};
colorblindSet = {[215,25,28];[253,174,97];[171,217,233]*.8;[44,123,182]*.8};

set3 = dim(set3,.93);

if colorblindFlag
    switch N
        %     sorry about this line folks. kind of legacy here because I used to
        %     use individual 1x3 cells instead of nx3 arrays
        case 4
            lineStyles = colorBrew2mat(colorblindSet);
        otherwise
            colorblindFlag = false;
            warning('sorry unsupported colorblind set for this number, using regular types');
    end
end
if ~colorblindFlag
    switch N
        case 1
            lineStyles = { [  55, 126, 184]/255};
        case {2, 3, 4, 5 }
            lineStyles = set1(1:N);
        case {6 , 7, 8, 9}
            lineStyles = set1JL(1:N)';
        case {10, 11, 12}
            if qualFlag % force qualitative graphs
                lineStyles = set3(1:N)';
            else % 10 is a good number to start with the sequential ones.
                lineStyles = cmap2linspecer(colorm(N));
            end
        otherwise % any old case where I need a quick job done.
            lineStyles = cmap2linspecer(colorm(N));
    end
end
lineStyles = cell2mat(lineStyles);

end

% extra functions
function varIn = colorBrew2mat(varIn)
for ii=1:length(varIn) % just divide by 255
    varIn{ii}=varIn{ii}/255;
end        
end

function varIn = brighten(varIn,varargin) % increase the brightness

if isempty(varargin),
    frac = .9; 
else
    frac = varargin{1}; 
end

for ii=1:length(varIn)
    varIn{ii}=varIn{ii}*frac+(1-frac);
end        
end

function varIn = dim(varIn,f)
    for ii=1:length(varIn)
        varIn{ii} = f*varIn{ii};
    end
end

function vOut = cmap2linspecer(vIn) % changes the format from a double array to a cell array with the right format
vOut = cell(size(vIn,1),1);
for ii=1:size(vIn,1)
    vOut{ii} = vIn(ii,:);
end
end
%%
% colorm returns a colormap which is really good for creating informative
% heatmap style figures.
% No particular color stands out and it doesn't do too badly for colorblind people either.
% It works by interpolating the data from the
% 'spectral' setting on http://colorbrewer2.org/ set to 11 colors
% It is modified a little to make the brightest yellow a little less bright.
function cmap = colorm(varargin)
n = 100;
if ~isempty(varargin)
    n = varargin{1};
end

if n==1
    cmap =  [0.2005    0.5593    0.7380];
    return;
end
if n==2
     cmap =  [0.2005    0.5593    0.7380;
              0.9684    0.4799    0.2723];
          return;
end

frac=.95; % Slight modification from colorbrewer here to make the yellows in the center just a bit darker
cmapp = [158, 1, 66; 213, 62, 79; 244, 109, 67; 253, 174, 97; 254, 224, 139; 255*frac, 255*frac, 191*frac; 230, 245, 152; 171, 221, 164; 102, 194, 165; 50, 136, 189; 94, 79, 162];
x = linspace(1,n,size(cmapp,1));
xi = 1:n;
cmap = zeros(n,3);
for ii=1:3
    cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
end
cmap = flipud(cmap/255);
end

function cmap = whiteFade(varargin)
n = 100;
if nargin>0
    n = varargin{1};
end

thisColor = 'blue';

if nargin>1
    thisColor = varargin{2};
end
switch thisColor
    case {'gray','grey'}
        cmapp = [255,255,255;240,240,240;217,217,217;189,189,189;150,150,150;115,115,115;82,82,82;37,37,37;0,0,0];
    case 'green'
        cmapp = [247,252,245;229,245,224;199,233,192;161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27];
    case 'blue'
        cmapp = [247,251,255;222,235,247;198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107];
    case 'red'
        cmapp = [255,245,240;254,224,210;252,187,161;252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13];
    otherwise
        warning(['sorry your color argument ' thisColor ' was not recognized']);
end

cmap = interpomap(n,cmapp);
end

% Eat a approximate colormap, then interpolate the rest of it up.
function cmap = interpomap(n,cmapp)
    x = linspace(1,n,size(cmapp,1));
    xi = 1:n;
    cmap = zeros(n,3);
    for ii=1:3
        cmap(:,ii) = pchip(x,cmapp(:,ii),xi);
    end
    cmap = (cmap/255); % flipud??
end




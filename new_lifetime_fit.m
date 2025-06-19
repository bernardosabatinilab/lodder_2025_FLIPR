%% This script assumes that there are two arrays in memory corresponding to
% intensity and lifetime measurements produced by FLIPR.
%   F_all : the intensity data as a function of time
%   T_all : the lifetime data as a function of time.
%
% The goal is to use an optimization procedure to fit F_all and T_all
% simultaneously. Athought not necessary, we typically invert the intensity
% measurements to make them positive and subtract off the electronic offset
% measured with the PMT blocked.

%% The fit variables and some initial conditions we used for dLight3.8

R=43; % ratio of flourscence of state2 to state1
f_offset=0; % offset in flourescence readout (no light gives this value). Typically 0
t1=1.55; % lifetime of state 1
t2=2; % lifetime of state 2
tb=1.55; % lifetime of background flourescence
fb=0.0; % the amoung of background flourescence in units of volts (same as F_all)
scale_factor=.1; % how to go from relative flourescence to the readout values in volts

params_0=[R f_offset t1 t2 tb fb scale_factor]; % a vector used to store the parameters


%% Do a test run of the model

f_state1=0:0.01:1; % the fraction of sensor in the unbound state

% Fl_model uses the equations in Lodder et al methods section plus
% additional parameters for background flourescence to calculate a measured
% lifetime and intensity for values of the fraction of sensor in state 1.
% For dLight3.8, state 1 corresponds to the unbound state and state 2 the
% bound state.  For easy interpretation, we plot version state2
% (1-f_unbound).
[F_calc, T_calc]=Fl_model(params_0, f_state1);

figure; plot(1-f_state1, F_calc); title('intensity vs frac bound')
figure; plot(1-f_state1, T_calc); title('lifetime vs frac bound')
figure; plot(1-f_state1, F_calc/scale_factor); title('dF/F vs frac bound')

%% set data boundaries for the fit
range_min=122800;
range_max=range_min+4000;

disp([range_min range_max]);

F_data=Fd_all(range_min:range_max); % the portion to fit
T_data=Td_all(range_min:range_max);

% rescale the intensity channel to have 1/2 the mean value of the lifetime.
% this ensures that both are used in the cost function but that the
% lifetime is more heavily weighted
scale_ratio=mean(T_data)/mean(F_data)/2;
F_data_scaled=F_data*scale_ratio; % rescale

% flag to control if we use the above guess for the params or continue from the last run
first_loop=true; 

%% set up initial guesses
f_guess=0.5*ones(1, length(F_data)); 

%params_guess=[43 0 1.55 1.92 1.55 0 0.4];
params_guess=[37.7944   0    1.5500    1.8620    1.5500    0    0.8035];

% this allows one to fix some parameters.  
% nan means fit it
% a value means hold it at that value
params_fix=[43 0 1.55 nan 1.55 0 nan];


%%
if first_loop
    p_est=[params_guess f_guess];
end


%%
p_init=p_est; % set up the initial parameters

% set upper and lower bounds
%    [R   f0  t1 t2 tb fb  sf   f1, f2, ...]
lb = [0     -Inf    0  0  0  0 0    zeros(1, length(F_data))];
ub = [50    Inf     3  3  3  1 Inf  ones(1, length(F_data))];

% set the fixed parameters
p_init(~isnan(params_fix))=params_fix(~isnan(params_fix));

disp('STARTING:')
disp(p_init(1:7))

% do the fit 
opts = optimoptions('lsqnonlin','Display','iter', 'TolFun',1e-8,'TolX',1e-8, 'MaxIterations', 300);
[p_est, resnorm] = lsqnonlin(@(p) fit_errors(p, F_data_scaled, T_data, params_fix) , p_init, lb, ub, opts);

%% analyze the returns
params_fit_scaled=p_est(1:7); % raw fit values 

% get the estimates of fraction in state 1
f_state1_fit=p_est(8:end);


% undo the scaling of F_data
p_est([2 6 7])=p_est([2 6 7])/scale_ratio;
params_fit=p_est(1:7);

first_loop=false;
% use the fit to calculate back F_data and T_data
[F_fit_scaled, T_fit_scaled]=Fl_model(params_fit_scaled, f_state1_fit); % with the scaling in place
[F_fit, T_fit]=Fl_model(params_fit, f_state1_fit);

%% plot the original data, scaled data, and the fits
figure; plot(F_fit); hold on; plot(F_data); title('intensity')
figure; plot(F_fit_scaled); hold on; plot(F_data_scaled); title('intensity scaled')
figure; plot(T_fit); hold on; plot(T_data); title('tau')
figure; plot(T_fit_scaled); hold on; plot(T_data); title('tau scaled')

% plot F vs. T
figure; plot(F_data, T_data, 'o'); hold on; plot(F_fit, T_fit, 'o'), title('T v F')

disp('FINAL FIT')
disp(params_fit)

%% Use the model to calculate the fraction in state 1 for all the data
[f_state1_from_F, f_state1_from_T]=Fl_model_invert(params_fit, Fd_all, Td_all);
figure; plot(f_state1_from_T); title('f state1 from tau')
figure; plot(f_state1_from_F); title('f state1 from intensity'); % !!! This one shows effects of bleaching as scale_factor is time dependent

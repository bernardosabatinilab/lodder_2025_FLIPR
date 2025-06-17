function res=fit_errors(p_est, F_data, T_data, fixed_params)
% calculate the residuals from the fit and the data 

% fix the parameters that need it
    if nargin==4
		p_est(~isnan(fixed_params))=fixed_params(~isnan(fixed_params));
    end

	[F_model, T_model]=Fl_model(p_est(1:7), p_est(8:end));

    res = [ F_model - F_data ; T_model - T_data];








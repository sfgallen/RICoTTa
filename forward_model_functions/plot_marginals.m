function plot_marginals(data,parameter_names,in_range,varargin)

%This function will generate a row plot of marginal posterior probability
% distribution functions using a 1D kernal denisity estimatior
%
% The data must be organized such that the number of rows is equal to the
% number of iterations and the number of columns is equal to the number of
% parameters
%
% The parameter_names input must be a cell of strings that is n-parameters
% long where each cell entry corresponds to parameter in each column of the
% data matrix
% 
% There is also an optional input for the MAP solution that will plot it as
% a dashed red line if provided.
%
% author: Sean F. Gallen
% data modified: 11/24/2019
% contact: sean.gallen[at]colostate.edu

% Parse Inputs
p = inputParser;
p.FunctionName = 'plot_marginals';

% required inputs
addRequired(p,'data', @(x) ismatrix(x));
addRequired(p,'parameter_names',@(x) iscell(x));
addRequired(p,'in_range', @(x) ismatrix(x));

% add optional inputs
addOptional(p,'map', [], @(x) isscalar(x)||isvector(x));

parse(p,data, parameter_names, in_range, varargin{:});
data  = p.Results.data;
parameter_names = p.Results.parameter_names;
in_range  = p.Results.in_range;
map = p.Results.map;


% make sure data and parameter names length matches
[~,nc] = size(data);
np = length(parameter_names);

if nc ~= np
    error('Error: number of data columns and parameter names is not equal...');
end

n_params = nc;

% find range of each parameter for x and y axis limit
bounds = in_range';

figure()
for i = 1:n_params
    % populate the figure by row
    data_x = data(:,i);
    x_bins = linspace(bounds(i,1),bounds(i,2),1000);
    
    subplot(1,n_params,i)
    
    [pdfx,xi] = ksdensity(data_x,x_bins);
    pdfx = pdfx./max(pdfx);
    plot(xi,pdfx,'k-'); hold on
    if ~isempty(map)
        plot([map(i),map(i)],[0 max(pdfx)],'r--')
    end
    
    xlim([bounds(i,:)])
    title(parameter_names{i});
    xlabel(parameter_names{i});

end
set(gcf,'Position',[100 400 2000 200])
end
% example script to show how to do the appraisal stage

load example_data.mat 
ndin=9;
rangein=zeros(2,ndin);
% allowed range of variables
rangein(1,:)=[10 0.1 0.1 0.01 0.1 1E-8 2/3 0.4 20];
rangein(2,:)=[30 3 2 0.5 2 1E-5 4/3 0.7 40];
scales(1:ndin)=-1;

% extract 1D and 2D marginals
nsobol=50000; % length of sobol sequence
dx1=20; % discretisation of 1D marginals
dx2=20; % discretisation of 2D marginals
[mean_1D,sd_1D,marg1D,marg2D]=sobol_resampling(models_all,misfit_all,ndin,rangein,scales,dx1,dx2,nsobol);
%[mean_1D,sd_1D,marg1Dthermo,marg2Dthermo]=sobol_resampling(models_all,logL_thermo_all,ndin,rangein,scales,dx1,dx2,nsobol)
%[mean_1D,sd_1D,marg1Dtcn,marg2Dtcn]=sobol_resampling(models_all,logL_tcn_all,ndin,rangein,scales,dx1,dx2,nsobol)
% combine 1D marginal probabilities
for i=1:ndin
    marg1Dsum=marg1D{i}.y;%.*marg1Dthermo{i}.y.*marg1Dtcn{i}.y;
    marg1Dsum=marg1Dsum/sum(marg1Dsum);
    x=marg1D{i}.x;
    figure
    bar(x,marg1Dsum,1,'k')
    xlabel(['P',num2str(i)])
    ylabel('Frequency')
    mean_parameter=sum(x.*marg1Dsum)
    sd_parameter=sqrt(sum(x.^2.*marg1Dsum)-mean_parameter^2)
end
% combine 2D marginal probabilities
c = combnk(1:ndin,2);
for i=1:size(c,1)
    Z=marg2D{i}.Z;
    Z=Z/sum(Z(:));
    x=marg2D{i}.x;
    y=marg2D{i}.y;
    figure
    imagesc(x,y,Z)
    set(gca,'YDir','normal')
    xlabel(['P',num2str(c(i,1))])
    ylabel(['P',num2str(c(i,2))])
end
function [mean_1D,sd_1D,marg1D,marg2D]=sobol_resampling(models,misfit,ndin,rangein,scales,dx1,dx2,nsobol)

% sobol resampling
sobol=NA_sobol(ndin,nsobol,0);
for i=1:length(misfit)
	models_sca(:,i)=transform2sca(models(:,i),ndin,rangein,scales);
end
IDX = knnsearch(models_sca',sobol');
prob = exp(-misfit(IDX));
sobol_models=models(:,IDX);

% plot 1D marginals and extract mean and sd
for j=1:ndin
%    figure
    x=linspace(rangein(1,j),rangein(2,j),dx1);
    for i=1:dx1-1
        y(i)=mean(prob(sobol_models(j,:)>x(i) & sobol_models(j,:)<=x(i+1)));
    end
    y(isnan(y))=0;
    y=y./sum(y);
    x2=x(1:end-1)+(x(2)-x(1))/2;
%     bar(x2,y,1,'k')
%     xlabel(['P',num2str(j)])
%     ylabel('Frequency')
    mean_1D(j)=sum(x2.*y);
    sd_1D(j)=sqrt(sum(x2.^2.*y)-mean_1D(j)^2);
    marg1D{j}.y=y;
    marg1D{j}.x=x2;
end


% plot all 2D marginals
c = combnk(1:ndin,2);
for o=1:size(c,1)
    p1=c(o,1);
    p2=c(o,2);
    x=linspace(rangein(1,p1),rangein(2,p1),dx2);
    x2=x(1:end-1)+(x(2)-x(1))/2;
    y=linspace(rangein(1,p2),rangein(2,p2),dx2);
    y2=y(1:end-1)+(y(2)-y(1))/2;
    [X,Y]=meshgrid(x,y);
    Z=zeros(size(X));
    for i=1:length(x)-1
        for j=1:length(y)-1
            Z(i,j)=mean(prob(sobol_models(p1,:)>X(i,j) & sobol_models(p1,:)<=X(i+1,j+1) & sobol_models(p2,:)>Y(i,j) & sobol_models(p2,:)<=Y(i+1,j+1)));
        end
    end
    Z(isnan(Z))=0;
    Z=Z./sum(Z(:));
%    figure
%    imagesc(x2,y2,Z)
%    set(gca,'YDir','normal')
%    xlabel(['P',num2str(c(o,1))])
%    ylabel(['P',num2str(c(o,2))])
    marg2D{o}.Z=Z;
    marg2D{o}.x=x2;
    marg2D{o}.y=y2;
end



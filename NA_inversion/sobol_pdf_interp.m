function [new_models_sca,ind]=sobol_pdf_interp(nsample,ncells,it,itmax,sobol,misfit,mfitord,I2)

n=ceil(nsample/ncells);
I=[];
for i=1:n
    I=[I mfitord(1:ncells)];
end
I=I(1:nsample);
nr=size(sobol,2)*(it/(itmax+1));
I3=find(misfit(1:nr)==0);

for i=1:nsample
    [IDX,D] = knnsearch(sobol(:,I3)',sobol(:,I2(I(i)))','K',1);
    ind(i)=I3(IDX);
    I3(IDX)=[];
    new_models_sca(:,i)=sobol(:,I3(IDX));
end


function corr_sample_size(sample,nu,nux,kflag)
% test correct size of free model parameters 

model_length=length(sample);
if kflag==0
    required_length=(nu-1)+nux^2*nu+3;
else
    required_length=(1+(2*(nu-1)))*nux^2+kflag*length(catchment.litho);
end

if model_length<required_length
    error_message=['Too few input parameters: ',int2str(required_length-model_length),' missing!'];
    error('ILM:Input',error_message)
elseif model_length>required_length
    error_message=['Too many input parameters: ',int2str(model_length-required_length),' remain!'];
    error('ILM:Input',error_message)
end


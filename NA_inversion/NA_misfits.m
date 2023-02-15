%----------------------------------------------------------------------------

%       NA_misfits - calculate performance statistics for NA algorithm.

%       Calls no other routines.

%					M. Sambridge, Oct. 1996

%----------------------------------------------------------------------------
function	[mfitmean,mfitminc,mfitmin,mfitord,mopt]=NA_misfits(misfit)

mfitmean=mean(misfit); % mean misfit
mfitminc = misfit(end); % misfit of last model
[mfitmin,mopt]=min(misfit); % overall lowest misfit and indice of model
% find models with lowest ncells misfit values
[C,mfitord]=sort(misfit);
end

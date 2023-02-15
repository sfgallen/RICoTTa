%----------------------------------------------------------------------------

%       NA_restart - resets NA walk to start from input model.

%       Calls no other routines.

%					M. Sambridge, Oct. 1996

%----------------------------------------------------------------------------

function [x,restartNA]=NA_restart(na_models,nd,mreset,x,debug)

if (debug)
    disp([' NA_restart: reset to model ',num2str(mreset)])
    disp(' current model on entry ')
    x(1:nd)
end

for i=1:nd
    x(i) = na_models(i,mreset);
end

restartNA = 1;

if (debug)
    disp(' current model on exit ')
    x(1:nd)
end

return
end
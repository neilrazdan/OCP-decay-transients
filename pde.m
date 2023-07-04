function [sol] = pde(X,T,R,eps,T_0,T_bulk,phi_bulk,alpha)

%% Solving stuff
m = 0;
t = T;
x = X;
% options = odeset('RelTol',1e-1);
sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);

% --------------------------------------------------------------
    function [c,f,s] = pdefun(x,t,u,dudx)
        c = [1; 1/R];
        f = [dudx(1); dudx(2)];
        s = [(1 - u(1).*u(2)); 1/R*(1 - u(1).*u(2))];
    end

% --------------------------------------------------------------
    function u_init = icfun(x)
        if x < alpha*eps^(1/2)
            u1 = alpha./((1-1/T_0).*(x*eps^(-1/2)) + alpha/T_0);
        else
            u1 = ((T_bulk-1)*(x*eps^(-1/2)) - (T_bulk*alpha-1))./(1-alpha);
        end
        u_init = [u1; 1./u1];
    end
% --------------------------------------------------------------
    function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
        pl = [0; 0];
        ql = [1; 1];
        pr = [ur(1) - T_bulk; ur(2) - phi_bulk];
        qr = [0; 0];
    end
end
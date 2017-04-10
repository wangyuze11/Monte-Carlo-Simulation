function [Nstep,Value] = BMBridge(M,T,b,K)
    h = 2^M;
    Nstep = 0;
    Path = zeros(h+1,1);
    Path(h+1) = normrnd(0,sqrt(T));
    Nstep = Nstep+1;
    Value = 1;
    if ((Path(h+1)>=b)||(Path(h+1)<K))
        Value = 0;
        return
    end
    for i = 1:M
        stp = 2^(M-i);
        t_l = 1;
        t_r = 2*stp+1;
        J_Max = 2^(i-1);
        for j = 1:J_Max
            mid = (t_l+t_r)/2;
            u = ((t_r-mid)*Path(t_l)+(mid-t_l)*Path(t_r))/(t_r-t_l);
            sigma = sqrt((t_r-mid)/h*(mid-t_l)/h/((t_r-t_l)/h));
            Path(mid) = u+sigma*normrnd(0,sqrt(T));
            Nstep = Nstep+1;
            if(Path(mid)>=b)
                Value = 0;
                return
            end
            t_l = t_l+2*stp;
            t_r = t_r+2*stp;
        end
    end
end
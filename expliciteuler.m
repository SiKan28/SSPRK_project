function [h,ts,ys] = expliciteuler(f, y0, t0, t_f,N)
    h=(t_f-t0)/N;
    ts = linspace(0, t_f, N+1);
    ys = zeros(length(y0),N+1);
    ys(:,1) = y0;
    for j=1:N
        ys(:,j+1)=ys(:,j)+h*f(ts(j), ys(:,j));
    end
end
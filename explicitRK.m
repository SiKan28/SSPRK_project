function [h, ts, ys] = explicitRK(f, y0, t0, t_f,N, A, b)
    c = sum(A,2);
    h=(t_f-t0)/N;
    ts = linspace(0, t_f, N+1);
    ys = zeros(length(y0),N+1);
    ys(:,1) = y0;
    for j=1:N
        %bk = zeros()
        K = zeros(length(y0), length(c));
        for i = 1:length(c)
            K(:,i) = f(ts(j)+h*c(i), ys(:,j)+h*K*A(i,:)');
        end
        ys(:,j+1)=ys(:,j)+h*K*b';
    end
end
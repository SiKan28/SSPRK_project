function [h, ts, ys] = RKsolver(f, y0, t_0, t_f, N, method)
    switch method
    % one-step methods
    % see explicit euler
    
    % two-step methods
        case 'genericTwoStep'
            A = [0, 0
                0.5, 0];
            b = [0, 1];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'heun'
            A = [0, 0
                 1, 0];
            b = [0.5, 0.5];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'ralston'
            A = [0, 0
                 2/3, 0];
            b = [1/4, 3/4];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
    % three-step methods
        case 'RK3'
            A = [0, 0, 0;
                 1/2, 0, 0;
                 -1, 2, 0];
            b = [1/6, 2/3, 1/6];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'heun3'
            A = [0, 0, 0;
                 1/3, 0, 0;
                 0, 2/3, 0];
            b = [1/4, 0, 3/4];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'ralston3'
            A = [0, 0, 0;
                 1/2, 0, 0;
                 0, 3/4, 0];
            b = [2/9, 2/3, 4/9];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'VDHW3'
            A = [0, 0, 0;
                 8/15, 0, 0;
                 1/4, 5/12, 0];
            b = [1/4, 0, 3/4];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'SSPRK3'
            A = [0, 0, 0;
                 1, 0, 0;
                 1/4, 1/4, 0];
            b = [1/6, 1/6, 2/3];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
    % four-step methods
        case 'RK4'
            A = [0, 0, 0, 0;
                 1/2, 0, 0, 0;
                 0, 1/2, 0, 0;
                 0, 0, 1, 0];
            b = [1/6, 1/3, 1/3, 1/6];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'threeEighths'
            A = [0, 0, 0, 0;
                 1/3, 0, 0, 0;
                 -1/3, 1, 0, 0;
                 1, -1, 1, 0];
            b = [1/8, 3/8, 3/8, 1/8];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'ralston4'
            A = [0, 0, 0, 0;
                 2/5, 0, 0, 0;
                 (-2889+1428*sqrt(5))/1024, (3785-1620*sqrt(5))/1024, 0, 0;
                 (-3365+2094*sqrt(5))/6040, (-975-3046*sqrt(5))/2552, (467040+203968*sqrt(5))/240845, 0];
            b = [(263+24*sqrt(5))/1812, (125-1000*sqrt(5))/3828, (5426304+1661952*sqrt(5))/5924787, (30-4*sqrt(5))/123];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        otherwise
            disp('That method is not supported yet')
    
    end
end
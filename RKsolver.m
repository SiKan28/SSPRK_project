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
    % five-step methods
    % six-step methods
    % seven-step methods
        case 'SSPRK75'      % CFL coefficient:1.178508348471858
            A = [0, 0, 0, 0, 0, 0, 0;
                 0.392382208054010, 0, 0, 0, 0, 0, 0;
                 0.310348765296963, 0.523846724909595, 0, 0, 0, 0, 0;
                 0.114817342432177, 0.248293597111781, 0, 0, 0, 0, 0;
                 0.136041285050893, 0.163250087363657, 0, 0.557898557725281, 0, 0, 0;
                 0.135252145083336, 0.207274083097540, -0.180995372278096, 0.326486467604174, 0.348595427190109, 0, 0;
                 0.082675687408986, 0.146472328858960, -0.160507707995237, 0.161924299217425, 0.028864227879979, 0.070259587451358, 0];
            b = [0.110184169931401, 0.122082833871843, -0.117309105328437, 0.169714358772186, 0.143346980044187, 0.348926696469455, 0.223054066239366];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
    % eigth-step methods
    % nine-step methods
    % ten-step methods
        case 'SSPRK104'
            A = zeros(10, 10)+tril(ones(10, 10)/6, -1);
            A(6:10, 1:5) = 1/15;
            b = ones(1, 10)/10;
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
        otherwise
            disp('That method is not supported yet')
    
    end
end

%TODO: https://nodepy.readthedocs.io/en/latest/modules/runge_kutta_method.html#nodepy.runge_kutta_method.shu_osher_to_butcher
%      https://ketch.github.io/numipedia/methods/SSPRK(10,4).html
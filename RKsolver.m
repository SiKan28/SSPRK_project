function [h, ts, ys] = RKsolver(f, y0, t_0, t_f, N, method)
    switch method
    % one-step methods
    % see explicit euler
    
    % two-step methods
        case 'midpoint'
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
        case 'SSPRK22'      % C = 1/2
            A = [0, 0
                 1, 0];
            b = [1/2, 1/2];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'ralston'
            A = [0, 0
                 2/3, 0];
            b = [1/4, 3/4];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
        case 'twoStepTest'
            A = [0, 0;
                 -20, 0];
            b = [41/40, -1/40];
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
        case 'SSPRK33'      % C = 1/3
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
        case 'SSPRK53'      % C = 0.53
            % alpha = [0, 0, 0, 0, 0; 1, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0.5666, 0, 0.4334, 0, 0; 0.0930, 0.00002, 0, 0.907, 0; 0.007, 0.2013, 0.0018, 0, 0.79];
            % beta = [0, 0, 0, 0, 0; 0.3773, 0, 0, 0, 0; 0, 0.3773, 0, 0, 0; 0, 0, 0.1635, 0, 0; 0.0007, 0, 0, 0.3422, 0; 0.0028, 0.00002, 0, 0, 0.2979];
            A = [0, 0, 0, 0, 0;
                 0.3773, 0, 0, 0, 0;
                 0.3773, 0.3773, 0, 0, 0;
                 0.1635, 0.1635, 0.1635, 0, 0;
                 0.1490, 0.1483, 0.1483, 0.3422, 0];
            b = [0.1972, 0.1179, 0.1172, 0.2703, 0.2979];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return

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
        case 'SSPRK85'      % CFL coefficient:1.8756849611641323
            A = [0, 0, 0, 0, 0, 0, 0, 0;
                 0.276409720937984, 0, 0, 0, 0, 0, 0, 0;
                 0.149896412080489, 0.289119929124728, 0, 0, 0, 0, 0, 0;
                 0.057048148321026, 0.110034365535150, 0.202903911101136, 0, 0, 0, 0, 0;
                 0.169059298369086, 0.326081269617717, 0.450795162456598, 0, 0, 0, 0, 0;
                 0.061792381825461, 0.119185034557281, 0.199236908877949, 0.521072746262762, -0.001094028365068, 0, 0, 0;
                 0.111048724765050, 0.214190579933444, 0.116299126401843, 0.223170535417453, -0.037093067908355, 0.228338214162494, 0, 0;
                 0.071096701602448, 0.137131189752988, 0.154859800527808, 0.043090968302309, -0.163751550364691, 0.044088771531945, 0.102941265156393, 0];
            b = [0.107263534301213, 0.148908166410810, 0.105268730914375, 0.124847526215373, -0.068303238298102, 0.127738462988848, 0.298251879839231, 0.156024937628252];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
    % nine-step methods
        case 'SSPRK95'      % CFL coefficient:2.695788289294857
            A = [0, 0, 0, 0, 0, 0, 0, 0, 0;
                 0.234806766829933, 0, 0, 0, 0, 0, 0, 0, 0;
                 0.110753442788106, 0.174968893063956, 0, 0, 0, 0, 0, 0, 0;
                 0.050146926953296, 0.079222388746543, 0.167958236726863, 0, 0, 0, 0, 0, 0;
                 0.143763164125647, 0.227117830897242, 0.240798769812556, 0, 0, 0, 0, 0, 0;
                 0.045536733856107, 0.071939180543530, 0.143881583463234, 0.298694357327376, -0.013308014505658, 0, 0, 0, 0;
                 0.058996301344129, 0.093202678681501, 0.109350748582257, 0.227009258480886, -0.010114159945349, 0.281923169534861, 0, 0, 0;
                 0.114111232336224, 0.180273547308430, 0.132484700103381, 0.107410821979346, -0.129172321959971, 0.133393675559324, 0.175516798122502, 0, 0;
                 0.096188287148324, 0.151958780732981, 0.111675915818310, 0.090540280530361, -0.108883798219725, 0.112442122530629, 0.147949153045843, 0.312685695043563, 0];
            b = [0.088934582057735, 0.102812792947845, 0.111137942621198, 0.158704526123705, -0.060510182639384, 0.197095410661808, 0.071489672566698, 0.151091084299943, 0.179244171360452];
            [h, ts, ys] = explicitRK(f, y0, t_0, t_f, N, A, b);
            return
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
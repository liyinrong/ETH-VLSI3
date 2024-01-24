%% Before running, set up the testbench cell name on line 8
% your decoder outputs should be labeled as Y0,Y1,Y2,Y3
% your decoder inputs should be labeled as A0,A1
%
% (c) 2022 Oscar Castaneda, Olalekan Afuye, Charles Jeon & Christoph Studer

% set up the name of your testbench cell
tb_name = 'decoder2to4_tb';

% set up cds_srr function
addpath('/usr/pack/cadence_mmsim-16.10.479-kgf/lnx86/tools.lnx86/spectre/matlab/64bit');

% directory that contains the simulation outputs
scratch_dir = strrep(getenv('HOME'),'home','scratch');
directory = sprintf('%s/simulation/%s/spectre/schematic/psf', scratch_dir, tb_name);

% set up basic parameters
Vdd = 1.2; % define vdd

% define how often the input signals change (in ps)
period_a1 =   5e3; % A1
period_a0 = 2.5e3; % A0

% get input signals
a1 = cds_srr(directory, 'tran-tran', 'A1', 0);
a0 = cds_srr(directory, 'tran-tran', 'A0', 0);

% convert time into ps
t_ps = a1.time*1e12;

% extract voltages of signals
a1 = a1.V;
a0 = a0.V;

% get output signals and put them together in a table where the i-th
% column corresponds to the 'Y(i-1)' output
y_mtx = [];
for i=1:4
    signal_name = ['Y',int2str(i-1)];
    y = cds_srr(directory, 'tran-tran', signal_name, 0);
    y_mtx = [y_mtx y.V];
end

exp_y_mtx = zeros(size(y_mtx));
sample_wvf = zeros(size(y_mtx));
mydecoder_output = zeros(4,4);
exp_decoder_output = zeros(4,4);

% during the fourth time the inputs repeat,
% we sample the inputs at 20% of the cycle's duration
t_ps_sample_in = 6*period_a1 + 0.2*period_a0 + (0:3)*period_a0;
% we sample the outputs at 80% of the cycle's duration
t_ps_sample_out = t_ps_sample_in + (0.8-0.2)*period_a0;

%% decoder output

% create base for expected output waveform
a_bits = (a1 > Vdd/2);
b_bits = (a0 > Vdd/2);
vec_bits = [a_bits b_bits];
exp_dec = bi2de(vec_bits,'left-msb');

%Check each one of the sampling points
err_flag = 0;
for i=1:4
    % find t_ps closest (from the right) to the t_ps_sample_in and _out
    t_ps_idx_in  = find(t_ps-t_ps_sample_in(i)>=0,1);
    t_ps_idx_out = find(t_ps-t_ps_sample_out(i)>=0,1);
    
    % measure the outputs and declare 1 if it is greater than Vdd/2    
    mydecoder_output(i,:) = y_mtx(t_ps_idx_out,:) > (Vdd/2);
    
    %create expected output waveform
    exp_y_mtx(:,i) = Vdd*(exp_dec == (i-1));
    
    %create sampling waveform
    sample_wvf(t_ps_idx_out,exp_dec(t_ps_idx_in)+1) = Vdd;
    
    % expected decoder output is given by:
    exp_decoder_output(i,exp_dec(t_ps_idx_in)+1) = 1;
    
    if (sum(exp_decoder_output(i,:) ~= mydecoder_output(i,:)))

        disp(['Expected output for input '...
            'A1=' num2str(vec_bits(t_ps_idx_in,1)) ...
            ' A0=' num2str(vec_bits(t_ps_idx_in,2)) ...        
            ' is y' num2str(4-i) ...
            ' but measured output is y' num2str(exp_dec(t_ps_idx_out))...
            ]) 
        err_flag  = err_flag + 1;
    end
end
sample_wvf = sum(sample_wvf,2);

if err_flag == 0
    disp('Your decoder circuit has no errors :)')
end
    %
%% plots
figure(1)
set(gcf,'units','pixel');
set(gcf,'position',[0,200,800,600]);
% we have 4 outputs so total time period elapsed is 2*period_a = 2000ps
subplot(2,1,1);
plot(t_ps,a1,'b',t_ps,a0,'r:','linewidth',2)
legend('A1','A0')
grid on
% note that our checks start for the fourth cycle
xlim(6*period_a1 + [0 2*period_a1])
title('inputs')
xlabel('time [ps]')
ylabel('input to decoder [V]')
set(gca,'fontsize',14)

% output from your circuit
subplot(2,1,2);
plot(t_ps,y_mtx,'linewidth',2)
legend('Y0','Y1','Y2','Y3')
grid on
xlim(6*period_a1 + [0 2*period_a1])
title('outputs')
xlabel('time [ps]')
ylabel('output [V]')
set(gca,'fontsize',14)

% print to SVG (vector graphics file)
print -dsvg Lab2_decoder_IO.svg

figure(2)
set(gcf,'units','pixel');
set(gcf,'position',[800,200,800,600]);
for i=1:4
    subplot(2,2,i)
    plot(t_ps,y_mtx(:,i),'k',t_ps,exp_y_mtx(:,i),'b:',t_ps,sample_wvf,'r--','linewidth',3);
    grid on
    legend('actual output','ideal output','sampling point','location','south')
    xlim(6*period_a1 + [0 2*period_a1])
    ylim([-1 1.5])
    xlabel('time [ps]')
    ylabel('output [V]')
    title(['Y' num2str(i-1)])
    set(gca,'fontsize',10)
end

% print to SVG (vector graphics file)
print -dsvg Lab2_decoder_expresp.svg
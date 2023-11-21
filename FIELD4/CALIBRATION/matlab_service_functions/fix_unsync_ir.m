%This function try to fix IR obtained trough Exponential Sine Sweep method 
%performed by unsynced play and record audio interfaces. 
%
%This result in spurious phase component proportional to f*(log(f)-1) 
%because the group delay, defined as the derivative of phase, 
%show a proportionality to log(f)
%
%The rx_speed estimation in ppm tell you the relative speed of the clock of
%the recording interface. If rx_speed>0 then recording interface clock is a
%little faster than players clock so the recorded sweep should be some 
%samples shorter.

function [irout,rx_speed] = fix_unsync_ir(ir,Fs,PlotName)
    ir = ir(:)';

    N = length(ir);

    %First of all move peak of IR in n=1 removing as much as possible IR pre-delay
    %This result in decreasing of linear phase component obtaining better from unwrap function
    [peak_v,peak_i] = max(abs(ir));    
    ir_shift = [ ir(peak_i:end) ir(1:peak_i-1) ];

    %compute unilateral FFT 
    H = fft(ir_shift);
    H = H(1:N/2+1)/N;
    H(2:N/2) = 2*H(2:N/2);
    f = linspace(0,Fs/2,N/2+1);

    %Compute dB magnitude
    mH = 20*log10(abs(H));
    mH = mH - max(mH);

    %Find useful frequency interval
    SNRth = 20;
    i1 =  find(mH>(-SNRth),1,'first')+2; 
    i2 =  find(mH>(-SNRth),1,'last')-round(N/16);

    %Compute unwrapped phase
    pH = unwrap(angle(H));
    pH = pH - pH(i1);



    % %Compute group delay
    % gd = [0 diff(pH)];
    % 
    % %Try to fit group delay
    % gd2 = gd(i1:i2);
    % f2 = f(i1:i2);
    % myfun = @(var) sum( ( (var(1)+var(2).*log(f2)) - gd2 ).^2 );
    % var = fminsearch(myfun, [0,0.01]);
    % gd_fit = (var(1)+var(2)*log(f));
    % 
    % %Plot fitting result
    % figure
    % plot(f2,gd2,'b-*',f,gd_fit,'r-')
    % title('Fitting of group delay dependency from clock desync ')
    % 
    % %Crea il filtro rifasatore
    % pR = [0 cumsum(gd_fit(2:end))];
    % pR(N/2+1) = 0;
    % HR = exp(-j*pR);
    % HR = [HR flipdim(conj(HR(2:end-1)),2)];




    %Try to fit directly unwarped phase
    pH2 = pH(i1:i2);
    f2 = f(i1:i2);
    myfun = @(k) sum( ( (k(1).*f2+k(2).*(f2.*(log(f2)-1))+k(3)) - pH2 ).^2 );
    options = optimset('TolFun',1e-9);
    [k,fval,exitflag,output] = fminsearch(myfun, [0 0 0], options);
    pH_fit = (k(1).*f+k(2).*(f.*(log(f)-1))+k(3));



    %Crea il filtro rifasatore utilizzando la sola componente di fase relativa
    %al termine logaritmico riconducibile allo skew dei clock
    %pR = (k(1).*f+k(2).*(f.*(log(f)-1)));
    pR = (k(2).*(f.*(log(f)-1)));
    pR(1) = 0;
    pR(N/2+1) = 0;
    HR = exp(-j*pR);
    HR = [HR flipdim(conj(HR(2:end-1)),2)];


    %Trasformalo nel dominio del tempo
    irR = ifft(HR);
    irR = [ irR(N/2:end) irR(1:N/2-1) ];

    %Rifasa il segnale
    irout = fd_conv(ir,irR);



    %Estimate RX relative clock speed difference
    rx_speed = k(2)*7.3288e+04;



    %Some plot to show result
    if exist('PlotName','var')
        figure('Name',PlotName,'NumberTitle','off')

        subplot(3,2,1)
        irplot(ir)
        title('Input IR')

        subplot(3,2,3)
        irplot(irR)
        title('IR of restoring filter')

        subplot(3,2,5)
        irplot(irout')
        title('Restored IR')

        ax(1)=subplot(3,2,2);
        [ax2,h1,h2] = plotyy(f,mH,f,angle(H),'semilogx');
        hold(ax2(1),'on');
        semilogx(f,-SNRth*ones(N/2+1,1),'r--')
        semilogx(f(i1)*ones(2,1),[min(mH),max(mH)],'r--',f(i2)*ones(2,1),[min(mH),max(mH)],'r--')
        ylim(ax2(1),[-SNRth*3 0])
        title('Magnitude and Phase of input IR')

        ax(2)=subplot(3,2,4);
        semilogx(f,pH)
        hold on;
        semilogx(f(i1)*ones(2,1),[min(pH),max(pH)],'r--',f(i2)*ones(2,1),[min(pH),max(pH)],'r--')
        ylim([min(pH(i1:i2)) max(pH(i1:i2))]);
        title('Unwrapped phase of input IR')

        linkaxes([ax(2) ax(1)],'x');

        subplot(3,2,6);
        plot(f2,pH2,'b-*',f,pH_fit,'r-')
        title('Fitting of phase dependency from clock desync ')
    end





end

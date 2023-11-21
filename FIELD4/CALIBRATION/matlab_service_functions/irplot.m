%
%irplot(ir_matrix[,fs[,xmode]])
%
%Plot a matrix of impulse response:
%ir_matrix : MxN impulse response matrix
%       fs : Sampling frequency (optional)
%    xunit : X axis unit 'n' sample, 's' second, 'm' meter (optional)
function irplot(ir_matrix,fs,xmode)

    M = size(ir_matrix,1);
    N = size(ir_matrix,2);
                          
    %Timebase conversion to meter is disabled as default setting sound speed to 1
    c = 1;
    %Timebase conversion to second is disabled as default setting sampling frequency to 1
    if not(exist('fs'))
        fs = 1;
    end
    %Manage timebase conversion from sample to second or meter
    if exist('xmode')
        switch xmode
            case 'n', fs = 1;  %Set fs = 1 to obtain timebase in sample
            case 'm', c = 343; %Set c = 1 to obtain timebase in meter
        end
    end
    
    
    %Upsample IR for smooth line visualization
    U = 50;
    ir_matrix_up = zeros(M,N*U);
    for m = 1:1:M
        ir_matrix_up(m,:) = interp(ir_matrix(m,:),U);
    end
    ir_timebase_up = (0:1:(N*U-1))/fs*c/U;
    plot(ir_timebase_up,ir_matrix_up')
    
    
    %Plot marker on original IR sample
    hold on
    ir_timebase = (0:1:(N-1))/fs*c;
    plot(ir_timebase,ir_matrix','.')
    hold off
    
    
%     %Frequency Domain Plot
%     FFTN = max( 2*ceil(N/2) , 1024 );
%     IR_matrix = fft(ir_matrix,2,FFTN);
%     IR_matrix = [2*IR_matrix(:,1) , IR_matrix(:,2:FFTN/2+1)];
%     df = fs / FFTN;
%     freqbase = df*(0:1:FFTN/2); 
%     figure
%     semilogx(freqbase,20*log10(abs(IR_matrix')))
%      hold on
%      semilogx(freqbase,20*log10(unwrap(angle(IR_matrix'))))
%      hold off
    
end


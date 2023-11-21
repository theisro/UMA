%flyingtime = flyingtime_estimation(ir_matrix[,n_start[,n_end]])
%
%Estimate flyingtime of each IR contained in an impulse response matrix
%Result is a vector of flyingtime expressed in sample with fractional part 
%
% ir_matrix : Matrix of IR (MxN) 

%flyingtime : Vector of M flyingtime expressed in sample

function flyingtime = fd_flyingtime_estimation(ir_matrix,fnl,fnh)
    
     %Reshape to obtain 2D matrix preserving last dimension
    ir_matrix_size = size(ir_matrix);
    ir_matrix = reshape (ir_matrix , prod(ir_matrix_size(1:end-1)) , ir_matrix_size(end) );
 
    M = size(ir_matrix,1);
    N = size(ir_matrix,2);
    N = 2*floor(N/2);
    
    if not(exist('fnl'))
        fnl = 0.005; %Fs=48kHz -> Fi=240Hz  
    end
    
    if not(exist('fnh'))
        fnh = 0.2; %Fs=48kHz -> Ff=9.6kHz  
    end
    
    %find approximate peak position
    [v,ipk] = max(abs(ir_matrix),[],2);
    ipk = round(mean(ipk));
    
    %remove initial delay
    ir_matrix = [ir_matrix(:,ipk+1:end) , ir_matrix(:,1:ipk)];
    
    
    %switch to frequency domain
    IR_matrix = fft(ir_matrix,N,2);
    
    %Calcola il rapporto complesso "riga per riga" tra gli spettri
    IRRatio = IR_matrix; %./IR2;

    %estrai modulo e fase unwarpata (solo primi N/2 campioni)
    %abs_IRRatio = abs(IRRatio(:,1:N/2));
    angle_IRRatio = zeros(M,N/2);
    for m = 1:1:M
        angle_IRRatio(m,:) = unwrap(angle(IRRatio(m,1:N/2)));
    end
    
    ki = ceil(N*fnl);
    kf = floor(N*fnh);
    K = kf-ki+1;
    angle_IRRatio = angle_IRRatio(:,ki:kf);

    %Ora bisogna calcolare la pendenza media del tratto lineare...ad esempio
    %facendo la media della derivata discreta
    %slope = -(mean(angle_IRRatio(:,2:K)-angle_IRRatio(:,1:K-1),2));
    slope = zeros(M,1);
    for m=1:1:M
        p = polyfit(1:1:K,angle_IRRatio(m,:),1);
        slope(m) = -p(1);
    end

    flyingtime = slope * N/(2*pi())+ipk;
 
     
    %Restore original matrix dimension
    if length(ir_matrix_size) > 2
        flyingtime = reshape(flyingtime,ir_matrix_size(1:end-1));
    end
    
end
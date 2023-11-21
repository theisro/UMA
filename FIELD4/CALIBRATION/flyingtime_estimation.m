%flyingtime = flyingtime_estimation(ir_matrix[,n_start[,n_end]])
%
%Estimate flyingtime of each IR contained in an impulse response matrix
%Result is a vector of flyingtime expressed in sample with fractional part 
%
% ir_matrix : Matrix of IR (...xN) 
%   n_start : Start of peak searching window (optional)
%     n_end : End of peak searching window (optional)
%
%flyingtime : Vector of M flyingtime expressed in sample
function flyingtime = flyingtime_estimation(ir_matrix,n_start,n_end)

    %Reshape to obtain 2D matrix preserving last dimension
    ir_matrix_size = size(ir_matrix);
    ir_matrix = reshape (ir_matrix , prod(ir_matrix_size(1:end-1)) , ir_matrix_size(end) );
    
    M = size(ir_matrix,1);
    N = size(ir_matrix,2);
    
    if not(exist('n_start'))
        n_start = 1;
    end
    
    if not(exist('n_end'))
        n_end = N;
    end
    
    W2 = 10;    %half search window
    U = 100;    %upscaling factor
    
    flyingtime = zeros(M,1);
    
    for m = 1:1:M
    
        %Rendi la polaritŕ della risposta all'impulso positiva
        if ( max(squeeze(ir_matrix(m,:))) < abs(min(squeeze(ir_matrix(m,:)))) )
            ir_matrix = -1.*ir_matrix;
        end
    
        %Find the impulse response aproximative peak position working on IR 
        [v,i] = max( abs( ir_matrix(m,n_start:n_end) ) );
        if (v~=0)
            %Mussi Fixing
            %[vmax,imax] = max( ( ir_matrix(m,n_start:n_end) ) );
            %[vmin,imin] = min( ( ir_matrix(m,n_start:n_end) ) );
            %i = min(imax,imin);
        
            i = i + n_start - 1;

            %Upsample a window around the peak
            ws = max( 1 , i-W2 );   %Upsampling window start
            we = min( i+W2 , N );   %Upsampling window end
            ir_upsampled = interp( double( ir_matrix(m,ws:we) ) , U );

            %Find the impulse response precise peak position working on upsampled IR
            [v,i] = max( abs( ir_upsampled ) );
            
            %Mussi Fixing
            %[vmax,imax] = max( ir_upsampled );
            %[vmin,imin] = min( ir_upsampled );
            %i = min(imax,imin);
            
            flyingtime(m) = ( i/U + ws-1 );
        end
    
    end
    
    %Restore original matrix dimension
    if length(ir_matrix_size) > 2
        flyingtime = reshape(flyingtime,ir_matrix_size(1:end-1));
    end

end
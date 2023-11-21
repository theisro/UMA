%flyingtime = flyingtime_estimation2(ir_matrix,[n_start,n_end],plot)
%
%This modified version keep in consideration first peak also if other peak
%are present
%
%Estimate flyingtime of each IR contained in an impulse response matrix
%Result is a vector of flyingtime expressed in sample with fractional part 
%
% ir_matrix : Matrix of IR (...xN) 
%   n_start : Start of peak searching window (optional)
%     n_end : End of peak searching window (optional)
%
%flyingtime : Vector of M flyingtime expressed in sample
function flyingtime = flyingtime_estimation2(ir_matrix,limit,PlotName)

    %Reshape to obtain 2D matrix preserving last dimension
    ir_matrix_size = size(ir_matrix);
    ir_matrix = reshape (ir_matrix , prod(ir_matrix_size(1:end-1)) , ir_matrix_size(end) );
    
    M = size(ir_matrix,1);
    N = size(ir_matrix,2);
    
    if (exist('limit','var') && (length(limit)==2))
        n_start = limit(1);
        n_end = limit(2);
    else
        n_start = 1;
        n_end = N;     
    end
    
    %Rendi la polaritŕ della matrice di IR positiva
    if ( sum(max(ir_matrix,[],2)) < abs(sum(min(ir_matrix,[],2))) )
        ir_matrix = -ir_matrix;
    end

    W2 = 4;    %half search window
    U = 100;    %upscaling factor
    th = 0.25;   %trigger treshold
    
    flyingtime = zeros(M,1);
    showmax = zeros(M,2);
    
    %For each IR in matrix
    for m = 1:1:M
        
        %Find the impulse response magnitude peak
        [v,i] = max( ir_matrix(m,n_start:n_end) );
        

        if (v~=0)
            %Find first trigger of 50% of peak magnitude
            i = find(ir_matrix(m,n_start:n_end)>(v*th),1);
            i2 = find(ir_matrix(m,i:n_end)<(v*th),1);
            
            %Mussi Fixing
            %[vmax,imax] = max( ( ir_matrix(m,n_start:n_end) ) );
            %[vmin,imin] = min( ( ir_matrix(m,n_start:n_end) ) );
            %i = min(imax,imin);
        
            i = i + n_start - 1;

            %Upsample a window around the peak
            ws = max( 1 , i-W2 );   %Upsampling window start
            we = min( i+W2 , N );   %Upsampling window end
            ir_upsampled = interp( double( ir_matrix(m,ws:we).*[0.25 0.5 1 1 1 1 1 0.5 0.25] ) , U );

            %Find the impulse response precise peak position working on upsampled IR
            [v,i] = max( ( ir_upsampled ) );
            
            %Mussi Fixing
            %[vmax,imax] = max( ir_upsampled );
            %[vmin,imin] = min( ir_upsampled );
            %i = min(imax,imin);
            
            flyingtime(m) = ( i/U + ws-1 );
            
            showmax(m,:) = [flyingtime(m);v];
        end
    
    end
    
    if (exist('PlotName','var'))
        if isnumeric(PlotName)
            set(0,'CurrentFigure',PlotName)
        else
            figure('Name',PlotName,'NumberTitle','off')
        end
        irplot(ir_matrix)
        hold on
        ColorArray = ['b','g','r','c','m','y','k'];
        for m = 1:1:M
            myColor = ColorArray(mod(m-1,7)+1);
            plot(showmax(m,1),showmax(m,2),'v','Color',myColor,'MarkerFaceColor',myColor,'markersize',10);
        end
        hold off
        xlim([min(flyingtime)*0.95,max(flyingtime)*1.05]);
    end
    
    %Restore original matrix dimension
    if length(ir_matrix_size) > 2
        flyingtime = reshape(flyingtime,ir_matrix_size(1:end-1));
    end

end
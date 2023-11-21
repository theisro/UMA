%out_data = TrimWinData(in_data,Skip,N[,WinType,NWin,FadeIn,FadeOut[,PlotName]])
%
%This function work on the last dimension of in_data matrix
%performing Trimming and Zero Padding of the data series and optionally
%apply fade in and fade out defined by standard windowing function like hann
%
%PlotName is an optional string becameing the name of a new plot showing windowing result
%If omitted plot is not generated
%
%Example:
%out_data = TrimWinData(in_data,128,1024,'hann',512,128,256,'MyWinPlot');
%Select 1024 samples starting from 128 then apply hann windowing
%with 128 samples of fade in and 256 of fade out.
%
function out_data = TrimWinData(in_data,Skip,N,WinType,NWin,FadeIn,FadeOut,PlotName)
   
    %Reshape to obtain 2D matrix preserving last dimension
    data_size = size(in_data);
    P = prod(data_size(1:end-1));
    L = data_size(end);
    in_data = reshape(in_data,P,L);
   
    %Create output matrix
    out_data = zeros(P,N);
    NCpy = min(L-Skip,N);
     
    if exist('WinType','var')
        %Create window function with selected profile, fadein and fadeout
        eval(sprintf('fadein_win = %s(%d);',WinType,2*FadeIn))
        fadein_win = fadein_win(1:FadeIn);
        eval(sprintf('fadeout_win = %s(%d);',WinType,2*FadeOut))
        fadeout_win = fadeout_win(FadeOut+1:end);    
        win = [ fadein_win ; ones(NWin-FadeIn-FadeOut,1) ; fadeout_win ; zeros(NCpy-NWin,1) ];
    else
        %No windowing
        win = ones(NCpy,1);
    end

     
    %Copy data applying windowing to the last dimension
    for p = 1:1:P
        out_data( p , 1:NCpy ) = in_data( p , (Skip+1:Skip+NCpy) ) .* win(1:NCpy)';
    end
    
    
     %Plot data and win
     if exist('PlotName','var')
        figure('Name',PlotName,'NumberTitle','off')
        subplot(2,1,1)
        hold on
        plot(in_data')
        plot([zeros(Skip,1) ; win]'*max(max(abs(in_data))))
        hold off
        subplot(2,1,2)
        plot(out_data')
     end
     
         
    %Restore original matrix dimension
    data_size(end) = N;
    out_data = reshape(out_data,data_size);
           
         
end

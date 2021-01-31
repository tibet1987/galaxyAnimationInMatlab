function [ gifFile ] = exportImagesToGif( figureHandle, mode, varargin )
% gifFile = exportImagesToGif( figureHandle, mode)
% gifFile = exportImagesToGif( figureHandle, mode, VarName, VarValue, ...)
%
%
% Thiis function collects images from figures and puts them together for an
% export to a gif file
%
% input:
%  figureHandle  - handle to figure, obtained with call 'get(gcf)'
%  mode          - options: 
%                     1) 'collect'   -> collect images 
%                     2) 'finalize'  -> create gif out of collected images
%  optional:        VarName      - Varvalue
%                     'Ts'       - delay time between frames e.g. 0.05
%                     'loopCnt'  - 'infinite', number of times
%
% output:
%  gifFile       - gif output file combining all the images

persistent im map

if nargin>2
    for k=1:nargin-2
        if strcmpi(varargin{k},'Ts')
            Ts = varargin{k+1};
        elseif strcmpi(varargin{k},'loopCnt') || strcmpi(varargin{k},'loopCount')
            loopCnt = varargin{k+1};
        elseif strcmpi(varargin{k},'gifName')
            gifname = varargin{k+1};
        end
    end
end

gifFile = 0; % default output

if isempty(im)
    % f = getframe(gcf);
    [im,map] = rgb2ind(figureHandle.cdata,256,'nodither'); %initialize image sequence
    if strcmpi(mode,'finalize')
        error('Function ''exportImagesToGif'': Mode ''finalize'' is only possible with more than one image')
    end
else
    if strcmpi(mode,'collect')
        im(:,:,1,end+1) = rgb2ind(figureHandle.cdata,map,'nodither');
    elseif strcmpi(mode,'finalize')
        imwrite(im,map,gifname,'DelayTime',Ts,'LoopCount',loopCnt) %g443800
        disp('gif-export complete')
    end
end


end


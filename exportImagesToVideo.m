function [ vidFile ] = exportImagesToVideo( figureHandle, mode, varargin )
% gifFile = exportImagesToVideo( figureHandle, mode)
% gifFile = exportImagesToVideo( figureHandle, mode, VarName, VarValue, ...)
%
%
% Thiis function collects images from figures and puts them together for an
% export to a video (mp4) file
%
% input:
%  figureHandle  - handle to figure, obtained with call 'get(gcf)'
%  mode          - options: 
%                     1) 'collect'   -> collect images 
%                     2) 'finalize'  -> create gif out of collected images
%  optional:        VarName      - Varvalue
%                     'Ts'       - delay time between frames e.g. 0.05
%                     'loopCnt'  - 'infinite', number of times
%                     'vidName'  - name of the video file
%                     'compRatio'- compression ratio (default is 10)
%
% output:
%  gifFile       - gif output file combining all the images

persistent writerObj

vidName = 'videoOut';
FrameRate = 30;
if nargin>2
    for k=1:nargin-2
        if strcmpi(varargin{k},'Ts')
            if varargin{k+1} > 0.01
                FrameRate = 1/varargin{k+1};
            else
                FrameRate = 100;
            end
        elseif strcmpi(varargin{k},'loopCnt') || strcmpi(varargin{k},'loopCount')
            loopCnt = varargin{k+1};
        elseif strcmpi(varargin{k},'vidName')
            if strcmpi(varargin{k+1}(end-3:end),'.mp4') || strcmpi(varargin{k+1}(1:end-4),'.avi')
                vidName = varargin{k+1}(1:end-4);
            else
                vidName = varargin{k+1};
            end
        end
    end
end

vidFile = 0; % default output

if isempty(writerObj)
    if strcmpi(mode,'finalize')
        error('Function ''exportImagesToVideo'': Mode ''finalize'' is only possible with more than one image')
    end
    writerObj = VideoWriter(vidName,'MPEG-4'); % create a new video object
    writerObj.FrameRate = FrameRate;
    writerObj.Quality = 100;
    open(writerObj); % open the newly created video object for writing
%     f = getframe(figureHandle); % get first frame
    writeVideo(writerObj,figureHandle); % write first frame into file
else
    if strcmpi(mode,'collect')
        writeVideo(writerObj,figureHandle); % write first frame into file
    elseif strcmpi(mode,'finalize')
        close(writerObj);
        disp('video export complete')
    end
end


end


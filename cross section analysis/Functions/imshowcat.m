function varargout = imshowcat(varargin)
%Combines and imshows multiple arrays in mixed fmt (uintX,double,logical)
%
%function imshowcat(varargin)
%
% This function takes up to three arrays and combines them together to make
% an 8-bit RGB array suitable for imshowing. If an array is passed, and the
% next argument is a scalar integer between zero and 255, then that is the
% intensity multiplier of that array.
%
% If a structure is passed, and that structure has the field "PerimXYList",
% then those pixels are also plotted on the image. Up to three such
% structures can be present and the colors that the curves will have will
% be c,m,y.
%
% At least one argument passed in must be an array (cannot have all
% structure variables), and the first argument must be an array. The only
% exception is that the first argument could be a figure handle. If so, it
% will be extracted and removed from varagin, and the rest of the function
% will proceed as normal, in which case, the original *second* argument
% must be an array.

%
% Super-special case in which the first argument is a figure handle. If so,
% then it will be extracted and removed from varargin, and we'll proceed as
% if it was not there previously. We needed to add several extra conditions
% here because apparently, ishandle will operate on arrays and return a
% logical array, and also because zero must be a handle (by default...the
% screen?).
%
if size(varargin{1},1) == 1 && size(varargin{1},2) == 1 && ...
		varargin{1} ~= 0 && ishandle(varargin{1})
	fh = varargin{1};
	varargin(1) = [];
end

%
% Unpacking varargin
%
i_array = 1;
i_struct = 1;
nArg = size(varargin,2); iArg = 1;
S = cell(nArg,1);
while iArg <= nArg
	I = varargin{iArg};
	
	%
	% Special case of iArg = 1, to pre-allocate and to make sure the first
	% input is an array.
	%
	if iArg == 1
		[H,W,D] = size(I);
		if H == 1 || W == 1
			error('First input must be image array')
		end
		if D > 1
			error('Image array must be 2D; cannot be 3D')
		end
		IM = uint8(zeros(H,W,3)); % set up RGB array		
	end
	
	if ~isstruct(I) 
		[H1,W1,D1] = size(I);
		if H1 ~= H || W1 ~= W || D1 ~= D
			error('image sizes must match')
		end
		
		if isa(I,'uint16') || isa(I,'uint8')
			I = mDU(double(I));
		elseif isa(I,'double')
			I = mDU(I);
		elseif islogical(I)
			I = double(I);
		end		
		
		%
		% Now check to see if subsequent argument is scalar multiplier
		%		
		if iArg < nArg && isnumeric(varargin{iArg+1}) && isscalar(varargin{iArg+1})
			chi = varargin{iArg+1};
			iArg = iArg + 1;
			
			if chi < 0 || chi > 255
				chi = 255;
			end
		else
			chi = 255;
		end
		
		%
		% Add the single image to the RGB image
		%
		IM(:,:,i_array) = uint8(chi*I);
		i_array = i_array + 1;
		
	elseif isfield(I,'PerimXYList') % special case of structure array
		S{i_struct} = I;
		i_struct = i_struct + 1;
	end	
	iArg = iArg + 1;
	
end
S(i_struct:end) = [];
n_struct = length(S);


if nargout == 0
	
	%
	% Focus on the given figure handle, if it was passed as the original
	% first argument of varargin. Otherwise, open a new figure.
	%
	if exist('fh','var')
		if any(findall(0,'Type','Figure') == fh) % check if figure is already open
			s = get(fh,'visib');
			
			%
			% If the figure is open and the visibility is off, then check
			% if desired figure is current. If so, then don't do anything.
			% If not, then we have to make it current, which will turn its
			% visibility on.
			%
			if strcmp(s,'off') && ~isequal(fh,gcf)
				figure(fh) % this flashes the figure
				set(fh,'visib','off')
			end
		else % open new version of the figure if it's not already open
			figure(fh)
		end
	else
		figure
	end
	
	%
	% imshow the RGB image array.
	%
	imshow(IM)
	hold on
	
	%
	% Now we plot the pixels in the any structure array we have obtained
	%
	COLR = {'c','y','g'};
	for i = 1:n_struct
		s = S{i};
		n = length(s);
		for j = 1:n
            
            if(COLR{i} == 'c')
                plot(s(j).PerimXYList(:,1),s(j).PerimXYList(:,2),COLR{i},'LineWidth',1)
            else
                plot(s(j).PerimXYList(:,1),s(j).PerimXYList(:,2),COLR{i})
            end
		end
	end
	
else
	varargout = {IM};
end







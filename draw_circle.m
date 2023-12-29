function [h,xdata,ydata] = draw_circle(varargin)
%
% [h,xdata,ydata] = DRAW_CIRCLE()
%  Draws a circle with the following default options
%   center = 0
%   radius = 1
%   N = 20 % Number of points
%   initialAngle = 0
%   finalAngle = 2*pi
%  Returns the handle h to the line and the (xdata,ydata) of each point in
%  the circle
%
% DRAW_CIRCLE('option_name',option_value)
%  The function accepts any number of pairs of option_name and option_value 
%  to change the default behavior
% 
% Author: Pedro Casau
%
    tti = 0;
    ttf = 2*pi;
    N = 20;
    c = zeros(2,1);
    r = 1;
    marker = 'none';
    arrow = 0;
    if mod(nargin,2) ~= 0 
        ME = MException('draw_circle:input','DRAW_CIRCLE requires pairwise inputs.');
        throw(ME);
    end
    for I = 1:2:nargin
        switch varargin{I}
            case 'center'
                c = varargin{I+1};
            case 'radius'
                r = varargin{I+1};
            case 'N'
                N = varargin{I+1};
            case 'initialAngle'
                tti = varargin{I+1};
            case 'finalAngle'
                ttf = varargin{I+1};
            case 'Marker'
                marker = varargin{I+1};
            case 'Arrow'
                arrow = varargin{I+1};
        end
    end
    tt = linspace(tti,ttf,N);
    xdata = c(1)+r*cos(tt);
    ydata = c(2)+r*sin(tt);
    h = plot(xdata,ydata,'Marker',marker);
    color = get(h,'color');
    if arrow 
        line([xdata(end-1)+0.1*r*(xdata(end-1)-c(1)),xdata(end-1)-0.1*r*(xdata(end-1)-c(1));
              xdata(end),                     xdata(end)],...
              [ydata(end-1)+0.1*r*(ydata(end-1)-c(2)),ydata(end-1)-0.1*r*(ydata(end-1)-c(2));
              ydata(end),                     ydata(end)],...
              'color',color);
    end
              
        
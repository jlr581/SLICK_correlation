function correlation = correlation(y1,y2,x1,x2,varargin)
%CORRELATION - calculate correlation of unevenly and differently sampled data
%              as per methods described in Roberts et al. (2019)
%
%   Inputs: 
%   y1      data of timeseries 1
%   y2      data of timeseries 2
%   x1      time stamps corresponding to y1
%   x2      time stamps corresponding to y2
%   optional: hc (default 0.4)
%
%   Usage:
%   [correlation] = correlation(series1,series2,time1,time2,0.3);

	n=length(y1);
	m=length(y2);

   if nargin<5
		disp(['No value for hc specified, using hc=0.4']);
      hc=0.4;
    else
      hc=varargin{1};
   end

	% find valid domain for both series, discard any data more than one point
	% outside range of other data series
	xmin=max(x1(1),x2(1));
	xmax=min(x1(n),x2(m));
	if xmin==x1(1)
		start_1=1;
		start_2=min(find(x2>xmin))-1;
	else
		start_2=1;
		start_1=min(find(x1>xmin))-1;
	end
	if xmax==x1(n)
		end_1=n;
		end_2=min(find(x2<xmax,1,'last'))+1;
	else
		end_2=m;
		end_1=min(find(x1<xmax,1,'last'))+1;
	end

	% find median and interquartile range of data spacing
	dx=x1(start_1+1:end_1)-x1(start_1:end_1-1);
	dx=sort(dx);
	if mod(end_1-start_1,2)==1 % odd number of points
		median1=dx(floor((end_1-start_1)/2)+1);
		i=floor((end_1-start_1)/4);
		if mod(end_1-start_1,4)==1
			temp1=(dx(i)+3*dx(i+1))/4;
			temp2=(3*dx(3*i+1)+dx(3*i+2))/4;
		else
			temp1=(3*dx(i+1)+dx(i+2))/4;
			temp2=(dx(3*i+2)+3*dx(3*i+3))/4;
		end
	else
		median1=(dx((end_1-start_1)/2)+dx((end_1-start_1)/2+1))/2;
		temp1=(dx(floor((end_1-start_1)/4))+dx(floor((end_1-start_1)/4)+1))/2;
		temp2=(dx(floor(3*(end_1-start_1)/4)+1)+dx(floor(3*(end_1-start_1)/4)+2))/2;
	end
	iqr1=temp2-temp1;

	dx(1:end_2-start_2)=x2(start_2+1:end_2)-x2(start_2:end_2-1);
	dx=sort(dx);
	if mod(end_2-start_2,2)==1 % odd number of points
		median2=dx(floor((end_2-start_2)/2)+1);
		i=floor((end_2-start_2)/4);
		if (mod(end_2-start_2,4)==1)
			temp1=(dx(i)+3*dx(i+1))/4;
			temp2=(3*dx(3*i+1)+dx(3*i+2))/4;
		else
			temp1=(3*dx(i+1)+dx(i+2))/4;
			temp2=(dx(3*i+2)+3*dx(3*i+3))/4;
		end
	else
		median2=(dx((end_2-start_2)/2)+dx((end_2-start_2)/2+1))/2; %something wrong here
		temp1=(dx(floor((end_2-start_2)/4))+dx(floor((end_2-start_2)/4)+1))/2;
		temp2=(dx(floor(3*(end_2-start_2)/4)+1)+dx(floor(3*(end_2-start_2)/4)+2))/2;
	end
	iqr2=temp2-temp1;

	% set interval around points
	h=hc*max([median1,median2,iqr1,iqr2]);

	% construct points where kernel starts, midpoint, stops
	j=1;
	for i=start_1:end_1
		if x1(i)-h > x1(start_1)
			xbase(j)=x1(i)-h;
			j=j+1;
		end
		xbase(j)=x1(i);
		j=j+1;
		if x1(i)+h < x1(end_1)
			xbase(j)=x1(i)+h;
			j=j+1;
		end
	end
	for i=start_2:end_2
		if (x2(i)-h > x2(start_2))
			xbase(j)=x2(i)-h;
			j=j+1;
		end
		xbase(j)=x2(i);
		j=j+1;
		if (x2(i)+h < x2(end_2))
			xbase(j)=x2(i)+h;
			j=j+1;
		end
	end
	n_xbase=j-1;

	xbase=sort(xbase(1:n_xbase));
	i=2;

	while i<=n_xbase
		if (xbase(i)==xbase(i-1))
			xbase(i:n_xbase-1)=xbase(i+1:n_xbase);
			n_xbase=n_xbase-1;
		else
			i=i+1;
		end
	end

	% generate points of valid data for integration
	idx1=2; 
	idx2=2; 
	j=1;
	for i=2:n_xbase
                if (xbase(i)<=max(x1(start_1),x2(start_2))) 
                  continue
                end
                if (xbase(i)>min(x1(end_1),x2(end_2))) 
                  continue
                end
		valid1=false;
		valid2=false;
		idx1=max(1,idx1-1);
		idx2=max(1,idx2-1);
		for k=idx1:n
			if (abs(x1(k)-(xbase(i-1)+xbase(i))/2.)<=h)
				valid1=true;
				idx1=max(1,k-1);
				break 
			end
		end
		for k=idx2:m
			if (abs(x2(k)-(xbase(i-1)+xbase(i))/2.)<=h)
				valid2=true;
				idx2=max(1,k-1);
				break
			end
		end

		if (valid1 && valid2)
			data_for_integration(j,1:2)=xbase(i-1:i);
			% linear interpolate
			for k=max(1,idx1-1):n
				if (x1(k)>=xbase(i-1))
					break;
				end
			end
			k=max(1,k-1);
			data_for_integration(j,3)=y1(k)+(y1(k+1)-y1(k))*(xbase(i-1)-x1(k))/(x1(k+1)-x1(k));
			for k=max(1,idx1-1):n
				if (x1(k)>=xbase(i))
					break
				end
			end
			k=max(1,k-1);
			data_for_integration(j,4)=(y1(k)+(y1(k+1)-y1(k))*(xbase(i)-x1(k))/(x1(k+1)-x1(k))-data_for_integration(j,3))/(data_for_integration(j,2)-data_for_integration(j,1));
			for k=max(1,idx2-1):m
				if (x2(k)>=xbase(i-1))
					break
				end
			end
			k=max(1,k-1);
			data_for_integration(j,5)=y2(k)+(y2(k+1)-y2(k))*(xbase(i-1)-x2(k))/(x2(k+1)-x2(k));
			for k=max(1,idx2-1):m
				if (x2(k)>=xbase(i))
					break
				end
			end
			k=max(1,k-1);
			data_for_integration(j,6)=(y2(k)+(y2(k+1)-y2(k))*(xbase(i)-x2(k))/(x2(k+1)-x2(k))-data_for_integration(j,5))/(data_for_integration(j,2)-data_for_integration(j,1));
			j=j+1;
		end
	end
	j=j-1;

	% calculate means by integration
	mean1=0.;
	mean2=0.;
	num=0.;
	for i=1:j
		delta_x=data_for_integration(i,2)-data_for_integration(i,1);
		a1=data_for_integration(i,3);
		b1=data_for_integration(i,4);
		mean1=mean1+delta_x*(a1+delta_x*b1/2.);
		a2=data_for_integration(i,5);
		b2=data_for_integration(i,6);
		mean2=mean2+delta_x*(a2+delta_x*b2/2.);
		num=num+delta_x;
	end
	mean1=mean1/num;
	mean2=mean2/num;

	% calculate correlation by integration
	data_for_integration(1:j,3)=data_for_integration(1:j,3)-mean1;
	data_for_integration(1:j,5)=data_for_integration(1:j,5)-mean2;
	num=0.;
	den1=0.;
	den2=0.;
	for i=1:j
		delta_x=data_for_integration(i,2)-data_for_integration(i,1);
		a1=data_for_integration(i,3);
		b1=data_for_integration(i,4);
		a2=data_for_integration(i,5);
		b2=data_for_integration(i,6);
		num=num+delta_x*(a1*a2+delta_x*(b1*a2/2+b2*a1/2+delta_x*b1*b2/3));
		den1=den1+delta_x*(a1^2+delta_x*(b1*a1+(delta_x*b1^2)/3));
		den2=den2+delta_x*(a2^2+delta_x*(b2*a2+(delta_x*b2^2)/3));
	end
	correlation=num/sqrt(den1*den2);

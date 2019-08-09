'''
Created on 24 May 2019

@author: Lenneke Jong (lenneke.jong@aad.gov.au)
'''
import numpy as np
import math


#default correlation method  
def correlation(y1,y2,x1,x2,**kwargs):
    """Correlate two data series. Calls the cslick correlation method by default. 
    
    (x1,y1), (x2,y2) are numpy arrays containing the two series to be correlated 
    hc: coefficient for tuning the data , defaults to 0.4 if not given as keyword argument
    """
    hc=0.4
    if "hc" in kwargs:
        hc=kwargs['hc']
    if "method" in kwargs:
        if kwargs['method']=="cslick":
            return correlation_cslick(y1,y2,x1,x2,hc)
        if kwargs['method']=="gaussian":
            return correlate_gaussian(y1, y2, x1, x2)
        print "Unknown method: {}. Defaulting to cslick".format(kwargs['method'])
        return correlation_cslick(y1,y2,x1,x2,hc)
    else:
        return correlation_cslick(y1,y2,x1,x2,hc)

def correlation_cslick(y1,y2,x1,x2,hc):
    """Calculate correlation for two series (x1,y1), (x2,y2) using the cslick method.
    """
    n=len(y1)
    m=len(y2)
    # find valid domain for both series, discard any data more than one point
    # outside range of other data series
    
    xmin=max(x1[0],x2[0])
    xmax=min(x1[-1],x2[-1])
    
    if xmin==x1[0]:
        start1=0
        start2=np.argwhere(x2>xmin)[0][0]-1
    else:
        start1=np.argwhere(x1>xmin)[0][0]-1
        start2=0
        
    if xmax==x1[-1]:
        end1=len(x1)-1
        end2=np.argwhere(x2<xmax)[-1][-1]+1
    else:
        end1=np.argwhere(x1<xmax)[-1][-1]+1
        end2=len(x2)-1
    
    x1base=x1[start1:end1+1]
    x2base=x2[start2:end2+1]
       
    dx1=np.sort(np.diff(x1base))
    median1=np.median(dx1)
    
    #  note that using the scipy function to calculate the interquartile range will
    #  give different results than the method used here.
    #  iqr1=scipy.stats.iqr(dx1)
    
    if (end1-start1)%2==1: # odd number of points
        i=int(math.floor((end1-start1)/4))
        if (end1-start1)%4==1:
            temp1=(dx1[i-1]+3*dx1[i])/4
            temp2=(3*dx1[3*i]+dx1[3*i+1])/4            
        else:
            temp1=(3*dx1[i]+dx1[i+1])/4
            temp2=(dx1[3*i+1]+3*dx1[3*i+2])/4
    else:
        temp1=(dx1[(end1-start1)/4-1]+dx1[(end1-start1)/4])/2.0
        temp2=(dx1[3*(end1-start1)/4]+dx1[3*(end1-start1)/4+1])/2.0
    iqr1=temp2-temp1
    
    dx2=np.sort(np.diff(x2base))
    median2=np.median(dx2)
    if (end2-start2)%2==1: # odd number of points
        i=int(math.floor((end2-start2)/4))
        if (end2-start2)%4==1:
            temp1=(dx2[i-1]+3*dx2[i])/4
            temp2=(3*dx2[3*i]+dx2[3*i+1])/4            
        else:
            temp1=(3*dx2[i]+dx2[i+1])/4
            temp2=(dx2[3*i+1]+3*dx2[3*i+2])/4
    else:
        temp1=(dx1[(end1-start1)/4-1]+dx1[(end1-start1)/4])/2.0
        temp2=(dx1[3*(end1-start1)/4]+dx1[3*(end1-start1)/4+1])/2.0
    
    iqr2=temp2-temp1
    
    #iqr2=scipy.stats.iqr(dx2)
    h=hc*max(median1,median2,iqr1,iqr2)
    xbase=np.hstack((x1base,x2base))
    

    xtemp=x1[(x1-h)>x1[start1]]
    xbase=np.hstack((xbase,xtemp-h))
    xtemp=x1[(x1+h)<x1[end1]]
    xbase=np.hstack((xbase,xtemp+h))
    
    xtemp=x2[(x2-h)>x2[start2]]
    xbase=np.hstack((xbase,xtemp-h))
    xtemp=x2[(x2+h)<x2[end2]]
    xbase=np.hstack((xbase,xtemp+h))
    
    xbase=np.unique(xbase)

    #generate points of valid data for integration    
    data_for_integration=np.zeros((3*(n+m),6),dtype=np.float64)
    
    idx1=0
    idx2=0
    j=0
    
    for i,xb in np.ndenumerate(xbase[1:]): #loop over all xbase points starting shifted over 1
        # find both series have a data point within valid interval
        i=i[0]+1
        if xb<xmin or xb>xmax:
            pass
            
        else:    
            valid1=False
            valid2=False
            xbm1=xbase[i-1] 
            
            tmp=np.argwhere(abs(x1-(xb+xbm1)/2.0)<h)
            if len(tmp)>0:
                valid1=True  
                idx1=tmp[0][0]
                
            tmp=np.argwhere(abs(x2-(xb+xbm1)/2.0)<h)
            if len(tmp)>0:
                valid2=True
                idx2=tmp[0][0]
                
            if valid1 and valid2:
                #print idx1,idx2
                data_for_integration[j,0]=xbm1
                data_for_integration[j,1]=xb
                # linear interpolate
                # search for index of first x1 greater than xbase[i-1]
                k = np.nonzero(x1[idx1:]>=xbm1)
                if k[0].size==0:
                    k=idx1
                else:
                    k=max(0,k[0][0])+idx1-1
                data_for_integration[j,2]=y1[k]+(y1[k+1]-y1[k])*(xbm1-x1[k])/(x1[k+1]-x1[k])
                # search for first x1 greater than xbase[i]
                k = np.nonzero(x1[idx1:]>=xb)
                if k[0].size==0:
                    k=idx1
                else:
                    k=max(0,k[0][0])+idx1-1
                data_for_integration[j,3]=(y1[k]+(y1[k+1]-y1[k])*(xb-x1[k])/(x1[k+1]-x1[k])-data_for_integration[j,2])/(data_for_integration[j,1]-data_for_integration[j,0])
    
                k=np.nonzero(x2[idx2:]>=xbm1)
                if k[0].size==0:
                    k=idx2
                else:
                    k=max(0,k[0][0])+idx2-1
                data_for_integration[j,4]=y2[k]+(y2[k+1]-y2[k])*(xbm1-x2[k])/(x2[k+1]-x2[k])
    
                k=np.nonzero(x2[idx2:]>=xb)
                if k[0].size==0:
                    k=idx2
                else:
                    k=max(0,k[0][0])+idx2-1
                data_for_integration[j,5]=(y2[k]+(y2[k+1]-y2[k])*(xb-x2[k])/(x2[k+1]-x2[k])-data_for_integration[j,4])/(data_for_integration[j,1]-data_for_integration[j,0])
                j=j+1
  
    # calculate means by integration
    mean1=0.0
    mean2=0.0
    num=0.0
    
    delta_x=data_for_integration[:j,1]-data_for_integration[:j,0]
    a1=data_for_integration[:j,2]
    b1=data_for_integration[:j,3]
   
    a2=data_for_integration[:j,4]
    b2=data_for_integration[:j,5]
    num=np.sum(delta_x)
    mean1=np.sum(delta_x*(a1+delta_x*b1/2.0))/num
    mean2=np.sum(delta_x*(a2+delta_x*b2/2.0))/num

    # calculate correlation by integration
    data_for_integration[:j,2]=data_for_integration[:j,2]-mean1
    data_for_integration[:j,4]=data_for_integration[:j,4]-mean2
    num=0.0
    den1=0.0
    den2=0.0

    num=np.sum(delta_x*(a1*a2+delta_x*(b1*a2/2.0+b2*a1/2.0+delta_x*b1*b2/3.0)))
    den1=np.sum(delta_x*(a1**2+delta_x*(b1*a1+(delta_x*b1**2)/3.0)))
    den2=np.sum(delta_x*(a2**2+delta_x*(b2*a2+(delta_x*b2**2)/3.0)))

    correlation=num/math.sqrt(den1*den2)
    return correlation

def correlate_gaussian(y1,y2,x1,x2):
    """ Calculate the correlation between two series (x1,y1), (x2,y2)
        using a Gaussian kernel method
        
    """
    y1mean=np.mean(y1)
    y2mean=np.mean(y2)
    
    
    x1mean=(x1[-1]-x1[0])/len(x1)
    x2mean=(x2[-1]-x2[0])/len(x2)

    delta_t=max(x1mean,x2mean)
    h=delta_t/4.0
    
    y1m,y2m=np.meshgrid(y1,y2)
    x1m,x2m=np.meshgrid(x1,x2)
    d=x2m-x1m
    b=np.exp(-np.power(d,2)/(2*h**2))/np.sqrt(2*np.pi*h)
    
    num=np.sum((y1m-y1mean)*(y2m-y2mean)*b)
    sdx=np.sum(b*np.power(y1m-y1mean,2))
    sdy=np.sum(b*np.power(y2m-y2mean,2)) 

    return num/np.sqrt(sdx*sdy)







    
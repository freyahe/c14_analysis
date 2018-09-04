# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:22:55 2017

@author: Freya
"""

import numpy as np
import matplotlib.pylab as plt
#import pylab as pl

from matplotlib.patches import Ellipse
from scipy.interpolate import interp1d
from scipy import misc

def Reservoir(UTh,intCal, decay,range_age):
    from scipy import optimize

    #Difference function
    def diff(UTh):
        return intCal(UTh) - decay(UTh)
#    print(diff(1.))
#    print(diff(49.))

    #Find intersection by finding root of difference
#    print(diff(0.))
#    print(diff(49000.))

#    UTh_intersection = optimize.brentq(diff,0.,50000.)
    UTh_intersection = optimize.brentq(diff,range_age[0],range_age[1])
    
    return UTh_intersection

def intCalUncertainty(slopeDecay, slopeIntCal, delta14CSE):
    return np.fabs(delta14CSE/(slopeDecay-slopeIntCal))



def calculate_Delta14C(C14_age,UTh_age):
    lambda_libby=8033 #[a]
    lambda_new=8266 #[a]

    deltaC14=(np.exp(-C14_age/lambda_libby)*np.exp(UTh_age/lambda_new)-1.)*1000.
    return deltaC14
    
def calculations_14c(Input_data,IntCal,n_random = 100, range_age=[0,50e3]):

    #Daten importieren, 2sigma-Fehler
    
    Data_UTh = Input_data['U/Th Age BP']
    Data_UTh2SD =  Input_data['U/Th Age 2SE']
    Data_14C =  Input_data['14C Age']
    Data_14C2SD=2*Input_data['14C Age SE']

    n_values=len(Data_UTh)
    
    CalAge = IntCal['CAL BP']
    deltaC14theo = IntCal['Delta 14C']
    sigma_deltaC14theo = IntCal['D14C SD']
    
    lambda_libby=8033 #[a]
    lambda_new=8266 #[a]
    
    
    delta14C_value=calculate_Delta14C(Data_14C,Data_UTh)

    
    #Delta14C ausrechnen, um Daten zusätzlich zu den Ellipsen zu plotten

    # seed for random variations of values:
    np.random.seed(3457625575)
    
    # create random variations of U-Th Age  
    x = np.linspace(range_age[0],range_age[1],1000)
    n_range_age_random=len(x)
    y_list = np.zeros(shape=(n_values,n_range_age_random))
    
    #lista for intersection
    intersectionx_list = np.zeros(n_values)
    intersectiony_list = np.zeros(n_values)
    Uncertainty_IntCal_list = np.zeros(n_values)
    DeltaR_list=np.zeros(n_random)
    DeltaDelta_list=np.zeros(n_random)
       
    #Arrazs to store output:

    Delta14C_mean=np.zeros(n_values)
    Delta14C_std=np.zeros(n_values)
    Delta14C_sigma_x_strich=np.zeros(n_values)
    Delta14C_sigma_y_strich=np.zeros(n_values)
    Delta14C_theta=np.zeros(n_values)
    Delta14C_theta2=np.zeros(n_values)
    
    DeltaR_mean=np.zeros(n_values)
    DeltaR_std=np.zeros(n_values)
    DeltaR_sigma_x_strich=np.zeros(n_values)
    DeltaR_sigma_y_strich=np.zeros(n_values)
    DeltaR_theta=np.zeros(n_values)
    DeltaR_theta2=np.zeros(n_values)

    DeltaDelta_mean=np.zeros(n_values)
    DeltaDelta_std=np.zeros(n_values)
    DeltaDelta_sigma_x_strich=np.zeros(n_values)
    DeltaDelta_sigma_y_strich=np.zeros(n_values)
    DeltaDelta_theta=np.zeros(n_values)
    DeltaDelta_theta2=np.zeros(n_values)

    xUTh_mean=np.zeros(n_values)
    xUTh_std=np.zeros(n_values)

    #Loop over data points:
    for index_UTh, value_UTh in enumerate(Data_UTh):
        print('data_point nr: ',index_UTh)

        for x_index, x_i in enumerate(x): 
            
            y_list[index_UTh,x_index]=((np.exp(x_i/lambda_new)*np.exp(-(Data_UTh[index_UTh]/lambda_libby)))-1)*1000.
        
        #interpolate IntCal between data points          
        intCal = interp1d(CalAge, deltaC14theo) 
        intCalSE = interp1d(CalAge, sigma_deltaC14theo)
        
        #decay curve through coral data point
        def decay(x):
            return (np.exp(x/lambda_new)*np.exp(-(Data_14C[index_UTh]/lambda_libby))-1)*1000.
        
        def slopeDecay(x):
            return np.exp(x/lambda_new)*np.exp(-(Data_14C[index_UTh])/lambda_libby)/lambda_new
        
        def slopeIntCal(x):
            return misc.derivative(intCal, x, dx=0.01)
        
        #use function reservoir_age which calculates the intersection
        intersection_x = Reservoir(Data_UTh[index_UTh],intCal, decay,range_age)
        intersection_y = decay(intersection_x)         
          
        intersectionx_list[index_UTh]=intersection_x
        intersectiony_list[index_UTh]=intersection_y
        
        #calculate uncertainty resulting from IntCal SE
        uncertaintyIntCal=intCalUncertainty(slopeDecay(intersection_x),slopeIntCal(intersection_x), intCalSE(intersection_x))                  
        
        #uncertainty_IntCal in a list
        Uncertainty_IntCal_list[index_UTh]=uncertaintyIntCal                 
                          
        xUTh = Data_UTh[index_UTh] + Data_UTh2SD[index_UTh] * np.random.randn(n_random) #make point cloud
        xC14 = Data_14C[index_UTh] + Data_14C2SD[index_UTh] * np.random.randn(n_random)

        #Filter our values below range:
    


        deltaC14=(np.exp(-xC14/lambda_libby)*np.exp(xUTh/lambda_new)-1.)*1000

    
        # for each point of cloud use the decay curve and calculate intersection with IntCal. write that into list intersction_point_list and calculate SD over those intersects
        for index_point in range(n_random):
            # Only perfrom calculations when xUTh is inside range_age, otherwise discard point from cloud
            if (xUTh[index_point]>range_age[0] and xUTh[index_point]<range_age[1]):
                def decaypoint(x):
                    return ((np.exp(x/lambda_new)*np.exp(-(xC14[index_point]/lambda_libby)))-1)*1000.
                intersection_point = Reservoir(xUTh[index_point],intCal, decaypoint,range_age) #preliminary intersection point wihtout uncertainty of IntCal
                intersection_point = Reservoir(xUTh[index_point],lambda x: intCal(x)+np.random.randn(1)*intCalSE(intersection_point), decaypoint,range_age) #intersection point wiht uncertainty of IntCal
                    
                DeltaR_list[index_point]=(intersection_point-xUTh[index_point])
                
                #DeltaDelta C Calculation
                DeltaDelta_list[index_point] = deltaC14[index_point]-intCal(xUTh[index_point])+np.random.randn(1)*intCalSE(xUTh[index_point])
        
        #Calculate parameters for plotting elypses (plots in ka!)

        (Delta14C_mean[index_UTh],Delta14C_std[index_UTh],Delta14C_sigma_x_strich[index_UTh],Delta14C_sigma_y_strich[index_UTh],Delta14C_theta[index_UTh],Delta14C_theta2[index_UTh],xUTh_mean[index_UTh],xUTh_std[index_UTh]) =  calculate_ellipses(xUTh/1000,deltaC14) 
        (DeltaR_mean[index_UTh],DeltaR_std[index_UTh],DeltaR_sigma_x_strich[index_UTh],DeltaR_sigma_y_strich[index_UTh],DeltaR_theta[index_UTh],DeltaR_theta2[index_UTh],xUTh_mean[index_UTh],xUTh_std[index_UTh]) =  calculate_ellipses(xUTh/1000,DeltaR_list/1000)    
        (DeltaDelta_mean[index_UTh],DeltaDelta_std[index_UTh],DeltaDelta_sigma_x_strich[index_UTh],DeltaDelta_sigma_y_strich[index_UTh],DeltaDelta_theta[index_UTh],DeltaDelta_theta2[index_UTh],xUTh_mean[index_UTh],xUTh_std[index_UTh]) =  calculate_ellipses(xUTh/1000,DeltaDelta_list)   


    # Write data into output DataFrame:
    Output_data=Input_data        
    Output_data['Delta14C']=delta14C_value
    Output_data['Delta14C_mean']=Delta14C_mean
    Output_data['Delta14C_std']=Delta14C_std
    Output_data['Delta14C_sigma_x_strich']=Delta14C_sigma_x_strich
    Output_data['Delta14C_sigma_y_strich']=Delta14C_sigma_y_strich
    Output_data['Delta14C_theta']=Delta14C_theta
    Output_data['Delta14C_theta2']=Delta14C_theta2
     
    Output_data['DeltaR_mean']=DeltaR_mean
    Output_data['DeltaR_std']=DeltaR_std
    Output_data['DeltaR_sigma_x_strich']=DeltaR_sigma_x_strich
    Output_data['DeltaR_sigma_y_strich']=DeltaR_sigma_y_strich
    Output_data['DeltaR_theta']=DeltaR_theta
    Output_data['DeltaR_theta2']=DeltaR_theta2
    
    Output_data['DeltaDelta_mean']=DeltaDelta_mean
    Output_data['DeltaDelta_std']=DeltaDelta_std
    Output_data['DeltaDelta_sigma_x_strich']=DeltaDelta_sigma_x_strich
    Output_data['DeltaDelta_sigma_y_strich']=DeltaDelta_sigma_y_strich
    Output_data['DeltaDelta_theta']=DeltaDelta_theta
    Output_data['DeltaDelta_theta2']=DeltaDelta_theta2

    Output_data['xUTh_mean']=xUTh_mean
    Output_data['xUTh_std']=xUTh_std


    return Output_data
     
        
def calculate_ellipses(xUTh,values):    
#Kovarianzellipsen für Delta berechnen
    xx = np.cov(xUTh , values )
    cov_delta = xx[0,1]

    values_mean = np.mean(values)
    values_std = np.sqrt(xx[1,1])
    xUTh_mean = np.mean(xUTh)
    xUTh_std = np.sqrt(xx[0,0])
    
    sigma_x_strich = np.sqrt( (values_std**2 + xUTh_std**2)/2 +\
        np.sqrt((values_std**2 - xUTh_std**2)**2/4 + cov_delta**2) )
    sigma_y_strich = np.sqrt( (values_std**2 + xUTh_std**2)/2 -\
        np.sqrt((values_std**2 - xUTh_std**2)**2/4 + cov_delta**2) )
    
    theta = 0.5 * np.arctan( 2*cov_delta / (xUTh_std**2 - values_std**2) )
    theta2 = theta*(180/np.pi)+270
    return values_mean,values_std,sigma_x_strich,sigma_y_strich,theta,theta2,xUTh_mean,xUTh_std
    

def plot_ellipses_Delta14C(Input,axes=plt.gca(),**kwargs):
    for index,row in Input.iterrows():
        plot_handle=plot_ellipses(xUTh_mean=row['xUTh_mean'],
                                  value_mean=row['Delta14C_mean'],
                                  sigma_x_strich=row['Delta14C_sigma_x_strich'],
                                  sigma_y_strich=row['Delta14C_sigma_y_strich'],
                                  theta2=row['Delta14C_theta2'],
                                  axes=axes,
                                  **kwargs)
    return plot_handle

def plot_ellipses_DeltaR(Input,axes=plt.gca(),**kwargs):
    for index,row in Input.iterrows():
        plot_handle=plot_ellipses(xUTh_mean=row['xUTh_mean'],
                                  value_mean=row['DeltaR_mean'],
                                  sigma_x_strich=row['DeltaR_sigma_x_strich'],
                                  sigma_y_strich=row['DeltaR_sigma_y_strich'],
                                  theta2=row['DeltaR_theta2'],
                                  axes=axes,
                                  **kwargs)
    return plot_handle

def plot_ellipses_DeltaDelta(Input,axes=plt.gca(),**kwargs):
    for index,row in Input.iterrows():
        plot_handle=plot_ellipses(xUTh_mean=row['xUTh_mean'],
                                  value_mean=row['DeltaDelta_mean'],
                                  sigma_x_strich=row['DeltaDelta_sigma_x_strich'],
                                  sigma_y_strich=row['DeltaDelta_sigma_y_strich'],
                                  theta2=row['DeltaDelta_theta2'],
                                  axes=axes,
                                  **kwargs)
    return plot_handle

def plot_ellipses(xUTh_mean,value_mean,sigma_x_strich,sigma_y_strich,theta2,color='k',alpha=0.15,lw=0.5,axes=plt.gca()):
#Plot sinle ellipsis:
    e = Ellipse((xUTh_mean, value_mean),\
        2*sigma_x_strich, 2*sigma_y_strich, theta2,
        color=color, alpha=alpha, lw=lw)
    axes.add_artist(e)
    e.set_edgecolor(color)
#    plot_marker=axes.plot(xUTh_mean, value_mean, color=color, marker='.')
    return e

# PLot intcal
def plot_int_cal(IntCal,axes=plt.gca(),color='k'):
    CalAge = IntCal['CAL BP']/1000
    deltaC14theo = IntCal['Delta 14C']
    sigma_deltaC14theo = IntCal['D14C SD']

    #1sigma Bereich der Intcal-Kurve berechnen
    deltaC14theoplus = deltaC14theo + sigma_deltaC14theo
    deltaC14theominus = deltaC14theo - sigma_deltaC14theo
    plot_intcal=axes.plot(CalAge,deltaC14theo, color=color, linestyle = '-')
    axes.plot(CalAge, deltaC14theoplus, color=color, linestyle = ':')
    axes.plot(CalAge, deltaC14theominus, color=color, linestyle = ':')
    return plot_intcal

 


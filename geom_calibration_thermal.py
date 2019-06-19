# Python code to plot the data from HELIX magnet test. Reads in .csv file downloaded form the HELIX wiki and plots versus time the quantities.  
# New code is to spline the levels on the magnet during the test to then take the derivative and plot as a function of time. 
# Update including taking the ratio of the flow per time to the derivative of the near level sensor to get a cross-sectional area of the dewar as a function of the level sensor. 
# Author Keith McBride 10/16/18

import matplotlib.pyplot as plt
import numpy
from scipy import interpolate
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular' #might be unneccessary



def load_flows_zero_to_near_lvl_sensor(seq):
# load lvlnear data
    datanear=numpy.genfromtxt('lvlSensorNear.csv', dtype=float, delimiter=',', names=True)
    timenear = datanear['time']
    lvlnear=datanear['StackSideLevelcm']
#load flow data
    data=numpy.genfromtxt('WhisperData_Cleaned.csv', dtype=float, delimiter=',', names=True)
    timeflow = data['time_s']
    pressure = data['pressure_PSI']
    temp = data['temperatureC']
    flow = data['volume_flow_LPM']
    mass = data['mass_flow']
    #calculate the mass flow from the volume flow
      #convert temp to kelvin
    kelvin=numpy.empty(len(temp))
    kelvin.fill(273.15)
    temp=temp+kelvin
    masscalc = 4.0/(1.20675)*numpy.true_divide(pressure,temp) #molar mass * Pressure/(R[PSIliters/(kelvinmol)]*Temp) has units of grams
    masscalc=numpy.multiply(masscalc,flow)
# convert unix time to seconds after tests started for all arrays of time, each array of time has different lengths since data was poorly recorded.
   # make array of same length as each time array
    a=numpy.empty(len(timeflow))
    b=numpy.empty(len(timenear))
   # fill each array of time with entry of earliest unix time (timenear[0])
    a.fill(timenear[0])
    b.fill(timenear[0])
  # subtract each start time array (a,b,c) from the corresponding time array.
    timeflow=timeflow-a
    timenear=timenear-b
# convert to minutes after test started
    timeflow= numpy.true_divide(timeflow,60.0)
    timenear= numpy.true_divide(timenear,60.0)
    data_list=numpy.array([timeflow,pressure,temp,flow,mass,masscalc])
    data_near=numpy.array([timenear,lvlnear])
    return data_list, data_near

def plot_mass_flow_and_level_der(seq):
    print('starting')
    der=spline_levels_near(1) # this loads the 2-d array of [0,:]= time (zeroed to near sensor start and in minutes) [1,:]= derivative of spline of level near sensor (cmpm)
    print('laoded derivative')
    flows=load_flows_zero_to_near_lvl_sensor(1)  # 5-d array with [0,:]=time (not zeroed to flow start and in minutes) [1,:]= pressure (PSI) [2,:]= temperature (Kelvin) [3,:]=volume flow (lpm) [4,:]= mass flow (? unknown units) [5,:]=mass flow calc from ideal gas law(kg per minute)
    print('loaded flows')
    datanear=load_levels_near(1) #[0,:]= time near from data taken (mins) [1,:]=lvl near cm reading (cm)
    print('abotu to make figure')
    # set flows time to zeroed 
    fig=plt.figure(figsize=(10, 8), dpi=800,)
    plt.title('Magnet Thermal Test, splined near sensor level derivative')
    ax_flow=plt.subplot(311)
    ax_flow.scatter(flows[0,:],flows[5,:],s=4)
    ax_flow.set_ylabel('Calculated massflow (kg per minute)')
    ax_flow.set_xlim([0,18000])
    ax_near=plt.subplot(312)
    ax_near.scatter(der[0,:],der[1,:],s=4)
    ax_near.set_xlim([0,18000])
    ax_near.set_ylabel('Level Near derivative(cm per min)')
    ax_lvl=plt.subplot(313)
    ax_lvl.scatter(datanear[0,:],datanear[1,:],s=4)
    ax_lvl.set_xlim([0,18000])
    ax_lvl.set_ylabel('Level Near (cm)')
    ax_lvl.set_xlabel('Time (mins)')
    sum=0.0
    i=0
    while i<len(flows[0,:]):
      if i==len(flows[0,:])-1: diff=0
      else: diff=flows[0,i+1]-flows[0,i]
      sum+=diff
      i+=1
    print (sum/len(flows[0,:]))
#    plt.scatter(der[0,:], der[1,:], c='b', marker='s', label='Near Sensor derivative')
#    plt.scatter(der[0,:], der[1,:], c='r', marker='s', label='Near Sensor derivative')
    fig.savefig('magnet_mass_flow_calc_and_lvl_spline_der_vs_time.png')
    
######################for 11/28/18
# functions that will help find sections

# moving average for finding flows around level sensor values
def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def integrate_time_range(lvltime,flowtime,flow):
    flow_at_lvl_time=[]
#   integrate_timesA=[] #time units to integrate
    l=0
    find_times=[]
    while l<len(lvltime):
 #      print('lvl time before checking is ' , lvl_time_spliceA[l])
 #      print(l)
       #match up flows
       m=0
       if l==len(lvltime)-1:
          delta_t=lvltime[l]-lvltime[l-1]
       if l!=len(lvltime)-1:
          delta_t=lvltime[l+1]-lvltime[l]
       while m<len(flowtime):
          if m!=0 or m!=len(flowtime)-1: #not the end points of the array
             if flowtime[m]>lvltime[l]-(delta_t/2.0) and flowtime[m]<lvltime[l]+(delta_t/2.0):
#            if (time_spliceA[m]<lvl_time_spliceA[l] and time_spliceA[m+1]>lvl_time_spliceA[l]) or time_spliceA[m]==lvl_time_spliceA[l]:
                find_times.append(m)
#               print('time flow is' , time_spliceA[m])
#               print('lvl time is ' , lvl_time_spliceA[l])
#integrate from one time of lvlnear to the next to get an average flow. 
#               average=(flows_spliceA[m-1]+flows_spliceA[m+1]+flows_spliceA[m])/3.0
          m+=1
       start=find_times[0]
       stop=find_times[-1]
       #print('start is ' , start)
       average=numpy.trapz(flow[start:stop], x=flowtime[start:stop])/(delta_t)
       #print (average)
       flow_at_lvl_time.append(average)
       find_times=[]
       l+=1
    return flow_at_lvl_time

def find_avg_errors(lvlA,lvlB,lvlC,ratA,ratB,ratC,num_bins):
    avg=[]
    errors=[]
    bins=[]
#stitch all the arrays together
    lvltots=[]
    rattots=[]
    l=0
    while l<len(lvlA):
       lvltots.append(lvlA[l])
       rattots.append(ratA[l])
       l+=1
    l=0
    while l<len(lvlB):
       lvltots.append(lvlB[l])
       rattots.append(ratB[l])
       l+=1
    l=0
    while l<len(lvlC):
       lvltots.append(lvlC[l])
       rattots.append(ratC[l])
       l+=1
#find the error and average in each bin
    #start with min and max of levels
    lvlmin=numpy.amin(lvltots)
    lvlmax=numpy.amax(lvltots)
    binwidth=(lvlmax-lvlmin)/num_bins
    print("{} is the binwidth in cm".format(binwidth))
    #outer loop for something about which bin
    lvliter=lvlmin-(binwidth/2.0)
    #array for storing the data points in the bin
    lvlbin_data=[]
    ratbin_data=[]
    while lvliter<lvlmax+(binwidth/2.0):
       #inner loop for checking all the elements
       l=0
       while l<len(lvltots):
          if lvltots[l]>=lvliter and lvltots[l]<=lvliter+(binwidth):
             lvlbin_data.append(lvltots[l]) # put that value in the array to be used for calculations
             ratbin_data.append(rattots[l])
             #find mean, reset values, declare values needed zero before this loop, and then calc standard deviation for those values somehow also
          l+=1
       print(len(ratbin_data))
#       if len(ratbin_data)==0:
#          ratbin_data.append(0)
#       avg.append(numpy.mean(ratbin_data))
#       errors.append(numpy.std(ratbin_data))
#       bins.append(lvliter)
       if len(ratbin_data)!=0:
          avg.append(numpy.mean(ratbin_data))
          errors.append(numpy.std(ratbin_data))
          bins.append(lvliter)
       lvliter+=binwidth
       lvlbin_data=[] # clear the bin data array
       ratbin_data=[]
    return bins, avg, errors

def plot_lvl_near_sensor_vs_time_section(seq):
###load everything###
    flows, lvlnear=load_flows_zero_to_near_lvl_sensor(1)
###declare all necessary arrays###
    lvl_spliceA=[]
    lvl_time_spliceA=[]
    lvl_spliceB=[]
    lvl_time_spliceB=[]
    lvl_spliceC=[]
    lvl_time_spliceC=[]
    flows_spliceA=[]
    time_spliceA=[]
    flows_spliceB=[]
    time_spliceB=[]
    flows_spliceC=[]
    time_spliceC=[]
###splice flows into sections###
    k=0
    while k<len(flows[0,:]):
#        if ((flows[0,k]>5138 and flows[0,k]<5470) or (flows[0,k]>5490 and flows[0,k]<6826.13)) and flows[3,k]>20 and flows[3,k]<28: # amended this with getting rid of the data points in the window of 5470 and 5490
        if (flows[0,k]>5138 and flows[0,k]<5470 and flows[3,k]>20 and flows[3,k]<28) or (flows[0,k]>5490 and flows[0,k]<6826.13 and flows[3,k]>20 and flows[3,k]<28): # amended this with getting rid of the data points in the window of 5470 and 5490
           flows_spliceA.append(flows[5,k])
           time_spliceA.append(flows[0,k])
           #how to splice simultaneously the lvl sensor? its a smaller array so just go through the elements?
        if (flows[0,k]>6900 and flows[0,k]<7040) or (flows[0,k]<8382.8 and flows[0,k]>7060): # amended to avoid window of [7040,7060]
#        if flows[0,k]>6900 and flows[0,k]<8382.8:
           flows_spliceB.append(flows[5,k])
           time_spliceB.append(flows[0,k])
        if flows[0,k]>10000 and flows[0,k]<16500:
           flows_spliceC.append(flows[5,k])
           time_spliceC.append(flows[0,k])
        k+=1
###splice levels into sections###
    j=0
    while j<len(lvlnear[0,:]):
#        if lvlnear[0,j] >time_spliceA[0] and lvlnear[0,j]<time_spliceA[-1]:
        if (lvlnear[0,j] >time_spliceA[0] and lvlnear[0,j] <5470) or (lvlnear[0,j]<time_spliceA[-1] and lvlnear[0,j]>5490): #amended to avoid the window of 5470 and 5490
           lvl_spliceA.append(lvlnear[1,j])
           lvl_time_spliceA.append(lvlnear[0,j])
#        if lvlnear[0,j] >time_spliceB[0] and lvlnear[0,j]<time_spliceB[-1]:
        if (lvlnear[0,j]>time_spliceB[0] and lvlnear[0,j]<7040) or (lvlnear[0,j]>7060 and lvlnear[0,j]<time_spliceB[-1]): #amended to avoid the window of 7040 and 7060
           lvl_spliceB.append(lvlnear[1,j])
           lvl_time_spliceB.append(lvlnear[0,j])
        if lvlnear[0,j] >time_spliceC[0] and lvlnear[0,j]<time_spliceC[-1]:
           lvl_spliceC.append(lvlnear[1,j])
           lvl_time_spliceC.append(lvlnear[0,j])
        j+=1
#    print(time_spliceC[-1])
###spline sections###
    #sectionA spline
    tckA, fpA, ierA, msgA= interpolate.splrep(lvl_time_spliceA, lvl_spliceA, s=0, full_output=True)
    splineA_times_more=numpy.arange(lvl_time_spliceA[0],lvl_time_spliceA[-1],0.1)
    splineA=interpolate.splev(lvl_time_spliceA, tckA, der=0)
    splineA_more=interpolate.splev(splineA_times_more, tckA, der=0)
    splineA_der=interpolate.splev(lvl_time_spliceA, tckA, der=1)# derivative only at lvl times
    splineA_der_more=interpolate.splev(splineA_times_more, tckA, der=1)# derivative at spline times
    splineA_der_flowtimes=interpolate.splev(time_spliceA, tckA, der=1)# derivative at flow
    #sectionB spline
    tckB, fpB, ierB, msgB= interpolate.splrep(lvl_time_spliceB, lvl_spliceB, s=0, full_output=True)
    splineB_times_more=numpy.arange(lvl_time_spliceB[0],lvl_time_spliceB[-1],0.1)
    splineB=interpolate.splev(lvl_time_spliceB, tckB, der=0)
    splineB_more=interpolate.splev(splineB_times_more, tckB, der=0)
    splineB_der=interpolate.splev(lvl_time_spliceB, tckB, der=1)# derivative only at lvl times
    splineB_der_more=interpolate.splev(splineB_times_more, tckB, der=1)# derivative only at lvl times
    splineB_der_flowtimes=interpolate.splev(time_spliceB, tckB, der=1)# derivative at flow
    #sectionC spline
    tckC, fpC, ierC, msgC= interpolate.splrep(lvl_time_spliceC, lvl_spliceC, s=0, full_output=True)
    splineC_times_more=numpy.arange(lvl_time_spliceC[0],lvl_time_spliceC[-1],0.1)
    splineC=interpolate.splev(lvl_time_spliceC, tckC, der=0)
    splineC_more=interpolate.splev(splineC_times_more, tckC, der=0)
    splineC_der=interpolate.splev(lvl_time_spliceC, tckC, der=1) # derivative only at lvl times
    splineC_der_more=interpolate.splev(splineC_times_more, tckC, der=1)# derivative only at lvl times
    splineC_der_flowtimes=interpolate.splev(time_spliceC, tckC, der=1)# derivative at flow
###############plotter section############
###plot all### 
    fig=plt.figure(figsize=(10, 8), dpi=800)
    if seq==1:
       fig.subplots_adjust(hspace=0)
       #need subplots, subplots(sharex=True)
       ax_flow=plt.subplot(211)
       plt.setp(ax_flow.get_xticklabels(), visible=False)
       ax_flow.scatter(flows[0,:],flows[3,:],s=3)
       # add in the section lines
       xposition=[5138,6826.13,6900,8382.8,10000,16500]
       labels=['section A', 'section B', 'section C']
       label_pos=[5300,7000,12000]
       i=0
       for xc in xposition:
          ax_flow.axvline(x=xc, color='k', linestyle='--') # this is for the lines to mark at what time notable things in the test occured.
       for xl in labels:
          ax_flow.text(label_pos[i], 40, labels[i], rotation=0, fontsize=8)
          i+=1
       ax_flow.set_ylabel('Flow (liters per minute)')
       ax_flow.set_xlim([0,18000])
       ax_near=plt.subplot(212,sharex=ax_flow)
       ax_near.scatter(lvlnear[0,:],lvlnear[1,:],s=3)
       for xi in xposition:
          ax_near.axvline(x=xi, color='k', linestyle='--') # this is for the lines to mark at what time notable things in the test occured.
       ax_near.set_ylabel('Level Sensor (cm)')
       ax_near.set_xlim([0,18000])
       ax_near.set_xlabel('Time (mins)')
       fig.savefig('ALL_flow_and_lvl_vs_time.png')
       plt.clf()
#################################################SECTIONS#############################################################
###plot sectionA###
    # plot flow
    if seq==2:
       plt.scatter(time_spliceA, flows_spliceA, c='b', s=3)
       plt.title('Flow vs Time, section A')
       plt.ylabel('Flow (liters per min)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionA_flow_vs_time.eps')
       plt.clf()
       # plot level data and spline
   #    plt.scatter(lvl_time_spliceA, splineA, c='r', s=4, marker='s', label='Spline')
       plt.scatter(splineA_times_more, splineA_more, c='r', s=4, marker='s', label='Spline')
       plt.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=3, label='Data')
       plt.legend(loc='upper right')
       plt.title('Level vs Time, section A')
       plt.ylabel('Level (cm)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionA_lvl_and_spline_vs_time.eps')
   #    fig.savefig('sectionA_lvl_vs_time.eps')
       plt.clf()
      #plot flow and level with subplots?
   #https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/shared_axis_demo.html#sphx-glr-gallery-subplots-axes-and-figures-shared-axis-demo-py
   #    ax_spline=plt.subplot(211,sharex
    if seq==1:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Aflow=plt.subplot(211)
       ax_Aflow.set_title('Section A')
       plt.setp(ax_Aflow.get_xticklabels(), visible=False)              
       ax_Aflow.scatter(time_spliceA, flows_spliceA, c='b', s=3)
       ax_Aflow.set_ylabel('Flow (liters per min)')
       ax_Aspline=plt.subplot(212)
       ax_Aspline.scatter(splineA_times_more, splineA_der_more, c='r', s=4, marker='s', label='Spline')
       ax_Aspline.scatter(lvl_time_spliceA, splineA_der, c='b', s=4, marker='s', label='Spline at data points')
#       ax_Aspline.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=3, label='Data')
       ax_Aspline.legend(loc='upper right')
       ax_Aspline.set_ylabel('Level derivative (cm per min)')
       ax_Aspline.set_xlabel('Time (mins)')
       fig.savefig('sectionA_lvl_der_spline_vs_time_subplots_both.eps')
#       fig.savefig('sectionA_lvl_der_spline_vs_time_subplots_more.eps')
       plt.clf()
    if seq==3:
    #plot flows and derivative at all flow points
       fig.subplots_adjust(hspace=0)
       ax_Aflow=plt.subplot(211)
       ax_Aflow.set_title('Section A')
       plt.setp(ax_Aflow.get_xticklabels(), visible=False)              
       ax_Aflow.scatter(time_spliceA, flows_spliceA, c='b', s=3)
       ax_Aflow.set_ylabel('Flow (liters per min)')
       ax_Aspline=plt.subplot(212)
       ax_Aspline.scatter(time_spliceA, splineA_der_flowtimes, c='r', s=4, marker='s', label='Spline')
#       ax_Aspline.scatter(lvl_time_spliceA, splineA_der, c='r', s=4, marker='s', label='Spline')
#       ax_Aspline.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=3, label='Data')
 #      ax_Aspline.legend(loc='upper right')
       ax_Aspline.set_ylabel('Level derivative (cm per min)')
       ax_Aspline.set_xlabel('Time (mins)')
       fig.savefig('sectionA_lvl_der_spline_vs_time_subplots_flowtimes.eps')
       plt.clf()
#    if seq==4:
       #make plot of ratios, which involves calculating the flows at the lvltimes?

###plot sectionB###
    if seq==2:
       plt.scatter(time_spliceB, flows_spliceB, c='b', s=3)
       plt.title('Flow vs Time, section B')
       plt.ylabel('Flow (liters per min)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionB_flow_vs_time.eps')
       plt.clf()
#    plt.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3)
       plt.scatter(splineB_times_more, splineB_more, c='r', s=4, marker='s', label='Spline')
       plt.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3, label='Data')
       plt.legend(loc='upper right')
       plt.title('Level vs Time, section B')
       plt.ylabel('Level (cm)')
       plt.xlabel('Time (mins)')
#       fig.savefig('sectionB_lvl_and_spline_vs_time.eps')
#    fig.savefig('sectionB_lvl_vs_time.eps')
       plt.clf()
    if seq==1:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Bflow=plt.subplot(211)
       ax_Bflow.set_title('Section B')
       plt.setp(ax_Bflow.get_xticklabels(), visible=False)              
       ax_Bflow.scatter(time_spliceB, flows_spliceB, c='b', s=3)
       ax_Bflow.set_ylabel('Flow (liters per min)')
       ax_Bspline=plt.subplot(212)
       ax_Bspline.scatter(splineB_times_more, splineB_der_more, c='r', s=4, marker='s', label='Spline')
       ax_Bspline.scatter(lvl_time_spliceB, splineB_der, c='b', s=4, marker='s', label='Spline at data points')
#       ax_Bspline.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3, label='Data')
       ax_Bspline.legend(loc='upper right')
       ax_Bspline.set_ylabel('Level derivative (cm per min)')
       ax_Bspline.set_xlabel('Time (mins)')
       fig.savefig('sectionB_lvl_der_spline_vs_time_subplots_both.eps')
#       fig.savefig('sectionB_lvl_der_spline_vs_time_subplots_more.eps')
       plt.clf()
    if seq==3:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Bflow=plt.subplot(211)
       ax_Bflow.set_title('Section B')
       plt.setp(ax_Bflow.get_xticklabels(), visible=False)              
       ax_Bflow.scatter(time_spliceB, flows_spliceB, c='b', s=3)
       ax_Bflow.set_ylabel('Flow (liters per min)')
       ax_Bspline=plt.subplot(212)
       ax_Bspline.scatter(time_spliceB, splineB_der_flowtimes, c='r', s=4, marker='s', label='Spline')
#       ax_Bspline.scatter(lvl_time_spliceB, splineB_der, c='r', s=4, marker='s', label='Spline')
#       ax_Bspline.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3, label='Data')
#       ax_Bspline.legend(loc='upper right')
       ax_Bspline.set_ylabel('Level derivative (cm per min)')
       ax_Bspline.set_xlabel('Time (mins)')
       fig.savefig('sectionB_lvl_der_spline_vs_time_subplots_flowtimes.eps')
       plt.clf()

###plot sectionC###
    if seq==2:
       plt.scatter(time_spliceC, flows_spliceC, c='b', s=3)
       plt.title('Flow vs Time, section C')
       plt.ylabel('Flow (liters per min)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionC_flow_vs_time.eps')
       plt.clf()
#       plt.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3)
       plt.scatter(splineC_times_more, splineC_more, c='r', s=4, marker='s', label='Spline')
       plt.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3, label='Data')
       plt.legend(loc='upper right')
       plt.title('Level vs Time, section C')
       plt.ylabel('Level (cm)')
       plt.xlabel('Time (mins)')
       fig.savefig('sectionC_lvl_and_spline_vs_time.eps')
       plt.clf()
#       fig.savefig('sectionC_lvl_vs_time.eps')
    if seq==1:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Cflow=plt.subplot(211)
       ax_Cflow.set_title('Section C')
       plt.setp(ax_Cflow.get_xticklabels(), visible=False)              
       ax_Cflow.scatter(time_spliceC, flows_spliceC, c='b', s=3)
       ax_Cflow.set_ylabel('Flow (liters per min)')
       ax_Cspline=plt.subplot(212)
       ax_Cspline.scatter(splineC_times_more, splineC_der_more, c='r', s=4, marker='s', label='Spline')
       ax_Cspline.scatter(lvl_time_spliceC, splineC_der, c='b', s=4, marker='s', label='Spline at data points')
#       ax_Cspline.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3, label='Data')
       ax_Cspline.legend(loc='upper right')
       ax_Cspline.set_ylabel('Level derivative (cm per min)')
       ax_Cspline.set_xlabel('Time (mins)')
       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots_both.eps')
#       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots_more.eps')
#       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots.eps')
       plt.clf()
    if seq==3:
# make subplots instead 
       fig.subplots_adjust(hspace=0)
       ax_Cflow=plt.subplot(211)
       ax_Cflow.set_title('Section C')
       plt.setp(ax_Cflow.get_xticklabels(), visible=False)              
       ax_Cflow.scatter(time_spliceC, flows_spliceC, c='b', s=3)
       ax_Cflow.set_ylabel('Flow (liters per min)')
       ax_Cspline=plt.subplot(212)
       ax_Cspline.scatter(time_spliceC, splineC_der_flowtimes, c='r', s=4, marker='s', label='Spline')
#       ax_Cspline.scatter(lvl_time_spliceC, splineC_der_more, c='r', s=4, marker='s', label='Spline')
#       ax_Cspline.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3, label='Data')
#       ax_Cspline.legend(loc='upper right')
       ax_Cspline.set_ylabel('Level derivative (cm per min)')
       ax_Cspline.set_xlabel('Time (mins)')
       fig.savefig('sectionC_lvl_der_spline_vs_time_subplots_flowtimes.eps')
       plt.clf()
###################PLOT RATIOS######################
# think about moving averages like this 
#def moving_average(a, n=3) :
#    ret = np.cumsum(a, dtype=float)
#    ret[n:] = ret[n:] - ret[:-n]
#    return ret[n - 1:] / n

##### calculate ratio in more than one scenario####
##moving average will not help as much as my own average since I need a different dimensional array for the times that the flow is then calculated for
#    flows_averageA=moving_average(flows_spliceA)
######################SECTION A RATIO#############################
    if seq==4:
#flip the derivative sign
       splineA_der=splineA_der*(-1.0)
       #calc the average at each time of lvl sensor data points for each section at once
       flow_at_lvlA=integrate_time_range(lvl_time_spliceA, time_spliceA, flows_spliceA)
       print(len(flow_at_lvlA))
       print(len(lvl_time_spliceA))
       fig.subplots_adjust(hspace=0.2)
#subplot 1
       ax_Alvl=plt.subplot(411)
       ax_Alvl.set_title('Section A')
       plt.setp(ax_Alvl.get_xticklabels(), visible=False)
       ax_Alvl.scatter(lvl_time_spliceA, lvl_spliceA, c='b', s=3)
#       step=(numpy.amax(lvl_spliceA)-numpy.amin(lvl_spliceA))/2.0
#       ax_Alvl.yaxis.set_ticks(numpy.arange(numpy.amin(lvl_spliceA)-step, numpy.amax(lvl_spliceA)+step, step))
       ax_Alvl.set_ylabel('Level (cm)')
#subplot 2
       ax_Aflow_lvl=plt.subplot(412, sharex=ax_Alvl)
#       ax_Aflow_lvl.set_title('Section A')
       plt.setp(ax_Aflow_lvl.get_xticklabels(), visible=False)              
       ax_Aflow_lvl.scatter(lvl_time_spliceA, flow_at_lvlA, c='b', s=3)
       ax_Aflow_lvl.set_ylabel('Flow (grams per min)')
#subplot 3
       ax_Ader=plt.subplot(413, sharex=ax_Alvl)
       plt.setp(ax_Ader.get_xticklabels(), visible=False)
       ax_Ader.scatter(lvl_time_spliceA, splineA_der, c='r', s=4, marker='s')
       ax_Ader.set_ylabel('der (cm per min)')
       # take the ratios?
       ratioA=numpy.true_divide(flow_at_lvlA,splineA_der)
       # convert from grams per cm to cm squared
       ratioA=ratioA/0.125 # density of the LHe is 0.125 g per cm cube. So divide by this
       # convert from grams per cm to cm squared
       ratioA=ratioA/10000.0 #convert to meters squared
       # find bad values of the ratio 
       n=0
       ratioA_refined=[]
       lvlA_refined=[]
       lvlA_time_refined=[]
       while n<len(ratioA):
          if ratioA[n]<1 and ratioA[n]>0:
             ratioA_refined.append(ratioA[n])
             lvlA_refined.append(lvl_spliceA[n])
             lvlA_time_refined.append(lvl_time_spliceA[n])
          n+=1
       #done finding bad values
#subplot 4
       ax_Aratio=plt.subplot(414, sharex=ax_Alvl)
       ax_Aratio.scatter(lvlA_time_refined, ratioA_refined, c='g', s=4)
       ax_Aratio.set_ylabel('Ratio (meters squared)')
       ax_Aratio.set_xlabel('Time (mins)')
       fig.savefig('sectionA_ratio_flows_level_and_der_vs_time_subplots_integration.png')
       plt.clf()
       plt.scatter(lvlA_refined,ratioA_refined,c='r', s=4)
       plt.title('Ratio vs Level, section A')
       plt.ylabel('Ratio (meters squared)')
       plt.xlabel('Level (cm)')
       fig.savefig('sectionA_ratio_vs_levels_refined_integration.png')
       plt.clf()
#subplots vs lvl
       A_flow_ax=plt.subplot(311)
       A_flow_ax.set_title('Section A')
       plt.setp(A_flow_ax.get_xticklabels(), visible=False)
       A_flow_ax.scatter(lvl_spliceA, flow_at_lvlA, c='b', s=3)
       A_flow_ax.set_ylabel('Flow (grams per min)')
       #subplot2
       A_der_ax=plt.subplot(312)
       plt.setp(A_der_ax.get_xticklabels(), visible=False)
       A_der_ax.scatter(lvl_spliceA, splineA_der, c='b', s=3)
       A_der_ax.set_ylabel('Der (cm per min)')
       #subplot3
       A_ratio_ax=plt.subplot(313)
#       plt.setp(A_ratio_ax.get_xticklabels(), visible=False)
       A_ratio_ax.scatter(lvlA_refined, ratioA_refined, c='b', s=3)
       A_ratio_ax.set_ylabel('Ratio (meters squared)')
       A_ratio_ax.set_xlabel('Level (cm)')
       fig.savefig('sectionA_der_flow_ratio_vs_levels_integration.png')
       plt.clf()
#subplot 2
######################SECTION B RATIO#############################
    if seq==4:
       #calc the average at each time of lvl sensor data points for each section at once
#flip the derivative sign
       splineB_der=splineB_der*(-1.0)
       #calc the average at each time of lvl sensor data points for each section at once
       flow_at_lvlB=integrate_time_range(lvl_time_spliceB, time_spliceB, flows_spliceB)
       print(len(flow_at_lvlB))
       print(len(lvl_time_spliceB))
       fig.subplots_adjust(hspace=0.2)
       ax_Blvl=plt.subplot(411)
       ax_Blvl.set_title('Section B')
       plt.setp(ax_Blvl.get_xticklabels(), visible=False)
       ax_Blvl.scatter(lvl_time_spliceB, lvl_spliceB, c='b', s=3)
       ax_Blvl.set_ylabel('Level (cm)')
       ax_Bflow_lvl=plt.subplot(412)
#       ax_Bflow_lvl.set_title('Section B')
       plt.setp(ax_Bflow_lvl.get_xticklabels(), visible=False)              
       ax_Bflow_lvl.scatter(lvl_time_spliceB, flow_at_lvlB, c='b', s=3)
       ax_Bflow_lvl.set_ylabel('Flow (grams per min)')
       ax_Bder=plt.subplot(413)
       plt.setp(ax_Bder.get_xticklabels(), visible=False)
       ax_Bder.scatter(lvl_time_spliceB, splineB_der, c='r', s=4, marker='s')
       ax_Bder.set_ylabel('der (cm per min)')
       # take the ratios?
       ratioB=numpy.true_divide(flow_at_lvlB,splineB_der)
       ratioB=ratioB/0.125 # density of the LHe is 0.125 g per cm cube. So divide by this
       # convert from grams per cm to cm squared
       ratioB=ratioB/10000.0 #convert to meters squared
       # find bad values of the ratio 
       n=0
       ratioB_refined=[]
       lvlB_refined=[]
       lvlB_time_refined=[]
       while n<len(ratioB):
          if ratioB[n]<1 and ratioB[n]>0:
             ratioB_refined.append(ratioB[n])
             lvlB_refined.append(lvl_spliceB[n])
             lvlB_time_refined.append(lvl_time_spliceB[n])
          n+=1
       #done finding bad values
       ax_Bratio=plt.subplot(414)
       ax_Bratio.scatter(lvlB_time_refined, ratioB_refined, c='g', s=4)
       ax_Bratio.set_ylabel('Ratio (meters squared)')
       ax_Bratio.set_xlabel('Time (mins)')
       fig.savefig('sectionB_ratio_flows_level_and_der_vs_time_subplots_integrate.png')
       plt.clf()
       plt.scatter(lvlB_refined,ratioB_refined,c='r', s=4)
       plt.title('Ratio vs Level, section B')
       plt.ylabel('Ratio (meters squared)')
       plt.xlabel('Level (cm)')
       fig.savefig('sectionB_ratio_vs_levels_refined_integrate.png')
       plt.clf()
#subplots vs lvl
       B_flow_ax=plt.subplot(311)
       B_flow_ax.set_title('Section B')
       plt.setp(B_flow_ax.get_xticklabels(), visible=False)
       B_flow_ax.scatter(lvl_spliceB, flow_at_lvlB, c='b', s=3)
       B_flow_ax.set_ylabel('Flow (grams per min)')
       #subplot2
       B_der_ax=plt.subplot(312)
       plt.setp(B_der_ax.get_xticklabels(), visible=False)
       B_der_ax.scatter(lvl_spliceB, splineB_der, c='b', s=3)
       B_der_ax.set_ylabel('Der (cm per min)')
       #subplot3
       B_ratio_ax=plt.subplot(313)
#       plt.setp(A_ratio_ax.get_xticklabels(), visible=False)
       B_ratio_ax.scatter(lvlB_refined, ratioB_refined, c='b', s=3)
       B_ratio_ax.set_ylabel('Ratio (meters squared)')
       B_ratio_ax.set_xlabel('Level (cm)')
       fig.savefig('sectionB_der_flow_ratio_vs_levels_integrate.png')
       plt.clf()
#subplot 2
######################SECTION C RATIO#############################
    if seq==4:
       #calc the average at each time of lvl sensor data points for each section at once
#flip the derivative sign
       splineC_der=splineC_der*(-1.0)
       #calc the average at each time of lvl sensor data points for each section at once
       flow_at_lvlC=integrate_time_range(lvl_time_spliceC, time_spliceC, flows_spliceC)
       print(len(flow_at_lvlC))
       print(len(lvl_time_spliceC))
       fig.subplots_adjust(hspace=0.2)
       ax_Clvl=plt.subplot(411)
       ax_Clvl.set_title('Section C')
       plt.setp(ax_Clvl.get_xticklabels(), visible=False)
       ax_Clvl.scatter(lvl_time_spliceC, lvl_spliceC, c='b', s=3)
       ax_Clvl.set_ylabel('Level (cm)')
       ax_Cflow_lvl=plt.subplot(412)
       plt.setp(ax_Cflow_lvl.get_xticklabels(), visible=False)
       ax_Cflow_lvl.scatter(lvl_time_spliceC, flow_at_lvlC, c='b', s=3)
       ax_Cflow_lvl.set_ylabel('Flow (grams per min)')
       ax_Cder=plt.subplot(413)
       plt.setp(ax_Cder.get_xticklabels(), visible=False)
       ax_Cder.scatter(lvl_time_spliceC, splineC_der, c='r', s=4, marker='s')
       ax_Cder.set_ylabel('der (cm per min)')
       # take the ratios?
       ratioC=numpy.true_divide(flow_at_lvlC,splineC_der)
       ratioC=ratioC/0.125 # density of the LHe is 0.125 g per cm cube. So divide by this
       # convert from grams per cm to cm squared
       ratioC=ratioC/10000.0 #convert to meters squared
       # find bad values of the ratio 
       n=0
       ratioC_refined=[]
       lvlC_refined=[]
       lvlC_time_refined=[]
       while n<len(ratioC):
          if ratioC[n]<1 and ratioC[n]>0:
             ratioC_refined.append(ratioC[n])
             lvlC_refined.append(lvl_spliceC[n])
             lvlC_time_refined.append(lvl_time_spliceC[n])
          n+=1
       #done finding bad values
       ax_Cratio=plt.subplot(414)
       ax_Cratio.scatter(lvlC_time_refined, ratioC_refined, c='g', s=4)
       ax_Cratio.set_ylabel('Ratio (meters squared)')
       ax_Cratio.set_xlabel('Time (mins)')
       fig.savefig('sectionC_ratio_flows_level_and_der_vs_time_subplots_integrate.png')
       plt.clf()
       plt.scatter(lvlC_refined,ratioC_refined,c='r', s=4)
       plt.title('Ratio vs Level, section C')
       plt.ylabel('Ratio (meters squared)')
       plt.xlabel('Level (cm)')
       fig.savefig('sectionC_ratio_vs_levels_refined_integrate.png')
       plt.clf()
#subplots vs lvl
       C_flow_ax=plt.subplot(311)
       C_flow_ax.set_title('Section C')
       plt.setp(C_flow_ax.get_xticklabels(), visible=False)
       C_flow_ax.scatter(lvl_spliceC, flow_at_lvlC, c='b', s=3)
       C_flow_ax.set_ylabel('Flow (grams per min)')
       #subplot2
       C_der_ax=plt.subplot(312)
       plt.setp(C_der_ax.get_xticklabels(), visible=False)
       C_der_ax.scatter(lvl_spliceC, splineC_der, c='b', s=3)
       C_der_ax.set_ylabel('Der (cm per min)')
       #subplot3
       C_ratio_ax=plt.subplot(313)
#       plt.setp(A_ratio_ax.get_xticklabels(), visible=False)
       C_ratio_ax.scatter(lvlC_refined, ratioC_refined, c='b', s=3)
       C_ratio_ax.set_ylabel('Ratio (meters squared)')
       C_ratio_ax.set_xlabel('Level (cm)')
       fig.savefig('sectionC_der_flow_ratio_vs_levels_integrate.png')
       plt.clf()
#subplot 2
#############PLOT ALL SECTIONS TOGETHER###################
       plt.scatter(lvlA_refined,ratioA_refined,c='r', label='Section A')
       plt.scatter(lvlB_refined,ratioB_refined,c='g', label='Section B')
       plt.scatter(lvlC_refined,ratioC_refined,c='b', label='Section C')
       plt.xlim([0,110])
       plt.legend(loc='upper right')
       plt.title('Ratio vs Level, ALL sections')
       plt.ylabel('Ratio (meters squared)')
       plt.xlabel('Level (cm)')
       fig.savefig('All_sections_ratio_vs_levels_refined.png')
       plt.clf()
#####################calc error bars and mean#################
       bins5, avg5, errors5=find_avg_errors(lvlA_refined,lvlB_refined,lvlC_refined,ratioA_refined,ratioB_refined,ratioC_refined,5)
       bins10, avg10, errors10=find_avg_errors(lvlA_refined,lvlB_refined,lvlC_refined,ratioA_refined,ratioB_refined,ratioC_refined,10)
       bins20, avg20, errors20=find_avg_errors(lvlA_refined,lvlB_refined,lvlC_refined,ratioA_refined,ratioB_refined,ratioC_refined,20)
       bins40, avg40, errors40=find_avg_errors(lvlA_refined,lvlB_refined,lvlC_refined,ratioA_refined,ratioB_refined,ratioC_refined,40)
# fitting section
#       print("{} is array bins".format(bins40))
#       print("{} is array avg".format(avg40))
       p40, res40, rank40, single40,rcond40=numpy.polyfit(bins40, avg40, 0, full=True)
#       pxvals40=numpy.linspace(bins40[0],bins40[-1],100)
       pvals40=numpy.poly1d(p40)
#20
       p20, res20, rank20, single20,rcond20=numpy.polyfit(bins20, avg20, 0, full=True)
#       pxvals20=numpy.linspace(bins20[0],bins20[-1],100)
       pvals20=numpy.poly1d(p20)
#10
       p10, res10, rank10, single10,rcond10=numpy.polyfit(bins10, avg10, 0, full=True)
#       pxvals10=numpy.linspace(bins10[0],bins10[-1],100)
       pvals10=numpy.poly1d(p10)
#5
       p5, res5, rank5, single5,rcond5=numpy.polyfit(bins5, avg5, 0, full=True)
#       pxvals5=numpy.linspace(bins5[0],bins5[-1],100)
       pvals5=numpy.poly1d(p5)
#old line for plotting evaluated points
#       ax_40.scatter(pxvals40, pvals40(pxvals40), linestyle='--', color='black',label='fit')

#put here my x array of lvl values       x = np.array([1, 2, 3, 4, 5])
# put here the mean in that bin of the ratio       y = np.power(x, 2) # Effectively y = x**2
# put here the standard deviations for each bin       e = np.array([1.5, 2.6, 3.7, 4.6, 5.5])
       plot_name40="{delta}{binwidth}".format(delta=r'$\Delta L=$',binwidth=bins40[1]-bins40[0])
       plot_name20="{delta}{binwidth}".format(delta=r'$\Delta L=$',binwidth=bins20[1]-bins20[0])
       plot_name10="{delta}{binwidth}".format(delta=r'$\Delta L=$',binwidth=bins10[1]-bins10[0])
       plot_name5="{delta}{binwidth}".format(delta=r'$\Delta L=$',binwidth=bins5[1]-bins5[0])

       fig.subplots_adjust(hspace=0.1)
       ax_40=plt.subplot(411)
       plt.setp(ax_40.get_xticklabels(), visible=False)
       ax_40.set_title('Ratio vs Level, Binned')
       ax_40.errorbar(bins40, avg40, errors40, linestyle='None', marker='^', capsize=3, color='r',label=plot_name40)
       ax_40.axhline(pvals40(50), color='black', lw=1)
       ax_40.legend(loc='upper right')
       ax_20=plt.subplot(412,sharex=ax_40, sharey=ax_40)
       plt.setp(ax_20.get_xticklabels(), visible=False)
       ax_20.errorbar(bins20, avg20, errors20, linestyle='None', marker='^', capsize=3, color='b',label=plot_name20)
       ax_20.axhline(pvals20(50), color='black', lw=1)
#       ax_40.set_ylabel('Ratio (meters squared)')
       ax_20.legend(loc='upper right')
       ax_10=plt.subplot(413,sharex=ax_40,sharey=ax_40)
       plt.setp(ax_10.get_xticklabels(), visible=False)
       ax_10.errorbar(bins10, avg10, errors10, linestyle='None', marker='^', capsize=3, color='g',label=plot_name10)
       ax_10.axhline(pvals10(50), color='black', lw=1)
       ax_10.set_ylabel('Ratio (meters squared)')
       ax_10.legend(loc='upper right')
       ax_5=plt.subplot(414,sharex=ax_40,sharey=ax_40)
#       plt.setp(ax_5.get_xticklabels(), visible=False)
       ax_5.errorbar(bins5, avg5, errors5, linestyle='None', marker='^', capsize=3, color='orange',label=plot_name5)
       ax_5.axhline(pvals5(50), color='black', lw=1)
       ax_5.set_xlabel('Level (cm)')
       ax_5.legend(loc='upper right')
#       plt.errorbar(bins20, avg20, errors20, linestyle='None', marker='^', capsize=3, color='b', label=plot_name20)
#       plt.errorbar(bins10, avg10, errors10, linestyle='None', marker='^', capsize=3, color='g', label=plot_name10)
#       plt.errorbar(bins5, avg5, errors5, linestyle='None', marker='^', capsize=3, color='orange', label=plot_name5)
#       plt.show()
#       plt.title('Ratio vs Level, Binned')
#       plt.legend(loc='upper left')
#       plt.ylabel('Ratio (meters squared)')
#       plt.xlabel('Level (cm)')
       filename="Ratio_vs_levels_binwidth_subplots.png"
       fig.savefig(filename)
       plt.clf()
       cryo_data=numpy.genfromtxt('vol_vs_lvl_cryo.csv', dtype=float, delimiter=',', names=True)
       lvl_cryo = cryo_data['lvl']
       vol_cryo = cryo_data['vol']
       
      
#       plt.scatter(lvl_cryo,vol_cryo,s=3, c='black')
       x=numpy.arange(lvl_cryo[2],100,0.01)
       #convert fits to liters (multiply by 10000 for cm^2 then divide by 1000 for the conversion to liters) and add in the value of cryo at around 30cm (since data stops there?)
       fit40=(pvals40(50)*10.0)*x+(vol_cryo[2]-(pvals40(50)*10.0*x[0]))
       fit20=(pvals20(50)*10.0)*x+(vol_cryo[2]-(pvals20(50)*10.0*x[0]))
       fit10=(pvals10(50)*10.0)*x+(vol_cryo[2]-(pvals10(50)*10.0*x[0]))
       fit5=(pvals5(50)*10.0)*x+(vol_cryo[2]-(pvals5(50)*10.0*x[0]))
       plt.plot(x,fit40,'-.',x,fit20,':',x,fit10,x,fit5,'--',lvl_cryo,vol_cryo,'o')
       plt.xlim(30,100)
       plt.ylim(0,300)
       
       nvol=[]
       nlvl=[]
       j = 0
       while j < len(lvl_cryo):
           if lvl_cryo[j]<=100 and lvl_cryo[j]>=30:
               nlvl.append(lvl_cryo[j])
               nvol.append(vol_cryo[j])
               print(nlvl[-1])
           j+=1
       
     
      
       
       
       x2 =numpy.arange(nlvl[0],nlvl[-1],0.01)
       y2=numpy.polyfit(nlvl,nvol,1)
       y2[1]=y2[1]+nvol[0]
       z2=numpy.polyval(y2,x2)
       plt.plot(x2,z2,c='red',dashes=[6,2])
       plt.xlim(30,100)
       plt.ylim(0,300)
       
       slope  = format(pvals40(50)*10.0,'4.3f')
       slope2 = format(pvals20(50)*10.0,'4.3f')
       slope3 = format(pvals10(50)*10.0,'4.3f')
       slope4 = format(pvals5(50)*10.0,'4.3f')
       slope5 = format(y2[0],'4.3f')
   
       
       plt.legend(("m="+slope+"  ∆L=1.5375","m="+slope2+"  ∆L=3.075","m="+slope3+"  ∆L=6.15","m="+slope4+"  ∆L=12.3",'Data',"m="+slope5+"  Fit data"),loc='upper left',prop={'size': 18})
       plt.title('Geometric calibration using gaseous flow',fontsize=20)
       plt.ylabel('Volume (liters)',fontsize=16)
       plt.xlabel('Sensor Level (cm)',fontsize=16)
       
       plt.tick_params(labelsize=16)

       filename="Volume_vs_levels_cryo_doc.png"
       fig.savefig(filename)
       ### take the constants from fits and that is the slope of the Vol vs level plot. 
####################### ~0.35meters *(level cm) ####################################### 1meter=1000liters
#       plt.clf()


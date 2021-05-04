import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import sys
import ode

#S-I-R Rate Equations
def rates(d,tn): #calculate rates, and return an array of rates
    
    deriv=np.zeros(3) #deriv[0]=S, deriv[1]=I, deriv[2]=R
    
    deriv[0]=-r_SI*d[0]*d[1]    #dS/dt, change in Susceptible Population
    deriv[2]=r_RI*d[1]          #dR/dt, change in Recovered Population
    deriv[1]=-deriv[0]-deriv[2] #dI/dt, change in Infected Population
    
    return deriv

def WSNSModel(CarrierInput,IntInput,InfRateInput,RecoveryInput,TimeInput):
    #Model for Wanek School of Natural Sciences
    global r_RI, r_SI
    Name="Wanek School of Natural Sciences"
    Shorthand="WSNS"
    TotalPop=300 #Total Population
    SusPop=240   #Initial Susceptible Population
    InfPop=2     #Initial Infected Population
    TotalInteractions=IntInput           #number of daily interactions an infected person has with the total population
    SusceptibleInteractions=InfRateInput #how many interactions between an infected and susceptible people it takes to result in infection
    RI=1/RecoveryInput   #Recovery Rate (1/length of time it takes a person to recover)
    Time=TimeInput       #Time, in days, the disease is being modelled over
    
    #Carrier Stuff
    #List of Potential Carriers of the Disease and how they Impact the Spread
    Carrier=float(CarrierInput)
    BirdMod=0.4
    RodentMod=0.3
    InsectMod=0.7
    BloodMod=0.7
    SalivaMod=0.3
    TouchMod=0.5
    AirMod=0.2
    WaterMod=0.3
    NoMod=0.0
    
    #If-statements determine which Carrier is in use and assigns the proper impact on the spread
    if Carrier==1:
        CarrierMod=BirdMod
    if Carrier==2:
        CarrierMod=RodentMod
    if Carrier==3:
        CarrierMod=InsectMod
    if Carrier==4:
        CarrierMod=BloodMod
    if Carrier==5:
        CarrierMod=SalivaMod
    if Carrier==6:
        CarrierMod=TouchMod
    if Carrier==7:
        CarrierMod=AirMod
    if Carrier==8:
        CarrierMod=WaterMod
    if Carrier==9:
        CarrierMod=NoMod
    
    #Assigns the impact a carrier as on the spread, as well as determines how much, as a percentage, the carrier increased the spread
    CarrierImpact=CarrierMod
    
    #Time Setup
    t=0     #initial time 
    dt=0.01 #change in time

    #Redefining Variables for S-I-R Formula
    N=TotalPop         #total population
    S=SusPop           #initial amount of susceptible people
    I=InfPop           #initial amount of infected people
    R=0                #initial amount of recovered people
    f=TotalInteractions              #number of daily interactions an infected person has with the total population
    Ptrans=1/SusceptibleInteractions #fraction of interactions between an infected person and a susceptible person that results in spread
    r_SI=((f*Ptrans)/N)*CarrierImpact   #transmission rate, how likely it is for an infected person to spread the disease
    r_RI=RI                             #recovery rate, time it takes for a person to recover
    
    #Creating lists for each variable in SIR formula
    tmodel=[]
    Smodel=[]
    Imodel=[]
    Rmodel=[]
    data=np.array([S,I,R])

    #Loop to actually run the Euler/RK method differential equations
    while t<Time:

        #Selecting which model to use
        #data=ode.Euler(rates,data,t,dt) #Uses Euler Model
        #data=ode.RK2(rates,data,t,dt)   #Uses Runge-Kutta 2 Model
        data=ode.RK4(rates,data,t,dt)    #Uses Runge-Kutta 4 Model

        #Advancing time
        t=t+dt
    
        #Appending each variable list
        tmodel.append(t)
        Smodel.append(data[0])
        Imodel.append(data[1])
        Rmodel.append(data[2])
  
    #Plotting Data
    fig=plt.figure(figsize=(10,6))
    plt.title("SIR Model of a Disease across Different Buildings on HPU's Campus: " +str(Shorthand))
    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.plot(tmodel,Smodel,'b-',label='Susceptible Population')
    plt.plot(tmodel,Imodel,'r-',label='Infected Population')
    plt.plot(tmodel,Rmodel,'g-',label='Recovered Population')
    plt.legend(loc="best")
    plt.grid()
    plt.show()
#-------------------------------------------------------------
def NQSCModel(CarrierInput,IntInput,InfRateInput,RecoveryInput,TimeInput):
    #Model for Nido R. Quebin School of Communications
    global r_RI, r_SI
    Name="Nido R. Quebin School of Communications"
    Shorthand="NQSC"
    TotalPop=500 #Total Population
    SusPop=499   #Initial Susceptible Population
    InfPop=1     #Initial Infected Population
    TotalInteractions=IntInput           #number of daily interactions an infected person has with the total population
    SusceptibleInteractions=InfRateInput #how many interactions between an infected and susceptible people it takes to result in infection
    RI=1/RecoveryInput   #Recovery Rate (1/length of time it takes a person to recover)
    Time=TimeInput       #Time, in days, the disease is being modelled over
    
    #Carrier Stuff
    #List of Potential Carriers of the Disease and how they Impact the Spread
    Carrier=float(CarrierInput)
    BirdMod=0.2
    RodentMod=0.6
    InsectMod=0.7
    BloodMod=0.3
    SalivaMod=0.2
    TouchMod=0.4
    AirMod=0.4
    WaterMod=0.3
    NoMod=0.0
    
    #If-statements determine which Carrier is in use and assigns the proper impact on the spread
    if Carrier==1:
        CarrierMod=BirdMod
    if Carrier==2:
        CarrierMod=RodentMod
    if Carrier==3:
        CarrierMod=InsectMod
    if Carrier==4:
        CarrierMod=BloodMod
    if Carrier==5:
        CarrierMod=SalivaMod
    if Carrier==6:
        CarrierMod=TouchMod
    if Carrier==7:
        CarrierMod=AirMod
    if Carrier==8:
        CarrierMod=WaterMod
    if Carrier==9:
        CarrierMod=NoMod
    
    #Assigns the impact a carrier as on the spread, as well as determines how much, as a percentage, the carrier increased the spread
    CarrierImpact=CarrierMod
    
    #Time Setup
    t=0     #initial time 
    dt=0.01 #change in time

    #Redefining Variables for S-I-R Formula
    N=TotalPop         #total population
    S=SusPop           #initial amount of susceptible people
    I=InfPop           #initial amount of infected people
    R=0                #initial amount of recovered people
    f=TotalInteractions              #number of daily interactions an infected person has with the total population
    Ptrans=1/SusceptibleInteractions #fraction of interactions between an infected person and a susceptible person that results in spread
    r_SI=((f*Ptrans)/N)*CarrierImpact   #transmission rate, how likely it is for an infected person to spread the disease
    r_RI=RI                             #recovery rate, time it takes for a person to recover
    
    #Creating lists for each variable in SIR formula
    tmodel=[]
    Smodel=[]
    Imodel=[]
    Rmodel=[]
    data=np.array([S,I,R])

    #Loop to actually run the Euler/RK method differential equations
    while t<Time:

        #Selecting which model to use
        #data=ode.Euler(rates,data,t,dt) #Uses Euler Model
        #data=ode.RK2(rates,data,t,dt)   #Uses Runge-Kutta 2 Model
        data=ode.RK4(rates,data,t,dt)    #Uses Runge-Kutta 4 Model

        #Advancing time
        t=t+dt
    
        #Appending each variable list
        tmodel.append(t)
        Smodel.append(data[0])
        Imodel.append(data[1])
        Rmodel.append(data[2])
  
    #Plotting Data
    fig=plt.figure(figsize=(10,6))
    plt.title("SIR Model of a Disease across Different Buildings on HPU's Campus: " +str(Shorthand))
    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.plot(tmodel,Smodel,'b-',label='Susceptible Population')
    plt.plot(tmodel,Imodel,'r-',label='Infected Population')
    plt.plot(tmodel,Rmodel,'g-',label='Recovered Population')
    plt.legend(loc="best")
    plt.grid()
    plt.show()
#------------------------------------------------------------------------
def MiDmModel(CarrierInput,IntInput,InfRateInput,RecoveryInput,TimeInput):
    #Model for Millis Dormatory
    global r_RI, r_SI
    Name="Millis Dormatory"
    Shorthand="Millis"
    TotalPop=100 #Total Population
    SusPop=99    #Initial Susceptible Population
    InfPop=1     #Initial Infected Population
    TotalInteractions=IntInput           #number of daily interactions an infected person has with the total population
    SusceptibleInteractions=InfRateInput #how many interactions between an infected and susceptible people it takes to result in infection
    RI=1/RecoveryInput   #Recovery Rate (1/length of time it takes a person to recover)
    Time=TimeInput       #Time, in days, the disease is being modelled over
    
    #Carrier Stuff
    #List of Potential Carriers of the Disease and how they Impact the Spread
    Carrier=float(CarrierInput)
    BirdMod=0.2
    RodentMod=0.8
    InsectMod=0.8
    BloodMod=0.7
    SalivaMod=0.7
    TouchMod=0.9
    AirMod=0.3
    WaterMod=0.4
    NoMod=0.0
    
    #If-statements determine which Carrier is in use and assigns the proper impact on the spread
    if Carrier==1:
        CarrierMod=BirdMod
    if Carrier==2:
        CarrierMod=RodentMod
    if Carrier==3:
        CarrierMod=InsectMod
    if Carrier==4:
        CarrierMod=BloodMod
    if Carrier==5:
        CarrierMod=SalivaMod
    if Carrier==6:
        CarrierMod=TouchMod
    if Carrier==7:
        CarrierMod=AirMod
    if Carrier==8:
        CarrierMod=WaterMod
    if Carrier==9:
        CarrierMod=NoMod
    
    #Assigns the impact a carrier as on the spread, as well as determines how much, as a percentage, the carrier increased the spread
    CarrierImpact=CarrierMod
    
    #Time Setup
    t=0     #initial time 
    dt=0.01 #change in time

    #Redefining Variables for S-I-R Formula
    N=TotalPop         #total population
    S=SusPop           #initial amount of susceptible people
    I=InfPop           #initial amount of infected people
    R=0                #initial amount of recovered people
    f=TotalInteractions              #number of daily interactions an infected person has with the total population
    Ptrans=1/SusceptibleInteractions #fraction of interactions between an infected person and a susceptible person that results in spread
    r_SI=((f*Ptrans)/N)*CarrierImpact   #transmission rate, how likely it is for an infected person to spread the disease
    r_RI=RI                             #recovery rate, time it takes for a person to recover
    
    #Creating lists for each variable in SIR formula
    tmodel=[]
    Smodel=[]
    Imodel=[]
    Rmodel=[]
    data=np.array([S,I,R])

    #Loop to actually run the Euler/RK method differential equations
    while t<Time:

        #Selecting which model to use
        #data=ode.Euler(rates,data,t,dt) #Uses Euler Model
        #data=ode.RK2(rates,data,t,dt)   #Uses Runge-Kutta 2 Model
        data=ode.RK4(rates,data,t,dt)    #Uses Runge-Kutta 4 Model

        #Advancing time
        t=t+dt
    
        #Appending each variable list
        tmodel.append(t)
        Smodel.append(data[0])
        Imodel.append(data[1])
        Rmodel.append(data[2])
  
    #Plotting Data
    fig=plt.figure(figsize=(10,6))
    plt.title("SIR Model of a Disease across Different Buildings on HPU's Campus: " +str(Shorthand))
    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.plot(tmodel,Smodel,'b-',label='Susceptible Population')
    plt.plot(tmodel,Imodel,'r-',label='Infected Population')
    plt.plot(tmodel,Rmodel,'g-',label='Recovered Population')
    plt.legend(loc="best")
    plt.grid()
    plt.show()
#----------------------------------------------------
def ViDmModel(CarrierInput,IntInput,InfRateInput,RecoveryInput,TimeInput):
    #Model for Aldridge Village Dormatory
    global r_RI, r_SI
    Name="Aldridge Village Dormatory"
    Shorthand="Village"
    TotalPop=250 #Total Population
    SusPop=249   #Initial Susceptible Population
    InfPop=1     #Initial Infected Population
    TotalInteractions=IntInput           #number of daily interactions an infected person has with the total population
    SusceptibleInteractions=InfRateInput #how many interactions between an infected and susceptible people it takes to result in infection
    RI=1/RecoveryInput   #Recovery Rate (1/length of time it takes a person to recover)
    Time=TimeInput       #Time, in days, the disease is being modelled over
    
    #Carrier Stuff
    #List of Potential Carriers of the Disease and how they Impact the Spread
    Carrier=float(CarrierInput)
    BirdMod=0.2
    RodentMod=0.4
    InsectMod=0.6
    BloodMod=0.5
    SalivaMod=0.2
    TouchMod=0.7
    AirMod=0.3
    WaterMod=0.4
    NoMod=0.0
    
    #If-statements determine which Carrier is in use and assigns the proper impact on the spread
    if Carrier==1:
        CarrierMod=BirdMod
    if Carrier==2:
        CarrierMod=RodentMod
    if Carrier==3:
        CarrierMod=InsectMod
    if Carrier==4:
        CarrierMod=BloodMod
    if Carrier==5:
        CarrierMod=SalivaMod
    if Carrier==6:
        CarrierMod=TouchMod
    if Carrier==7:
        CarrierMod=AirMod
    if Carrier==8:
        CarrierMod=WaterMod
    if Carrier==9:
        CarrierMod=NoMod
    
    #Assigns the impact a carrier as on the spread, as well as determines how much, as a percentage, the carrier increased the spread
    CarrierImpact=CarrierMod
    
    #Time Setup
    t=0     #initial time 
    dt=0.01 #change in time

    #Redefining Variables for S-I-R Formula
    N=TotalPop         #total population
    S=SusPop           #initial amount of susceptible people
    I=InfPop           #initial amount of infected people
    R=0                #initial amount of recovered people
    f=TotalInteractions              #number of daily interactions an infected person has with the total population
    Ptrans=1/SusceptibleInteractions #fraction of interactions between an infected person and a susceptible person that results in spread
    r_SI=((f*Ptrans)/N)*CarrierImpact   #transmission rate, how likely it is for an infected person to spread the disease
    r_RI=RI                             #recovery rate, time it takes for a person to recover
    
    #Creating lists for each variable in SIR formula
    tmodel=[]
    Smodel=[]
    Imodel=[]
    Rmodel=[]
    data=np.array([S,I,R])

    #Loop to actually run the Euler/RK method differential equations
    while t<Time:

        #Selecting which model to use
        #data=ode.Euler(rates,data,t,dt) #Uses Euler Model
        #data=ode.RK2(rates,data,t,dt)   #Uses Runge-Kutta 2 Model
        data=ode.RK4(rates,data,t,dt)    #Uses Runge-Kutta 4 Model

        #Advancing time
        t=t+dt
    
        #Appending each variable list
        tmodel.append(t)
        Smodel.append(data[0])
        Imodel.append(data[1])
        Rmodel.append(data[2])
  
    #Plotting Data
    fig=plt.figure(figsize=(10,6))
    plt.title("SIR Model of a Disease across Different Buildings on HPU's Campus: " +str(Shorthand))
    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.plot(tmodel,Smodel,'b-',label='Susceptible Population')
    plt.plot(tmodel,Imodel,'r-',label='Infected Population')
    plt.plot(tmodel,Rmodel,'g-',label='Recovered Population')
    plt.legend(loc="best")
    plt.grid()
    plt.show()
#-----------------------------------------------------------
def CustomModel(CustomName,IntInput,InfRateInput,RecoveryInput,TimeInput):
    #Custom Model
    global r_RI, r_SI
    Name=CustomName
    print()
    PopInput=int(input("Enter an amount for the total population: "))
    print()
    SusInput=int(input("Enter an amount for the initial susceptible population: "))
    print()
    InfInput=int(input("Enter an amount for the initial infected population: "))
    print()
    print("As of right now, this custom model does not support any carrier modifications. Maybe that will change in the future when a more advanced version of this code is created...")
    print()
    
    #Initial Population Counts
    TotalPop=PopInput   #Total Population
    SusPop=SusInput     #Initial Susceptible Population
    InfPop=InfInput     #Initial Infected Population
    TotalInteractions=IntInput           #number of daily interactions an infected person has with the total population
    SusceptibleInteractions=InfRateInput #how many interactions between an infected and susceptible people it takes to result in infection
    RI=1/RecoveryInput   #Recovery Rate (1/length of time it takes a person to recover)
    Time=TimeInput       #Time, in days, the disease is being modelled over
    
    #Time Setup
    t=0     #initial time 
    dt=0.01 #change in time

    #Redefining Variables for S-I-R Formula
    N=TotalPop         #total population
    S=SusPop           #initial amount of susceptible people
    I=InfPop           #initial amount of infected people
    R=0                #initial amount of recovered people
    f=TotalInteractions              #number of daily interactions an infected person has with the total population
    Ptrans=1/SusceptibleInteractions #fraction of interactions between an infected person and a susceptible person that results in spread
    r_SI=((f*Ptrans)/N)   #transmission rate, how likely it is for an infected person to spread the disease
    r_RI=RI               #recovery rate, time it takes for a person to recover
    
    #Creating lists for each variable in SIR formula
    tmodel=[]
    Smodel=[]
    Imodel=[]
    Rmodel=[]
    data=np.array([S,I,R])

    #Loop to actually run the Euler/RK method differential equations
    while t<Time:

        #Selecting which model to use
        #data=ode.Euler(rates,data,t,dt) #Uses Euler Model
        #data=ode.RK2(rates,data,t,dt)   #Uses Runge-Kutta 2 Model
        data=ode.RK4(rates,data,t,dt)    #Uses Runge-Kutta 4 Model

        #Advancing time
        t=t+dt
    
        #Appending each variable list
        tmodel.append(t)
        Smodel.append(data[0])
        Imodel.append(data[1])
        Rmodel.append(data[2])
  
    #Plotting Data
    fig=plt.figure(figsize=(10,6))
    plt.title("SIR Model of Disease Spread for "  +str(Name))
    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.plot(tmodel,Smodel,'b-',label='Susceptible Population')
    plt.plot(tmodel,Imodel,'r-',label='Infected Population')
    plt.plot(tmodel,Rmodel,'g-',label='Recovered Population')
    plt.legend(loc="best")
    plt.grid()
    plt.show()
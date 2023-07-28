#Begin Post Processing
#Open the Output Data Base for the current Job
from visualization import *
odb = openOdb(path='strain3.odb', readOnly=True);#avoid warning
myAssembly = odb.rootAssembly;

#Creating a temporary variable to hold the frame repository provides the same functionality and speeds up the process
frameRepository = odb.steps['Step-1'].frames;
frameS=[];
frameIVOL=[];
# Create a Coordinate System in the Compsite direction (Global) 
coordSys=odb.rootAssembly.DatumCsysByThreePoints(name='NEW', coordSysType=CARTESIAN,origin=(0,0,0),point1=(1,0,0),point2=(0,1,0))
# Transform stresses from Lamina Coordinates 
# to Laminate Coordinate System defined in NEW 
stressTrans=frameRepository[-1].fieldOutputs['S'].getTransformedField(datumCsys=coordSys)

#Get only the last frame [-1]
frameS.insert(0,stressTrans.getSubset(position=INTEGRATION_POINT));
frameIVOL.insert(0,frameRepository[-1].fieldOutputs['IVOL'].getSubset(position=INTEGRATION_POINT));
#Total Volume
Tot_Vol=0;
#Stress Sum
Tot_Stress=0;
#
for II in range(0,len(frameS[-1].values)):
     Tot_Vol=Tot_Vol+frameIVOL[0].values[II].data;
     Tot_Stress=Tot_Stress+frameS[0].values[II].data * frameIVOL[0].values[II].data;

#Calculate Average
Avg_Stress = Tot_Stress/Tot_Vol;
#print 'Abaqus/Standard Stress Tensor Order:'
#From Abaqus Analysis User's Manual - 1.2.2 Conventions - Convention used for stress and strain components
#print 'Average stresses Global CSYS: 11-22-33-12-13-23';
#print Avg_Stress;
C11 = Avg_Stress[2]#z-component,1-direction
C21 = Avg_Stress[0]#x-component,2-direction
C31 = Avg_Stress[1]#y-component,3-direction 

from visualization import *
odb = openOdb(path='strain1.odb', readOnly=True);#avoid warning
myAssembly = odb.rootAssembly;
#Creating a temporary variable to hold the frame repository provides the same functionality and speeds up the process
frameRepository = odb.steps['Step-1'].frames;
frameS=[];
frameIVOL=[];
#Get only the last frame [-1]
stressTrans=frameRepository[-1].fieldOutputs['S'].getTransformedField(datumCsys=coordSys)
frameS.insert(0,stressTrans.getSubset(position=INTEGRATION_POINT));
frameIVOL.insert(0,frameRepository[-1].fieldOutputs['IVOL'].getSubset(position=INTEGRATION_POINT));
#Total Volume
Tot_Vol=0;
#Stress Sum
Tot_Stress=0;
#
for II in range(0,len(frameS[-1].values)):
     Tot_Vol=Tot_Vol+frameIVOL[0].values[II].data;
     Tot_Stress=Tot_Stress+frameS[0].values[II].data * frameIVOL[0].values[II].data;

#Calculate Average
Avg_Stress = Tot_Stress/Tot_Vol;
#print 'Abaqus/Standard Stress Tensor Order:'
#From Abaqus Analysis User's Manual - 1.2.2 Conventions - Convention used for stress and strain components
#print 'Average stresses Global CSYS: 11-22-33-12-13-23';
#print Avg_Stress;
C12 = Avg_Stress[2]#z-component,1-direction
C22 = Avg_Stress[0]#x-component,2-direction
C32 = Avg_Stress[1]#y-component,3-direction 

from visualization import *
odb = openOdb(path='strain2.odb', readOnly=True);#avoid warning
myAssembly = odb.rootAssembly;
#Creating a temporary variable to hold the frame repository provides the same functionality and speeds up the process
frameRepository = odb.steps['Step-1'].frames;
frameS=[];
frameIVOL=[];
#Get only the last frame [-1]
stressTrans=frameRepository[-1].fieldOutputs['S'].getTransformedField(datumCsys=coordSys)
frameS.insert(0,stressTrans.getSubset(position=INTEGRATION_POINT));
frameIVOL.insert(0,frameRepository[-1].fieldOutputs['IVOL'].getSubset(position=INTEGRATION_POINT));
#Total Volume
Tot_Vol=0;
#Stress Sum
Tot_Stress=0;
#
for II in range(0,len(frameS[-1].values)):
     Tot_Vol=Tot_Vol+frameIVOL[0].values[II].data;
     Tot_Stress=Tot_Stress+frameS[0].values[II].data * frameIVOL[0].values[II].data;

#Calculate Average
Avg_Stress = Tot_Stress/Tot_Vol;
#print 'Abaqus/Standard Stress Tensor Order:'
#From Abaqus Analysis User's Manual - 1.2.2 Conventions - Convention used for stress and strain components 
#print 'Average stresses in Global CSYS: 11-22-33-12-13-23';
#print Avg_Stress;
C13 = Avg_Stress[2]#z-component,1-direction
C23 = Avg_Stress[0]#x-component,2-direction
C33 = Avg_Stress[1]#y-component,3-direction 

from visualization import *
odb = openOdb(path='strain4.odb', readOnly=True);#avoid warning
myAssembly = odb.rootAssembly;

#Creating a temporary variable to hold the frame repository provides the same functionality and speeds up the process
frameRepository = odb.steps['Step-1'].frames;
frameS=[];
frameIVOL=[];
# Create a Coordinate System in the Compsite direction (Global) 
coordSys=odb.rootAssembly.DatumCsysByThreePoints(name='NEW', coordSysType=CARTESIAN,origin=(0,0,0),point1=(1,0,0),point2=(0,1,0))
# Transform stresses from Lamina Coordinates 
# to Laminate Coordinate System defined in NEW 
stressTrans=frameRepository[-1].fieldOutputs['S'].getTransformedField(datumCsys=coordSys)

#Get only the last frame [-1]
frameS.insert(0,stressTrans.getSubset(position=INTEGRATION_POINT));
frameIVOL.insert(0,frameRepository[-1].fieldOutputs['IVOL'].getSubset(position=INTEGRATION_POINT));
#Total Volume
Tot_Vol=0;
#Stress Sum
Tot_Stress=0;
#
for II in range(0,len(frameS[-1].values)):
     Tot_Vol=Tot_Vol+frameIVOL[0].values[II].data;
     Tot_Stress=Tot_Stress+frameS[0].values[II].data * frameIVOL[0].values[II].data;

#Calculate Average
Avg_Stress = Tot_Stress/Tot_Vol;
#print 'Abaqus/Standard Stress Tensor Order:'
#From Abaqus Analysis User's Manual - 1.2.2 Conventions - Convention used for stress and strain components
#print 'Average stresses Global CSYS: 11-22-33-12-13-23';
#print Avg_Stress;
C44 = Avg_Stress[3]#xy-component,2-3-direction

from visualization import *
odb = openOdb(path='strain5.odb', readOnly=True);#avoid warning
myAssembly = odb.rootAssembly;

#Creating a temporary variable to hold the frame repository provides the same functionality and speeds up the process
frameRepository = odb.steps['Step-1'].frames;
frameS=[];
frameIVOL=[];
# Create a Coordinate System in the Compsite direction (Global) 
coordSys=odb.rootAssembly.DatumCsysByThreePoints(name='NEW', coordSysType=CARTESIAN,origin=(0,0,0),point1=(1,0,0),point2=(0,1,0))
# Transform stresses from Lamina Coordinates 
# to Laminate Coordinate System defined in NEW 
stressTrans=frameRepository[-1].fieldOutputs['S'].getTransformedField(datumCsys=coordSys)

#Get only the last frame [-1]
frameS.insert(0,stressTrans.getSubset(position=INTEGRATION_POINT));
frameIVOL.insert(0,frameRepository[-1].fieldOutputs['IVOL'].getSubset(position=INTEGRATION_POINT));
#Total Volume
Tot_Vol=0;
#Stress Sum
Tot_Stress=0;
#
for II in range(0,len(frameS[-1].values)):
     Tot_Vol=Tot_Vol+frameIVOL[0].values[II].data;
     Tot_Stress=Tot_Stress+frameS[0].values[II].data * frameIVOL[0].values[II].data;

#Calculate Average
Avg_Stress = Tot_Stress/Tot_Vol;
#print 'Abaqus/Standard Stress Tensor Order:'
#From Abaqus Analysis User's Manual - 1.2.2 Conventions - Convention used for stress and strain components
#print 'Average stresses Global CSYS: 11-22-33-12-13-23';
#print Avg_Stress;
C66 = Avg_Stress[4]#xz-component,2-1-direction

from visualization import *
odb = openOdb(path='strain6.odb', readOnly=True);#avoid warning
myAssembly = odb.rootAssembly;

#Creating a temporary variable to hold the frame repository provides the same functionality and speeds up the process
frameRepository = odb.steps['Step-1'].frames;
frameS=[];
frameIVOL=[];
# Create a Coordinate System in the Compsite direction (Global) 
coordSys=odb.rootAssembly.DatumCsysByThreePoints(name='NEW', coordSysType=CARTESIAN,origin=(0,0,0),point1=(1,0,0),point2=(0,1,0))
# Transform stresses from Lamina Coordinates 
# to Laminate Coordinate System defined in NEW 
stressTrans=frameRepository[-1].fieldOutputs['S'].getTransformedField(datumCsys=coordSys)

#Get only the last frame [-1]
frameS.insert(0,stressTrans.getSubset(position=INTEGRATION_POINT));
frameIVOL.insert(0,frameRepository[-1].fieldOutputs['IVOL'].getSubset(position=INTEGRATION_POINT));
#Total Volume
Tot_Vol=0;
#Stress Sum
Tot_Stress=0;
#
for II in range(0,len(frameS[-1].values)):
     Tot_Vol=Tot_Vol+frameIVOL[0].values[II].data;
     Tot_Stress=Tot_Stress+frameS[0].values[II].data * frameIVOL[0].values[II].data;

#Calculate Average
Avg_Stress = Tot_Stress/Tot_Vol;
#print 'Abaqus/Standard Stress Tensor Order:'
#From Abaqus Analysis User's Manual - 1.2.2 Conventions - Convention used for stress and strain components
#print 'Average stresses Global CSYS: 11-22-33-12-13-23';
#print Avg_Stress;
C55 = Avg_Stress[5]#zy-component,1-3-direction


 


EL=C11-2*C12*C21/(C22+C23)              # Longitudinal E1 modulus
nuL=C12/(C22+C23)                       # 12 Poisson coefficient
ET=(C11*(C22+C23)-2*C12*C12)*(C22-C23)/(C11*C22-C12*C21)
                                        # Transversal E2 modulus
nuT=(C11*C23-C12*C21)/(C11*C22-C12*C21) # 23 Poisson coefficient
GT=(C22-C23)/2 # or GT=ET/2/(1+nuT)     # 23 Shear stiffness



print "If Moduli are entered in TPa and dimensions in microns, results are in TPa"
print "E1=",EL,"TPa"
print "E2=",ET,"TPa"
print "PR12=",nuL
print "PR23=",nuT
print "G23=",GT,"TPa"
print "G12=",C55,"TPa"



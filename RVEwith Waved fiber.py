#3D RVE
#Isotropic Fiber
#Elastic Matrix
#Whole model
#Square Array
# import modules
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

# Geometry & Material input data
E_fiber, nu_fiber, v_fiber, d_fiber = 237000.0, 0.3, 0.6, 7.0
E_matrix, nu_matrix = 2500.0, 0.35

a_2 = d_fiber *sqrt(pi/(v_fiber))/4.0 # half RVE along z
a_3 = a_2 # half RVE along y 
r_f = 3.5 # fiber radius
amplitude=10.0
wavelength=100.0
partname='WavyRVE'
MeshingSize=0.5
strain = [[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]    # epsilon_xx, epsilon_yy, epsilon__zz
a_1=wavelength/2.0

for jj in range(len(strain)):
    mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-'+str(jj+1))
    mdb.models['Model-'+str(jj+1)].ConstrainedSketch(name='__sweep__', sheetSize=30.0)
    mdb.models['Model-'+str(jj+1)].sketches['__sweep__'].Spline(points=((-wavelength/2-5.0, 0.0), (-wavelength/2, 0.0), (
        -wavelength/4, amplitude), (0.0, 2*amplitude), (wavelength/4, amplitude), (wavelength/2, 0.0), (wavelength/2+5, 0.0)))
    mdb.models['Model-'+str(jj+1)].ConstrainedSketch(name='__profile__', sheetSize=200.0, 
        transform=(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 
        -0.0, -0.0, -0.0, -wavelength/4, 0.0, 0.0))
    mdb.models['Model-'+str(jj+1)].sketches['__profile__'].ConstructionLine(point1=(-15.0, 
        0.0), point2=(15.0, 0.0))
    mdb.models['Model-'+str(jj+1)].sketches['__profile__'].ConstructionLine(point1=(0.0, 
        -15.0), point2=(0.0, 15.0))    
    mdb.models['Model-'+str(jj+1)].sketches['__profile__'].rectangle(point1=(-a_3, 
        -a_2), point2=(a_3,a_2))
    mdb.models['Model-'+str(jj+1)].Part(dimensionality=THREE_D, name=partname, type=
        DEFORMABLE_BODY)
    mdb.models['Model-'+str(jj+1)].parts[partname].BaseSolidSweep(path=
        mdb.models['Model-'+str(jj+1)].sketches['__sweep__'], sketch=
        mdb.models['Model-'+str(jj+1)].sketches['__profile__'])
    del mdb.models['Model-'+str(jj+1)].sketches['__profile__']
    del mdb.models['Model-'+str(jj+1)].sketches['__sweep__']
    p=mdb.models['Model-'+str(jj+1)].parts[partname]
    mdb.models['Model-'+str(jj+1)].ConstrainedSketch(gridSpacing=24.25, name='__profile__', 
        sheetSize=970.29, transform=
        mdb.models['Model-'+str(jj+1)].parts['WavyRVE'].MakeSketchTransform(
        sketchPlane=p.faces[4], 
        sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-1'].parts['WavyRVE'].edges[10], 
        sketchOrientation=RIGHT, origin=(wavelength/2+5, 0.0, 0.0)))
    mdb.models['Model-'+str(jj+1)].parts['WavyRVE'].projectReferencesOntoSketch(filter=
        COPLANAR_EDGES, sketch=mdb.models['Model-'+str(jj+1)].sketches['__profile__'])
    mdb.models['Model-'+str(jj+1)].sketches['__profile__'].CircleByCenterPerimeter(center=(
        0.0, 0.0), point1=(r_f, 0.0))
    mdb.models['Model-'+str(jj+1)].parts['WavyRVE'].PartitionFaceBySketch(faces= p.faces[4], sketch=mdb.models['Model-'+str(jj+1)].sketches['__profile__'], sketchUpEdge=
        mdb.models['Model-'+str(jj+1)].parts['WavyRVE'].edges[10])
    del mdb.models['Model-'+str(jj+1)].sketches['__profile__']
    mdb.models['Model-'+str(jj+1)].parts['WavyRVE'].PartitionCellBySweepEdge(cells= p.cells [0], edges=(mdb.models['Model-'+str(jj+1)].parts['WavyRVE'].edges[4], ), sweepPath=
        mdb.models['Model-'+str(jj+1)].parts['WavyRVE'].edges[8])

    #Material 
    mdb.models['Model-'+str(jj+1)].Material(name='Fiber')
    mdb.models['Model-'+str(jj+1)].materials['Fiber'].Elastic(table=((E_fiber,nu_fiber), )) 
    mdb.models['Model-'+str(jj+1)].Material(name='Matrix')
    mdb.models['Model-'+str(jj+1)].materials['Matrix'].Elastic(table=((E_matrix, nu_matrix), ))

    #section
    mdb.models['Model-'+str(jj+1)].HomogeneousSolidSection(material='Fiber', name='Fiber', 
        thickness=None)
    mdb.models['Model-'+str(jj+1)].HomogeneousSolidSection(material='Matrix', name='Matrix', thickness=None)  
    #Assign section
    e=p.edges
    f=p.faces
    c=p.cells 
    v=p.vertices
    p.Set(cells= c.findAt(((wavelength/2,a_3,a_2),),), name='Set-Matrix') 
    p.Set(cells= c.findAt(((wavelength/2+5,0,0),),), name='Set-Fiber') 

    p.SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=p.sets['Set-Matrix'], sectionName='Matrix', 
        thicknessAssignment=FROM_SECTION)  
    p.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=p.sets['Set-Fiber'], sectionName='Fiber', thicknessAssignment=FROM_SECTION) 
    #Assembly
    mdb.models['Model-'+str(jj+1)].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-'+str(jj+1)].rootAssembly.Instance(dependent=OFF, name='RVE', part=p)
    mdb.models['Model-'+str(jj+1)].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 10.0, 
        0.0), axisPoint=(0.0, 0.0, 0.0), instanceList=('RVE', ))
    aa=mdb.models['Model-'+str(jj+1)].rootAssembly.instances['RVE']
    a=mdb.models['Model-'+str(jj+1)].rootAssembly
    c=aa.cells
    f=aa.faces
    e=aa.edges
    v=aa.vertices 
    #Mesh
    a.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges= e[7:8], minSizeFactor=0.1, size=MeshingSize)
    a.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges= e[2:3], minSizeFactor=0.1, size=MeshingSize)
    a.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges= e[8:9], minSizeFactor=0.1, size=MeshingSize)
    a.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges= e[9:10], minSizeFactor=0.1, size=MeshingSize)
    a.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges= e[10:11], minSizeFactor=0.1, size=MeshingSize)
    a.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges= e[14:15], minSizeFactor=0.1, size=MeshingSize*5)
    a.generateMesh(regions=( mdb.models['Model-1'].rootAssembly.instances['RVE'], ))  
    # Step
    mdb.models['Model-'+str(jj+1)].StaticStep(name='Step-1', nlgeom=OFF, previous='Initial')
    #Field Output Request
    mdb.models['Model-'+str(jj+1)].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'LE', 'U', 'RF', 'IVOL'))    
    # Constraint equations
    f=aa.faces
    e=aa.edges
    v=aa.vertices
    a.Set(faces= (f[1:2],f[8:9]), name='back')
    a.Set(faces= f[4:5], name='left')
    a.Set(faces= f[5:6], name='top')
    a.Set(faces= f[6:7], name='right')
    a.Set(faces= f[3:4], name='bottom')
    a.Set(faces= (f[7:8],f[2:3]), name='front')
    a.Set(edges= e[3:4], name='cg')
    a.Set(edges= e[4:5], name='hg')
    a.Set(edges= e[5:6], name='dh')
    a.Set(edges= e[6:7], name='dc')
    a.Set(edges= e[7:8], name='ef')
    a.Set(edges= e[8:9], name='bf')
    a.Set(edges= e[9:10], name='ab')
    a.Set(edges= e[10:11], name='ae')
    a.Set(edges= e[11:12], name='eh')
    a.Set(edges= e[12:13], name='fg')
    a.Set(edges= e[13:14], name='ad')
    a.Set(edges= e[14:15], name='bc')
    a.Set(name='c', vertices= v[2:3])
    a.Set(name='g', vertices= v[3:4])
    a.Set(name='h', vertices= v[4:5])
    a.Set(name='d', vertices= v[5:6])
    a.Set(name='e', vertices= v[6:7])
    a.Set(name='f', vertices= v[7:8])
    a.Set(name='b', vertices= v[8:9])
    a.Set(name='a', vertices= v[9:10])    
    front=a.sets['front']
    back=a.sets['back']
    right=a.sets['right']
    left=a.sets['left']
    top=a.sets['top']
    bottom=a.sets['bottom']
    edgebf=a.sets['bf']
    edgeae=a.sets['ae']
    edgeab=a.sets['ab']
    edgeef=a.sets['ef']
    edgecg=a.sets['cg']
    edgedh=a.sets['dh']
    edgedc=a.sets['dc']
    edgehg=a.sets['hg']
    edgebc=a.sets['bc']
    edgead=a.sets['ad']
    edgeeh=a.sets['eh']
    edgefg=a.sets['fg']
    VertexA=a.sets['a']
    VertexB=a.sets['b']
    VertexC=a.sets['c']
    VertexD=a.sets['d']
    VertexE=a.sets['e']
    VertexF=a.sets['f']
    VertexG=a.sets['g']
    VertexH=a.sets['h']
    def TakeVertexOut(edge):
        edge.pop(0)
        edge.pop(-1)
        return edge
    def SortListOfNodes(edge,coordinate):
        newlist = []
        oldlist = []
        for ii in range(len(edge.nodes)):
            oldlist.append( edge.nodes[ii].coordinates[coordinate])
        orderedlist = sorted(oldlist)
        for ii in range(len(oldlist)):
            vecindex = oldlist.index(orderedlist[ii])
            #newlist.append(oldlist[vecindex])
            newlist.append(edge.nodes[vecindex].label-1)
        return newlist
    # Construction Pairs of nodes for edges
    # Edge Front right and back left #soriting based on y position which is the coordinate number1 BF DH
    ParingFacesbfdh = []
    ParingFacesbfdh.append(TakeVertexOut(SortListOfNodes(edgebf,1)))
    ParingFacesbfdh.append(TakeVertexOut(SortListOfNodes(edgedh,1)))
    #soriting based on z position which is the coordinate number0 AB HG
    ParingFacesabhg = []
    ParingFacesabhg.append(TakeVertexOut(SortListOfNodes(edgeab,0)))
    ParingFacesabhg.append(TakeVertexOut(SortListOfNodes(edgehg,0)))
    #soriting based on x position which is the coordinate number2 BC EH
    ParingFacesbceh = []
    ParingFacesbceh.append(TakeVertexOut(SortListOfNodes(edgebc,2)))
    ParingFacesbceh.append(TakeVertexOut(SortListOfNodes(edgeeh,2)))
    #soriting based on y position which is the coordinate number1 AE CG
    ParingFacesaecg = []
    ParingFacesaecg.append(TakeVertexOut(SortListOfNodes(edgeae,1)))
    ParingFacesaecg.append(TakeVertexOut(SortListOfNodes(edgecg,1)))
    #soriting based on z position which is the coordinate number0 EF DC
    ParingFacesefdc = []
    ParingFacesefdc.append(TakeVertexOut(SortListOfNodes(edgeef,0)))
    ParingFacesefdc.append(TakeVertexOut(SortListOfNodes(edgedc,0)))
    #soriting based on z position which is the coordinate number2 GF DA
    ParingFacesgfda = []
    ParingFacesgfda.append(TakeVertexOut(SortListOfNodes(edgefg,2)))
    ParingFacesgfda.append(TakeVertexOut(SortListOfNodes(edgead,2)))   
    a.ReferencePoint(point=(0.0, 0.0,wavelength/2+20.0))
    MasterIndex=a.features['RP-1']
    a.Set(name='MasterNode', referencePoints=(a.referencePoints[MasterIndex.id], ))
    # U1=1, to write the CE's
    mdb.models['Model-'+str(jj+1)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(referencePoints=(
        a.referencePoints[MasterIndex.id], 
        )), u1=1.0, u2=UNSET, ur3=UNSET)
    #Rigidbodymotion
    Centralfornt=aa.nodes.getClosest([0.0,0.0,wavelength/2+5],1,).label-1
    Centralback=aa.nodes.getClosest([0.0,0.0,-wavelength/2-5],1,).label-1
    a.Set(name='centerfront', nodes=(aa.nodes[Centralfornt:Centralfornt+1],))
    a.Set(name='centerback', nodes=(aa.nodes[Centralback:Centralback+1],))
    Centralbody=aa.nodes.getClosest([0, 2*amplitude,0],1,).label-1
    a.Set(name='Centralbody', nodes=(aa.nodes[Centralbody:Centralbody+1],))
    mdb.models['Model-'+str(jj+1)].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
        region=
        a.sets['centerfront'], 
        u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models['Model-'+str(jj+1)].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-3', 
        region=
        a.sets['centerback'], u1=
        SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mdb.models['Model-'+str(jj+1)].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-4', 
        region=
        a.sets['Centralbody'], 
        u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    #Writting Constrain equations 
    #For the Vertices
    #Master V1 (Contained in 1 and 3)   
    #z component=1 direction  
    #x component = 2 dicertion
    #y component = 3 direction
    ezz=strain[jj][2]
    exx=strain[jj][0]
    eyy=strain[jj][1]
    ezx=strain[jj][4]/2.0
    ezy=strain[jj][5]/2.0
    exy=strain[jj][3]/2.0 
    mdb.models['Model-'+str(jj+1)].Equation(name='B-H-3', terms=((1.0, 'b', 3),  (-1.0, 'h', 3),(-ezz*2*a_1-ezx*2*a_2-ezy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='B-H-1', terms=((1.0, 'b', 1),  (-1.0, 'h', 1),(-ezx*2*a_1-exx*2*a_2-exy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='B-H-2', terms=((1.0, 'b', 2),  (-1.0, 'h', 2),(-ezy*2*a_1-exy*2*a_2-eyy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='C-E-3', terms=((1.0, 'c', 3),  (-1.0, 'e', 3),(ezz*2*a_1-ezx*2*a_2-ezy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='C-E-1', terms=((1.0, 'c', 1),  (-1.0, 'e', 1),(ezx*2*a_1-exx*2*a_2-exy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='C-E-2', terms=((1.0, 'c', 2),  (-1.0, 'e', 2),(ezy*2*a_1-exy*2*a_2-eyy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='A-G-3', terms=((1.0, 'a', 3),  (-1.0, 'g', 3),(-ezz*2*a_1+ezx*2*a_2-ezy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='A-G-1', terms=((1.0, 'a', 1),  (-1.0, 'g', 1),(-ezx*2*a_1+exx*2*a_2-exy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='A-G-2', terms=((1.0, 'a', 2),  (-1.0, 'g', 2),(-ezy*2*a_1+exy*2*a_2-eyy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='F-D-3', terms=((1.0, 'f', 3),  (-1.0, 'd', 3),(-ezz*2*a_1-ezx*2*a_2+ezy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='F-D-1', terms=((1.0, 'f', 1),  (-1.0, 'd', 1),(-ezx*2*a_1-exx*2*a_2+exy*2*a_3,'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].Equation(name='F-D-2', terms=((1.0, 'f', 2),  (-1.0, 'd', 2),(-ezy*2*a_1-exy*2*a_2+eyy*2*a_3,'MasterNode', 1)))
    #masternode is a dummy node in step one apply a unite displacement to it
    # For the Edges (Between pairs of nodes)
    for ii in range(len(ParingFacesbfdh[0])):
        a.Set(name='MasterFaceBF-'+str(ii), nodes=(
            aa.nodes[ParingFacesbfdh[0][ii]:ParingFacesbfdh[0][ii]+1],))
        a.Set(name='SlaveFaceDH-'+str(ii), nodes=(
            aa.nodes[ParingFacesbfdh[1][ii]:ParingFacesbfdh[1][ii]+1],))   
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintfbdh-3-'+str(ii), terms=((1.0, 
            'MasterFaceBF-'+str(ii), 3), (-1.0, 'SlaveFaceDH-'+str(ii), 3), (-ezz*2*a_1-ezx*2*a_2, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintfbdh-1-'+str(ii), terms=((1.0, 
            'MasterFaceBF-'+str(ii), 1), (-1.0, 'SlaveFaceDH-'+str(ii), 1), (-ezx*2*a_1-exx*2*a_2, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintfbdh-2-'+str(ii), terms=((1.0, 
            'MasterFaceBF-'+str(ii), 2), (-1.0, 'SlaveFaceDH-'+str(ii), 2), (-ezy*2*a_1-exy*2*a_2, 
            'MasterNode', 1)))
    for ii in range(len(ParingFacesabhg[0])):
        a.Set(name='MasterFaceAB-'+str(ii), nodes=(
            aa.nodes[ParingFacesabhg[0][ii]:ParingFacesabhg[0][ii]+1],))
        a.Set(name='SlaveFaceHG-'+str(ii), nodes=(
            aa.nodes[ParingFacesabhg[1][ii]:ParingFacesabhg[1][ii]+1],))   
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintabhg-3-'+str(ii), terms=((1.0, 
            'MasterFaceAB-'+str(ii), 3), (-1.0, 'SlaveFaceHG-'+str(ii), 3), (-ezz*2*a_1-ezy*2*a_3, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintabhg-1-'+str(ii), terms=((1.0, 
            'MasterFaceAB-'+str(ii), 1), (-1.0, 'SlaveFaceHG-'+str(ii), 1), (-ezx*2*a_1-exy*2*a_3, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintabhg-2-'+str(ii), terms=((1.0, 
            'MasterFaceAB-'+str(ii), 2), (-1.0, 'SlaveFaceHG-'+str(ii), 2), (-ezy*2*a_1-eyy*2*a_3, 
            'MasterNode', 1)))
    for ii in range(len(ParingFacesbceh[0])):
        a.Set(name='MasterFaceBC-'+str(ii), nodes=(
            aa.nodes[ParingFacesbceh[0][ii]:ParingFacesbceh[0][ii]+1],))
        a.Set(name='SlaveFaceEH-'+str(ii), nodes=(
            aa.nodes[ParingFacesbceh[1][ii]:ParingFacesbceh[1][ii]+1],))   
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintbceh-3-'+str(ii), terms=((1.0, 
            'MasterFaceBC-'+str(ii), 3), (-1.0, 'SlaveFaceEH-'+str(ii), 3), (-ezx*2*a_2-ezy*2*a_3, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintbceh-1-'+str(ii), terms=((1.0, 
            'MasterFaceBC-'+str(ii), 1), (-1.0, 'SlaveFaceEH-'+str(ii), 1), (-exx*2*a_2-exy*2*a_3, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintbceh-2-'+str(ii), terms=((1.0, 
            'MasterFaceBC-'+str(ii), 2), (-1.0, 'SlaveFaceEH-'+str(ii), 2), (-ezy*2*a_2-eyy*2*a_3, 
            'MasterNode', 1)))
    for ii in range(len(ParingFacesaecg[0])):
        a.Set(name='MasterFaceAE-'+str(ii), nodes=(
            aa.nodes[ParingFacesaecg[0][ii]:ParingFacesaecg[0][ii]+1],))
        a.Set(name='SlaveFaceCG-'+str(ii), nodes=(
            aa.nodes[ParingFacesaecg[1][ii]:ParingFacesaecg[1][ii]+1],))   
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintaecg-3-'+str(ii), terms=((1.0, 
            'MasterFaceAE-'+str(ii), 3), (-1.0, 'SlaveFaceCG-'+str(ii), 3), (-ezz*2*a_1+ezx*2*a_2, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintaecg-1-'+str(ii), terms=((1.0, 
            'MasterFaceAE-'+str(ii), 1), (-1.0, 'SlaveFaceCG-'+str(ii), 1), (-ezx*2*a_1+exx*2*a_2, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintaecg-2-'+str(ii), terms=((1.0, 
            'MasterFaceAE-'+str(ii), 2), (-1.0, 'SlaveFaceCG-'+str(ii), 2), (-ezy*2*a_1+exy*2*a_2, 
            'MasterNode', 1)))
    for ii in range(len(ParingFacesefdc[0])):
        a.Set(name='MasterFaceEF-'+str(ii), nodes=(
            aa.nodes[ParingFacesefdc[0][ii]:ParingFacesefdc[0][ii]+1],))
        a.Set(name='SlaveFaceCD-'+str(ii), nodes=(
            aa.nodes[ParingFacesefdc[1][ii]:ParingFacesefdc[1][ii]+1],))   
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintefcd-3-'+str(ii), terms=((1.0, 
            'MasterFaceEF-'+str(ii), 3), (-1.0, 'SlaveFaceCD-'+str(ii), 3), (-ezz*2*a_1+ezy*2*a_3, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintefcd-1-'+str(ii), terms=((1.0, 
            'MasterFaceEF-'+str(ii), 1), (-1.0, 'SlaveFaceCD-'+str(ii), 1), (-ezx*2*a_1+exy*2*a_3, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintefcd-2-'+str(ii), terms=((1.0, 
            'MasterFaceEF-'+str(ii), 2), (-1.0, 'SlaveFaceCD-'+str(ii), 2), (-ezy*2*a_1+eyy*2*a_3, 
            'MasterNode', 1)))    
    for ii in range(len(ParingFacesgfda[0])):
        a.Set(name='MasterFaceGF-'+str(ii), nodes=(
            aa.nodes[ParingFacesgfda[0][ii]:ParingFacesgfda[0][ii]+1],))
        a.Set(name='SlaveFaceDA-'+str(ii), nodes=(
            aa.nodes[ParingFacesgfda[1][ii]:ParingFacesgfda[1][ii]+1],))   
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintgfda-3-'+str(ii), terms=((1.0, 
            'MasterFaceGF-'+str(ii), 3), (-1.0, 'SlaveFaceDA-'+str(ii), 3), (-ezx*2*a_2+ezy*2*a_3, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintgfda-1-'+str(ii), terms=((1.0, 
            'MasterFaceGF-'+str(ii), 1), (-1.0, 'SlaveFaceDA-'+str(ii), 1), (-exx*2*a_2+exy*2*a_3, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintgfda-2-'+str(ii), terms=((1.0, 
            'MasterFaceGF-'+str(ii), 2), (-1.0, 'SlaveFaceDA-'+str(ii), 2), (-ezy*2*a_2+eyy*2*a_3, 
            'MasterNode', 1)))
    ##nodes at vertices    
    VertexAlabel=aa.nodes.getClosest([-a_2, a_3, wavelength/2+5],1,).label-1
    VertexBlabel=aa.nodes.getClosest([a_2, a_3, wavelength/2+5],1,).label-1
    VertexClabel=aa.nodes.getClosest([a_2, a_3, -wavelength/2-5],1,).label-1
    VertexDlabel=aa.nodes.getClosest([-a_2, a_3, -wavelength/2-5],1,).label-1
    VertexElabel=aa.nodes.getClosest([-a_2, -a_3, wavelength/2+5],1,).label-1
    VertexFlabel=aa.nodes.getClosest([a_2,- a_3, wavelength/2+5],1,).label-1
    VertexGlabel=aa.nodes.getClosest([a_2, -a_3, -wavelength/2-5],1,).label-1
    VertexHlabel=aa.nodes.getClosest([-a_2, -a_3, -wavelength/2-5],1,).label-1
    toplabel=[]
    for ii in range(len(top.nodes)):
        toplabel.append(top.nodes[ii].label-1)                
    for ii in range(len(ParingFacesbceh[0])) :     #bc
        toplabel.remove(ParingFacesbceh[0][ii])        
    for ii in range(len(ParingFacesabhg[0])) :        #ab
        toplabel.remove(ParingFacesabhg[0][ii])         
    for ii in range(len(ParingFacesefdc[1])) :          #dc
        toplabel.remove(ParingFacesefdc[1][ii])       
    for ii in range(len(ParingFacesgfda[1])) :                #da
         toplabel.remove(ParingFacesgfda[1][ii])   
    toplabel.remove(VertexAlabel)
    toplabel.remove(VertexBlabel)
    toplabel.remove(VertexClabel)
    toplabel.remove(VertexDlabel)              
    for ii in range(len(toplabel)):
        a.Set(name='MasterFaceTop-'+str(ii), nodes= (aa.nodes[toplabel[ii]:toplabel[ii]+1],))
        slavelabel=aa.nodes.getClosest([(aa.nodes[toplabel[ii]].coordinates[0]),(aa.nodes[toplabel[ii]].coordinates[1]-2*a_3),(aa.nodes[toplabel[ii]].coordinates[2])],1,).label-1
        a.Set(name='SlaveFaceBottom-'+str(ii), nodes=(aa.nodes[slavelabel:slavelabel+1],))
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constrainttopbottom-3-'+str(ii), terms=((1.0, 
            'MasterFaceTop-'+str(ii), 3), (-1.0, 'SlaveFaceBottom-'+str(ii), 3), (-ezy*2*a_3, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constrainttopbottom-1-'+str(ii), terms=((1.0, 
            'MasterFaceTop-'+str(ii), 1), (-1.0, 'SlaveFaceBottom-'+str(ii), 1), (-exy*2*a_3, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constrainttopbottom-2-'+str(ii), terms=((1.0, 
            'MasterFaceTop-'+str(ii), 2), (-1.0, 'SlaveFaceBottom-'+str(ii), 2), (-eyy*2*a_3, 
            'MasterNode', 1)))
    rightlabel=[]
    for ii in range(len(right.nodes)):
        rightlabel.append(right.nodes[ii].label-1)
                
    for ii in range(len(ParingFacesbceh[0])) :     #bc
        rightlabel.remove(ParingFacesbceh[0][ii])        
    for ii in range(len(ParingFacesaecg[1])) :        #cg
        rightlabel.remove(ParingFacesaecg[1][ii])         
    for ii in range(len(ParingFacesgfda[0])) :          #gf
        rightlabel.remove(ParingFacesgfda[0][ii])        
    for ii in range(len(ParingFacesbfdh[0])) :                #bf
        rightlabel.remove(ParingFacesbfdh[0][ii])  
    rightlabel.remove(VertexBlabel)
    rightlabel.remove(VertexClabel)
    rightlabel.remove(VertexFlabel)
    rightlabel.remove(VertexGlabel)               
    for ii in range(len(rightlabel)):
        a.Set(name='MasterFaceright-'+str(ii), nodes= (aa.nodes[rightlabel[ii]:rightlabel[ii]+1],))
        slavelabel=aa.nodes.getClosest([-(aa.nodes[rightlabel[ii]].coordinates[0]),(aa.nodes[rightlabel[ii]].coordinates[1]),(aa.nodes[rightlabel[ii]].coordinates[2])],1,).label-1
        a.Set(name='SlaveFaceleft-'+str(ii), nodes=(aa.nodes[slavelabel:slavelabel+1],))
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintrightleft-3-'+str(ii), terms=((1.0, 
            'MasterFaceright-'+str(ii), 3), (-1.0, 'SlaveFaceleft-'+str(ii), 3), (-ezx*2*a_2, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintrightleft-1-'+str(ii), terms=((1.0, 
            'MasterFaceright-'+str(ii), 1), (-1.0, 'SlaveFaceleft-'+str(ii), 1), (-exx*2*a_2, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintrightleft-2-'+str(ii), terms=((1.0, 
            'MasterFaceright-'+str(ii), 2), (-1.0, 'SlaveFaceleft-'+str(ii), 2), (-exy*2*a_2, 
            'MasterNode', 1)))       
    frontlabel=[]
    for ii in range(len(front.nodes)):
        frontlabel.append(front.nodes[ii].label-1)                    
    for ii in range(len(ParingFacesabhg[0])) :     #ab
        frontlabel.remove(ParingFacesabhg[0][ii])        
    for ii in range(len(ParingFacesbfdh[0])) :        #bf
        frontlabel.remove(ParingFacesbfdh[0][ii])        
    for ii in range(len(ParingFacesefdc[0])) :          #ef
        frontlabel.remove(ParingFacesefdc[0][ii])       
    for ii in range(len(ParingFacesaecg[0])) :                #ae
        frontlabel.remove(ParingFacesaecg[0][ii])          
    frontlabel.remove(VertexAlabel)
    frontlabel.remove(VertexBlabel)
    frontlabel.remove(VertexElabel)
    frontlabel.remove(VertexFlabel)            
    for ii in range(len(frontlabel)):
        a.Set(name='MasterFacefront-'+str(ii), nodes= (aa.nodes[frontlabel[ii]:frontlabel[ii]+1],))
        slavelabel=aa.nodes.getClosest([(aa.nodes[frontlabel[ii]].coordinates[0]),(aa.nodes[frontlabel[ii]].coordinates[1]),-(aa.nodes[frontlabel[ii]].coordinates[2])],1,).label-1
        a.Set(name='SlaveFaceback-'+str(ii), nodes=(aa.nodes[slavelabel:slavelabel+1],))
        # ( i=1)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintfrontback-3-'+str(ii), terms=((1.0, 
            'MasterFacefront-'+str(ii), 3), (-1.0, 'SlaveFaceback-'+str(ii), 3), (-ezz*2*a_1, 
            'MasterNode', 1)))
        # ( i=2)
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintfrontback-1-'+str(ii), terms=((1.0, 
            'MasterFacefront-'+str(ii), 1), (-1.0, 'SlaveFaceback-'+str(ii), 1), (-ezx*2*a_1, 
            'MasterNode', 1)))
        # (i=3)    
        mdb.models['Model-'+str(jj+1)].Equation(name='Constraintfrontback-2-'+str(ii), terms=((1.0, 
            'MasterFacefront-'+str(ii), 2), (-1.0, 'SlaveFaceback-'+str(ii), 2), (-ezy*2*a_1, 
            'MasterNode', 1)))
    mdb.models['Model-'+str(jj+1)].rootAssembly.regenerate()
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-'+str(jj+1), modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name='strain'+str(jj+1), nodalOutputPrecision=SINGLE, 
        numCpus=4, numDomains=4, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
        ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    mdb.jobs['strain'+str(jj+1)].submit(consistencyChecking=OFF)
    mdb.jobs['strain'+str(jj+1)].waitForCompletion()

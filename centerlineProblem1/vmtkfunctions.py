from math import *
from vmtk import pypes, vmtkscripts
from vtk import *
import csv


def vmtkbifurcationreferencesystems(centerline):
    computer = vmtkscripts.vmtkBifurcationReferenceSystems()
    computer.Centerlines = centerline
    computer.RadiusArrayName = 'MaximumInscribedSphereRadius'
    computer.BlankingArrayName = 'Blanking'
    computer.GroupIdsArrayName = 'GroupIds'
    computer.Execute()
    return computer.ReferenceSystems


def vmtkbifurcationsections(surface, centerline, distance=1):
    sectioner = vmtkscripts.vmtkBifurcationSections()
    sectioner.Surface = surface
    sectioner.Centerlines = centerline
    sectioner.NumberOfDistanceSpheres = distance
    sectioner.RadiusArrayName = 'MaximumInscribedSphereRadius'
    sectioner.GroupIdsArrayName = 'GroupIds'
    sectioner.CenterlineIdsArrayName = 'CenterlineIds'
    sectioner.TractIdsArrayName = 'TractIds'
    sectioner.BlankingArrayName = 'Blanking'
    sectioner.Execute()
    return sectioner.BifurcationSections


def vmtkbifurcationvectors(centerline, referencesystem):
    computer = vmtkscripts.vmtkBifurcationVectors()
    computer.Centerlines = centerline
    computer.ReferenceSystems = referencesystem
    computer.RadiusArrayName = 'MaximumInscribedSphereRadius'
    computer.GroupIdsArrayName = 'GroupIds'
    computer.CenterlineIdsArrayName = 'CenterlineIds'
    computer.TractIdsArrayName = 'TractIds'
    computer.BlankingArrayName = 'Blanking'
    computer.ReferenceSystemsNormalArrayName = 'Normal'
    computer.ReferenceSystemsUpNormalArrayName = 'UpNormal'
    computer.Execute()
    return computer.BifurcationVectors

def vmtkbranchclipper(surface, centerline, groupidlist=[0], insideout=0):
    clipper = vmtkscripts.vmtkBranchClipper()
    clipper.Surface = surface
    clipper.Centerlines = centerline
    # clipper.InsideOut = insideout
    # clipper.GroupIds = groupidlist
    clipper.GroupIdsArrayName = 'GroupIds'
    clipper.RadiusArrayName = 'MaximumInscribedSphereRadius'
    clipper.BlankingArrayName = 'Blanking'
    clipper.Execute()
    return clipper.Surface


def vmtkbranchextractor(centerline):
    extractor = vmtkscripts.vmtkBranchExtractor()
    extractor.Centerlines = centerline
    extractor.RadiusArrayName = 'MaximumInscribedSphereRadius'
    extractor.Execute()
    return extractor.Centerlines


def vmtkbranchmapping(surface, centerline, referencesystem):
    mapper = vmtkscripts.vmtkBranchMapping()
    mapper.Surface = surface
    mapper.Centerlines = centerline
    mapper.ReferenceSystems = referencesystem
    mapper.AbscissasArrayName = 'Abscissas'
    mapper.NormalsArrayName = 'ParallelTransportNormals'
    mapper.GroupIdsArrayName = 'GroupIds'
    mapper.CenterlineIdsArrayName = 'CenterlineIds'
    mapper.TractIdsArrayName = 'TractIds'
    mapper.ReferenceSystemsNormalArrayName = 'Normal'
    mapper.RadiusArrayName = 'MaximumInscribedSphereRadius'
    mapper.BlankingArrayName = 'Blanking'
    mapper.AngularMetricArrayName = 'AngularMetric'
    mapper.AbscissaMetricArrayName = 'AbscissaMetric'
    mapper.Execute()
    return mapper.Surface


def vmtkbranchmetrics(surface, centerline):
    computer = vmtkscripts.vmtkBranchMetrics()
    computer.Surface = surface
    computer.Centerlines = centerline
    computer.AbscissasArrayName = 'Abscissas'
    computer.NormalsArrayName = 'ParallelTransportNormals'
    computer.RadiusArrayName = 'MaximumInscribedSphereRadius'
    computer.GroupIdsArrayName = 'GroupIds'
    computer.CenterlineIdsArrayName = 'CenterlineIds'
    computer.TractIdsArrayName = 'TractIds'
    computer.BlankingArrayName = 'Blanking'
    computer.Execute()
    return computer.Surface


def vmtkbranchpatching(surface, longitudinalpatchsize=1.0,
                       circularnumberofpatches=12):
    patcher = vmtkscripts.vmtkBranchPatching()
    patcher.Surface = surface
    patcher.LongitudinalPatchSize = longitudinalpatchsize
    patcher.CircularNumberOfPatches = circularnumberofpatches
    patcher.UseConnectivity = 1
    patcher.GroupIdsArrayName = 'GroupIds'
    patcher.LongitudinalMappingArrayName = 'StretchedMapping'
    patcher.CircularMappingArrayName = 'AngularMetric'
    patcher.Execute()
    return (patcher.Surface, patcher.PatchedData)


def vmtkcenterlineattributes(centerline):
    computer = vmtkscripts.vmtkCenterlineAttributes()
    computer.Centerlines = centerline
    computer.Execute()
    return computer.Centerlines


def vmtkcenterlinegeometry(centerline):
    computer = vmtkscripts.vmtkCenterlineGeometry()
    computer.Centerlines = centerline
    computer.Execute()
    return computer.Centerlines

def vmtkcenterlinemerge(centerline, length=.1):
    merger = vmtkscripts.vmtkCenterlineMerge()
    merger.Centerlines = centerline
    merger.Length = length
    merger.RadiusArrayName = 'MaximumInscribedSphereRadius'
    merger.GroupIdsArrayName = 'GroupIds'
    merger.CenterlineIdsArrayName = 'CenterlineIds'
    merger.BlankingArrayName = 'Blanking'
    merger.TractIdsArrayName = 'TractIds'
    merger.Execute()
    return merger.Centerlines


def vmtkcenterlinemodeller(centerline, size=[64, 64, 64]):
    modeller = vmtkscripts.vmtkCenterlineModeller()
    modeller.Centerlines = centerline
    modeller.RadiusArrayName = 'MaximumInscribedSphereRadius'
    modeller.SampleDimensions = size
    modeller.Execute()
    return modeller.Image


def vmtkcenterlineoffsetattributes(centerline, referencesystems,
                                   refgroupid=1):
    computer = vmtkscripts.vmtkCenterlineOffsetAttributes()
    computer.Centerlines = centerline
    computer.ReferenceSystems = referencesystems
    computer.ReferenceGroupId = refgroupid
    computer.AbscissasArrayName = 'Abscissas'
    computer.NormalsArrayName = 'ParallelTransportNormals'
    computer.GroupIdsArrayName = 'GroupIds'
    computer.CenterlineIdsArrayName = 'CenterlineIds'
    computer.ReferenceSystemsNormalArrayName = 'Normal'
    computer.Execute()
    return computer.Centerlines


def vmtkcenterlineresampling(centerline, length=.1):
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = centerline
    resampler.Length = length
    resampler.Execute()
    return resampler.Centerlines


def vmtkcenterlines(surface, sourcepoints, targetpoints,endpoints=0):
    computer = vmtkscripts.vmtkCenterlines()
    computer.Surface = surface
    computer.SeedSelectorName = 'pointlist'
    computer.SourcePoints = sourcepoints
    computer.TargetPoints = targetpoints
    computer.AppendEndPoints = endpoints

    # computer.DelaunayTessellationOutputFileName = 'delaunay.vtp'
    # computer.VoronoiDiagramOutputFileName = 'voronoi.vtp'
    computer.Execute()
    return computer.Centerlines


def vmtkcenterlinesinteractive(surface, endpoints=0):
    computer = vmtkscripts.vmtkCenterlines()
    computer.Surface = surface
    computer.AppendEndPoints = endpoints
    computer.SeedSelectorName = 'pickpoint'
    computer.Execute()
    return computer.Centerlines


def vmtkcenterlinesections(surface, centerline):
    sectioner = vmtkscripts.vmtkCenterlineSections()
    sectioner.Surface = surface
    sectioner.Centerlines = centerline
    sectioner.Execute()
    return sectioner.CenterlineSections


def vmtkcenterlinesmoothing(centerline, iterations=100, factor=0.1):
    smoother = vmtkscripts.vmtkCenterlineSmoothing()
    smoother.Centerlines = centerline
    smoother.NumberOfSmoothingIterations = iterations
    smoother.SmoothingFactor = factor
    smoother.Execute()
    return smoother.Centerlines

def vmtkcenterlinesvoronoi(surface,sourcepoints, targetpoints):
    computer = vmtkscripts.vmtkCenterlines()
    computer.Surface = surface
    computer.AppendEndPoints = 1
    computer.SeedSelectorName = 'pointlist'
    computer.SourcePoints = sourcepoints
    computer.TargetPoints = targetpoints
    computer.Execute()
    return computer.Centerlines, computer.VoronoiDiagram


def vmtkdelaunayvoronoi(surface, removesubresolution=1):
    computer = vmtkscripts.vmtkDelaunayVoronoi()
    computer.Surface = surface
    computer.RemoveSubresolutionTetrahedra = removesubresolution
    computer.Execute()
    return computer.VoronoiDiagram


def vmtkdistancetocenterlines(surface, centerline):
    distance = vmtkscripts.vmtkDistanceToCenterlines()
    distance.Surface = surface
    distance.Centerlines = centerline
    distance.RadiusArrayName = 'MaximumInscribedSphereRadius'
    distance.Execute()
    return distance.Surface

def vmtkflowextensions(surface, interactive=1, extensionlength=10):
    extender = vmtkscripts.vmtkFlowExtensions()
    extender.Surface = surface
    extender.ExtensionMode = 'boundarynormal'
    extender.TransitionRatio = .25
    extender.Interactive = interactive
    extender.ExtensionLength = extensionlength
    extender.Execute()
    return extender.Surface

def vmtkflowextensionscenterline(surface, cl, interactive=1, extensionlength=10):
    extender = vmtkscripts.vmtkFlowExtensions()
    extender.Surface = surface
    extender.ExtensionMode = 'centerlinedirection'
    extender.Centerlines = cl
    extender.TransitionRatio = .25
    extender.Interactive = interactive
    extender.ExtensionLength = extensionlength
    extender.Execute()
    return extender.Surface


def vmtkicpregistration(surface, referencesurface,iters=100,nlmks=1000):
    registration = vmtkscripts.vmtkICPRegistration()
    registration.Surface = surface
    registration.ReferenceSurface = referencesurface
    registration.DistanceArrayName = 'original_distance'
    registration.MaximumNumberOfIterations = iters
    registration.MaximumNumberOfLandmarks = nlmks
    registration.Execute()
    return registration.Surface


def vmtkimagereader(filename):
    reader = vmtkscripts.vmtkImageReader()
    reader.InputFileName = filename
    reader.Execute()
    return reader.Image


def vmtkimagewriter(image, filename):
    writer = vmtkscripts.vmtkImageWriter()
    writer.Image = image
    writer.OutputFileName = filename
    writer.Execute()


def vmtkmarchingcubes(image):
    marcher = vmtkscripts.vmtkMarchingCubes()
    marcher.Image = image
    marcher.Execute()
    return marcher.Surface


def vmtkmeshreader(filename):
    reader = vmtkscripts.vmtkMeshReader()
    reader.InputFileName = filename
    reader.Execute()
    return reader.Mesh


def vmtkmeshvectorfromcomponents(mesh):
    computer = vmtkscripts.vmtkMeshVectorFromComponents()
    computer.Mesh = mesh
    computer.VectorArrayName = 'Velocity'
    computer.ComponentsArrayNames = ['VelocityX', 'VelocityX', 'VelocityX']
    computer.Execute()
    return computer.Mesh


def vmtkmeshwriter(mesh, filename):
    writer = vmtkscripts.vmtkMeshWriter()
    writer.Mesh = mesh
    writer.OutputFileName = filename
    writer.Execute()


def vmtknetworkextraction(surface):
    extractor = vmtkscripts.vmtkNetworkExtraction()
    extractor.Surface = surface
    extractor.Execute()
    return extractor.Network


def vmtkpointsplitextractor(centerline, splitpoint, gap=1.0, tolerance=0.0001):
    extractor = vmtkscripts.vmtkPointSplitExtractor()
    extractor.Centerlines = centerline
    extractor.RadiusArrayName = 'MaximumInscribedSphereRadius'
    extractor.GroupIdsArrayName = 'GroupIds'
    extractor.SplitPoint = splitpoint
    extractor.Execute()
    return extractor.Centerlines

def vmtkpolyballmodeller(voronoi, size=[64, 64, 64]):
    modeller = vmtkscripts.vmtkPolyBallModeller()
    modeller.Surface = voronoi
    modeller.RadiusArrayName = 'MaximumInscribedSphereRadius'
    modeller.SampleDimensions = size
    modeller.Execute()
    return modeller.Image

def vmtksurfacebooleanoperation(surface,surface2,operation='intersection'):
    boolean = vmtkscripts.vmtkSurfaceBooleanOperation()
    boolean.Surface(surface)
    boolean.surface2(surface2)
    boolean.Operation(operation)
    boolean.Execute()
    return boolean.Surface

def vmtksurfacedecimation(surface):
    decimator = vmtkscripts.vmtkSurfaceDecimation()
    decimator.Surface = surface
    decimator.Execute()
    return decimator.Surface


def vmtksurfacesubdivision(surface):
    divider = vmtkscripts.vmtkSurfaceSubdivision()
    divider.Surface = surface
    divider.Execute()
    return divider.Surface


def vmtksurfacecapper(surface,method='centerpoint',nrings=4,const=0.1):
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = surface
    capper.Method = method
    capper.NumberOfRings = nrings
    capper.ConstraintFactor = const
    capper.Interactive = 0
    capper.Execute()
    return capper.Surface


def vmtksurfacecenterlineprojection(surface, centerline):
    projection = vmtkscripts.vmtkSurfaceCenterlineProjection()
    projection.Surface = surface
    projection.Centerlines = centerline
    projection.RadiusArrayName = 'MaximumInscribedSphereRadius'
    projection.Execute()
    return projection.Surface


def vmtksurfacedistance(surface, referencesurface,
                        distancearrayname='Distance'):
    computer = vmtkscripts.vmtkSurfaceDistance()
    computer.Surface = surface
    computer.ReferenceSurface = referencesurface
    computer.DistanceArrayName = distancearrayname
    computer.DistanceVectorsArrayName = ''
    computer.SignedDistanceArrayName = ''
    computer.Execute()
    return computer.Surface

# def vmtksurfacecurvature(surface, curvtype='gaussian', abos= False, median=False):
#     computer = vmtkscripts.vmtkSurfaceCurvature()
#     computer.Surface = surface
#     computer.CurvatureType = curvtype
#     computer.AbsoluteCurvature = abso
#     computer.MedianFiltering = median
#     computer.Execute()
#     return computer.Surface


def vmtksurfacekiteremoval(surface):
    computer = vmtkscripts.vmtkSurfaceKiteRemoval()
    computer.Surface = surface
    computer.Execute()
    return computer.Surface

def vmtksurfmesh(surface):
    foo = vmtkscripts.vmtkSurfMesh()
    foo.Surface(surface)
    foo.Execute()
    return foo.Surface


def vmtksurfacenormals(surface):
    computer = vmtkscripts.vmtkSurfaceNormals()
    computer.Surface = surface
    computer.Execute()
    return computer.Surface


def vmtksurfaceprojection(surface, refsurface):
    """Project pointdata of refsurface onto surface"""
    projector = vmtkscripts.vmtkSurfaceProjection()
    projector.Surface = surface
    projector.ReferenceSurface = refsurface
    projector.Execute()
    return projector.Surface


def vmtksurfacereader(filename):
    reader = vmtkscripts.vmtkSurfaceReader()
    reader.InputFileName = filename
    reader.Execute()
    return reader.Surface


def vmtksurfaceremeshing(surface, edgelength=1.0, iterations=10):
    remesher = vmtkscripts.vmtkSurfaceRemeshing()
    remesher.Surface = surface
    remesher.NumberOfIterations = iterations
    remesher.ElementSizeMode = 'edgelength'
    remesher.TargetEdgeLength = edgelength
    remesher.Execute()
    return remesher.Surface

def vmtksurfaceremeshingarea(surface, iterations=10):
    remesher = vmtkscripts.vmtkSurfaceRemeshing()
    remesher.Surface = surface
    remesher.NumberOfIterations = iterations
    remesher.Execute()
    return remesher.Surface


def vmtksurfacescaling(polydata, scalefactor):
    scaler = vmtkscripts.vmtkSurfaceScaling()
    scaler.Surface = polydata
    scaler.ScaleFactor = scalefactor
    scaler.Execute()
    return scaler.Surface


def vmtksurfacesmoothing(surface, iterations=100, method='taubin'):
    smoother = vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.Method = method
    smoother.NumberOfIterations = iterations
    smoother.PassBand = 0.1
    smoother.Execute()
    return smoother.Surface


def vmtksurfacewriter(polydata, filename):
    print "Writing",filename
    writer = vmtkscripts.vmtkSurfaceWriter()
    writer.Surface = polydata
    writer.OutputFileName = filename
    writer.Execute()

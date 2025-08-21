from matplotlib import path
import vispy
from vispy.scene import SceneCanvas, visuals
from vispy.app import use_app
from vispy.visuals.filters import ShadingFilter
from scipy.spatial.transform import Rotation
from vispy.visuals.transforms.linear import *
from vispy.visuals.transforms.chain import *
from vispy.scene.visuals import InstancedMesh
from vispy import app, scene, use

ElementSymbols = ["H","He","Li","Be","B","C","N","O","F",
          "Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K",
          "Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",
          "Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y",
          "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In",
          "Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr",
          "Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
          "Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au",
          "Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",
          "Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
          "Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt",
          "Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og" ]


class MolView2:
    def __init__(self, scene):
        self.scene = scene
        self.principalAxisMesh = None
        self.bondMesh = None
        self.atomMesh = None
        self.selectedAtoms = []
        self.labelsVisual = None
        self.bonds = []
        self.atomPositions =  []
        self.atomsElement = []
        self.ellipsoids = []
        self.labelColor = (0.0,0.0,0.0)
        self.labelFontSize = 16
        self.labelFont = 'Arial Black'
    
    def setAtomsAndBonds(self, atomPositions, atomsElement, bonds = [], ellipsoids = []):
        self.atomPositions = atomPositions
        self.atomsElement = atomsElement
        self.ellipsoids = ellipsoids
        self.bonds = bonds
        self.updateMeshes()
    
    def showLabels(self, labels, atomIndices):
        radious = 100
        pos = []
        for idx in atomIndices:
            r = self.atomPositions[idx]
            pos.append([r[0]*radious,r[1]*radious,r[2]*radious])
        
        self.labelsVisual = vispy.scene.visuals.Text(labels, pos=pos, face = self.labelFont, font_size=self.labelFontSize, color=self.labelColor, parent=self.scene)

    def unsetLabels(self):
        self.labelsVisual.parent = None                
        
    def lassoSelection(self, lasso_points):
        atom_positions = []
        #for atom_view_position_3d in self.atomMesh.instance_positions:
        for atom_xyz in self.atomPositions:
            atom_xyz_view = [atom_xyz[0]*100,atom_xyz[1]*100,atom_xyz[2]*100]
            atom_position_2d = self.atomMesh.get_transform(map_from="visual", map_to="canvas").map(atom_xyz_view)
            atom_position_2d /= atom_position_2d[3:]
            atom_positions.append(atom_position_2d[:2])
        # Select vertices inside the polygon.
        polygon = path.Path(lasso_points, closed = True)
        pts = np.ndarray(shape=(len(atom_positions),2),dtype=float)
        for i in range(0,len(atom_positions)):
            pts[i,0]=atom_positions[i][0]
            pts[i,1]=atom_positions[i][1]
        polygon_mask = polygon.contains_points(pts)
        self.selectedAtoms = []
        for i in range(0,len(polygon_mask)):
            if polygon_mask[i]:
                self.selectedAtoms.append(i)
        self.updateAtomMesh()
    
    def clear(self):
        if self.principalAxisMesh is not None:
            self.principalAxisMesh.parent = None
            self.principalAxisMesh = None
        if self.bondMesh is not None:
            self.bondMesh.parent = None
            self.bondMesh = None
        if self.atomMesh is not None:
            self.atomMesh.parent = None
            self.atomMesh = None

        
    def make_shading_filter(self):
        return vispy.visuals.filters.mesh.ShadingFilter(
             shading='smooth', 
             ambient_coefficient=(1, 1, 1, 1), 
             diffuse_coefficient=(1, 1, 1, 1), 
             specular_coefficient=(1, 1, 1, 1),
             shininess=100, #light_dir=(10, 5, -5),
             ambient_light=(1, 1, 1, 0.7), 
             diffuse_light=(1, 1, 1, 0.7), 
             specular_light=(1, 1, 1, 0.25), 
             enabled=True)
    
    def vec100Rot(self, target):
        v_t = np.array(target)
        v_t_n = v_t/np.linalg.norm(v_t)
        v_x = np.array([1,0,0])
        v_y = np.array([0,1,0])
        v_1 = v_x
        if abs(np.inner(v_t, v_x))>abs(np.inner(v_t, v_y)):
            v_1 = v_y
        ort_1 = np.add(v_1,-1*np.inner(v_1, v_t_n)*v_t_n)
        ort_1 = ort_1/np.linalg.norm(ort_1)
        ort_2 = np.cross(ort_1,v_t_n)
        return [[target[0], ort_1[0], ort_2[0]],
                [target[1], ort_1[1], ort_2[1]],
                [target[2], ort_1[2], ort_2[2]]]
    
    def updateBondMeshe(self):
        if self.bondMesh is not None:
            self.bondMesh.parent = None
            self.bondMesh = None
        if not self.bonds:
            return
        #bonds
        tube = vispy.scene.visuals.Tube(radius=5, color=[1,0,0],points=[[0,0,0],[100,0,0]],tube_points=32)
        shading_filter_2 = vispy.visuals.filters.mesh.ShadingFilter(
             shading='smooth', 
             ambient_coefficient=(1, 1, 1, 1), 
             diffuse_coefficient=(1, 1, 1, 1), 
             specular_coefficient=(1, 1, 1, 1),
             shininess=100, #light_dir=(10, 5, -5),
             ambient_light=(1, 1, 1, 0.9), 
             diffuse_light=(1, 1, 1, 0.7), 
             specular_light=(1, 1, 1, 0.2), 
             enabled=True)

        tube.attach(shading_filter_2)
        
        
        vertices_t = tube.mesh_data.get_vertices()
        faces_t = tube.mesh_data.get_faces()
        n_instances = len(self.bonds)
        instance_colors_t = [[0.5,0.5,0.5]]*n_instances
        #instance_colors_t = np.random.rand(n_instances, 3).astype(np.float32)
        radious = 100
        instance_positions_t = [[ self.atomPositions[bond[0]][0]*radious,self.atomPositions[bond[0]][1]*radious,self.atomPositions[bond[0]][2]*radious] for bond in self.bonds]
        instance_transforms_t = [] 
        for bond in self.bonds:
            target = [self.atomPositions[bond[1]][0]-self.atomPositions[bond[0]][0],
                      self.atomPositions[bond[1]][1]-self.atomPositions[bond[0]][1],
                      self.atomPositions[bond[1]][2]-self.atomPositions[bond[0]][2]]
            rot   = self.vec100Rot(target)
            instance_transforms_t.append(rot)
        self.bondMesh = InstancedMesh(
                 vertices_t,
                 faces_t,
                 instance_colors=instance_colors_t,
                 instance_positions=instance_positions_t,
                 instance_transforms=instance_transforms_t,
                 parent=self.scene)
        self.bondMesh.interactive = True
        self.bondMesh.attach(shading_filter_2)

    def updateAtomMesh(self, colors=None):
        if self.atomMesh is not None:
            self.atomMesh.parent = None
            self.atomMesh = None
        element_color = {
          "H": [0.7,0.7,0.7],
          "C": [0.4,0.4,0.4],
          "N": [0.0,0.0,1.0],
          "O": [1.0,0.0,0.0]
        }
        selected_colour = [0.5,1.0,0.0]
        sph = vispy.scene.visuals.Sphere(radius=1, color=[1,0,0])
        shading_filter = self.make_shading_filter()
        
        vertices = sph.mesh.mesh_data.get_vertices()
        faces = sph.mesh.mesh_data.get_faces()
        nAtoms = len(self.atomPositions)
        
        if nAtoms == 0:
            return
        
        atom_colors = []
        if colors is None:
            for atomicNumber in self.atomsElement:
                atom_colors.append(element_color[ElementSymbols[atomicNumber-1]])
            for selectedIdx in self.selectedAtoms:
                atom_colors[selectedIdx] = selected_colour
        else:
            atom_colors = colors
        
        
        radious = 100
        atom_positions = [[r[0]*radious,r[1]*radious,r[2]*radious] for r in self.atomPositions]
        atom_transforms = [] 
        
        for u in self.ellipsoids:
            d=radious
            if len(u)>1:
                atom_transforms.append([[u[0]*u[9]*d,u[3]*u[10]*d,u[6]*u[11]*d],
                                            [u[1]*u[9]*d,u[4]*u[10]*d,u[7]*u[11]*d],
                                            [u[2]*u[9]*d,u[5]*u[10]*d,u[8]*u[11]*d]])
            elif len(u)==1:
                atom_transforms.append([[u[0]*d,0.0,0.0],
                                            [0.0,u[0]*d,0.0],
                                            [0.0,0.0,u[0]*d]])
            else:
                atom_transforms.append([[0.1*d,0.0,0.0],
                                            [0.0,0.1*d,0.0],
                                            [0.0,0.0,0.1*d]])
        self.atomMesh = InstancedMesh(
                 vertices,
                 faces,
                 instance_colors=atom_colors,
                 instance_positions=atom_positions,
                 instance_transforms=atom_transforms,
                 parent=self.scene)
        self.atomMesh.interactive = True
        self.atomMesh.attach(shading_filter)
        if self.labelsVisual is not None:
            if self.labelsVisual.parent is not None:
                self.setLabels()
    
    #def updatePrincipalAxesMesh(self, scene):
    def updatePrincipalAxesMesh(self):
        radious = 100
        instance_positions = [[r[0]*radious,r[1]*radious,r[2]*radious] for r in self.atomPositions]
        if self.principalAxisMesh is not None:
            self.principalAxisMesh.parent = None
            self.principalAxisMesh = None
        principal_axis_mesh = vispy.scene.visuals.Tube(radius=1.0, color=[0,0,0],points=[[-1.5,0,0],[1.5,0,0]],tube_points=32)
        principal_axis_vertices = principal_axis_mesh.mesh_data.get_vertices()
        principal_axis_faces = principal_axis_mesh.mesh_data.get_faces()
        
        principal_axes_positions = []
        principal_axes_transforms = []
        principal_axes_colors = [[0.1,0.1,0.1]] * 3 * len(self.ellipsoids)
        idx=0
        shading_filter_pa = vispy.visuals.filters.mesh.ShadingFilter(
             shading='smooth', 
             ambient_coefficient=(1, 1, 1, 1), 
             diffuse_coefficient=(1, 1, 1, 1), 
             specular_coefficient=(1, 1, 1, 1),
             shininess=100, #light_dir=(10, 5, -5),
             ambient_light=(1, 1, 1, 0.9), 
             diffuse_light=(1, 1, 1, 0.7), 
             specular_light=(1, 1, 1, 0.2), 
             enabled=True)

        for u in self.ellipsoids:
            d=radious*1.02
            principal_axes_transforms.append([[u[0],u[3]*u[10]*d,u[6]*u[11]*d],
                                              [u[1],u[4]*u[10]*d,u[7]*u[11]*d],
                                              [u[2],u[5]*u[10]*d,u[8]*u[11]*d]])
            principal_axes_transforms.append([[u[3],u[0]*u[9]*d,u[6]*u[11]*d],
                                              [u[4],u[1]*u[9]*d,u[7]*u[11]*d],
                                              [u[5],u[2]*u[9]*d,u[8]*u[11]*d]])
            principal_axes_transforms.append([[u[6],u[0]*u[9]*d,u[3]*u[10]*d],
                                              [u[7],u[1]*u[9]*d,u[4]*u[10]*d],
                                              [u[8],u[2]*u[9]*d,u[5]*u[10]*d]])
            principal_axes_positions.append(instance_positions[idx])
            principal_axes_positions.append(instance_positions[idx])
            principal_axes_positions.append(instance_positions[idx])
            idx += 1
        self.principalAxisMesh = InstancedMesh(
                 principal_axis_vertices,
                 principal_axis_faces,
                 instance_colors=principal_axes_colors,
                 instance_positions=principal_axes_positions,
                 instance_transforms=principal_axes_transforms,
                 parent=self.scene)
        self.principalAxisMesh.interactive = True
        self.principalAxisMesh.attach(shading_filter_pa)
        
    def updateMeshes(self):
        self.updateAtomMesh()
        self.updatePrincipalAxesMesh()
        self.updateBondMeshe()
        

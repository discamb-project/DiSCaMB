from enum import IntEnum
import copy

from matplotlib import path
import numpy as np
#from ElementData import vdvRadious

import vispy
from vispy.scene import SceneCanvas
from vispy.visuals.filters import ShadingFilter
from vispy.visuals.transforms.linear import *
from vispy.scene.visuals import InstancedMesh
from vispy import app, scene, use
#from CrystalStructurePresenter import CrystalStructurePresenter

CANVAS_SIZE = (800, 600)  # (width, height)

def defaultSelectionCallback(molCanvas3d, atoms):
    if atoms:
        molCanvas3d.selectedAtoms = atoms
        molCanvas3d.updateAtomMesh()
    else:
        if molCanvas3d.selectedAtoms:
            molCanvas3d.selectedAtoms = []
            molCanvas3d.updateAtomMesh()



class MolCanvas3d:
    def __init__(self):
        # Create canvas and view
        self.atomSelectionCallback = defaultSelectionCallback
        self.keyPressCallback = None
        self.canvas = SceneCanvas(keys='interactive', size=CANVAS_SIZE, show=True)#, bgcolor='white')
        self.view = self.canvas.central_widget.add_view()
        self.view.camera = vispy.scene.cameras.ArcballCamera(fov=0)
        self.view.camera.scale_factor = 500
        self.view.camera.interactive = True
        self.principalAxisMesh = None
        self.bondMesh = None
        self.stickBondMesh = None
        self.wireframeBondMesh = None
        self.atomMesh = None
        self.wireAtomMesh = None
        self.ellipsoidalAtomMesh = None
        self.ellipsoidAtomToAtomIdx = []
        self.wireAtomToAtomIdx = []
        
        self.atomRadious = 0.1
        self.bondRadious = 0.05
        #used only for spacefill atom style
        self.spaceFillAtomRadious = []
        
        self.atomStyle = []
        # self.allAtomStyle can survive if you change 'molecule'
        #self.allAtomStyle = AtomDisplayStyle.none
        
        self.unitCellVisuals = []
        self.selectedAtoms = []
        self.labelsVisual = None
        self.bonds = []
        self.wire_bonds = []
        self.nonBondedAtoms = []
        self.atomPositions =  np.empty(shape = [0,3])
        self.atomPositionsScaled = np.empty(shape = [0,3])
        self.atomColors = []
        self.ellipsoids = []
        self.__labels = []
        self.__labeledAtomsIndices = []
        self.selected_color = [0.0,1.0,0.0]
        self.labelColor = (0.0,1.0,1.0)
        self.labelFontSize = 16
        self.labelFont = 'Arial'
        self.scaleFactor = 100
        self.labelShift = [30, -10, 0]

        #self.crystalStructurePresenter = CrystalStructurePresenter(self)
        #self.pos=[]      
        
        @self.canvas.events.mouse_press.connect
        def on_mouse_press(event):
            clicked_meshes = self.canvas.visuals_at(event.pos)
            clicked_mesh = None
            if self.atomMesh in clicked_meshes:
                clicked_mesh = self.atomMesh
            selectedAtomIdx = -1
            if isinstance(clicked_mesh, InstancedMesh):
                pos1, min, selectedItemIdx = self.get_view_axis_in_scene_coordinates(
                    event.pos, clicked_mesh)
                if clicked_mesh == self.atomMesh:
                    selectedAtomIdx = selectedItemIdx
                    self.updateAtomMesh()
            if selectedAtomIdx > -1:
                self.onAtomPicking([selectedAtomIdx])
            else:
                self.onAtomPicking([])

        @self.canvas.events.key_press.connect
        def on_key_press(event):
            if self.keyPressCallback is not None:
                self.keyPressCallback(event.native)

        @self.canvas.events.key_release.connect
        def on_key_release(event):
            if self.keyReleaseCallback is not None:
                self.keyReleaseCallback(event.native)

    def set_canvas_color(self, rgb_int):
        self.canvas.bgcolor = vispy.color.Color([rgb_int[0]/255.0, rgb_int[1]/255.0, rgb_int[2]/255.0])
    def onResize(self):
        labels = self.__labels
        labeledAtomsIndices = self.__labeledAtomsIndices
        self.unsetLabels()
        self.showLabels(labels, labeledAtomsIndices)
        
    
    def onAtomPicking(self, atoms):
        if self.atomSelectionCallback is not None:
            self.atomSelectionCallback(self, atoms)

    def attach_headlight(self, shading_filter):
        view = self.view
        light_dir = (0, 1, 0, 0)
        shading_filter.light_dir = light_dir[:3]
        initial_light_dir = view.camera.transform.imap(light_dir)
    
        @view.scene.transform.changed.connect
        def on_transform_change(event):
            transform = view.camera.transform
            shading_filter.light_dir = transform.map(initial_light_dir)[:3]
    
    def centerCamera(self):
        newCenter = [0.0,0.0,0.0]
        if self.atomColors:
            newCenter = np.mean(self.atomPositionsScaled, 0)
        self.view.camera.center = newCenter
                    
    def get_view_axis_in_scene_coordinates(self,pos, mesh):
        event_pos = np.array([pos[0], pos[1], 0, 1])
        instances_on_canvas = []
        # Translate each position to corresponding 2d canvas coordinates
        for instance in mesh.instance_positions:
            on_canvas = mesh.get_transform(map_from="visual", map_to="canvas").map(instance)
            on_canvas /= on_canvas[3:]
            instances_on_canvas.append(on_canvas)
        
        min = 10000
        min_pos = None
        # Find the closest position to the clicked position
        for i, instance_pos in enumerate(instances_on_canvas):
            # Not minding z axis
            temp_min = np.linalg.norm(
                np.array(event_pos[:2]) - np.array(instance_pos[:2])
            )
            if temp_min < min:
                min = temp_min
                min_pos = i
        
        return instances_on_canvas, min, min_pos
    
    def removeMolecule(self):
        self.clear()
            
    def switchInteractive(self):
        self.view.camera.interactive = not self.view.camera.interactive   
    

    def __findNonbondedAtoms(self):
        self.nonBondedAtoms = []        
        bonded = []
        for bond in self.bonds:
            bonded.append(bond[0])
            bonded.append(bond[1])
        for idx in range(0, len(self.atomColors)):
            if idx not in bonded:
                self.nonBondedAtoms.append(idx)
    
    def setAtomsAndBonds(
            self, 
            atomPositions, 
            atomColors, 
            bonds = [],
            wire_bonds = None,
            wire_bond_width = 1,
            tube_bond_r = 0.1,
            ellipsoids = [],
            atomicNumbers = []):
        # remove existing labels
        self.labelsVisual = None
        self.__labels = []
        self.__labeledAtomsIndices = []
        #
        self.wire_bonds = wire_bonds
        self.atomicNumbers = atomicNumbers
        self.atomPositions = np.array(atomPositions)
        self.atomPositionsScaled = self.atomPositions * self.scaleFactor
        self.atomColors = atomColors
        self.ellipsoids = ellipsoids
        self.bonds = bonds
        self.centerCamera()
        self.__findNonbondedAtoms()
        #self.setAtomStyle(atomStyle)
        self.updateMeshes()
        
    def showUnitCell(self, show, translationVectors=[], untCellRange = [[0,1],[0,1],[0,1]]):
        
        if self.unitCellVisuals:
            for visual in self.unitCellVisuals:
                visual.parent = None
            unitCellVisuals = []
        
        if show is not True:
            return
        
        if not translationVectors:
            return
        
        a = np.array(translationVectors[0])*self.scaleFactor
        b = np.array(translationVectors[1])*self.scaleFactor
        c = np.array(translationVectors[2])*self.scaleFactor
        zero = np.array([0.0,0.0,0.0])
        #self.scaleFactor = 100
        pos_down = [zero, a,  a+b, b, zero]  
        pos_up   = [r+c for r in pos_down]
        pos_00   = [zero,c]
        pos_01   = [a, a+c]
        pos_11   = [a+b, a+b+c]
        pos_10   = [b, b+c]
        yellow = [1.0,1.0,0.0]
        color5 = [yellow,yellow,yellow,yellow,yellow]
        color2 = [yellow,yellow]
        self.unitCellVisuals = []
        self.unitCellVisuals.append(scene.visuals.Line(pos=pos_down, color=color5, width=3, parent=self.view.scene))
        self.unitCellVisuals.append(scene.visuals.Line(pos=pos_up, color=color5, width=3,parent=self.view.scene))
        self.unitCellVisuals.append(scene.visuals.Line(pos=pos_00, color=color2, width=3,parent=self.view.scene))
        self.unitCellVisuals.append(scene.visuals.Line(pos=pos_10, color=color2, width=3,parent=self.view.scene))
        self.unitCellVisuals.append(scene.visuals.Line(pos=pos_11, color=color2, width=3,parent=self.view.scene))
        self.unitCellVisuals.append(scene.visuals.Line(pos=pos_01, color=color2, width=3,parent=self.view.scene))
        

    def hideUnitCell(self, unitCellParameters):
        pass

    def redrawLabels(self):
        if self.labelsVisual is not None:
            self.labelsVisual.parent = None
        pos = []
        if self.__labeledAtomsIndices:
            pos = [self.atomPositionsScaled[i] for i in self.__labeledAtomsIndices] 
        else:
            #pos = self.atomPositionsScaled
            return
        self.labelsVisual = vispy.scene.visuals.Text(self.__labels, pos=pos, face = self.labelFont, font_size=self.labelFontSize, color=self.labelColor, parent=self.view.scene)
        chainTransform = vispy.visuals.transforms.chain.ChainTransform([copy.deepcopy(self.labelsVisual.transforms.canvas_transform),STTransform(translate=self.labelShift)])
        self.labelsVisual.transforms.canvas_transform = chainTransform
        #self.labelsVisual.transform = STTransform(translate=[20,0,0])
        
        
    def showLabels(self, labels, atomIndices=[]):
        if not labels:
            self.unsetLabels()
        self.__labels = labels
        if not atomIndices:
            self.__labeledAtomsIndices = [idx for idx in range(0, len(self.atomPositionsScaled))]
        else:
            self.__labeledAtomsIndices = atomIndices
        pos = []
        if atomIndices:
            pos = np.array([self.atomPositionsScaled[i] for i in atomIndices]) 
        else:
            pos = self.atomPositionsScaled
        if not pos.any():
            return
        self.labelsVisual = vispy.scene.visuals.Text(labels, pos=pos, face = self.labelFont, font_size=self.labelFontSize, color=self.labelColor, parent=self.view.scene)
        chainTransform = vispy.visuals.transforms.chain.ChainTransform([copy.deepcopy(self.labelsVisual.transforms.canvas_transform),STTransform(translate=self.labelShift)])
        self.labelsVisual.transforms.canvas_transform = chainTransform

    def unsetLabels(self):
        if self.labelsVisual is not None:
            self.labelsVisual.parent = None
            self.labelsVisual = None
        self.__labels = []
        self.__labeledAtomsIndices = []
        
        
    def lassoSelection(self, lasso_points):
        atom_positions = []
        
        #for atom_view_position_3d in self.atomMesh.instance_positions:
        for atom_xyz in self.atomPositionsScaled:
            atom_position_2d = self.atomMesh.get_transform(map_from="visual", map_to="canvas").map(atom_xyz)
            atom_position_2d /= atom_position_2d[3:]
            atom_positions.append(atom_position_2d[:2])
        
        # Select vertices inside the polygon.
        polygon = path.Path(lasso_points, closed = True)
        pts = np.ndarray(shape=(len(atom_positions),2),dtype=float)
        for i in range(0,len(atom_positions)):
            pts[i,0]=atom_positions[i][0]
            pts[i,1]=atom_positions[i][1]
        polygon_mask = polygon.contains_points(pts)
        selected = []
        for i in range(0,len(polygon_mask)):
            if polygon_mask[i]:
                selected.append(i)
        self.atomSelectionCallback(self, selected)
    
    def clear(self):
        self.unsetLabels()
        if self.principalAxisMesh is not None:
            self.principalAxisMesh.parent = None
            self.principalAxisMesh = None
        if self.bondMesh is not None:
            self.bondMesh.parent = None
            self.bondMesh = None
        if self.atomMesh is not None:
            self.atomMesh.parent = None
            self.atomMesh = None
        if self.stickBondMesh is not None:
            self.stickBondMesh.parent = None
            self.stickBondMesh = None
        if self.wireframeBondMesh is not None:
            self.wireframeBondMesh.parent = None
            self.wireframeBondMesh = None

        
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
    
    def updateWireframeBondMesh(self, bonds):
        if self.wireframeBondMesh is not None:
            self.wireframeBondMesh.parent = None
        if not bonds:
            return
        pos = []
        colors = []
        for bond in bonds:
            midpoint = (self.atomPositionsScaled[bond[0]]+self.atomPositionsScaled[bond[1]])*0.5
            pos.append(self.atomPositionsScaled[bond[0]])
            pos.append(midpoint)
            pos.append(midpoint)
            pos.append(self.atomPositionsScaled[bond[1]])
            color = self.selected_color if bond[0] in self.selectedAtoms else self.atomColors[bond[0]]
            colors.append(color)
            colors.append(color)
            #colors.append(self.atomColors[bond[0]])
            #colors.append(self.atomColors[bond[0]])
            color = self.selected_color if bond[1] in self.selectedAtoms else self.atomColors[bond[1]]
            colors.append(color)
            colors.append(color)
            #colors.append(self.atomColors[bond[1]])
            #colors.append(self.atomColors[bond[1]])
        self.wireframeBondMesh = vispy.scene.visuals.Line(pos=pos, color=colors, connect="segments",parent=self.view.scene)


    def updateStickBondMesh(self, bonds):
        if self.stickBondMesh is not None:
            self.stickBondMesh.parent = None
        if not bonds:
            if self.stickBondMesh is not None:
                self.stickBondMesh.parent = None
                self.stickBondMesh = None
            return
        tube = vispy.scene.visuals.Tube(radius=self.bondRadious*self.scaleFactor, color=[1,0,0],points=[[0,0,0],[self.scaleFactor,0,0]],tube_points=32)
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
        n_instances = len(bonds)
        colors = [] #[[0.5,0.5,0.5]]*n_instances
        instance_positions_t = [] #[ self.atomPositionsScaled[bond[0]] for bond in bonds]
        instance_transforms_t = [] 
        for bond in bonds:
            colors.append(self.atomColors[bond[0]] if bond[0] not in self.selectedAtoms else self.selected_color)
            colors.append(self.atomColors[bond[1]] if bond[1] not in self.selectedAtoms else self.selected_color)
            bond_vec = self.atomPositions[bond[1]]-self.atomPositions[bond[0]]
            midpoint = self.atomPositions[bond[0]] + 0.5*bond_vec
            instance_positions_t.append(self.atomPositionsScaled[bond[0]])
            instance_positions_t.append(midpoint*self.scaleFactor)
            target = self.atomPositions[bond[1]]-self.atomPositions[bond[0]]
            rot   = self.vec100Rot(bond_vec*0.5)#target)
            instance_transforms_t.append(rot)
            instance_transforms_t.append(rot)
        self.stickBondMesh = InstancedMesh(
                 vertices_t,
                 faces_t,
                 instance_colors=colors,
                 instance_positions=instance_positions_t,
                 instance_transforms=instance_transforms_t,
                 parent=self.view.scene)
        self.stickBondMesh.interactive = True
        self.stickBondMesh.attach(shading_filter_2)

    
    def updateBondMesh(self):
        if self.bondMesh is not None:
            self.bondMesh.parent = None
            self.bondMesh = None
        if not self.bonds:
            #self.updateWireframeBondMesh([])
            self.updateStickBondMesh([])
            return
        #self.updateWireframeBondMesh(self.bonds)
        #self.updateStickBondMesh(self.bonds)
        #return
        
        wireframe_bonds = []
        stick_bonds = []
        n_bonds = len(self.bonds)
        for i in range(0,n_bonds):
            if i in self.wire_bonds:
                wireframe_bonds.append(self.bonds[i])
            else:
                stick_bonds.append(self.bonds[i])
        self.updateWireframeBondMesh(wireframe_bonds)
        self.updateStickBondMesh(stick_bonds)
        
        return


    def updateWireAtomMesh(self, wireAtoms):
        if self.wireAtomMesh is not None:
            self.wireAtomMesh.parent = None
        if not wireAtoms:
            return
        colors = []
        positions = []
        edges = np.array([[0.2,0,0],[-0.2,0,0],[0,0.2,0],[0,-0.2,0],[0,0,0.2],[0,0,-0.2]])
        self.wireAtomToAtomIdx = wireAtoms
        for atomIdx in wireAtoms:
            for edge in edges:
                positions.append((self.atomPositions[atomIdx]+edge)*self.scaleFactor)
                colors.append(self.atomColors[atomIdx])
        for selected in self.selectedAtoms:
            if selected in wireAtoms:
                idx = wireAtoms.index(selected)
                for i in range(0,6):
                    colors[i+idx] = self.selected_color
        self.wireAtomMesh = vispy.scene.visuals.Line(pos=positions, color=colors, connect="segments",
                                                     parent=self.view.scene)
        
        

            
    def updateEllispoidalAtomMesh(self, sphericalAtoms, radious):
        if self.sphericalAtomMesh is not None:
            self.sphericalAtomMesh.parent = None
        
        self.ellipsoidAtomToAtomIdx = []
        ellipsoidalAtoms = []
        sphericalAtoms = []
        radious = []
        wireAtoms = []
        idx=0
            
        
        #selected_colour = [0.5,1.0,0.0]
        sph = vispy.scene.visuals.Sphere(radius=1, color=[1,0,0])
        shading_filter = self.make_shading_filter()
        self.attach_headlight(shading_filter)
        vertices = sph.mesh.mesh_data.get_vertices()
        faces = sph.mesh.mesh_data.get_faces()
        nAtoms = len(sphericalAtoms)
        if nAtoms == 0:
            return
        
        atom_colors = []#copy.deepcopy(self.atomColors)
        atom_transforms = [] 
        
        for idx,r in zip(sphericalAtoms,radious):
            atom_colors.append(self.atomColors[idx])
            atom_transforms.append( np.eye(3)*r*self.scaleFactor)
                               
        self.atomMesh = InstancedMesh(
                 vertices,
                 faces,
                 instance_colors=atom_colors,
                 instance_positions= self.atomPositionsScaled, 
                 instance_transforms=atom_transforms,
                 parent=self.view.scene)
        self.sphericalAtomMesh.interactive = True
        self.sphericalAtomMesh.attach(shading_filter)

    def updateEllipsoidalAtomMesh(self, spherical_atoms, radius, ellipsoidal_atoms):
        if not self.atomPositionsScaled.any():
            return
        if self.atomMesh is not None:
            self.atomMesh.parent = None
            self.atomMesh = None
        self.ellipsoidAtomToAtomIdx = spherical_atoms + ellipsoidal_atoms
        atom_colors = copy.deepcopy(self.atomColors)
        for selected in self.selectedAtoms:
            atom_colors[selected] = self.selected_color
        atom_transforms = []
        n_atoms = len(self.atomPositionsScaled)
        for idx in range(0, n_atoms):
            if idx in spherical_atoms:
                r = radius[spherical_atoms.index(idx)]
                atom_transforms.append(np.eye(3) * r * self.scaleFactor)
            elif idx in ellipsoidal_atoms:
                u = [self.atomRadious]
                if self.ellipsoids:
                    u = self.ellipsoids[idx]
                if len(u) > 1:
                    transform = np.array([[u[0], u[3], u[6]],
                                          [u[1], u[4], u[7]],
                                          [u[2], u[5], u[8]]]) * self.scaleFactor
                    transform[:, 0] *= u[9]
                    transform[:, 1] *= u[10]
                    transform[:, 2] *= u[11]
                    atom_transforms.append(transform)
                elif len(u) == 1:
                    atom_transforms.append(np.eye(3) * (self.scaleFactor * u[0]))
                else:
                    atom_transforms.append(self.atomRadious * np.eye(3) * self.scaleFactor)
            else:
                atom_transforms.append(np.eye(3) * 0.01 * self.scaleFactor)

        # selected_colour = [0.5,1.0,0.0]
        sph = vispy.scene.visuals.Sphere(radius=1, color=[1, 0, 0])
        shading_filter = self.make_shading_filter()
        self.attach_headlight(shading_filter)
        vertices = sph.mesh.mesh_data.get_vertices()
        faces = sph.mesh.mesh_data.get_faces()
        self.atomMesh = InstancedMesh(
                 vertices,
                 faces,
                 instance_colors=atom_colors,
                 instance_positions=self.atomPositionsScaled,
                 instance_transforms=atom_transforms,
                 parent=self.view.scene)
        self.atomMesh.interactive = True
        self.atomMesh.attach(shading_filter)

    def updateEllipsoidalAtomMesh_bck(self, sphericalAtoms, radious, ellipsoidalAtoms):
        if self.atomMesh is not None:
            self.atomMesh.parent = None
            self.atomMesh = None
        self.ellipsoidAtomToAtomIdx = sphericalAtoms + ellipsoidalAtoms
        atom_colors = []#copy.deepcopy(self.atomColors)
        atom_transforms = []
        positions = []
        
        for idx,r in zip(sphericalAtoms,radious):
            atom_colors.append(self.atomColors[idx])
            atom_transforms.append( np.eye(3)*r*self.scaleFactor)
            positions.append(self.atomPositionsScaled[idx])
        
        #selected_colour = [0.5,1.0,0.0]
        sph = vispy.scene.visuals.Sphere(radius=1, color=[1,0,0])
        shading_filter = self.make_shading_filter()
        self.attach_headlight(shading_filter)
        vertices = sph.mesh.mesh_data.get_vertices()
        faces = sph.mesh.mesh_data.get_faces()
        nAtoms = len(sphericalAtoms)+len(ellipsoidalAtoms)
        if nAtoms == 0:
            return
        
        for idx in ellipsoidalAtoms:
            positions.append(self.atomPositionsScaled[idx])
            u = self.atomRadious
            if self.ellipsoids:
                u = self.ellipsoids[idx]
            if len(u)>1:
                transform = np.array([[u[0],u[3],u[6]],
                                     [u[1],u[4],u[7]],
                                     [u[2],u[5],u[8]]])*self.scaleFactor
                transform[:,0] *= u[9]
                transform[:,1] *= u[10]
                transform[:,2] *= u[11]
                atom_transforms.append(transform)
            elif len(u)==1:
                atom_transforms.append(np.eye(3) * (self.scaleFactor * u[0]))
            else:
                atom_transforms.append(self.atomRadious * np.eye(3) * self.scaleFactor)
            atom_colors.append(self.atomColors[idx])


        self.atomMesh = InstancedMesh(
                 vertices,
                 faces,
                 instance_colors=atom_colors,
                 instance_positions= positions,#self.atomPositionsScaled,
                 instance_transforms=atom_transforms,
                 parent=self.view.scene)
        self.atomMesh.interactive = True
        self.atomMesh.attach(shading_filter)
        
    def updateAtomMesh(self, selected_colour = [0.5,1.0,0.0]):
        if self.atomMesh is not None:
            self.atomMesh.parent = None
            self.atomMesh = None
        
        ellipsoidalAtoms = []
        sphericalAtoms = []
        radious = []
        wireAtoms = []
        idx=0
        for ellipsoid in self.ellipsoids:
            if ellipsoid:
                if len(ellipsoid)>1:
                    ellipsoidalAtoms.append(idx)
                else:
                    sphericalAtoms.append(idx)
                    radious.append(ellipsoid[0])
            else:
                if idx in self.nonBondedAtoms:
                    wireAtoms.append(idx)
            idx += 1
        self.updateWireAtomMesh(wireAtoms)
        #self.sphericalAtomMesh(sphericalAtoms, radious)
        self.updateEllipsoidalAtomMesh(sphericalAtoms, radious, ellipsoidalAtoms)
                
        if self.labelsVisual is not None:
            if self.labelsVisual.parent is not None:
                self.redrawLabels()
                #self.setLabels()

    def updateAtomMeshBck(self, selected_colour = [0.5,1.0,0.0]):
        if self.atomMesh is not None:
            self.atomMesh.parent = None
            self.atomMesh = None
        
        ellipsoidalAtoms = []
        sphericalAtoms = []
        radious = []
        wireAtoms = []
        idx=0
        for style in self.atomStyle:
            if style != AtomDisplayStyle.line:
                spherical = False
                if style == AtomDisplayStyle.none or style == AtomDisplayStyle.ball_and_stick:
                    radious.append(self.atomRadious)
                if style == AtomDisplayStyle.stick:
                    radious.append(self.bondRadious)
                if style == AtomDisplayStyle.cpk:
                    if self.spaceFillAtomRadious:
                        radious.append(self.spaceFillAtomRadious[idx])
                    else:
                        radious.append(self.atomRadious)
                if style == AtomDisplayStyle.ellipsoid_and_line or style == AtomDisplayStyle.ellipsoid_and_stick:
                    if ellipsoids:
                        if ellispoid[idx]:
                            if len(ellispoid[idx])==1:
                                radious.append(ellispoid[idx][0])
                            if len(ellispoid[idx])==6:
                                nonSpherical = True
                if spherical:
                    sphericalAtoms.append(idx)
                else:
                    ellipsoidalAtoms.append(idx)
            else:
                if idx in self.nonBondedAtoms:
                    wireAtoms.append(idx)
            idx += 1
        self.updateWireAtomMesh(wireAtoms)
        self.sphericalAtomMesh(sphericalAtoms, radious)
        self.ellipsoidalAtomMesh(ellipsoidalAtoms)
            
        
        #selected_colour = [0.5,1.0,0.0]
        sph = vispy.scene.visuals.Sphere(radius=1, color=[1,0,0])
        shading_filter = self.make_shading_filter()
        self.attach_headlight(shading_filter)
        vertices = sph.mesh.mesh_data.get_vertices()
        faces = sph.mesh.mesh_data.get_faces()
        nAtoms = len(self.atomColors)
        if nAtoms == 0:
            return
        
        atom_colors = copy.deepcopy(self.atomColors)
        if selected_colour:
            for selectedIdx in self.selectedAtoms:
                atom_colors[selectedIdx] = selected_colour
        
        
        radious = 100
        atom_transforms = [] 
        
        for u in self.ellipsoids:
            d=radious
            if len(u)>1:
                transform = np.array([[u[0],u[3],u[6]],
                                     [u[1],u[4],u[7]],
                                     [u[2],u[5],u[8]]])*self.scaleFactor
                transform[:,0] *= u[9]
                transform[:,1] *= u[10]
                transform[:,2] *= u[11]
                atom_transforms.append(transform)
            elif len(u)==1:
                atom_transforms.append(np.eye(3) * (self.scaleFactor * u[0]))
            else:
                atom_transforms.append(np.eye(3) * self.scaleFactor)
                
        self.atomMesh = InstancedMesh(
                 vertices,
                 faces,
                 instance_colors=atom_colors,
                 instance_positions= self.atomPositionsScaled, 
                 instance_transforms=atom_transforms,
                 parent=self.view.scene)
        self.atomMesh.interactive = True
        self.atomMesh.attach(shading_filter)
        if self.labelsVisual is not None:
            if self.labelsVisual.parent is not None:
                self.redrawLabels()
                #self.setLabels()
    
    #def updatePrincipalAxesMesh(self, scene):
    def updatePrincipalAxesMesh(self):
        #radious = 100
        #instance_positions = [[r[0]*radious,r[1]*radious,r[2]*radious] for r in self.atomPositions]
        if self.principalAxisMesh is not None:
            self.principalAxisMesh.parent = None
            self.principalAxisMesh = None
        principal_axis_mesh = vispy.scene.visuals.Tube(radius=1.0, color=[0,0,0],points=[[-0.015*self.scaleFactor,0,0],[0.015*self.scaleFactor,0,0]],tube_points=32)
        principal_axis_vertices = principal_axis_mesh.mesh_data.get_vertices()
        principal_axis_faces = principal_axis_mesh.mesh_data.get_faces()
        
        principal_axes_positions = []
        principal_axes_transforms = []
        
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
            d=self.scaleFactor*1.02
            if len(u)>1:
                principal_axes_transforms.append([[u[0],u[3]*u[10]*d,u[6]*u[11]*d],
                                                [u[1],u[4]*u[10]*d,u[7]*u[11]*d],
                                                [u[2],u[5]*u[10]*d,u[8]*u[11]*d]])
                principal_axes_transforms.append([[u[3],u[0]*u[9]*d,u[6]*u[11]*d],
                                                [u[4],u[1]*u[9]*d,u[7]*u[11]*d],
                                                [u[5],u[2]*u[9]*d,u[8]*u[11]*d]])
                principal_axes_transforms.append([[u[6],u[3]*u[10]*d,u[0]*u[9]*d],
                                                [u[7],u[4]*u[10]*d,u[1]*u[9]*d],
                                                [u[8],u[5]*u[10]*d,u[2]*u[9]*d]])
                principal_axes_positions.append(self.atomPositionsScaled[idx])
                principal_axes_positions.append(self.atomPositionsScaled[idx])
                principal_axes_positions.append(self.atomPositionsScaled[idx])
            idx += 1
        if not principal_axes_positions:
            return
        principal_axes_colors = [[0.1,0.1,0.1]] * len(principal_axes_positions)
        self.principalAxisMesh = InstancedMesh(
                 principal_axis_vertices,
                 principal_axis_faces,
                 instance_colors=principal_axes_colors,
                 instance_positions=principal_axes_positions,
                 instance_transforms=principal_axes_transforms,
                 parent=self.view.scene)
        self.principalAxisMesh.interactive = True
        self.principalAxisMesh.attach(shading_filter_pa)
        
    def updateMeshes(self):
        self.updateAtomMesh()
        self.updateBondMesh()
        self.updatePrincipalAxesMesh()

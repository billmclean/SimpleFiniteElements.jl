import Gmsh: gmsh

gmsh.initialize()
gmsh.open("rect.geo")
gmsh.option.set_number("Mesh.ElementOrder", 2)
gmsh.model.mesh.generate(2)
gmsh.finalize()

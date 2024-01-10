import Gmsh: gmsh

gmsh.initialize()
gmsh.open("roundL.geo")
#gmsh.option.set_number("General.FontSize", 24)
gmsh.option.set_number("Geometry.Points", 0)
gmsh.option.set_number("Geometry.LabelType", 4)
gmsh.option.set_number("Geometry.CurveLabels", 1)
gmsh.option.set_number("Geometry.Surfaces", 1)
gmsh.option.set_number("Geometry.SurfaceLabels", 1)
#gmsh.model.mesh.generate(2)
if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()


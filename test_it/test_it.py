import math
import EasyFsi
from EasyFsi import*

# create FieldInfo
field_info_s = FieldInfo("TestField","m",1,FieldLocation.NodeCentered,FieldIO.OutgoingDofs)
field_info_t = FieldInfo("TestField","m",1,FieldLocation.NodeCentered,FieldIO.IncomingDofs)

# create boundary
bound_s = Boundary()
bound_t = Boundary()

bound_s.load("fe.msh")
bound_t.load("fv.msh")
bound_s.register_field(field_info_s)
bound_t.register_field(field_info_t)

# create interpolator
interp = Interpolator()
interp.add_source_boundary(bound_s)
interp.add_target_boundary(bound_t)
interp.compute_interp_coeff(InterpolationMethod.LocalXPS,20)
interp.save_coefficients("coeff.txt")

# setup field value
field_s = bound_s.get_field("TestField")
for i in range(0,bound_s.nnode):
    coord=bound_s.node_coords(i)
    field_s[i,0]=math.sin(coord.x*math.pi)*math.sin(coord.y*math.pi)

# do interpolating
interp.interp_dofs_s2t("TestField")

# save results
bound_s.save("fe.plt")
bound_t.save("fv.plt")

#exit()

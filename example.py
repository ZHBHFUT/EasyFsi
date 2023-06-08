#
# example for using EasyFsi library
#

import EasyFsi
from EasyFsi import*

#-------------------------------------------------------
# define helper functions for field reading and writing
#-------------------------------------------------------

# function used to read outgoing field of boundary.
#    app: application object
#    bd:  boundary object
#    fieldname: name of the field
#    location: location of the field, see FieldLocation
#    values: field data, type=MatView object
#    user_obj: object passed from app.exchange_solution
def get_bound_field(app,bd,fieldname,ncomp,location,values,user_obj):
    # TODO: update values
    for i in range(bd.nnode):
        values[i,0]=???; # update component-0
        values[i,1]=???; # update component-1
        #...
    # END of get_bound_field

# function used to write incoming field of boundary.
#    app: application object
#    bd:  boundary object
#    fieldname: name of the field
#    location: location of the field, see FieldLocation
#    values: field data, type=MatView object
#    user_obj: object passed from app.exchange_solution
def set_bound_field(app,bd,fieldname,ncomp,location,values,user_obj):
    # TODO: update values
    for i in range(bd.nnode):
        # ... = values[i,0]; # update component-0
        # ... = values[i,1]; # update component-1
        #...
    # END of set_bound_field

#-------------------------------------------------------
# preprocessing
#-------------------------------------------------------

# define application
app = Application("PythonSolver");

# define boundary
bd0 = app.add_coupled_boundary()
bd0.name = "bd0"
# create boundary manually:
#bd0.reserve(200,100,400)
#bd0.add_node(x,y,z,unique_id) # define node
#...
#bd0.add_face(type,nodes)  # define face
#...
# create boundary from file:
#bd0.load("???.msh")  # read Gmsh file
#bd0.compute_metics(5.0)

# define field
#    arg0  The name of this field
#    arg1  The component number of this field, >=1
#    arg2  Location of field, NodeCentered or FaceCentered
#    arg3  Input/Output type, see FieldIO
#    arg4  Units of this field
app.register_field("displacement",3,FieldLocation.NodeCentered,FieldIO.OutgoingDofs,"m")
app.register_field("force",3,FieldLocation.NodeCentered,FieldIO.IncomingLoads,"N")
app.set_field_function(get_bound_field,set_bound_field)

# define communicator between applications
#    arg0  True/False, this application is master?
#    arg1  Number of applications for this coupling problem, >=2
#    arg2  IP address of machine that master application is running.
#    arg3  Port number
#    arg4  Timeout value in second.
inter_comm = SocketCommunicator(False,2,"127.0.0.1",1234,60)

# start coupling
app.start_coupling(inter_comm)

#-------------------------------------------------------
# solving
#-------------------------------------------------------

# solving
dt = 0.001  # timestep size
nt = 1000   # timestep number
for it in range(nt)
    # TODO: solve one step
    #...
    
    # exchange solution
    app.exchange_solution(dt*(it+1), None)
    
    # other post operations
    # ...

# stop coupling when finished
app.stop_coupling()

#-------------------------------------------------------
# postprocessing
#-------------------------------------------------------

# save results
app.save_tecplot("pysolver.res.plt")

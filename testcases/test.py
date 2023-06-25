import EasyLib
from EasyLib import*

#def MyComm(Communicator):
#    def set_constant(self,name,value):
#        self.name = value
#    
#    def set_function(self,name,func):
#        self.name = func
#    
#    def rank(self):
#        return self.rank_
#    
#    def size(self):
#        return self.size_
#    
#    def send(self,data,count,dest_rank,tag):
#        
#        
#    def recv(self,data,count,src_rank,tag):
#        
#
#    def disconnect(self):
#        # do nothing

def get_bound_field(app,bound,fieldname,location,data,user_data):
    for i=0:bound.nnode
        xyz=bound.node_coords(i)
        data[i,0]=0
        data[i,1]=0
        data[i,2]=0

def set_bound_field(app,bound,fieldname,location,data,user_data):
    

inter_comm = SocketCommunicator(True,2,"127.0.0.1",1234,60)

app1 = Application("app1")
bd1 = app1.add_coupled_boundary()
app1.set_field_function(get_bound_field,set_bound_field)
app1.register_field("displacement",3,NodeCentered,OutgoingDofs,"m")
app1.register_field("force",3,NodeCentered,IncomingLoads,"N")
app1.start_coupling(inter_comm)

app1.exchange_solution(1.0,)
app1.stop_coupling()
app1.save_tecplot("app1.res.plt", True)
